// main_BEeffects.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
using namespace std;

//==========================================================================


class MyImpactParameterGenerator : public ImpactParameterGenerator
{

	public:

		/// The default constructor takes a general width (in femtometers) as
		/// argument.
		MyImpactParameterGenerator()
			: bMinSave(0.0), bMaxSave(0.0) {}

		/// Virtual destructor.
		virtual ~MyImpactParameterGenerator() {}

		/// Virtual init method.
		virtual bool init()
		{
			bMinSave = 0.0;		// the minimum impact parameter you want
			bMaxSave = 16.0;	// the maximum impact parameter you want

			return true;

		}

		/// Set the bMin and bMax (in femtometers).
		void bMin(double bMinIn) { bMinSave = bMinIn; }
		void bMax(double bMaxIn) { bMaxSave = bMaxIn; }

		/// Get bMin and bMax.
		double bMin() const { return bMinSave; }
		double bMax() const { return bMaxSave; }

		/// Return a new impact parameter (sampled uniformly)
		// and set the corresponding weight provided.
		virtual Vec4 generate(double & weight) const
		{
			double b = bMin() + ( bMax() - bMin() ) * rndPtr->flat();
			double phi = 2.0*M_PI*rndPtr->flat();
			weight = M_PI*( bMax()*bMax() - bMin()*bMin() );
			return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
		}

	private:

		/// The max and min of the distribution.
		double bMinSave;
		double bMaxSave;

};



class MyHIUserHooks : public HIUserHooks
{

	public:

		// Constructor creates impact parameter generator.
		MyHIUserHooks() { myImpactParameterGeneratorPtr = new MyImpactParameterGenerator(); }

		// Destructor deletes impact parameter generator.
		~MyHIUserHooks() { delete myImpactParameterGeneratorPtr; }

		// need this method to tell Pythia that you have defined
		// your own generator to replace the default
		virtual bool hasImpactParameterGenerator() const { cout << "Evaluated this function" << endl; return true; }

		// need this method to return pointer to your generator
		virtual ImpactParameterGenerator * impactParameterGenerator() const { return myImpactParameterGeneratorPtr; }

	private:

		// your private generator pointer
		ImpactParameterGenerator * myImpactParameterGeneratorPtr;

};

//==========================================================================


int main(int argc, char *argv[])
{
	// Check number of command-line arguments.
	if (argc != 2)
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: ./main_template [Number of events]" << endl;
		exit(8);
	}

	//=====================================================
	// A little dictionary to convert particle names into ID #s
	std::unordered_map<string, string> particle_IDs
		= {
			{ "p"  , "2212"       },
			{ "d"  , "1000010020" },
			{ "t"  , "1000010030" },
			{ "He3", "1000020030" },
			{ "He4", "1000020040" },
			{ "Li" , "1000030060" },
			{ "C"  , "1000060120" },
			{ "O"  , "1000080160" },
			{ "Cu" , "1000290630" },
			{ "Xe" , "1000541290" },
			{ "Au" , "1000791970" },
			{ "Pb" , "1000822080" },
			{ "U"  , "1000922380" }
		  };

	//=====================================================
	// Sanity check
	cout << "Running with number of OpenMP threads = "
			<< omp_get_max_threads() << endl;

	//=====================================================
	// Need this to hold a separate Pythia instance
	// for computations in each separate thread
	vector<Pythia> pythiaVector( omp_get_max_threads() );

	//=====================================================
	// Loop serially over OpenMP threads to initialize.
	for (int iThread = 0; iThread < omp_get_max_threads(); iThread++)
	{

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// N.B. - THIS FOR LOOP JUST INITIALIZES THE PYTHIA
		//        VECTOR NEEDED FOR THE PARALLEL COMPUTATIONS.
		//        WE ARE NOT YET IN A PARALLEL REGION
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		//=====================================================
		// Turn off printing
		pythiaVector[iThread].readString("Print:quiet = off");

		// Seed RNG different for each thread to avoid redundant events
		pythiaVector[iThread].readString("Random:setSeed = on");
		pythiaVector[iThread].readString("Random:seed = " + to_string(iThread-1));

		//=====================================================
		// Setup the beams.
		pythiaVector[iThread].readString("Beams:idA = " + particle_IDs["Pb"]);
		pythiaVector[iThread].readString("Beams:idB = " + particle_IDs["Pb"]);
		pythiaVector[iThread].readString("Beams:eCM = 2760");
		pythiaVector[iThread].readString("SoftQCD:all = on");
		pythiaVector[iThread].readString("Beams:frameType = 1");

		//=====================================================
		// Initialize the Angantyr model to fit the total and semi-inclusive
		// cross sections in Pythia within some tolerance.
		pythiaVector[iThread].readString("HeavyIon:SigFitErr = "
							"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
		// These parameters are typically suitable for sqrt(S_NN)=5 TeV
		pythiaVector[iThread].readString("HeavyIon:SigFitDefPar = "
		                "13.91,1.78,0.22,0.0,0.0,0.0,0.0,0.0");
		// A simple genetic algorithm is run for 20 generations to fit the
		// parameters.
		pythiaVector[iThread].readString("HeavyIon:SigFitNGen = 0");


		//=====================================================
		// Initialize impact parameter selection over finite, user-defined range
		MyHIUserHooks* myHIUserHooks = new MyHIUserHooks();
		pythiaVector[iThread].setHIHooks( myHIUserHooks );

		//=====================================================
		// Initialise Pythia here.
		pythiaVector[iThread].init();

		cout << "Completed initialization of Pythia for thread# = " << iThread << endl;


	}	//end serial initialization

	//======================================
	// Total number of events we want to generate
	const int total_number_of_events = atoi(argv[1]);

	cout << "Output: [EventID] [N^{ch}] [N^{\\pi^+}]" << endl;

	//======================================
	// The event we're currently processing
	int iEvent = 0;
	//======================================
	// Loop over events and OpenMP threads.
	// N.B. - NOW WE ARE REALLY IN A PARALLEL REGION
	#pragma omp parallel for
	for (int iThread = 0; iThread < omp_get_max_threads(); iThread++)
	{

		//======================================
		// break the loop when enough events have been generated
		bool need_to_continue = true;

		do
		{
			//======================================
			// Get the next event in this thread.
			if ( not pythiaVector[iThread].next() )
				continue;

			int event_multiplicity = 0;
			int pion_multiplicity = 0;

			for (int i = 0; i < pythiaVector[iThread].event.size(); ++i)
			{
				Particle & p = pythiaVector[iThread].event[i];
				if ( p.isFinal() and p.isHadron() )
				{
					//count all final hadrons in multiplicity
					event_multiplicity++;

					//count all final pi^+s
			 		if ( p.id() == 211 )
						pion_multiplicity++;
				}
			}

			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// N.B. - NO OTHER THREADS CAN ENTER THIS BLOCK
			//        UNTIL THE CURRENT THREAD FINISHES
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			#pragma omp critical
			{
				// Check if this thread can stop generating events (once we've gotten enough)
				need_to_continue = ( iEvent < total_number_of_events );

				// only one thread in this region at a time;
				// IF need_to_continue has evaluated to false in
				// a different thread, then DO NOT print and just
				// allow this thread to terminate
				if ( need_to_continue )
				{
					//========================================
					// print output and increment the event index
					cout << "Output: "
						 << iEvent++ << "   "
						 << event_multiplicity << "   "
						 << pion_multiplicity << endl;
				}
			}

		} while ( need_to_continue );	// continue until we've hit the desired number of events
	}

	// And we're done!
	return 0;
}


// End of file
