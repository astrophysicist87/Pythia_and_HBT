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

		/// The default constructor.
		MyImpactParameterGenerator()
			: bMinSave(0.0), bMaxSave(0.0) {}

		/// Virtual destructor.
		virtual ~MyImpactParameterGenerator() {}

		/// Virtual init method.
		virtual bool init()
		{
			bMinSave = 0.0;
			bMaxSave = 20.0;

			return true;

		}

		/// Set the bMin and bMax (in femtometers).
		void bMin(double bMinIn) { bMinSave = bMinIn; }
		void bMax(double bMaxIn) { bMaxSave = bMaxIn; }

		/// Get bMin and bMax.
		double bMin() const { return bMinSave; }
		double bMax() const { return bMaxSave; }

		/// Return a new impact parameter and set the corresponding weight
		/// provided.
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

		virtual bool hasImpactParameterGenerator() const { cout << "Evaluated this function" << endl; return true; }

		virtual ImpactParameterGenerator * impactParameterGenerator() const { return myImpactParameterGeneratorPtr; }

	private:

		ImpactParameterGenerator * myImpactParameterGeneratorPtr;

};

//==========================================================================




vector<int> get_centrality_limits(
			const double centrality_class_lower_limit,
			const double centrality_class_upper_limit,
			const int n_events_to_use, Pythia & pythia );

void print_particle_record(
		int iEvent,
		vector<Particle> & particles_to_output,
		vector<int> & particle_is_thermal_or_decay,
		ofstream & record_stream );

int main(int argc, char *argv[])
{
	// Check number of command-line arguments.
	if (argc != 8)
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: ./main_BEeffects_arbtryParticle_OpenMP [Projectile nucleus] [Target nucleus] [Beam energy in GeV]"
				<< " [Number of events] [Results directory]"
				<< " [Lower centrality %] [Upper centrality %]" << endl;
		exit(8);
	}

	// Can easily add more later
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

	// particles to consider for HBT effects
	std::unordered_map<int, int> HBT_particle_IDs
		= { 
			{  211 , 0 },	// pion(+)
			{ -211 , 1 },	// pion(-)
			{  321 , 2 },	// Kaon(+)
			{ -321 , 3 }	// Kaon(-)
		  };
	/*std::unordered_map<int, int> HBT_particle_IDs
		= { 
			{  211 , 0 },	// pion(+)
			{ -211 , 1 }	// pion(-)
		  };*/
	/*std::unordered_map<int, int> HBT_particle_IDs
		= { 
			{  211 , 0 }	// pion(+)
		  };*/

	// thermal particles only or resonance decays included
	bool thermal_only = false;
	bool track_unshifted_particles = true;
	//if ( string(argv[-1]) == "false" )
	//	track_unshifted_particles = false;

	//if ( momentum_space_modifications )
	//	cout << "Using momentum space modifications!" << endl;
	//else
	//	cout << "Not using momentum space modifications!" << endl;

	if ( thermal_only )
		cout << "Using thermal particles only!" << endl;
	else
		cout << "Using all particles!" << endl;

	cout << "Set " << "Beams:idA = " + particle_IDs[string(argv[1])]
		<< " and " << "Beams:idB = " + particle_IDs[string(argv[2])]
		<< " and " << "Beams:eCM = " + string(argv[3]) << endl;

	//const int total_number_of_events = 100000;
	const int total_number_of_events = atoi(argv[4]);
	const int max_events_per_file = 10000;
	int current_file_index = 0;
	string file_index_string = "";
	if (total_number_of_events > max_events_per_file)
	{
		file_index_string = "_" + to_string(current_file_index);
	}

	//string systemSpecs = string(argv[1]) + string(argv[2]) + "_" + string(argv[3]) + "GeV_"
	//						+ "C" + string(argv[4]) + "_" + string(argv[5]) + "_Nev" + string(argv[4]);
	string systemSpecs = string(argv[1]) + string(argv[2]) + "_" + string(argv[3]) + "GeV_Nev" + string(argv[4]);

	string path = string(argv[5]) + "/";

	//====================================
	// files to hold actual output
	ostringstream filename_stream, mult_fn_stream;
	filename_stream //<< path
					<< systemSpecs
					<< file_index_string
					<< ".dat";
	ofstream outmain( ( path + filename_stream.str()).c_str());
	mult_fn_stream //<< path
					<< systemSpecs << "_mult"
					<< ".dat";
	ofstream outMultiplicities( (path + mult_fn_stream.str()).c_str());

	//====================================
	// in case we want separate files for unshifted particles
	// (probably for debugging purposes)
	//ofstream outmain_noShift;
	//ofstream outfilenames_noShift;
	//if ( track_unshifted_particles )
	//{
		ostringstream filename_noShift_stream;
		filename_noShift_stream //<< path
						<< systemSpecs
						<< "_unshifted"
						<< file_index_string
						<< ".dat";
		//outmain_noShift = ofstream( ( path + filename_noShift_stream.str()).c_str());
		ofstream outmain_noShift( ( path + filename_noShift_stream.str()).c_str());

		//outfilenames_noShift = ofstream(path + systemSpecs + "_unshifted_S_x_p_filenames.dat");
		ofstream outfilenames_noShift(path + systemSpecs + "_unshifted_S_x_p_filenames.dat");
		outfilenames_noShift << filename_noShift_stream.str() << endl;
	//}

	//====================================
	// meta-files to record output files and locations
	ofstream outfilenames(path + systemSpecs + "_S_x_p_filenames.dat");
	outfilenames << filename_stream.str() << endl;

	ofstream outmult_filenames(path + systemSpecs + "_total_N_filename.dat");
	outmult_filenames << mult_fn_stream.str() << endl;
	outmult_filenames.close();


	bool printing_particle_records = true;
	if ( not printing_particle_records )
	{
		outmain << "Not printing particle records!" << endl;
	}


	//int count = 0;

	// Estimate centrality class limits
	//const int n_events_to_use = 10000;
	const double centrality_class_lower_limit = atof( argv[6] );
	const double centrality_class_upper_limit = atof( argv[7] );

	cout << "Read in these centrality limits: "
			<< centrality_class_lower_limit << " to "
			<< centrality_class_upper_limit << endl;

	/*vector<int> centrality_limits
		= get_centrality_limits(
			centrality_class_lower_limit,
			centrality_class_upper_limit,
			n_events_to_use, pythia );*/
	cout << "Only minimum bias supported at the moment!  Using these centrality limits: "
			<< 0 << " to " << 1e+9 << endl;
	
	vector<int> centrality_limits = { 0, 1000000000 };
	//vector<int> centrality_limits = { 90, 110 };

	const int multiplicity_lower_limit = centrality_limits[0];
	const int multiplicity_upper_limit = centrality_limits[1];

	cout << "Accepted multiplicity range: "
			<< multiplicity_lower_limit << " to "
			<< multiplicity_upper_limit << endl;

	cout << "Running with number of OpenMP threads = "
			<< omp_get_max_threads() << endl;

	string fittedHI_SigFitDefPar = "";

	vector<Pythia> pythiaVector( omp_get_max_threads() );

	// Loop over OpenMP threads.
	for (int iThread = 0; iThread < omp_get_max_threads(); iThread++)
	{
		// Local copy of pythia object for each thread
		//Pythia pythia;

		// Turn off printing inside parallel region
		pythiaVector[iThread].readString("Print:quiet = off");

		// Seed RNG different for each thread to avoid redundant events
 		if ( omp_get_max_threads() > 1 )
		{
			pythiaVector[iThread].readString("Random:setSeed = on");
 			pythiaVector[iThread].readString("Random:seed = " + to_string(iThread+1));
		}
 		//else
 		//	pythiaVector[iThread].readString("Random:seed = -1");


		//========================================
		// Read in any standard Pythia options
		pythiaVector[iThread].readFile( "main_BEeffects.cmnd" );

		// use this to turn off energy-momentum
		// conservation, etc. for debugging purposes
		//pythia.readString("Check:event = off");

		// if true, consider BE effects only for directly (i.e., thermall) produced particles
		bool momentum_space_modifications = pythiaVector[iThread].settings.flag("HadronLevel:BoseEinstein");
		if ( thermal_only and momentum_space_modifications )
			pythiaVector[iThread].readString("BoseEinstein:widthSep = 1.0");	// allows to shift only "thermal" particles

		// Setup the beams.
		pythiaVector[iThread].readString("Beams:idA = " + particle_IDs[string(argv[1])]);
		pythiaVector[iThread].readString("Beams:idB = " + particle_IDs[string(argv[2])]);
		pythiaVector[iThread].readString("Beams:eCM = " + string(argv[3]));
		pythiaVector[iThread].readString("SoftQCD:all = on");
		//pythiaVector[iThread].readString("HardQCD:all = on");
		pythiaVector[iThread].readString("Beams:frameType = 1");

		if ( iThread == 0 )
		{
			// Initialize the Angantyr model to fit the total and semi-inclusive
			// cross sections in Pythia within some tolerance.
			pythiaVector[iThread].readString("HeavyIon:SigFitErr = "
								"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
			// These parameters are typically suitable for sqrt(S_NN)=5 TeV
			pythiaVector[iThread].readString("HeavyIon:SigFitDefPar = "
				            "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
			// A simple genetic algorithm is run for 20 generations to fit the
			// parameters.
			pythiaVector[iThread].readString("HeavyIon:SigFitNGen = 20");
		}
		else
		{
			// use results of previous fit
			pythiaVector[iThread].readString("HeavyIon:SigFitNGen = 0");
			pythiaVector[iThread].readString("HeavyIon:SigFitErr = "
								"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
			pythiaVector[iThread].readString("HeavyIon:SigFitDefPar = " + fittedHI_SigFitDefPar );
		}

		// Initialize impact parameter selection over finite, user-defined range
		MyHIUserHooks* myHIUserHooks = new MyHIUserHooks();
		pythiaVector[iThread].setHIHooks( myHIUserHooks );

		// Initialise Pythia.
		pythiaVector[iThread].init();

		vector<double> print_this = pythiaVector[iThread].settings.pvec("HeavyIon:SigFitDefPar");
		for (int iprint = 0; iprint < (int)print_this.size(); iprint++)
			cout << "CHECK: " << print_this[iprint] << endl;

		// Update fittedHI_SigFitDefPar if on first iteration.
		if ( iThread == 0 )
		{
			vector<double> result_fittedHI_SigFitDefPar
							= pythiaVector[iThread].settings.pvec("HeavyIon:SigFitDefPar");

			fittedHI_SigFitDefPar
				 = to_string( result_fittedHI_SigFitDefPar[0] );

			for (int iResult = 1; iResult < (int)result_fittedHI_SigFitDefPar.size(); iResult++)
				fittedHI_SigFitDefPar
					= fittedHI_SigFitDefPar + ","
						+ to_string( result_fittedHI_SigFitDefPar[iResult] );
			cout << "Using these settings for subsequent threads:" << endl;
			cout << "HeavyIon:SigFitNGen = 0" << endl;
			cout << "HeavyIon:SigFitDefPar = " + fittedHI_SigFitDefPar << endl;

		}

		cout << "Completed initialization of Pythia in thread# = " << iThread << endl;

	}	//end serial initialization



	// With all threads initialized,
	// set particle information for
	// HBT particles and output to files
	for (const auto & itPID : HBT_particle_IDs)
	{
		int PID = itPID.first;		// 211,-211,321,...
		int Pindex = itPID.second;	// 0,1,2,...

		ofstream out(path + "HBT_particle_" + to_string( Pindex ) + ".dat");
		out << "name = " << pythiaVector[0].particleData.name( PID )
			<< "\t\t# Particle name" << endl
			<< "monval = " << PID
			<< "\t\t# Monte-Carlo number" << endl
			<< "mass = " << pythiaVector[0].particleData.m0( PID )
			<< "\t\t# mass" << endl
			<< "charge = " << pythiaVector[0].particleData.charge( PID )
			<< "\t\t# charge" << endl
			<< "spinType = " << pythiaVector[0].particleData.spinType( PID )
			<< "\t\t# spin type" << endl
			<< "chargeType = " << pythiaVector[0].particleData.chargeType( PID )
			<< "\t\t# charge type" << endl;
		out.close();
	}


	// Loop over events and OpenMP threads.
	int iEvent = 0;
	#pragma omp parallel for
	for (int iThread = 0; iThread < omp_get_max_threads(); iThread++)
	{

		// needed below
		bool momentum_space_modifications = pythiaVector[iThread].settings.flag("HadronLevel:BoseEinstein");

		// break the loop when enough events have been generated
		bool need_to_continue = true;
		//int this_thread_count = 0;

		do
		{
			// Generate next Pythia event in this thread.
			if ( not pythiaVector[iThread].next() )
				continue;

			int event_multiplicity = 0, charged_multiplicity = 0,
				hadron_multiplicity = 0, charged_hadron_multiplicity = 0;

			// vector to hold number of each HBT particle species
			vector<int> HBT_particle_multiplicities(HBT_particle_IDs.size(), 0);

			vector<Particle> particles_to_output, unshifted_particles_to_output;
			vector<int> particle_is_thermal_or_decay;

			for (int i = 0; i < pythiaVector[iThread].event.size(); ++i)
			{
				Particle & p = pythiaVector[iThread].event[i];
				if ( p.isFinal() )
				{
					event_multiplicity++;
					if ( p.isCharged() )
						charged_multiplicity++;

					if ( p.isHadron() )
					{
						//count all final hadrons in multiplicity
						hadron_multiplicity++;
						if ( p.isCharged() )
							charged_hadron_multiplicity++;

				 		if ( HBT_particle_IDs.count( p.id() ) > 0 )	// i.e., is pion(+) or another HBT particle in the HBT_particle_IDs map
						{

							const int pmother1 = p.mother1();
							const int pmother2 = p.mother2();

							bool particle_is_decay = ( 	momentum_space_modifications
														and p.status() != 99 )	// particle is decay not affected by modifications
													or ( ( not momentum_space_modifications )
														and p.status() >= 90 );	// particle is just a normal decay product

							// if only looking at thermal (i.e., non-decay) particles
							if ( thermal_only and particle_is_decay )
								continue;

							//=================================
							// if also recording unshifted particles
							if ( momentum_space_modifications
									and track_unshifted_particles )
							{
								if ( pmother1 == pmother2 )
									unshifted_particles_to_output.push_back( pythiaVector[iThread].event[pmother1] );
								else
								{
									cerr << "Encountered a problem!  Exiting..." << endl;
									exit(8);
								}
							}

							// save the final (shifted) particles, no matter what
							particles_to_output.push_back( p );

							// Save whether the particle is thermal or decay
							// push back 0 for thermal
							// push back 1 for decay
							particle_is_thermal_or_decay.push_back( static_cast<int>( particle_is_decay ) );

							// Increment this particle count by 1
							//pion_multiplicity++;
							HBT_particle_multiplicities[ HBT_particle_IDs[ p.id() ] ]++;

						}
					}
				}
			}

			bool event_in_chosen_centrality_class
				= ( event_multiplicity >= multiplicity_lower_limit)
					and ( event_multiplicity <= multiplicity_upper_limit);

			// only save this event if N^{ch} is correct range
			if ( not event_in_chosen_centrality_class )
				continue;

			//just for now
			//if ( pion_multiplicity < 50 or pion_multiplicity > 100 )
			//	continue;

			#pragma omp critical
			{
				need_to_continue = ( iEvent < total_number_of_events );

				// only one thread in this region at a time;
				// IF need_to_continue has evaluated to false in
				// a different thread, then DO NOT print and just
				// allow this thread to terminate
				if ( need_to_continue )
				{
					//========================================
					// output physical particles here
					if ( printing_particle_records )
						print_particle_record( iEvent, particles_to_output, particle_is_thermal_or_decay, outmain );
	
					//========================================
					// output unshifted particles here in case
					// also tracking these
					if ( momentum_space_modifications
							and track_unshifted_particles
							and printing_particle_records )
						print_particle_record( iEvent, unshifted_particles_to_output, particle_is_thermal_or_decay, outmain_noShift );
	
					// If too many events for single file, set-up new file here
					if ( (iEvent + 1) % max_events_per_file == 0
						and total_number_of_events > max_events_per_file
						and iEvent + 1 < total_number_of_events )
					{
						outmain.close();
						file_index_string = "_" + to_string(++current_file_index);
						filename_stream.str("");
						filename_stream //<< path
							<< systemSpecs
							<< file_index_string
							<< ".dat";
						outmain.open(( path + filename_stream.str()).c_str());
						outfilenames << filename_stream.str() << endl;
	
					}
	
					bool verbose = true;

					//outMultiplicities
					//			<< iEvent << "   "
					//			<< event_multiplicity << "   "
					//			<< pion_multiplicity;
					outMultiplicities
								<< iEvent << "   "
								<< event_multiplicity << "   "
								<< charged_multiplicity << "   "
								<< hadron_multiplicity << "   "
								<< charged_hadron_multiplicity;

					// output particle multiplicities (order hardcoded for now)
					for ( int iHBTParticle = 0; iHBTParticle < (int)HBT_particle_IDs.size(); ++iHBTParticle )
						outMultiplicities
								<< "   " << HBT_particle_multiplicities[ iHBTParticle ];

					if ( verbose )
						outMultiplicities
								<< "   "
								<< pythiaVector[iThread].info.hiinfo->b() << "   "
								<< pythiaVector[iThread].info.hiinfo->nPartProj() << "   "
								<< pythiaVector[iThread].info.hiinfo->nPartTarg() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollTot() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollND() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollNDTot() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollSDP() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollSDT() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollDD() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollCD() << "   "
								<< pythiaVector[iThread].info.hiinfo->nCollEL() << "   "
								<< pythiaVector[iThread].info.hiinfo->nAbsProj() << "   "
								<< pythiaVector[iThread].info.hiinfo->nDiffProj() << "   "
								<< pythiaVector[iThread].info.hiinfo->nElProj() << "   "
								<< pythiaVector[iThread].info.hiinfo->nAbsTarg() << "   "
								<< pythiaVector[iThread].info.hiinfo->nDiffTarg() << "   "
								<< pythiaVector[iThread].info.hiinfo->nElTarg();

					outMultiplicities
								<< endl;

					//cout << "CHECK: " << iEvent << "   " << total_number_of_events << "   " << (iEvent < total_number_of_events) << "   " ;
	
					//need_to_continue = ( ++iEvent < total_number_of_events );
					iEvent++;

					//cout << iEvent << "   " << (iEvent < total_number_of_events) << endl;
				}
			}

		} while ( need_to_continue );	// continue until we've hit the desired number of events
	}

	outmain.close();
	outMultiplicities.close();
	outfilenames.close(); 

	// And we're done!
	return 0;
}









bool decreasing (int i,int j) { return (i>j); }


vector<int> get_centrality_limits(
			const double centrality_class_lower_limit,
			const double centrality_class_upper_limit,
			const int n_events_to_use, Pythia & pythia )
{
	int iEvent = 0;
	vector<int> event_multiplicities;

	vector<int> results(2);

	// if we're just doing all the events, no matter what
	if ( centrality_class_lower_limit < 1.e-6
			and 100.0 - centrality_class_upper_limit < 1.e-6 )
	{
		cout << "main_BEeffects(): setting default centralities" << endl;
		results[0] = 0;
		results[1] = 1e+9;
	}
	else
	{
		do
		{
			if ( !pythia.next() )
				continue;

			int event_multiplicity = 0;

			for (int i = 0; i < pythia.event.size(); ++i)
			{
				Particle & p = pythia.event[i];
				if ( p.isFinal() and p.isHadron() )
				{
					//count all final hadrons in multiplicity
					//to estimate centrality
					event_multiplicity++;
				}
			}

			event_multiplicities.push_back( event_multiplicity );

			++iEvent;

		} while ( iEvent < n_events_to_use );

		sort( event_multiplicities.begin(), event_multiplicities.end(), decreasing );

		int lower_index = max( 0, (int)floor( 0.01*centrality_class_lower_limit*n_events_to_use+0.5 ) );
		int upper_index = min( n_events_to_use - 1, (int)floor( 0.01*centrality_class_upper_limit*n_events_to_use+0.5 ) );
		results[0] = ( 100.0 - centrality_class_upper_limit < 1.e-6 ) ?
						0 : 
						event_multiplicities[upper_index];	//smaller multiplicity limit first
		results[1] = ( centrality_class_lower_limit < 1.e-6 ) ?
						1e+9 : 
						event_multiplicities[lower_index];	//larger multiplicity limit second
	}

	return ( results );
}




void print_particle_record(
		int iEvent,
		vector<Particle> & particles_to_output,
		vector<int> & particle_is_thermal_or_decay,
		ofstream & record_stream )
{

	bool OSCARformat = true;

	if ( OSCARformat )
	{
		// output this event header
		record_stream
			<< iEvent << "   "
			<< (int)particles_to_output.size()
			<< endl;

		for (int i = 0; i < (int)particles_to_output.size(); ++i)
		{
			Particle & p = particles_to_output[i];

			int thermal_or_decay = particle_is_thermal_or_decay[i];

			record_stream
				//<< iEvent << "   "				// deprecated
				//<< i << "   "						// deprecated
				<< p.id() << "   "					// Monte-Carlo ID (need not all be same particle species!)
				<< thermal_or_decay << "   "		// 0 for thermal particle; 1 for resonance decay particle 
				<< p.e() << "   "
				<< p.px() << "   "
				<< p.py() << "   "
				<< p.pz() << "   "
				<< p.tProd() << "   "
				<< p.xProd() << "   "
				<< p.yProd() << "   "
				<< p.zProd()
				<< endl;
		}
	}
	else
	{
		for (int i = 0; i < (int)particles_to_output.size(); ++i)
		{
			Particle & p = particles_to_output[i];

			record_stream
				<< iEvent << "   " << i
				<< "   " << p.e()
				<< "   " << p.px()
				<< "   " << p.py()
				<< "   " << p.pz()
				<< "   " << p.tProd()
				<< "   " << p.xProd()
				<< "   " << p.yProd()
				<< "   " << p.zProd()
				<< endl;
		}
	}

	return;
}

// End of file
