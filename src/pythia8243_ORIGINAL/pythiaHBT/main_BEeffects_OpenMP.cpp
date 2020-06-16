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

#include "./include/ParameterReader.h"

using namespace Pythia8;
using namespace std;

//==========================================================================

constexpr double bMinDefault = 0.0;
constexpr double bMaxDefault = 20.0;

class MyImpactParameterGenerator : public ImpactParameterGenerator
{

	public:

		/// The default constructor.
		MyImpactParameterGenerator()
			: bMinSave(bMinDefault), bMaxSave(bMaxDefault) {}

		/// A different constructor.
		MyImpactParameterGenerator(double bMinIn, double bMaxIn)
			: bMinSave(bMinIn), bMaxSave(bMaxIn) {}

		/// Virtual destructor.
		virtual ~MyImpactParameterGenerator() {}

		/// Virtual init method.
		virtual bool init()
		{
			//bMinSave = 10.0;
			//bMaxSave = 10.001;

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

		// Constructors create impact parameter generator.
		MyHIUserHooks()
		{
			myImpactParameterGeneratorPtr = new MyImpactParameterGenerator();
		}
		MyHIUserHooks(double bMinIn, double bMaxIn)
		{
			myImpactParameterGeneratorPtr = new MyImpactParameterGenerator(bMinIn, bMaxIn);
		}

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
		ofstream & record_stream,
		double ev_b_phi,
		bool include_tau_and_eta );

int main(int argc, char *argv[])
{
	// Check number of command-line arguments.
	/*if (argc != 9 and argc != 12)
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: " << endl
				<< "    ./main_BEeffects_OpenMP"
				<< " [Projectile nucleus] [Target nucleus]"
				<< " [Beam energy in GeV] [Number of events]"
				<< " [HBT particle ID] [Results directory]"
				<< " [Lower centrality %] [Upper centrality %]" << endl
				<< "  <<< OR >>> " << endl
				<< "    ./main_BEeffects_OpenMP"
				<< " [Projectile nucleus] [Target nucleus]"
				<< " [Beam energy in GeV] [Number of events]"
				<< " [HBT particle ID] [Results directory]"
				<< " [Lower centrality %] [Upper centrality %]" 
				<< " [Minimum b] [Maximum b] [Output Bjorken variables]" << endl;
				
		exit(8);
	}*/
	if ( argc < 4 )
	{
		cerr << "Incorrect number of arguments!" << endl;
		cerr << "Usage: " << endl
				<< "    ./main_BEeffects_OpenMP"
				<< " [Projectile nucleus] [Target nucleus]"
				<< " [Results directory] [[Miscellaneous options]]"
				<< endl;
		exit(8);
	}

	ParameterReader * paraRdr = new ParameterReader;

	paraRdr->readFromFile( "./parameters.dat" );
	paraRdr->readFromArguments( argc, argv, (string)("#"), 4 );
	paraRdr->echo();

	// Initialize needed command-line arguments and parameters here
	string projectile_name 				= string(argv[1]);
	string target_name 					= string(argv[2]);
	string results_directory			= string(argv[3]);

	const int dataset_seed 				= paraRdr->getVal("pythiaHBT::seed");

	int chosen_HBT_particle_ID 			= paraRdr->getVal("pythiaHBT::chosenParticleID");
	const int total_number_of_events 	= paraRdr->getVal("pythiaHBT::number_of_events");
	const int max_events_per_file 		= paraRdr->getVal("pythiaHBT::max_events_per_file");
	double beam_energy 					= paraRdr->getVal("pythiaHBT::beam_energy");
	int event_selection_mode			= paraRdr->getVal("pythiaHBT::event_selection_mode");
	double event_class_lowerLimit		= paraRdr->getVal("pythiaHBT::lowerLimit");
	double event_class_upperLimit		= paraRdr->getVal("pythiaHBT::upperLimit");
	double b_min						= paraRdr->getVal("pythiaHBT::bmin");
	double b_max						= paraRdr->getVal("pythiaHBT::bmax");

	// thermal particles only or resonance decays included
	bool thermal_only 					= static_cast<bool>( paraRdr->getVal("pythiaHBT::thermal_only") );
	bool track_unshifted_particles 		= static_cast<bool>( paraRdr->getVal("pythiaHBT::track_unshifted_particles") );
	bool store_Bjorken_coordinates 		= static_cast<bool>( paraRdr->getVal("pythiaHBT::output_Bjorken_variables") );
	bool printing_particle_records 		= static_cast<bool>( paraRdr->getVal("pythiaHBT::printing_particle_records") );


	// Use this to convert boolean variables
	// to appropriate flag toggles in Pythia.
	std::unordered_map<bool, string> boolean_toggle
		= {
			{ true,  "on"  },
			{ false, "off" }
		  };

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

	// particle(s) to consider for HBT effects
	/*std::unordered_map<int, int> HBT_particle_IDs
		= { 
			{  211 , 0 },	// pion(+)
			{ -211 , 1 },	// pion(-)
			{  321 , 2 },	// Kaon(+)
			{ -321 , 3 }	// Kaon(-)
		  };*/
	std::unordered_map<int, int> HBT_particle_IDs;
	HBT_particle_IDs.insert ( { { chosen_HBT_particle_ID, 0 } } );

	// Some Bose-Einstein options (add to *.cmnd file eventually)
	bool useInvariantSourceSize             = false;     // Lorentz-invariant size vs. spatial size only
	bool useDistribution                    = false;     // Estimate QRef vs. take as input parameter
	bool useRelativeDistance                = true;     // Use relative distances or absolute sizes
	bool useRestFrame 						= true;		// Use rest frame vs. lab frame
	bool includePhaseSpace					= true;		// Include phase-space factor
	bool linearInterpolateCDF				= false;     // Estimate pair density via linear interpolation
	bool computeBEEnhancementExactly		= false;     // Whether to evaluate BE enhancement approximately or exactly

	int shiftingSet      = 1;
	int compensationSet  = 0;
	int compensationMode = 1;

	//if ( momentum_space_modifications )
	//	cout << "Using momentum space modifications!" << endl;
	//else
	//	cout << "Not using momentum space modifications!" << endl;

	if ( thermal_only )
		cout << "Using thermal particles only!" << endl;
	else
		cout << "Using all particles!" << endl;

	cout << "Set " << "Beams:idA = " + particle_IDs[projectile_name]
		<< " and " << "Beams:idB = " + particle_IDs[target_name]
		<< " and " << "Beams:eCM = " << beam_energy << endl;

	int current_file_index = 0;
	string file_index_string = "";
	if (total_number_of_events > max_events_per_file)
	{
		file_index_string = "_" + to_string(current_file_index);
	}

	string systemSpecs = projectile_name + target_name + "_" + to_string(int(beam_energy)) + "GeV_Nev" + to_string(total_number_of_events);
	string path = results_directory + "/";


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


	if ( not printing_particle_records )
	{
		outmain << "Not printing particle records!" << endl;
	}


	//int count = 0;

	int multiplicity_lower_limit = -1;
	int multiplicity_upper_limit = -1;

	// select by centrality class, convert to multiplicity range
	if ( event_selection_mode == 0 )
	{
		// Estimate centrality class limits
		//const int n_events_to_use = 10000;
		const double centrality_class_lower_limit = event_class_lowerLimit;
		const double centrality_class_upper_limit = event_class_upperLimit;
	
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
	
		multiplicity_lower_limit = centrality_limits[0];
		multiplicity_upper_limit = centrality_limits[1];
	}
	// Set multiplicity range directly
	else if ( event_selection_mode == 1)
	{
		multiplicity_lower_limit = event_class_lowerLimit;
		multiplicity_upper_limit = event_class_upperLimit;
	}
	else
	{
		cerr << "event_selection_mode == " << event_selection_mode << " not supported!  Exiting..." << endl;
		exit(8);
	}

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
		auto & pythia = pythiaVector[iThread];

		// ==============================================
		// Turn off printing inside parallel region
		if ( omp_get_max_threads() > 1 ) pythia.readString("Print:quiet = off");


		// ==============================================
		// Seed RNG different for each thread to avoid redundant events
		bool use_random_time = false;
 		if ( omp_get_max_threads() > 1 or dataset_seed > 0 )
		{
			const int seed_for_this_thread = dataset_seed*omp_get_max_threads()
												+ iThread + 1;
			pythia.readString("Random:setSeed = on");
 			pythia.readString("Random:seed = " + to_string(seed_for_this_thread));
		}
 		else if ( use_random_time )
		{
			pythia.readString("Random:setSeed = on");
 			pythia.readString("Random:seed = 0");
		}
 		/*else
		{
			pythia.readString("Random:setSeed = on");
 			pythia.readString("Random:seed = -1");
		}*/


		//----------------------------------------------
		// N.B.: THESE OPTIONS MAY BE OVERWRITTEN BY
		//       OPTIONS PASSED IN THROUGH *.CMND FILES.
		/*pythia.readString("BoseEinstein:useInvariantSourceSize = "      + boolean_toggle[useInvariantSourceSize]);
		pythia.readString("BoseEinstein:useDistribution = "             + boolean_toggle[useDistribution]);
		pythia.readString("BoseEinstein:useRelativeDistance = "         + boolean_toggle[useRelativeDistance]);
		pythia.readString("BoseEinstein:useRestFrame = "                + boolean_toggle[useRestFrame]);
		pythia.readString("BoseEinstein:includePhaseSpace = "           + boolean_toggle[includePhaseSpace]);
		pythia.readString("BoseEinstein:linearInterpolateCDF = "        + boolean_toggle[linearInterpolateCDF]);
		pythia.readString("BoseEinstein:computeBEEnhancementExactly = " + boolean_toggle[computeBEEnhancementExactly]);

		pythia.readString("BoseEinstein:shiftingSet = "      + std::to_string(shiftingSet));
		pythia.readString("BoseEinstein:compensationSet = "  + std::to_string(compensationSet));
		pythia.readString("BoseEinstein:compensationMode = " + std::to_string(compensationMode));*/


		// For now.
		//pythia.readString("BoseEinstein:Kaon = off");
		//pythia.readString("BoseEinstein:Eta = off");


		//========================================
		// Read in any standard Pythia options
		//pythia.readFile( "main_BEeffects.cmnd" );
		pythia.readFile( path + "main_BEeffects.cmnd" );
		cout << "main_BEeffects_OpenMP: loading " << path + "main_BEeffects.cmnd..." << endl;


		// ==============================================
		// use this to turn off energy-momentum
		// conservation, etc. for debugging purposes
		//pythia.readString("Check:event = off");


		// ==============================================
		// Further BE options
		// if true, consider BE effects only for directly
		// (i.e., thermally) produced particles
		bool momentum_space_modifications = pythia.settings.flag("HadronLevel:BoseEinstein");
		if ( thermal_only and momentum_space_modifications )
			pythia.readString("BoseEinstein:widthSep = 1.0");	// allows to shift only "thermal" particles


		// ==============================================
		// Setup the beams.
		pythia.readString("Beams:idA = " + particle_IDs[ projectile_name ]);
		pythia.readString("Beams:idB = " + particle_IDs[ target_name ]);
		pythia.readString("Beams:eCM = " + to_string( beam_energy ));
		pythia.readString("Beams:frameType = 1");


		// ==============================================
		// Setup the processes.
		pythia.readString("SoftQCD:all = on");
		//pythia.readString("SoftQCD:nonDiffractive = on");
		//pythia.readString("PhaseSpace:pTHatMax = 10.0");
		//pythia.readString("HardQCD:all = on");
		//pythia.readString("Angantyr:CollisionModel = 0");


		// ==============================================
		if ( iThread == 0 )
		{
			// just for the moment

			// Initialize the Angantyr model to fit the total and semi-inclusive
			// cross sections in Pythia within some tolerance.
			pythia.readString("HeavyIon:SigFitErr = "
								"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
			// These parameters are typically suitable for sqrt(S_NN)=5 TeV
			pythia.readString("HeavyIon:SigFitDefPar = "
				            "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
			// A simple genetic algorithm is run for 20 generations to fit the
			// parameters.
			pythia.readString("HeavyIon:SigFitNGen = 20");

			/*pythia.readString("HeavyIon:SigFitDefPar = "
							"13.91,1.78,0.22,0.0,0.0,0.0,0.0,0.0" );
			pythia.readString("HeavyIon:SigFitErr = "
								"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
			pythia.readString("HeavyIon:SigFitNGen = 0");*/
		}
		else
		{
			// use results of previous fit
			pythia.readString("HeavyIon:SigFitNGen = 0");
			pythia.readString("HeavyIon:SigFitErr = "
								"0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
			pythia.readString("HeavyIon:SigFitDefPar = "
												+ fittedHI_SigFitDefPar );
		}

		// ==============================================
		// Cut out excessively long-lived particles.
		bool omit_longlived_decays = false;
		if ( omit_longlived_decays )
		{
			///*
			//const double longest_permitted_proper_lifetime = 100.0; //fm/c
			//const double longest_permitted_radius = 100.0; //fm
			pythia.readString("ParticleDecays:limitTau0 = on");
			if ( thermal_only )
				pythia.readString("ParticleDecays:tau0Max = 0.0");
			else
				pythia.readString("ParticleDecays:tau0Max = 0.000001");
			//pythia.readString("331:mayDecay = off");		//eta'
			//pythia.readString("ParticleDecays:limitTau = on");
			//pythia.readString("ParticleDecays:tauMax = 1.0");
			//pythia.readString("ParticleDecays:limitRadius = on");
			//pythia.readString("ParticleDecays:rMax = " + to_string( 1.e-12*longest_permitted_radius ) );
			//pythia.readString("ParticleDecays:limitCylinder = off");
			//pythia.readString("ParticleDecays:xyMax");
			//pythia.readString("ParticleDecays:zMax");
			//*/
			/*
			// try turning off specific decays instead
			pythia.readString("310:mayDecay = off");		//K_S0
			pythia.readString("3112:mayDecay = off");	//Sigma-
			pythia.readString("3122:mayDecay = off");	//Lambda0
			pythia.readString("3222:mayDecay = off");	//Sigma+
			pythia.readString("3312:mayDecay = off");	//Xi-
			//pythia.readString("3322:mayDecay = off");	//Xi0
			pythia.readString("3334:mayDecay = off");	//Omega-
			// next tier down in decay length
			pythia.readString("411:mayDecay = off");		//D+
			pythia.readString("421:mayDecay = off");		//D0
			pythia.readString("431:mayDecay = off");		//D_s+
			pythia.readString("511:mayDecay = off");		//B0
			pythia.readString("521:mayDecay = off");		//B+
			pythia.readString("531:mayDecay = off");		//B_s0
			pythia.readString("541:mayDecay = off");		//B_c+
			*/
		}


		// ==============================================
		// Initialize impact parameter selection over
		// finite, user-defined range.
		MyHIUserHooks* myHIUserHooks;
		myHIUserHooks = new MyHIUserHooks( b_min, b_max );
		pythia.setHIHooks( myHIUserHooks );


		// ==============================================
		// Initialise Pythia.
		pythia.init();


		// ======================================================
		// Update fittedHI_SigFitDefPar if on first iteration.
		// NOTE: only constructs new fit string, does NOT execute fit!!!
		if ( iThread == 0 )
		{
			vector<double> result_fittedHI_SigFitDefPar
							= pythia.settings.pvec("HeavyIon:SigFitDefPar");

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

		// ==============================================
		// Print progress.
		cout << "Completed initialization of Pythia in thread# = " << iThread << endl;

	}	//end serial initialization



	// ==============================================
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
		auto & pythia = pythiaVector[iThread];

		// needed below
		bool momentum_space_modifications = pythia.settings.flag("HadronLevel:BoseEinstein");

		// break the loop when enough events have been generated
		bool need_to_continue = true;
		const int seed_for_this_thread = dataset_seed*omp_get_max_threads()
												+ iThread + 1;
		//int this_thread_count = 0;

		do
		{
			// Generate next Pythia event in this thread.
			if ( not pythia.next() )
			{
				cout << "Event generation failed for some reason!  Continuing!" << endl;
exit(8);
				continue;
			}
			//else
			//	cout << "Event generation was successful!" << endl;


			// If using BE shifts, check that this event
			// was shifted successfully.
			/*if ( momentum_space_modifications
					and not pythia.info.hasBECShifts() )
			{
				cout << "WARNING: event not shifted correctly!  Trying again..." << endl;
				cout << "WARNING: momentum_space_modifications = " << momentum_space_modifications << endl;
				cout << "WARNING: pythia.info.hasBECShifts() = " << pythia.info.hasBECShifts() << endl;
				continue;
			}*/

if (false)
{
	cout << "START EVENT LISTING" << endl;
	pythia.event.list(true, true, 8);
	cout << "END EVENT LISTING" << endl;
}

			int event_multiplicity = 0;
			int charged_multiplicity = 0;
			int hadron_multiplicity = 0;
			int charged_hadron_multiplicity = 0;    // N_{ch}
			int Nch_absEta_lt_0_5 = 0;              // N_{ch}|_{|\eta|<=0.5}
			int Nch_absEta_lt_1_2 = 0;              // N_{ch}|_{|\eta|<=1.2}

			// vector to hold number of each HBT particle species
			vector<int> HBT_particle_multiplicities(HBT_particle_IDs.size(), 0);
			vector<int> HBT_particle_multiplicities_absEta_lt_0_5(HBT_particle_IDs.size(), 0);
			vector<int> HBT_particle_multiplicities_absEta_lt_1_2(HBT_particle_IDs.size(), 0);

			vector<Particle> particles_to_output, unshifted_particles_to_output;
			vector<int> particle_is_thermal_or_decay;

			for (int i = 0; i < pythia.event.size(); ++i)
			{
				Particle & p = pythia.event[i];
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
						{
							charged_hadron_multiplicity++;
							const double abs_p_eta = abs(p.eta());
							if ( abs_p_eta < 1.2 )
							{
								Nch_absEta_lt_1_2++;
								if ( abs_p_eta < 0.5 )
									Nch_absEta_lt_0_5++;
							}
						}

						// i.e., is pion(+) or another HBT particle in the HBT_particle_IDs map
				 		if ( HBT_particle_IDs.count( p.id() ) > 0 )
						{

							const int pmother1 = p.mother1();
							const int pmother2 = p.mother2();

							// particle is decay not affected by modifications
							bool particle_is_decay_withMSM
                                 = momentum_space_modifications
                                   and p.status() != 99;

							// particle is just a normal decay product
							bool particle_is_decay_withoutMSM
                                 = ( not momentum_space_modifications )
                                   and p.status() >= 90;

							bool particle_is_decay
                                 = particle_is_decay_withMSM
                                   or particle_is_decay_withoutMSM;	

							// if only looking at thermal (i.e., non-decay) particles
							if ( thermal_only and particle_is_decay )
								continue;

							//=================================
							// if also recording unshifted particles
							if ( momentum_space_modifications
									and track_unshifted_particles )
							{
								if ( pmother1 == pmother2 )
									unshifted_particles_to_output.push_back( pythia.event[pmother1] );
								/*else
								{
									cerr << "Encountered a problem!  Exiting..." << endl;
									exit(8);
								}*/
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
							const double abs_p_eta = abs(p.eta());
							if ( abs_p_eta < 1.2 )
							{
								HBT_particle_multiplicities_absEta_lt_1_2[ HBT_particle_IDs[ p.id() ] ]++;
								if ( abs_p_eta < 0.5 )
									HBT_particle_multiplicities_absEta_lt_0_5[ HBT_particle_IDs[ p.id() ] ]++;
							}
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

			// ignore events which don't generate particles!
			if ( event_multiplicity < 3 )
				continue;

			//just for now
			//const int pion_multiplicity = HBT_particle_multiplicities[ HBT_particle_IDs[ 211 ] ];
			//if ( pion_multiplicity < 50 or pion_multiplicity > 100 )
			//	continue;
			// just pick something to guarantee large multiplicity
			//if ( event_multiplicity < 70000 )
			//	continue;

			// Apply any additional cuts here
			const int pion_multiplicity_absEta_lt_1_2
                      = HBT_particle_multiplicities_absEta_lt_1_2[ HBT_particle_IDs[ 211 ] ];
			const double dNch_deta_absEta_lt_0_5 = Nch_absEta_lt_0_5 / 1.0;   // dN_{ch}/d\eta|_{|\eta|<=0.5}
			const double dNch_deta_absEta_lt_1_2 = Nch_absEta_lt_1_2 / 2.4;   // dN_{ch}/d\eta|_{|\eta|<=1.2}

			// only keep events with at least two charged pions with |eta|<=1.2
			// (only one will not contribute to final results)
			if ( pion_multiplicity_absEta_lt_1_2 < 2 ) continue;


			#pragma omp critical
			{
				need_to_continue = ( iEvent < total_number_of_events );

				// only one thread in this region at a time;
				// IF need_to_continue has evaluated to false in
				// a different thread, then DO NOT print and just
				// allow this thread to terminate
				if ( need_to_continue )
				{
					//const double phi_to_use = ( particle_IDs[ projectile_name ] == "p"
					//							and particle_IDs[ target_name ] == "p" ) ?
					//							0.0 : pythia.info.hiinfo->phi();
					const double phi_to_use = 0.0;

					//========================================
					// output physical particles here
					if ( printing_particle_records )
						print_particle_record( dataset_seed * total_number_of_events + iEvent,
												particles_to_output,
												particle_is_thermal_or_decay, outmain,
												phi_to_use, store_Bjorken_coordinates );
	
					//========================================
					// output unshifted particles here in case
					// also tracking these
					if ( momentum_space_modifications
							and track_unshifted_particles
							and printing_particle_records )
						print_particle_record( dataset_seed * total_number_of_events + iEvent,
												unshifted_particles_to_output,
												particle_is_thermal_or_decay, outmain_noShift,
												phi_to_use, store_Bjorken_coordinates );
	
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
	
					bool verbose = false;

					outMultiplicities
								<< dataset_seed * total_number_of_events + iEvent << "   "
								<< dataset_seed << "   "
								<< iEvent << "   "
								<< charged_hadron_multiplicity;

					// output particle multiplicities (order hardcoded for now)
					for ( int iHBTParticle = 0; iHBTParticle < (int)HBT_particle_IDs.size(); ++iHBTParticle )
						outMultiplicities
								<< "   " << HBT_particle_multiplicities[ iHBTParticle ];

					outMultiplicities
								<< "   " << Nch_absEta_lt_1_2
								<< "   " << dNch_deta_absEta_lt_0_5;

					// output particle multiplicities with |eta|<=1.2
					for ( int iHBTParticle = 0; iHBTParticle < (int)HBT_particle_IDs.size(); ++iHBTParticle )
						outMultiplicities
								<< "   " << HBT_particle_multiplicities_absEta_lt_1_2[ iHBTParticle ];

					if ( verbose )
						outMultiplicities
								<< "   "
								<< pythia.info.hiinfo->b() << "   "
								<< phi_to_use;
								//<< pythia.info.hiinfo->phi();// << "   "
								/*<< pythia.info.hiinfo->nPartProj() << "   "
								<< pythia.info.hiinfo->nPartTarg() << "   "*/
								//<< pythia.info.hiinfo->nCollTot() << "   "
								/*<< pythia.info.hiinfo->nCollND() << "   "
								<< pythia.info.hiinfo->nCollNDTot() << "   "
								<< pythia.info.hiinfo->nCollSDP() << "   "
								<< pythia.info.hiinfo->nCollSDT() << "   "
								<< pythia.info.hiinfo->nCollDD() << "   "
								<< pythia.info.hiinfo->nCollCD() << "   "
								<< pythia.info.hiinfo->nCollEL() << "   "
								<< pythia.info.hiinfo->nAbsProj() << "   "
								<< pythia.info.hiinfo->nDiffProj() << "   "
								<< pythia.info.hiinfo->nElProj() << "   "
								<< pythia.info.hiinfo->nAbsTarg() << "   "
								<< pythia.info.hiinfo->nDiffTarg() << "   "
								<< pythia.info.hiinfo->nElTarg() << "   "*/
								//<< pythia.info.nMPI() << "   "
								//<< pythia.info.hiinfo->sigmaTotThisEvent() << "   "
								//<< pythia.info.hiinfo->sigmaNDThisEvent() << "   "
								//<< pythia.info.hiinfo->sigmaTot() << "   "
								//<< pythia.info.hiinfo->sigmaTotErr() << "   "
								//<< pythia.info.hiinfo->sigmaND() << "   "
								//<< pythia.info.hiinfo->sigmaNDErr();
								//<< pythia.info.hiinfo->sigTot() << "   "
								//<< pythia.info.hiinfo->sigEl() << "   "
								//<< pythia.info.hiinfo->sigCDE() << "   "
								//<< pythia.info.hiinfo->sigSDE() << "   "
								//<< pythia.info.hiinfo->sigSDEP() << "   "
								//<< pythia.info.hiinfo->sigSDET() << "   "
								//<< pythia.info.hiinfo->sigDDE() << "   "
								//<< pythia.info.hiinfo->sigND();

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

	delete paraRdr;

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
	if ( centrality_class_lower_limit < 1e-6
			and 100.0 - centrality_class_upper_limit < 1e-6 )
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
		results[0] = ( 100.0 - centrality_class_upper_limit < 1e-6 ) ?
						0 : 
						event_multiplicities[upper_index];	//smaller multiplicity limit first
		results[1] = ( centrality_class_lower_limit < 1e-6 ) ?
						1e+9 : 
						event_multiplicities[lower_index];	//larger multiplicity limit second
	}

	return ( results );
}




void print_particle_record(
		int iEvent,
		vector<Particle> & particles_to_output,
		vector<int> & particle_is_thermal_or_decay,
		ofstream & record_stream,
		double ev_b_phi,
		bool include_tau_and_eta )
{

	double cphi = cos(ev_b_phi), sphi = sin(ev_b_phi);

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

			double tau2 = p.vProd().pPos()*p.vProd().pNeg();
			double tau_loc = (tau2 < 0.0) ? -sqrt(-tau2) : sqrt(tau2);

			record_stream
				//<< iEvent << "   "				// deprecated
				//<< i << "   "						// deprecated
				<< p.id() << "   "					// Monte-Carlo ID (need not all be same particle species!)
				<< thermal_or_decay << "   "		// 0 for thermal particle; 1 for resonance decay particle 
				<< p.e() << "   "
				<< p.px() << "   "
				<< p.py() << "   "
				<< p.pz() << "   ";
			if ( include_tau_and_eta )
				record_stream
					<< p.p().mCalc() << "   " 
					<< p.p().rap() << "   ";
			record_stream
				<< p.tProd() << "   "
				<<  cphi*p.xProd() + sphi*p.yProd() << "   "
				<< -sphi*p.xProd() + cphi*p.yProd() << "   "
				<< p.zProd();
			if ( include_tau_and_eta )
				record_stream
					<< "   " << tau_loc
					<< "   " << p.vProd().rap();
			record_stream
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
				//<< "   " << p.vProd().mCalc()
				//<< "   " << p.vProd().rap()
				<< endl;
		}
	}

	return;
}

// End of file
