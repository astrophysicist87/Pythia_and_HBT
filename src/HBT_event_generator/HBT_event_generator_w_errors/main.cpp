#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <sys/time.h>
#include <fenv.h>

#include "src/Stopwatch.h"
#include "src/HBT_event_generator.h"
#include "src/ParameterReader.h"
//#include "/home/blixen/plumberg/src/ArsenalAndParameterReaderSource/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "src/random_events.h"
#include "main.h"

// these are included remotely
/*#include "FourVector.h"
#include "ParticleRecord.h"
#include "ParameterReader.h"
#include "shifter.h"*/

using namespace std;

// By default, events are unshifted or shifted externally
constexpr bool shift_events = false;


/*void convert_event_to_shifter_format(
		const EventRecord & event,
		vector<shift_lib::ParticleRecord> & event_to_shift,
		shift_lib::ParameterReader * converted_paraRdr );

void convert_shifter_format_to_event( 
		const vector<shift_lib::ParticleRecord> & event_to_shift,
		EventRecord & event );*/


int main(int argc, char *argv[])
{
	// Display intro
	cout << endl
			<< "              HBT event generator              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 04/2020" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo(2); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
   
	// Set parameter and catalogue filenames and results directory name from command-line
	string path                      = string(argv[1]);	// results directory
	string parametersFilename        = string(argv[2]);	
	string particleCatalogueFilename = string(argv[3]);	
	string catalogueFilename         = string(argv[4]);	
	string ensembleCatalogueFilename = string(argv[5]);	

	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	//paraRdr->readFromFile("./parameters.dat");
	paraRdr->readFromFile(parametersFilename);

	// Read-in particle and ensemble information
	vector<string> particle_info_filename, ensemble_info_filename;
	//read_file_catalogue("./particle_catalogue.dat", particle_info_filename);
	read_file_catalogue(particleCatalogueFilename, particle_info_filename);
	paraRdr->readFromFile(particle_info_filename[0]);

	// Read-in command-line arguments
	//paraRdr->readFromArguments(argc, argv);
	paraRdr->readFromArguments( argc, argv, (string)("#"), 6 );
	paraRdr->echo();

	// Check OpenMP
	cout << "Check OpenMP: omp_get_max_threads() = " << omp_get_max_threads() << endl;

	// Start timing
	Stopwatch sw;
	sw.Start();

	// Throws exception if NaNs are encountered
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	int file_mode = paraRdr->getVal("file_mode");


	// Set-up output files
	//string path = "./results/";	// make sure this directory exists
	string chosen_particle_name = "pi";
	ostringstream out_filename_stream, err_filename_stream;
	out_filename_stream << path << "/HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.out";
	err_filename_stream << path << "/HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.err";
	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());




	if ( file_mode == 1 )	// default==1: read-in particles from file
	{

		// Specify files containing all position-momentum information
		// from which to construct HBT correlation function
		vector<string> all_file_names;
		//read_file_catalogue("./catalogue.dat", all_file_names);
		read_file_catalogue(catalogueFilename, all_file_names);

		// Process multiplicity and ensemble information
		vector<string> ensemble_info;
		//read_file_catalogue("./ensemble_catalogue.dat", ensemble_info);
		read_file_catalogue(ensembleCatalogueFilename, ensemble_info);

		// allows to give files appropriate names
		string collision_system_info = ensemble_info[0];
		string target_name, projectile_name, beam_energy;
		int Nevents;
		istringstream iss(collision_system_info);
		iss >> target_name
			>> projectile_name
			>> beam_energy
			//>> centrality_minimum
			//>> centrality_maximum
			>> Nevents;

		double centrality_minimum = paraRdr->getVal("centrality_minimum");
		double centrality_maximum = paraRdr->getVal("centrality_maximum");

		cout << "run_HBT_event_generator(): "
				<< "Using centrality class: "
				<< centrality_minimum << "-"
				<< centrality_maximum << "%!" << endl;
	

		// select only those events falling into specificed centrality range
		string multiplicity_filename = ensemble_info[1];
		cout << "Reading in " << multiplicity_filename << endl;
		get_events_in_centrality_class(
					multiplicity_filename, ensemble_multiplicites,
					centrality_minimum, centrality_maximum );

		cout << "run_HBT_event_generator(): "
				<< "Using " << ensemble_multiplicites.size()
				<< " events in centrality class "
				<< centrality_minimum << "-"
				<< centrality_maximum << "%!" << endl;

		cout << "Check events in this centrality class: " << endl;
		for (int iEvent = 0; iEvent < ensemble_multiplicites.size(); ++iEvent)
			cout << ensemble_multiplicites[iEvent].eventID << "   "
					<< ensemble_multiplicites[iEvent].total_multiplicity << "   "
					<< ensemble_multiplicites[iEvent].particle_multiplicity << endl;


		// Proceed with HBT calculations
		string mode = "read stream";

		if ( mode == "default" or mode == "read stream" )
		{

			// Vector to hold all event information
			vector<EventRecord> allEvents;


			// Read in the first file
			int iFile = 0;
			cout << "Processing " << all_file_names[iFile] << "..." << endl;

			get_all_events(all_file_names[iFile], allEvents, paraRdr);


			// Create HBT_event_generator object here
			// note: numerator and denominator computed automatically
			HBT_event_generator
				HBT_event_ensemble( paraRdr, allEvents,
									outmain, errmain );


			// Loop over the rest of the files
			for (iFile = 1; iFile < all_file_names.size(); ++iFile)
			{

				cout << "Processing " << all_file_names[iFile] << "..." << endl;

				// Read in the next file
				get_all_events(all_file_names[iFile], allEvents, paraRdr);


				// - for each file, update numerator and denominator
				HBT_event_ensemble.Update_records( allEvents );

			}

			// Compute correlation function itself (after
			// all events have been read in)
			HBT_event_ensemble.Compute_correlation_function();


			// Output results
			HBT_event_ensemble.Output_correlation_function( path + "/HBT_pipiCF.dat" );

		}
		else if ( mode == "read all" )
		{

			// Vector to hold all event information
			vector<EventRecord> allEvents;


			// Read in the files
			get_all_events(all_file_names, allEvents, paraRdr);


			// Create HBT_event_generator object from allEvents
			HBT_event_generator
				HBT_event_ensemble( paraRdr, allEvents,
									outmain, errmain );


			// Compute correlation function itself
			HBT_event_ensemble.Compute_correlation_function();


			// Output correlation function
			HBT_event_ensemble.Output_correlation_function( path + "/HBT_pipiCF.dat" );

		}
	
	}
	else if ( file_mode == 0 )	//use random number generation
	{

		// Vector to hold all event information
		vector<EventRecord> allEvents;


		// Read in the files
		//generate_events(allEvents, paraRdr);
		generate_events_v2(allEvents, paraRdr);


		// Shift events here.
		/*if ( shift_events )
		{
			paraRdr->setVal("BE_mode", 1);

			shift_lib::ParameterReader * converted_paraRdr = new shift_lib::ParameterReader;
			converted_paraRdr->readFromFile("./parameters.dat");
			converted_paraRdr->readFromFile(particle_info_filename[0]);
			converted_paraRdr->readFromArguments(argc, argv);
			converted_paraRdr->setVal("BE_mode", 1);

			for ( auto & event: allEvents )
			{
				vector<shift_lib::ParticleRecord> event_to_shift;
				convert_event_to_shifter_format( event, event_to_shift, converted_paraRdr );

				shift_lib::shifter shifted_event( converted_paraRdr, event_to_shift, cout, cerr );

				convert_shifter_format_to_event( event_to_shift, event );
			}
		}*/


		// Create HBT_event_generator object from allEvents
		HBT_event_generator
			HBT_event_ensemble( paraRdr, allEvents,
								outmain, errmain );


		// Loop a few more times to build up statistics
		// N.B. - 2147483647 is max value for int
		//const int nLoops = 100;  //say
		const int nLoops = paraRdr->getVal("RNG_nLoops");
		for (int iLoop = 1; iLoop < nLoops; ++iLoop)
		{

			if ( iLoop >= nLoops )
				break;

			if ( iLoop % 1 == 0 )
				cout << "Starting iLoop = " << iLoop << endl;

			// Read in the next file
			//generate_events(allEvents, paraRdr);
			generate_events_v2(allEvents, paraRdr);


			// Shift events also here.
			/*if ( shift_events )
			{
				paraRdr->setVal("BE_mode", 1);

				shift_lib::ParameterReader * converted_paraRdr = new shift_lib::ParameterReader;
				converted_paraRdr->readFromFile("./parameters.dat");
				converted_paraRdr->readFromFile(particle_info_filename[0]);
				converted_paraRdr->readFromArguments(argc, argv);
				converted_paraRdr->setVal("BE_mode", 1);

				for ( auto & event: allEvents )
				{
					vector<shift_lib::ParticleRecord> event_to_shift;
					convert_event_to_shifter_format( event, event_to_shift, converted_paraRdr );
					shift_lib::shifter shifted_event( converted_paraRdr, event_to_shift, cout, cerr );
					convert_shifter_format_to_event( event_to_shift, event );
				}
			}*/


			// - for each file, update numerator and denominator
			HBT_event_ensemble.Update_records( allEvents );

		}


		// Compute correlation function itself
		HBT_event_ensemble.Compute_correlation_function();


		// Output correlation function
		HBT_event_ensemble.Output_correlation_function( path + "/HBT_pipiCF.dat" );

	}


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}


/*void convert_event_to_shifter_format(
		const EventRecord & event,
		vector<shift_lib::ParticleRecord> & event_to_shift,
		shift_lib::ParameterReader * converted_paraRdr )
{
	event_to_shift.clear();

	int particleIndex = 0;
	double pion_mass = converted_paraRdr->getVal("mass");
	for ( auto & particle : event.particles )
	{
		shift_lib::ParticleRecord particle_to_shift;

		particle_to_shift.particleID 	= particleIndex++;
		particle_to_shift.m 			= pion_mass;
		particle_to_shift.m2 			= pion_mass*pion_mass;
		particle_to_shift.x 			= shift_lib::Vec4( particle.x, particle.y, particle.z, particle.t );
		particle_to_shift.p 			= shift_lib::Vec4( particle.px, particle.py, particle.pz, particle.E );
		particle_to_shift.pShift 		= shift_lib::Vec4( 0.0, 0.0, 0.0, 0.0 );
		particle_to_shift.pComp 		= shift_lib::Vec4( 0.0, 0.0, 0.0, 0.0 );

		event_to_shift.push_back( particle_to_shift );
	}
	return;
}



void convert_shifter_format_to_event( 
		const vector<shift_lib::ParticleRecord> & event_to_shift,
		EventRecord & event )
{
	int particleIndex = 0;
	for ( auto & particle : event_to_shift )
	{
		event.particles.at(particleIndex).x = particle.x.px();
		event.particles.at(particleIndex).y = particle.x.py();
		event.particles.at(particleIndex).z = particle.x.pz();
		event.particles.at(particleIndex).t = particle.x.e();
		event.particles.at(particleIndex).px = particle.p.px();
		event.particles.at(particleIndex).py = particle.p.py();
		event.particles.at(particleIndex).pz = particle.p.pz();
		event.particles.at(particleIndex).E = particle.p.e();

		particleIndex++;
	}
	return;
}*/


//End of file
