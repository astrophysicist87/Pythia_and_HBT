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
#include "src/BalanceFunction.h"
#include "src/ParameterReader.h"
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "src/random_events.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Display intro
	cout << endl
			<< "              HBT event generator              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 10/2018" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo(2); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
   

	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	// Read-in particle and ensemble information
	vector<string> particle_info_filename, ensemble_info_filename;
	read_file_catalogue("./particle_catalogue.dat", particle_info_filename);
	paraRdr->readFromFile(particle_info_filename[0]);

	// Read-in command-line arguments
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();

	// Start timing
	Stopwatch sw;
	sw.Start();

	// Throws exception if NaNs are encountered
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	
	// Set-up output files
	string path = "./results/";	// make sure this directory exists
	string chosen_particle_name = "pi";
	ostringstream out_filename_stream, err_filename_stream;
	out_filename_stream << path << "BF_"
						<< chosen_particle_name << chosen_particle_name
						<< ".out";
	err_filename_stream << path << "BF_"
						<< chosen_particle_name << chosen_particle_name
						<< ".err";
	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());
	


	// Specify files containing all position-momentum information
	// from which to construct HBT correlation function
	vector<string> all_file_names;
	read_file_catalogue("./catalogue.dat", all_file_names);

	// Process multiplicity and ensemble information
	vector<string> ensemble_info;
	read_file_catalogue("./ensemble_catalogue.dat", ensemble_info);

	// allows to give files appropriate names
	string collision_system_info = ensemble_info[0];
	string target_name, projectile_name, beam_energy;
	int Nevents;
	double centrality_minimum, centrality_maximum;
	istringstream iss(collision_system_info);
	iss >> target_name
		>> projectile_name
		>> beam_energy
		>> centrality_minimum
		>> centrality_maximum
		>> Nevents;


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


	// Vector to hold all event information
	vector<EventRecord> allEvents;


	// Read in the first file
	int iFile = 0;
	cout << "Processing " << all_file_names[iFile] << "..." << endl;

	get_all_events(all_file_names[iFile], allEvents, paraRdr);


	// Create BalanceFunction object here
	BalanceFunction
		balance_function( 211, 211,
							paraRdr, allEvents,
							outmain, errmain );

	/*
	// Loop over the rest of the files
	for (iFile = 1; iFile < all_file_names.size(); ++iFile)
	{

		cout << "Processing " << all_file_names[iFile] << "..." << endl;

		// Read in the next file
		get_all_events(all_file_names[iFile], allEvents, paraRdr);


		// - for each file, update numerator and denominator
		//balance_function.Update_distributions( allEvents );

	}*/

	// Compute correlation function itself (after
	// all events have been read in)
	//balance_function.Compute_balance_function();
	balance_function.Compute_1p_spectra(0);

	balance_function.Compute_2p_spectra(0, 0);

	balance_function.Compute_rho1(0);

	balance_function.Compute_rho2(0, 0);

	balance_function.Check_normalizations(0, 0, 0);

	// Output results
	//balance_function.Output_balance_function( "./results/HBT_pipiCF.dat" );
	balance_function.Output_1p_spectra( 0, "./results/pi_1p_spectra.dat" );


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}

//End of file
