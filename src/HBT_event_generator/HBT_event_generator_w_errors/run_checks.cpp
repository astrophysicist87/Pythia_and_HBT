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
#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "src/random_events.h"
#include "run_checks.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Display intro
	cout << endl
			<< "              HBT event generator ( RUNNING CHECKS )             " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 2/2020" << endl;
	cout << endl << "**********************************************************" << endl;
	//display_logo(2); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
   

	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	// Read-in particle and ensemble information
	vector<string> particle_info_filename;
	read_file_catalogue("./particle_catalogue.dat", particle_info_filename);
	paraRdr->readFromFile(particle_info_filename[0]);

	// Read-in command-line arguments
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();

	// Check OpenMP
	cout << "Check OpenMP: omp_get_max_threads() = " << omp_get_max_threads() << endl;

	// Start timing
	Stopwatch sw;
	sw.Start();

	// Throws exception if NaNs are encountered
	feenableexcept(FE_INVALID | FE_OVERFLOW);

	// Set-up output files
	string path = "./results/";	// make sure this directory exists
	string chosen_particle_name = "pi";
	ostringstream out_filename_stream, err_filename_stream;
	out_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.out";
	err_filename_stream << path << "HBT_"
						<< chosen_particle_name << chosen_particle_name
						<< "CF.err";
	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());

	//-----------------------------
	//Use random number generation
	//-----------------------------

	// Vector to hold all event information
	vector<EventRecord> allEvents;


	// Read in the files
	//generate_events(allEvents, paraRdr);
	generate_events_v2(allEvents, paraRdr);


	// Shift events here.
	shiftEvents( allEvents );


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


		// Shift events here.
		shiftEvents();


		// - for each file, update numerator and denominator
		HBT_event_ensemble.Update_records( allEvents );

	}


	// Compute correlation function itself
	HBT_event_ensemble.Compute_correlation_function();


	// Output correlation function
	HBT_event_ensemble.Output_correlation_function( "./results/HBT_pipiCF.dat" );



	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}

//End of file
