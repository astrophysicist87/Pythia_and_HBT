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
			<< "              Balance Function              " << endl
			<< endl
			<< "  Ver 1.0   ----- Christopher Plumberg, 9/2019" << endl;
	cout << endl << "**********************************************************" << endl;
	display_logo(2); // Hail to the king~
	cout << endl << "**********************************************************" << endl << endl;
   

	if (true)
	{
		cerr << "You still need to fix not eventID == ... error!!!" << endl;
		exit(8);
	}

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
	ostringstream out_filename_stream, err_filename_stream;
	out_filename_stream << path << "BF.out";
	err_filename_stream << path << "BF.err";
	ofstream outmain(out_filename_stream.str().c_str());
	ofstream errmain(err_filename_stream.str().c_str());
	


	// Specify files containing all position-momentum information
	// from which to construct HBT correlation function
	vector<string> all_file_names;
	read_file_catalogue("./catalogue.dat", all_file_names);

	// Process multiplicity and ensemble information
	vector<string> ensemble_info;
	read_file_catalogue("./ensemble_catalogue.dat", ensemble_info);

	// ensemble_info contains chosen centrality class information, etc.
	vector<EventMultiplicity> ensemble;
	get_events_in_centrality_class( ensemble_info, ensemble );

	const int pion = 211;
	const int kaon = 321;

	// Create BalanceFunction object here (all calculations done automatically)
	BalanceFunction
		balance_function( pion, pion,
				paraRdr, all_file_names,
				ensemble, outmain, errmain );

	// Output results
	//balance_function.Output_1p_spectra( 0, "./results/pi_1p_spectra.dat" );
	//balance_function.Output_2p_spectra( 0, 0, "./results/pi_2p_spectra.dat" );
	balance_function.Output_integrated_BF( "./results/integrated_BF_Dely_Delphi.dat" );
	balance_function.Output_integrated_BF_Dely( "./results/integrated_BF_Dely.dat" );
	balance_function.Output_integrated_BF_Delphi( "./results/integrated_BF_Delphi.dat" );
	//balance_function.Output_balance_function( "./results/BF_pipiCF.dat" );


	// Print out run-time
	sw.Stop();
	cout 	<< "Finished everything in "
			<< sw.printTime() << " seconds." << endl;


	// Wrap it up!
	return (0);
}

//End of file
