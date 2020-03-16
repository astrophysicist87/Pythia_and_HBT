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

#include "/scratch/blixen/plumberg/test_shifter/include/FourVector.h"
#include "/scratch/blixen/plumberg/test_shifter/include/ParticleRecord.h"
#include "/scratch/blixen/plumberg/test_shifter/include/shifter.h"

using namespace std;

void convert_event_to_shifter_format(
		const EventRecord & event,
		vector<shift_lib::ParticleRecord> & event_to_shift );
void convert_shifter_format_to_event( 
		const vector<shift_lib::ParticleRecord> & event_to_shift,
		EventRecord & event );




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
	constexpr bool shift_events = true;
	if ( shift_events )
	{
		paraRdr->setVal("BE_mode", 1);
		for ( auto & event: allEvents )
		{
			vector<shift_lib::ParticleRecord> event_to_shift;
			convert_event_to_shifter_format( event, event_to_shift );
			shift_lib::shifter shifted_event( paraRdr, event_to_shift, cout, cerr );
			convert_shifter_format_to_event( event_to_shift, event );
		}
	}


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
		if ( shift_events )
			for ( auto & event: allEvents )
			{
				vector<shift_lib::ParticleRecord> event_to_shift;
				convert_event_to_shifter_format( event, event_to_shift );
				shift_lib::shifter shifted_event( paraRdr, event_to_shift, cout, cerr );
				convert_shifter_format_to_event( event_to_shift, event );
			}


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




void convert_event_to_shifter_format(
		const EventRecord & event,
		vector<shift_lib::ParticleRecord> & event_to_shift )
{
	event_to_shift.clear();
	event_to_shift.resize(event.particles.size());

	int particleIndex = 0;
	double pion_mass = 0.13957;
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
}


//End of file
