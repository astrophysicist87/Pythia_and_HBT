#ifndef READ_IN_DATA_H
#define READ_IN_DATA_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "../main.h"
#include "EventRecord.h"
#include "ParticleRecord.h"
#include "ParameterReader.h"
#include "ensemble.h"

using namespace std;

// function to read in catalogue of event files and return number of files to read
int read_file_catalogue(string catalogue_name, vector<string> & allLines)
{
	ifstream catalogue_in(catalogue_name.c_str());

	string line;
	while (getline(catalogue_in, line))
		allLines.push_back(line);

	return ( allLines.size() );
}

inline void complete_particle(ParticleRecord & p)
{
	double E = p.E, px = p.px, py = p.py, pz = p.pz;
	double t = p.t, x = p.x, y = p.y, z = p.z;

	p.pT 		= sqrt(px*px+py*py);
	p.pMag 		= sqrt(px*px+py*py+pz*pz);
	p.pphi 		= atan2(py, px);
	p.pY 		= 0.5*log(abs((E+pz)/(E-pz+1.e-100)));
	//p.pY = 0.0;
	p.ps_eta 	= 0.5*log(abs((p.pMag+pz)/(p.pMag-pz+1.e-100)));

	p.rT 		= sqrt(x*x+y*y);
	p.r 		= sqrt(x*x+y*y+z*z);
	p.phi 		= atan2(y, x);

	return;
}


void read_in_file_OSCAR(string filename, vector<EventRecord> & eventsInFile, ParameterReader * paraRdr)
{
	// open file
	ifstream infile(filename.c_str());

	// set the particle we want to do HBT on
	int chosen_MCID 		= paraRdr->getVal("chosen_MCID");
	int include_Bjorken 	= paraRdr->getVal("store_Bjorken_coordinates");
	bool include_thermal 	= (bool)paraRdr->getVal("include_thermal");
	bool include_decays 	= (bool)paraRdr->getVal("include_decays");

	// this vector contains events to include (for specific centrality class)
	int nextEventIndex = 0;
	int nextEventID = ensemble_multiplicites[nextEventIndex].eventID;
	int n_events_read_from_this_file = 0;

	do
	{
		int eventID = -1, nParticles = -1;
		int particleCount = 0;

		// read in first line: contains event index and number of particles stored from this event
		infile >> eventID;
		//cout << "Reading in eventID = " << eventID;
		if ( infile.eof() ) break;
		infile >> nParticles;
		//cout << " containing " << nParticles << " particles." << endl;

		// if this event not in chosen centrality class,
		// skip all the particles in it automatically
		if ( not ( eventID == nextEventID ) )
		{
			string line;
			for ( int iParticle = 0; iParticle < nParticles; ++iParticle )
				getline(infile, line);
		}
		else
		{
			//set next event to include
			// (negative means we're done reading in selected events)
			++nextEventIndex;
			nextEventID = ( nextEventIndex == ensemble_multiplicites.size() ) ?
							-1 : ensemble_multiplicites[nextEventIndex].eventID;
		}

		// to hold relevant particles for this event
		EventRecord event;

		int MCID, thermal_or_decay;
		double E, px, py, pz;
		double t, x, y, z;
		string dummy;

		// loop over all particles
		for ( int iParticle = 0; iParticle < nParticles; ++iParticle )
		{
			infile  //>> eventID		// these are now redundant
					//>> particleID		// and unnecessary
					>> MCID >> thermal_or_decay
					>> E >> px >> py >> pz;
			if ( include_Bjorken )
				infile >> dummy >> dummy;	//m and y
			infile  >> t >> x >> y >> z;
			if ( include_Bjorken )
				infile >> dummy >> dummy;	// tau and eta

			// skip the particles we don't care about
			if ( 	( not ( MCID == chosen_MCID ) )							// isn't the particle we're doing HBT on
					or ( not include_thermal and thermal_or_decay == 0 )	// or is thermal when we're excluding thermal
					or ( not include_decays  and thermal_or_decay == 1 )	// or is a decay when we're excluding decays
				)
			{
				//cout << "Skipping particle we don't care about: " << MCID << "   " << chosen_MCID << endl;
				continue;
			}

			// if we care, create a particle entry
			ParticleRecord particle;

			//cout << "\t - setting particle info..." << endl;
			double m = paraRdr->getVal("mass");
			particle.eventID 	= eventID;
			particle.particleID = iParticle;
			//particle.E 			= E;
			particle.E 			= sqrt( m*m + px*px + py*py + pz*pz );
			particle.px 		= px;
			particle.py 		= py;
			particle.pz 		= pz;
			particle.t 			= t / MmPerFm;
			particle.x 			= x / MmPerFm;
			particle.y 			= y / MmPerFm;
			particle.z 			= z / MmPerFm;

			//cout << "\t - completing particle information..." << endl;
			complete_particle(particle);

			// save particle
			event.particles.push_back(particle);

			particleCount++;
		}

		// this event completed;
		// store results in events vector
		if ( particleCount > 0 )
			eventsInFile.push_back(event);

		++n_events_read_from_this_file;

		// break if done with this centrality class
		if ( nextEventID < 0 )
		{
			//cout << "nextEventID < 0!  Exiting..." << endl;
			//cout << nextEventIndex << "   " << nextEventID << endl;
			break;
		}

	} while ( not infile.eof() );

	infile.close();

	// discard number of events read in from this file
	if ( n_events_read_from_this_file > 0 )
		ensemble_multiplicites.erase( ensemble_multiplicites.begin(),
									ensemble_multiplicites.begin()
									+ n_events_read_from_this_file );

	// debugging
	bool verbose = false;
	if (verbose)
	{
		cout << "=============================================================================================" << endl;
		for (int iEvent = 0; iEvent < (int)eventsInFile.size(); ++iEvent)
		{
			cout << "DEBUG: " << iEvent << endl;

			EventRecord event = eventsInFile[iEvent];
			for (int iParticle = 0; iParticle < (int)event.particles.size(); ++iParticle)
			{
				ParticleRecord particle = event.particles[iParticle];
				cout << "DEBUG: "
						<< particle.eventID << "   "
						<< particle.particleID << "   "
						<< particle.E << "   "
						<< particle.px << "   "
						<< particle.py << "   "
						<< particle.pz << "   "
						<< particle.t << "   "
						<< particle.x << "   "
						<< particle.y << "   "
						<< particle.z << endl;
			}
		}
		cout << "=============================================================================================" << endl;
	}

	return;
}


// function to read in a file containing some number of events
void read_in_file(string filename, vector<EventRecord> & eventsInFile, ParameterReader * paraRdr)
{
	bool OSCARformat = static_cast<bool>( paraRdr->getVal( "use_OSCAR_format" ) );
	if ( OSCARformat )
	{
		read_in_file_OSCAR(filename, eventsInFile, paraRdr);
		return;
	}

	bool check_particle_range = false;

	ifstream infile(filename.c_str());

	int count = 0;
	string line;
	int previous_eventID = -1, current_eventID = -1;
	
	EventRecord event;

	if ( ensemble_multiplicites.size() < 1 )
	{
		cout << "Finished reading in all necessary files!" << endl;
		return;
	}

	// filter out particles produced unreasonably far from collision point
	double max_accepted_tau = paraRdr->getVal("max_accepted_tau");

	// this vector contains events to include (for specific centrality class)
	int nextEventIndex = 0;
	int nextEventID = ensemble_multiplicites[nextEventIndex].eventID;
	cout << "nextEventID = " << nextEventID << endl;
	int countthis = 0;
	int n_events_read_from_this_file = 0;

	while (getline(infile, line))
	{
		istringstream iss(line);

		//cout << "Made it to line#" << count << endl;
		//cout << "Check event size: " << __LINE__ << "   " << event.particles.size() << endl;

		ParticleRecord particle;
		int eventID, particleID;
		double E, px, py, pz;
		double t, x, y, z;

		//cout << "\t - splitting up input line..." << endl;
		if ( !( iss >> eventID
					>> particleID
					>> E >> px >> py >> pz
					>> t >> x >> y >> z
			 ) ) { /*cout << "no success reading in!" << endl;*/ break; }

		//cout << "This particle: " << eventID << "   " << particleID << "   "
		//		<< E << "   " << px << "   " << py << "   " << pz << "   "
		//		<< t << "   " << x << "   " << y << "   " << z << endl;

		if ( check_particle_range )
		{
			const double tconv = t / MmPerFm;
			const double xconv = x / MmPerFm;
			const double yconv = y / MmPerFm;
			const double zconv = z / MmPerFm;

			const double current_tau2 = tconv*tconv - xconv*xconv - yconv*yconv - zconv*zconv;

			const double upper = 10.0;

			// if this particle is out of range, go to next particle
			if ( max_accepted_tau*max_accepted_tau < current_tau2
					or tconv*tconv - zconv*zconv > upper*upper
					or xconv*xconv + yconv*yconv > upper*upper
					or tconv*tconv > upper*upper
					or xconv*xconv > upper*upper
					or yconv*yconv > upper*upper
					or zconv*zconv > upper*upper
				)	//just some rough upper limits
				continue;
		}

		//cout << "eventID = " << eventID << endl;
		//cout << "nextEventID = " << nextEventID << endl;

		if ( eventID < nextEventID )
		{
			//if (countthis == 0) cout << eventID << " < " << nextEventID << endl;
			countthis++;
			continue;
		}
		countthis = 0;

		// apply momentum-space cuts, if any
		/*bool apply_momentum_space_cuts = false;
		if ( apply_momentum_space_cuts
				and ( abs(px) > max_pT
				or abs(py) > max_pT
				or abs(pz) > max_pz ) )
		{
			continue;
		}*/

		// apply position-space cuts, if any
		//if ()
		//{
		//	;
		//}

		//cout << "\t - setting particle info..." << endl;
		double m = paraRdr->getVal("mass");
		particle.eventID 	= eventID;
		particle.particleID = particleID;
		//particle.E 			= E;
		particle.E 			= sqrt( m*m + px*px + py*py + pz*pz );
		particle.px 		= px;
		particle.py 		= py;
		particle.pz 		= pz;
		particle.t 			= t / MmPerFm;
		particle.x 			= x / MmPerFm;
		particle.y 			= y / MmPerFm;
		particle.z 			= z / MmPerFm;

		//cout << "\t - completing particle information..." << endl;
		complete_particle(particle);

		//cout << "\t - storing particle information..." << endl;
		// Decide what to do with new particle
		// if on first iteration
		if (count == 0)
		{
			// initialize previous eventID
			// and current eventID
			previous_eventID = eventID;
			current_eventID = eventID;

			// push particle to event
			event.particles.push_back(particle);
		}
		// otherwise...
		else
		{
			current_eventID = eventID;

			// if newest particle does not
			// correspond to a new event
			if (current_eventID == previous_eventID)
			{
				event.particles.push_back(particle);
			}
			// if newest particle corresponds
			// to a new event
			else
			{
				// push event to eventsInFile
				event.eventID = previous_eventID;
				//cout << "current_eventID = " << current_eventID << endl;
				//cout << "Pushing previous_eventID = " << previous_eventID << " while reading in " << filename << endl;
				//cout << "event.particles.size() = " << event.particles.size() << endl;
				eventsInFile.push_back(event);
				++n_events_read_from_this_file;

				// reset event
				event = EventRecord();

				//set next event to include
				// (negative means we're done reading in selected events)
				++nextEventIndex;
				nextEventID = ( nextEventIndex == ensemble_multiplicites.size() ) ?
								-1 : ensemble_multiplicites[nextEventIndex].eventID;

				// break if done with this centrality class
				if ( nextEventID < 0 )
					goto finish;

				// skip if next event not included
				// in this centrality class
				if ( current_eventID < nextEventID )
				{
					count = 0;
					continue;
				}

				// otherwise, push new particle to new event
				event.particles.push_back(particle);

				//cout << "Check event size: " << __LINE__ << "   "
				//		<< event.particles.size() << endl;
			}
		}
		//cout << "\t - finished this loop!" << endl;

		previous_eventID = current_eventID;
		++count;
	}

	// push final event to eventsInFile
	if ( current_eventID > -1 and event.particles.size() > 0 )
	{
		event.eventID = current_eventID;
		//cout << "Pushing current_eventID = " << current_eventID << " after reading in " << filename << endl;
		//cout << "event.particles.size() = " << event.particles.size() << endl;
		eventsInFile.push_back(event);
		++n_events_read_from_this_file;
	}

	// reading in events terminates to here
	// if all events from specified centrality
	// class have been read in
	finish:

	infile.close();
	if ( n_events_read_from_this_file > 0 )
	{
		int n_current_events = ensemble_multiplicites.size();

		/*cout << "Deleting " << n_events_read_from_this_file
				<< " computed events from list of "
				<< n_current_events << " events..."
				<< endl;

		cout << "Before: "
				<< ensemble_multiplicites[0].eventID << " "
				<< ensemble_multiplicites[1].eventID << " "
				<< ensemble_multiplicites[2].eventID << "..."
				<< ensemble_multiplicites[n_current_events-3].eventID << " "
				<< ensemble_multiplicites[n_current_events-2].eventID << " "
				<< ensemble_multiplicites[n_current_events-1].eventID << endl;*/

		ensemble_multiplicites.erase( ensemble_multiplicites.begin(),
									ensemble_multiplicites.begin()
									+ n_events_read_from_this_file );

		/*n_current_events = ensemble_multiplicites.size();

		cout << "After: "
				<< ensemble_multiplicites[0].eventID << " "
				<< ensemble_multiplicites[1].eventID << " "
				<< ensemble_multiplicites[2].eventID << "..."
				<< ensemble_multiplicites[n_current_events-3].eventID << " "
				<< ensemble_multiplicites[n_current_events-2].eventID << " "
				<< ensemble_multiplicites[n_current_events-1].eventID << endl;

		cout << "n_current_events = " << n_current_events << endl;*/
	}

	//if (1) exit(8);

	return;
}




void get_all_events(vector<string> & all_file_names, vector<EventRecord> & allEvents, ParameterReader * paraRdr)
{
	// Read in the files
	vector<EventRecord> eventsInFile;

	allEvents.clear();
	for (int iFile = 0; iFile < all_file_names.size(); ++iFile)
	{
		// Reset
		eventsInFile.clear();

		// Read in this file
		read_in_file(all_file_names[iFile], eventsInFile, paraRdr);

		// Append these events to allEvents vector
		allEvents.insert( allEvents.end(),
							eventsInFile.begin(),
							eventsInFile.end() );
	}

	return;
}



void get_all_events(string file_name, vector<EventRecord> & allEvents, ParameterReader * paraRdr)
{
	// Read in the files
	vector<EventRecord> eventsInFile;

	// Reset
	allEvents.clear();
	eventsInFile.clear();

	// Read in this file
	read_in_file(file_name, eventsInFile, paraRdr);

	// Append these events to allEvents vector
	allEvents.insert( allEvents.end(),
						eventsInFile.begin(),
						eventsInFile.end() );

	return;
}


#endif
