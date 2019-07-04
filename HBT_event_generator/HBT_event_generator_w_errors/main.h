#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "src/EventRecord.h"
#include "src/ParticleRecord.h"
#include "src/ParameterReader.h"
#include "src/ensemble.h"

using namespace std;

//this is just to give this file a reason to exist for the moment...
const double plumbergtest = 0.;
const bool CONVERT_MM_TO_FM = true;	// needs to be true if running on Pythia output, false for Vishnu output
const double MmPerFm = ( CONVERT_MM_TO_FM ) ? 1.e-12 : 1.0;	//mm-to-fm conversion

vector<EventMultiplicity> ensemble_multiplicites;

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


// function to read in a file containing some number of events
void read_in_file(string filename, vector<EventRecord> & eventsInFile, ParameterReader * paraRdr)
{
	ifstream infile(filename.c_str());

	int count = 0;
	string line;
	int previous_eventID = -1, current_eventID = -1;
	
	EventRecord event;

	//=============================================
	// Set momentum ranges to try to reduce number
	// of unnecessary particles

	// pair momentum
	double KT_max 			= paraRdr->getVal("KTmax");
	double KL_max 			= paraRdr->getVal("KLmax");

	// relative momentum
	double n_qo_pts 		= paraRdr->getVal("n_qo_pts");
	double n_qs_pts 		= paraRdr->getVal("n_qs_pts");
	double n_ql_pts 		= paraRdr->getVal("n_ql_pts");
	double delta_qo 		= paraRdr->getVal("delta_qo");
	double delta_qs 		= paraRdr->getVal("delta_qs");
	double delta_ql 		= paraRdr->getVal("delta_ql");
	double max_qo 			= 0.5*double(n_qo_pts-1)*delta_qo;
	double max_qs 			= 0.5*double(n_qs_pts-1)*delta_qs;
	double max_ql 			= 0.5*double(n_ql_pts-1)*delta_ql;

	// now get limits
	double max_pT = 1.01*(KT_max + max(max_qo, max_qs));	//cut applies to both px and py
	double max_pz = 1.01*(KL_max + max_ql);
	//=============================================

	if ( ensemble_multiplicites.size() < 1 )
	{
		cout << "Finished reading in all necessary files!" << endl;
		return;
	}

	// this vector contains events to include (for specific centrality class)
	int nextEventIndex = 0;
	int nextEventID = ensemble_multiplicites[nextEventIndex].eventID;
	//cout << "nextEventID = " << nextEventID << endl;
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
		bool apply_momentum_space_cuts = false;
		if ( apply_momentum_space_cuts
				and ( abs(px) > max_pT
				or abs(py) > max_pT
				or abs(pz) > max_pz ) )
		{
			continue;
		}

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
