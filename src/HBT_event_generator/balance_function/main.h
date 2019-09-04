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


void read_in_file(string filename, vector<EventRecord> & eventsInFile, ParameterReader * paraRdr)
{
	// open file
	ifstream infile(filename.c_str());

	// set the particle we want to do HBT on
	int chosen_MCID 		= paraRdr->getVal("chosen_MCID");
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
		if ( not eventID == nextEventID )
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

		// loop over all particles
		for ( int iParticle = 0; iParticle < nParticles; ++iParticle )
		{
			infile  //>> eventID		// these are now redundant
					//>> particleID		// and unnecessary
					>> MCID >> thermal_or_decay
					>> E >> px >> py >> pz
					>> t >> x >> y >> z;

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
