#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "EventMultiplicity.h"

using namespace std;

bool compareByTotalMultiplicity(const EventMultiplicity &a, const EventMultiplicity &b)
{
	// a should come first if it has larger multiplicity
    return ( a.total_multiplicity > b.total_multiplicity );
}

bool compareByEventID(const EventMultiplicity &a, const EventMultiplicity &b)
{
	// a should come first if it has smaller eventID
    return ( a.eventID < b.eventID );
}

void get_sorted_multiplicity(string filename, vector<EventMultiplicity> & ensemble)
{
	// load multiplicities from filename
	ifstream infile(filename.c_str());

	ensemble.clear();

	int eventID, total_multiplicity, particle_multiplicity;

	// read in data here
	string line;
	while (getline(infile, line))
	{
		istringstream iss(line);

		if ( !( iss >> eventID
					>> total_multiplicity
					>> particle_multiplicity
			 ) ) { break; }

		// if this line has three valid entries, store and move on
		EventMultiplicity thisEvent;
		thisEvent.eventID = eventID;
		thisEvent.total_multiplicity = total_multiplicity;
		thisEvent.particle_multiplicity = particle_multiplicity;

		ensemble.push_back( thisEvent );
	}

	// now sort the results
	sort( ensemble.begin(), ensemble.end(), compareByTotalMultiplicity );

	infile.close();

	return;
}

//void get_events_in_centrality_class(
//			string filename,
//			vector<EventMultiplicity> & ensemble,
//			double centrality_minimum, double centrality_maximum, bool verbose = false )
void get_events_in_centrality_class(
			vector<string> & ensemble_info,
			vector<EventMultiplicity> & ensemble,
			bool verbose )
{

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

	// contains EbE multiplicities
	string filename = ensemble_info[1];

	if (verbose)
		cout << "get_events_in_centrality_class(): "
				<< "Using centrality class: "
				<< centrality_minimum << "-"
				<< centrality_maximum << "%!" << endl
				<< "Reading in " << filename << endl;

	get_sorted_multiplicity(filename, ensemble);

	double x1 = 0.01*centrality_minimum*static_cast<double>(ensemble.size());
	double x2 = 0.01*centrality_maximum*static_cast<double>(ensemble.size());

	int firstEvent = static_cast<int>(x1);
	int lastEvent = static_cast<int>(x2);

	vector<EventMultiplicity>::const_iterator first = ensemble.begin() + firstEvent;
	vector<EventMultiplicity>::const_iterator last = ensemble.begin() + lastEvent;
	vector<EventMultiplicity> centralityClass(first, last);

	ensemble = centralityClass;
	sort( ensemble.begin(), ensemble.end(), compareByEventID );


	if (verbose)
	{
		cout << "get_events_in_centrality_class(): "
				<< "Using " << ensemble.size()
				<< " events in centrality class "
				<< centrality_minimum << "-"
				<< centrality_maximum << "%!" << endl;

		cout << "Check events in this centrality class: " << endl;
		for (int iEvent = 0; iEvent < ensemble.size(); ++iEvent)
			cout << ensemble[iEvent].eventID << "   "
					<< ensemble[iEvent].total_multiplicity << "   "
					<< ensemble[iEvent].particle_multiplicity << endl;
	}


	return;
}

