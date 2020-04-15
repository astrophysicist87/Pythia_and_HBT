#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

typedef struct
{
	int eventID;
	int chosen_multiplicity;	// the multiplicity used for sorting
	vector<string> fields;		// event ID, total multiplicity, dNch/deta, etc...
} EventMultiplicity;


bool compareByMultiplicity(const EventMultiplicity &a, const EventMultiplicity &b)
{
	// a should come first if it has larger multiplicity
    return ( a.chosen_multiplicity > b.chosen_multiplicity );
}

bool compareByEventID(const EventMultiplicity &a, const EventMultiplicity &b)
{
	// a should come first if it has smaller eventID
    return ( a.eventID < b.eventID );
}

void read_in_multiplicities(string filename, int column_to_read, vector<EventMultiplicity> & ensemble)
{
	// load multiplicities from filename
	ifstream infile(filename.c_str());

	ensemble.clear();

	// read in data here
	string line;
	while ( getline(infile, line) )
	{
		EventMultiplicity thisEvent;
		istringstream iss( line );

		string field = "";
		while ( iss >> field ) thisEvent.fields.push_back( field );

		// if this line has three valid entries, store and move on
		thisEvent.eventID = stoi( thisEvent.fields[0] );
		thisEvent.chosen_multiplicity = stoi( thisEvent.fields[ column_to_read ] );
		ensemble.push_back( thisEvent );
	}

	infile.close();

	return;
}

void get_events_in_event_class(
			string filename,
			vector<EventMultiplicity> & ensemble,
			double minimum, double maximum,
			int column_to_use = 1,	// default column for total multiplicity
			string selection_mode = "centrality" )
{
	// Load ensemble information
	read_in_multiplicities(filename, column_to_use, ensemble);

	if ( selection_mode == "centrality" )
	{
		// sort by multiplicity
		sort( ensemble.begin(), ensemble.end(), compareByMultiplicity );
	
		// centrality selection
		double x1 = 0.01*minimum*static_cast<double>(ensemble.size());
		double x2 = 0.01*maximum*static_cast<double>(ensemble.size());
	
		int firstEvent = static_cast<int>(x1);
		int lastEvent = static_cast<int>(x2);
		//cout << ensemble.size() << "   " << x1 << "   " << x2 << "   " << firstEvent << "   " << lastEvent << endl;
	
		vector<EventMultiplicity>::const_iterator first = ensemble.begin() + firstEvent;
		vector<EventMultiplicity>::const_iterator last = ensemble.begin() + lastEvent;
		vector<EventMultiplicity> eventClass(first, last);
	
		ensemble = eventClass;
	}
	else if ( selection_mode == "multiplicity" )
	{
		vector<EventMultiplicity> eventClass;
		for ( const auto & event : ensemble )
			if ( event.chosen_multiplicity >= (int)minimum and event.chosen_multiplicity <= (int)maximum )
				eventClass.push_back( event );

		ensemble = eventClass;
	}
	else
	{
		cerr << "This selection mode not supported!  Crashing..." << endl;
		exit(8);
	}

	sort( ensemble.begin(), ensemble.end(), compareByEventID );

	return;
}

#endif
