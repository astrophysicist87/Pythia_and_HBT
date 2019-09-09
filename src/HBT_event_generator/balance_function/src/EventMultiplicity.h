#ifndef EVENT_MULTIPLICITY_H
#define EVENT_MULTIPLICITY_H

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
	int total_multiplicity;		// all particles from event
	int particle_multiplicity;	// all particles used in HBT for this event
} EventMultiplicity;


bool compareByTotalMultiplicity(const EventMultiplicity &a, const EventMultiplicity &b);
bool compareByEventID(const EventMultiplicity &a, const EventMultiplicity &b);
void get_sorted_multiplicity(string filename, vector<EventMultiplicity> & ensemble);
void get_events_in_centrality_class( vector<string> & ensemble_info, vector<EventMultiplicity> & ensemble, bool verbose = false );

#endif
