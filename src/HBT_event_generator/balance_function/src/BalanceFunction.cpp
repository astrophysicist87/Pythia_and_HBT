#include <iostream>
#include <fstream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <random>
#include <algorithm>

#include "BalanceFunction.h"
#include "estimate_error.h"
#include "Arsenal.h"
#include "Stopwatch.h"

using namespace std;


void BalanceFunction::initialize_all(
	int particle_MCID1_in,
	int particle_MCID2_in,
	ParameterReader * paraRdr_in,
	const vector<EventRecord> & allEvents_in )
{
	// set particles to study
	particle_species_map
		= {
			{ particle_MCID1_in, 0 },
			{ particle_MCID2_in, 1 }
		  };

	// Load parameters
	paraRdr			= paraRdr_in;

	// Copy in records of all events
	allEvents		= allEvents_in;
	total_N_events	= allEvents.size();
	number_of_completed_events
					= 0;

	dN_pTdpTdpphidpY
		= vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) )
	N = vector(2, 0.0);

	return;
}

BalanceFunction::~BalanceFunction()
{
	//clear everything

	return;
}



//End of file
