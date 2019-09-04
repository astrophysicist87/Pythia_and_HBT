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
	particle_MCID1 = particle_MCID1_in;
	particle_MCID2 = particle_MCID2_in;
	particle_species_map
		= {
			{ particle_MCID1, 0 },
			{ particle_MCID2, 1 }
		  };

	// Load parameters
	paraRdr			= paraRdr_in;

	// Copy in records of all events
	allEvents		= allEvents_in;
	total_N_events	= allEvents.size();
	number_of_completed_events
					= 0;

	//Define various grid sizes
	// - pair momenta points at which to evaluate correlation function
	n_pT_pts 		= 11;
	pT_min 			= 0.0;
	pT_max 			= 1.0;
	n_pphi_pts 		= 11;
	pphi_min 		= -M_PI;
	pphi_max 		= M_PI;
	n_pY_pts 		= 11;
	pY_min 			= -5.0;
	pY_max 			= 5.0;

	n_pT_bins 		= n_pT_pts - 1;
	n_pphi_bins 	= n_pphi_pts - 1;
	n_pY_bins 		= n_pY_pts - 1;

	pT_pts 			= vector<double> (n_pT_pts);
	pphi_pts 		= vector<double> (n_pphi_pts);
	pY_pts 			= vector<double> (n_pY_pts);

	linspace(pT_pts, pT_min, pT_max);
	linspace(pphi_pts, pphi_min, pphi_max);
	linspace(pY_pts, pY_min, pY_max);

	pT_bin_width = pT_pts[1] - pT_pts[0];
	pphi_bin_width = pphi_pts[1] - pphi_pts[0];
	pY_bin_width = pY_pts[1] - pY_pts[0];

	cout << "CHECK: " << total_N_events << "   " << pT_bin_width << "   " << pphi_bin_width << "   " << pY_bin_width << endl;

	dN_pTdpTdpphidpY
		= vector<vector<double> >( 1,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );
	dN_dpphidpY
		= vector<vector<double> >( 1,
			vector<double>( n_pphi_bins*n_pY_bins, 0.0 ) );
	N = vector<double>(1, 0.0);



	rho1_pT_pphi_pY
		= vector<vector<double> >( 1,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );
	rho1_pphi_pY
		= vector<vector<double> >( 1,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );

	dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y
		= vector<vector<vector<double> > >( 1,
			vector<vector<double> >( 1,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins
							*n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	dN2_dp1phidp1Y_dp2phidp2Y
		= vector<vector<vector<double> > >( 1,
			vector<vector<double> >( 1,
			vector<double>( n_pphi_bins*n_pY_bins
							*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	rho2_p1Tp1phip1Y_p2Tp2phip2Y
		= vector<vector<vector<double> > >( 1,
			vector<vector<double> >( 1,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins
							*n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	rho2_p1phip1Y_p2phip2Y
		= vector<vector<vector<double> > >( 1,
			vector<vector<double> >( 1,
			vector<double>( n_pphi_bins*n_pY_bins
							*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	N2 = vector<vector<double> >( 1, vector<double>( 1, 0.0 ) );


	// Initialize all spectra here
	Get_spectra();


	return;
}


void BalanceFunction::initialize_all(
	int particle_MCID1_in,
	int particle_MCID2_in,
	ParameterReader * paraRdr_in,
	const vector<string> & allEvents_filenames_in,
	const vector<EventMultiplicity> & ensemble_in )
{
	int iFile = 0;
	all_file_names = allEvents_filenames_in;
	ensemble = ensemble_in;

	// use this to pass into initialization routine
	vector<EventRecord> allEvents_in;
	get_all_events(all_file_names[iFile], allEvents_in, paraRdr_in);

	// only call this once
	initialize_all( particle_MCID1_in, particle_MCID2_in,
					paraRdr_in, allEvents_in );

	// Loop over the rest of the files
	for (iFile = 1; iFile < all_file_names.size(); ++iFile)
	{

		// Read in the next file
		get_all_events(all_file_names[iFile], allEvents, paraRdr);	//N.B. - everything now initialized


		// - for each file, update spectra
		Get_spectra();

	}

}


BalanceFunction::~BalanceFunction()
{
	//clear everything

	return;
}



void BalanceFunction::Get_spectra()
{
	Compute_1p_spectra(0);

	Compute_2p_spectra(0, 0);

	return;
}



void BalanceFunction::Check_normalizations(int particle_index, int ip1, int ip2)
{
	//====================================
	// Check rho1 normalizations
	double normalization = 0.0;
	int idx3D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
		double pphi_bin_center = 0.5*(pphi_pts[ipphi]+pphi_pts[ipphi+1]);
		double pY_bin_center = 0.5*(pY_pts[ipY]+pY_pts[ipY+1]);
		normalization += pT_bin_center * rho1_pT_pphi_pY[particle_index][idx3D++]
							* pT_bin_width * pphi_bin_width * pY_bin_width;
	}
	out << "rho1_pT_pphi_pY normalization = "
		<< normalization << endl;

	//===================
	normalization = 0.0;
	int idx2D = 0;
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pphi_bin_center = 0.5*(pphi_pts[ipphi]+pphi_pts[ipphi+1]);
		double pY_bin_center = 0.5*(pY_pts[ipY]+pY_pts[ipY+1]);
		normalization += rho1_pphi_pY[particle_index][idx2D++]
							* pphi_bin_width * pY_bin_width;
	}
	out << "rho1_pphi_pY normalization = "
		<< normalization << endl;
	//====================================

	//====================================
	// Check rho2 normalizations
	normalization = 0.0;
	int idx6D = 0;
	for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
	for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
	for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
	{
		double p1T_bin_center = 0.5*(pT_pts[ip1T]+pT_pts[ip1T+1]);
		double p1phi_bin_center = 0.5*(pphi_pts[ip1phi]+pphi_pts[ip1phi+1]);
		double p1Y_bin_center = 0.5*(pY_pts[ip1Y]+pY_pts[ip1Y+1]);
		double p2T_bin_center = 0.5*(pT_pts[ip2T]+pT_pts[ip2T+1]);
		double p2phi_bin_center = 0.5*(pphi_pts[ip2phi]+pphi_pts[ip2phi+1]);
		double p2Y_bin_center = 0.5*(pY_pts[ip2Y]+pY_pts[ip2Y+1]);
		normalization += p1T_bin_center * p2T_bin_center
							* rho2_p1Tp1phip1Y_p2Tp2phip2Y[ip1][ip2][idx6D++]
							* pT_bin_width * pphi_bin_width * pY_bin_width
							* pT_bin_width * pphi_bin_width * pY_bin_width;
	}
	out << "rho2_p1Tp1phip1Y_p2Tp2phip2Y normalization = "
		<< normalization << endl;

	//===================
	normalization = 0.0;
	int idx4D = 0;
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
	for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
	{
		double p1phi_bin_center = 0.5*(pphi_pts[ip1phi]+pphi_pts[ip1phi+1]);
		double p1Y_bin_center = 0.5*(pY_pts[ip1Y]+pY_pts[ip1Y+1]);
		double p2phi_bin_center = 0.5*(pphi_pts[ip2phi]+pphi_pts[ip2phi+1]);
		double p2Y_bin_center = 0.5*(pY_pts[ip2Y]+pY_pts[ip2Y+1]);
		normalization += rho2_p1phip1Y_p2phip2Y[ip1][ip2][idx4D++]
							* pphi_bin_width * pY_bin_width
							* pphi_bin_width * pY_bin_width;
	}
	out << "rho2_p1phip1Y_p2phip2Y normalization = "
		<< normalization << endl;
	//====================================

	return;
}



void BalanceFunction::Output_1p_spectra( int particle_index, string filename )
{

	ofstream ofs( filename.c_str() );

	int idx3D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
		double pphi_bin_center = 0.5*(pphi_pts[ipphi]+pphi_pts[ipphi+1]);
		double pY_bin_center = 0.5*(pY_pts[ipY]+pY_pts[ipY+1]);
		ofs << setprecision(18) << setw(25)
			<< pT_bin_center   << "   "
			<< pphi_bin_center << "   "
			<< pY_bin_center   << "   "
			<< dN_pTdpTdpphidpY[particle_index][idx3D++]
			<< endl;
	}

	ofs.close();

	return;
}

void BalanceFunction::Output_2p_spectra( int ip1, int ip2, string filename )
{

	ofstream ofs( filename.c_str() );

	int idx6D = 0;
	for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
	for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
	for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
	{
		double p1T_bin_center = 0.5*(pT_pts[ip1T]+pT_pts[ip1T+1]);
		double p1phi_bin_center = 0.5*(pphi_pts[ip1phi]+pphi_pts[ip1phi+1]);
		double p1Y_bin_center = 0.5*(pY_pts[ip1Y]+pY_pts[ip1Y+1]);
		double p2T_bin_center = 0.5*(pT_pts[ip2T]+pT_pts[ip2T+1]);
		double p2phi_bin_center = 0.5*(pphi_pts[ip2phi]+pphi_pts[ip2phi+1]);
		double p2Y_bin_center = 0.5*(pY_pts[ip2Y]+pY_pts[ip2Y+1]);

		ofs << setprecision(18) << setw(25)
			<< p1T_bin_center   << "   "
			<< p1phi_bin_center << "   "
			<< p1Y_bin_center   << "   "
			<< p2T_bin_center   << "   "
			<< p2phi_bin_center << "   "
			<< p2Y_bin_center   << "   "
			<< dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[ip1][ip2][idx6D++]
			<< endl;
	}

	ofs.close();

	return;
}


//End of file
