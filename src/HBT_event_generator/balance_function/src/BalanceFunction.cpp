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
	int reference_MCID_in,
	int associate_MCID_in,
	ParameterReader * paraRdr_in,
	const vector<EventRecord> & allEvents_in )
{
	// set particles to study
	reference_MCID 	= reference_MCID_in;
	associate_MCID 	= associate_MCID_in;
	MCID_indices
		= {
			{ reference_MCID, 0 },
			{ associate_MCID, 1 }
		  };
	cout << "Have " << reference_MCID << " and " << associate_MCID << endl;
	reference_ID 	= 0;
	associate_ID 	= 1;
	particle_MCIDs	= { reference_MCID, associate_MCID };

	// Load parameters
	paraRdr			= paraRdr_in;

	// Copy in records of all events
	allEvents		= allEvents_in;
	total_N_events	= allEvents.size();
	number_of_completed_events
					= 0;

	//Define various grid sizes
	// - pair momenta points at which to evaluate correlation function
	n_pT_pts 		= 2;
	pT_min 			= 0.2;
	pT_max 			= 2.0;
	n_pphi_pts 		= 21;
	pphi_min 		= -M_PI;
	pphi_max 		= M_PI;
	n_pY_pts 		= 21;
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

	// set Delta pY and Delta pphi grids
	//n_Delta_pphi_pts
	//				= 2*n_pphi_bins;
	n_Delta_pphi_pts
					= n_pphi_pts;
	//n_Delta_pY_pts	= 2*n_pY_bins;
	n_Delta_pY_pts	= 20;
	//n_Delta_pphi_bins
	//				= n_Delta_pphi_pts - 1;
	n_Delta_pphi_bins
					= n_pphi_bins;
	n_Delta_pY_bins	= n_Delta_pY_pts - 1;
	Delta_pphi_min	= pphi_min;
	//Delta_pY_min	= 2.0*pY_min + 0.5*pY_bin_width;
	Delta_pY_min	= -9.5;
	Delta_pphi_max	= pphi_max;
	//Delta_pY_max	= 2.0*pY_max - 0.5*pY_bin_width;
	Delta_pY_max	= 9.5;

	Delta_pphi_pts	= vector<double> (n_Delta_pphi_pts);
	Delta_pY_pts	= vector<double> (n_Delta_pY_pts);

	//linspace(Delta_pphi_pts, Delta_pphi_min, Delta_pphi_max);
	Delta_pphi_pts = pphi_pts;
	linspace(Delta_pY_pts, Delta_pY_min, Delta_pY_max);

	Delta_pphi_binwidth = Delta_pphi_pts[1] - Delta_pphi_pts[0];
	Delta_pY_binwidth = Delta_pY_pts[1] - Delta_pY_pts[0];

	cout << "CHECK: " << total_N_events << "   " << pT_bin_width << "   "
			<< pphi_bin_width << "   " << pY_bin_width << "   "
			<< Delta_pphi_binwidth << "   " << Delta_pY_binwidth << endl;

	dN_pTdpTdpphidpY
		= vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );
	dN_dpphidpY
		= vector<vector<double> >( 2,
			vector<double>( n_pphi_bins*n_pY_bins, 0.0 ) );
	N = vector<double>(1, 0.0);



	rho1_pT_pphi_pY
		= vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );
	rho1_pphi_pY
		= vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) );

	dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y
		= vector<vector<vector<double> > >( 2,
			vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins
							*n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	dN2_dp1phidp1Y_dp2phidp2Y
		= vector<vector<vector<double> > >( 2,
			vector<vector<double> >( 2,
			vector<double>( n_pphi_bins*n_pY_bins
							*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	rho2_p1Tp1phip1Y_p2Tp2phip2Y
		= vector<vector<vector<double> > >( 2,
			vector<vector<double> >( 2,
			vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins
							*n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	rho2_p1phip1Y_p2phip2Y
		= vector<vector<vector<double> > >( 2,
			vector<vector<double> >( 2,
			vector<double>( n_pphi_bins*n_pY_bins
							*n_pphi_bins*n_pY_bins, 0.0 ) ) );

	N2 = vector<vector<double> >( 2, vector<double>( 2, 0.0 ) );


	differential3D_bf
		= vector<double>( n_pT_bins*n_pphi_bins*n_pY_bins
							*n_pT_bins*n_pphi_bins*n_pY_bins, 0.0 );
	differential2D_bf
		= vector<double>( n_pphi_bins*n_pY_bins*n_pphi_bins*n_pY_bins, 0.0 );
	integrated_bf
		= vector<double>( n_Delta_pphi_bins*n_Delta_pY_bins, 0.0 );
	integrated_bf_Dely
		= vector<double>( n_Delta_pY_bins, 0.0 );
	integrated_bf_Delphi
		= vector<double>( n_Delta_pphi_bins, 0.0 );


	// Initialize all spectra here
	Get_spectra();


	return;
}


void BalanceFunction::initialize_all(
	int reference_MCID_in,
	int associate_MCID_in,
	ParameterReader * paraRdr_in,
	const vector<string> & allEvents_filenames_in,
	const vector<EventMultiplicity> & ensemble_in )
{
	int iFile = 0;
	all_file_names = allEvents_filenames_in;
	ensemble = ensemble_in;

	// use this to pass into initialization routine
	vector<EventRecord> allEvents_in;
	get_all_events(all_file_names[iFile], allEvents_in, paraRdr_in, true);

	// only call this once
	initialize_all( reference_MCID_in, associate_MCID_in,
					paraRdr_in, allEvents_in );

	// Loop over the rest of the files
	for (iFile = 1; iFile < all_file_names.size(); ++iFile)
	{

		// Read in the next file
		get_all_events(all_file_names[iFile], allEvents, paraRdr, true);	//N.B. - everything now initialized


		// keep running total
		total_N_events += allEvents.size();


		// - for each file, update spectra
		Get_spectra();

	}

	// everything else done automatically here
	Get_balance_functions();

	return;
}


BalanceFunction::~BalanceFunction()
{
	//clear everything

	return;
}



void BalanceFunction::Get_spectra()
{
	Compute_1p_spectra(reference_MCID);

	Compute_2p_spectra(reference_MCID, associate_MCID);

	return;
}




void BalanceFunction::Check_normalizations()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
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
			normalization += pT_bin_center * rho1_pT_pphi_pY[iRefID][idx3D++]
								* pT_bin_width * pphi_bin_width * pY_bin_width;
		}
		out << "rho1_pT_pphi_pY["
			<< iRefID << "] normalization = "
			<< setw(10) << normalization << endl;

		//===================
		normalization = 0.0;
		int idx2D = 0;
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			double pphi_bin_center = 0.5*(pphi_pts[ipphi]+pphi_pts[ipphi+1]);
			double pY_bin_center = 0.5*(pY_pts[ipY]+pY_pts[ipY+1]);
			normalization += rho1_pphi_pY[iRefID][idx2D++]
								* pphi_bin_width * pY_bin_width;
		}
		out << "rho1_pphi_pY["
			<< iRefID << "] normalization = "
			<< setw(10) << normalization << endl;
		//====================================

		for (int iAssocID = 0; iAssocID < 2; iAssocID++)
		{

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
									* rho2_p1Tp1phip1Y_p2Tp2phip2Y[iRefID][iAssocID][idx6D++]
									* pT_bin_width * pphi_bin_width * pY_bin_width
									* pT_bin_width * pphi_bin_width * pY_bin_width;
			}
			out << "rho2_p1Tp1phip1Y_p2Tp2phip2Y["
				<< iRefID << "][" << iAssocID << "] normalization = "
				<< setw(10) << normalization << endl;

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
				normalization += rho2_p1phip1Y_p2phip2Y[iRefID][iAssocID][idx4D++]
									* pphi_bin_width * pY_bin_width
									* pphi_bin_width * pY_bin_width;
			}
			out << "rho2_p1phip1Y_p2phip2Y["
				<< iRefID << "][" << iAssocID << "] normalization = "
				<< setw(10) << normalization << endl;
			//====================================
		}
	}


	//=========================================
	// Check balance function normalizations
	//=========================================
	out << "=== Checking balance function normalizations ===" << endl;
	{
		double normalization = 0.0, normalization_pT = 0.0;
		int idx6D = 0;
		for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
		for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
		for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
		for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			double p1T_bin_center = 0.5*(pT_pts[ip1T]+pT_pts[ip1T+1]);
			double p2T_bin_center = 0.5*(pT_pts[ip2T]+pT_pts[ip2T+1]);

			normalization
				+= pT_bin_width * pphi_bin_width * pY_bin_width
					* pT_bin_width * pphi_bin_width * pY_bin_width
					* differential3D_bf[idx6D];
			normalization_pT
				+= pT_bin_width * pphi_bin_width * pY_bin_width
					* pT_bin_width * pphi_bin_width * pY_bin_width
					* p1T_bin_center * p2T_bin_center * differential3D_bf[idx6D++];
		}

		out << "normalization = " << setw(10) << normalization << endl;
		out << "normalization_pT = " << setw(10) << normalization_pT << endl;

	}

	

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



void BalanceFunction::Output_integrated_BF( string filename )
{

	ofstream ofs( filename.c_str() );

	for (int iDelta_pphi = 0; iDelta_pphi < n_Delta_pphi_bins; iDelta_pphi++)
	for (int iDelta_pY = 0; iDelta_pY < n_Delta_pY_bins; iDelta_pY++)
	{
		const double Delta_pphi_loc = 0.5*(Delta_pphi_pts[iDelta_pphi]+Delta_pphi_pts[iDelta_pphi+1]);
		const double Delta_pY_loc = 0.5*(Delta_pY_pts[iDelta_pY]+Delta_pY_pts[iDelta_pY+1]);

		ofs << Delta_pphi_loc << "   "
			<< Delta_pY_loc << "   "
			<< integrated_bf[iDelta_pphi * n_Delta_pY_bins + iDelta_pY] << endl;
	}

	ofs.close();

	return;
}

void BalanceFunction::Output_integrated_BF_Dely( string filename )
{

	ofstream ofs( filename.c_str() );

	for (int iDelta_pY = 0; iDelta_pY < n_Delta_pY_bins; iDelta_pY++)
	{
		const double Delta_pY_loc = 0.5*(Delta_pY_pts[iDelta_pY]+Delta_pY_pts[iDelta_pY+1]);
		ofs << Delta_pY_loc << "   " << integrated_bf_Dely[iDelta_pY] << endl;
	}

	ofs.close();

	return;
}

void BalanceFunction::Output_integrated_BF_Delphi( string filename )
{

	ofstream ofs( filename.c_str() );

	for (int iDelta_pphi = 0; iDelta_pphi < n_Delta_pphi_bins; iDelta_pphi++)
	{
		const double Delta_pphi_loc = 0.5*(Delta_pphi_pts[iDelta_pphi]+Delta_pphi_pts[iDelta_pphi+1]);
		ofs << Delta_pphi_loc << "   " << integrated_bf_Delphi[iDelta_pphi] << endl;
	}

	ofs.close();

	return;
}



//End of file
