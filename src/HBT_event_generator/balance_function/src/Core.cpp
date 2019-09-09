#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "BalanceFunction.h"
#include "Arsenal.h"
#include "Stopwatch.h"

void BalanceFunction::Get_balance_functions()
{
	Compute_event_averages();
	
	Compute_rho1();
	Compute_rho2();

	Check_normalizations();

	Compute_balance_functions();

}

void BalanceFunction::Compute_event_averages()
{
	double nev = static_cast<double>(total_N_events);
	cout << "CHECK: " << total_N_events << "   " << nev << endl;

	for (int iRefID = 0; iRefID < 2; iRefID++)
	{
		// average single-particle spectra
		int idx = 0;
		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
			dN_pTdpTdpphidpY[iRefID][idx++] /= nev;

		idx = 0;
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
			dN_dpphidpY[iRefID][idx++] /= nev;

		N[iRefID] /= nev;

		for (int iAssocID = 0; iAssocID < 2; iAssocID++)
		{
			// average two-particle spectra
			idx = 0;
			for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
			for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
			for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
			for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
			for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
			for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
				dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[iRefID][iAssocID][idx++] /= nev;

			idx = 0;
			for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
			for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
			for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
			for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
				dN2_dp1phidp1Y_dp2phidp2Y[iRefID][iAssocID][idx++] /= nev;

			N2[iRefID][iAssocID] /= nev;
		}
	}

	return;
}



void BalanceFunction::Compute_balance_functions()
{
	Compute_differential3D_balance_function();

	Compute_differential2D_balance_function();

	Compute_integrated_balance_function();

	return;
}


void BalanceFunction::Compute_differential3D_balance_function()
{
	int idx3D = 0;
	int idx6D = 0;
	for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	{
		for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			double term_pp  = rho2_p1Tp1phip1Y_p2Tp2phip2Y[0][0][idx6D]
								/ (rho1_pT_pphi_pY[0][idx3D] + 1.e-100);
			double term_pm  = rho2_p1Tp1phip1Y_p2Tp2phip2Y[0][1][idx6D]
								/ (rho1_pT_pphi_pY[0][idx3D] + 1.e-100);
			double term_mp  = rho2_p1Tp1phip1Y_p2Tp2phip2Y[1][0][idx6D]
								/ (rho1_pT_pphi_pY[1][idx3D] + 1.e-100);
			double term_mm  = rho2_p1Tp1phip1Y_p2Tp2phip2Y[1][1][idx6D]
								/ (rho1_pT_pphi_pY[1][idx3D] + 1.e-100);

			differential3D_bf[idx6D]
							= 0.5 * ( term_pm + term_mp - term_pp - term_mm );
			idx6D++;
		}
		idx3D++;
	}
	return;
}


void BalanceFunction::Compute_differential2D_balance_function()
{
	int idx2D = 0;
	int idx4D = 0;
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	{
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			double term_pp  = rho2_p1phip1Y_p2phip2Y[0][0][idx4D]
								/ (rho1_pphi_pY[0][idx2D] + 1.e-100);
			double term_pm  = rho2_p1phip1Y_p2phip2Y[0][1][idx4D]
								/ (rho1_pphi_pY[0][idx2D] + 1.e-100);
			double term_mp  = rho2_p1phip1Y_p2phip2Y[1][0][idx4D]
								/ (rho1_pphi_pY[1][idx2D] + 1.e-100);
			double term_mm  = rho2_p1phip1Y_p2phip2Y[1][1][idx4D]
								/ (rho1_pphi_pY[1][idx2D] + 1.e-100);

			differential2D_bf[idx4D]
							= 0.5 * ( term_pm + term_mp - term_pp - term_mm );
			idx4D++;
		}
		idx2D++;
	}
	return;
}


/*double interpolate_2D(const vector<double> & grid, const vector<double> & xpts, const vector<double> & ypts, double x0, double y0)
{
	
}*/



void BalanceFunction::Compute_integrated_balance_function()
{
	// just sum available cells?
	bool naive = true;
	if ( naive )
	{
		// Loop over reference and associate particle together
		for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
		for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			const double Delta_pphi = place_in_range(pphi_pts[ip1phi] - pphi_pts[ip2phi], -M_PI, M_PI);
			const double Delta_pY 	= pY_pts[ip1Y] - pY_pts[ip2Y];
			const int iDelta_pphi 	= floor( ( Delta_pphi - Delta_pphi_min ) / Delta_pphi_binwidth );
			const int iDelta_pY 	= floor( ( Delta_pY - Delta_pY_min ) / Delta_pY_binwidth );

			// skip if we're out of range
			if ( iDelta_pphi < 0 or iDelta_pphi >= n_Delta_pphi_bins )
				continue;
			if ( iDelta_pY   < 0 or iDelta_pY   >= n_Delta_pY_pts )
				continue;

			const double bin_volume = pphi_bin_width * pphi_bin_width * pY_bin_width * pY_bin_width;

			// otherwise, just lump into appropriate cell
			integrated_bf[ iDelta_pphi * n_Delta_pY_bins + iDelta_pY ]
				+= bin_volume * differential2D_bf[indexer(ip1phi, ip1Y, ip2phi, ip2Y)];
		}
	}
	/*else	// or try interpolating?
	{
		// Loop over relative separation
		for (int iDelta_phi = 0; iDelta_phi < 20; iDelta_phi++)
		for (int iDelta_Y = 0; iDelta_Y < 20; iDelta_Y++)
		{
			const double Delta_Y = Delta_Y_pts[iDelta_Y];
			const double Delta_phi = Delta_phi_pts[iDelta_phi];

			// Loop over reference particle
			for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
			for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
			{
				const double p2phi = pphi_pts[ip1phi] + Delta_phi;
				const double p2Y = pY_pts[ip1Y] + Delta_Y;

				// Interpolate differential2D_bf[indexer(ip1phi, ip1Y, :, :)] to (p2phi, p2Y)
				// ...
				// ...
				// ...
			}
		}
	}*/

	// now get BF's vs. Delta_pY and Delta_pphi
	for (int iDelta_pphi = 0; iDelta_pphi < n_Delta_pphi_bins; iDelta_pphi++)
	for (int iDelta_pY = 0; iDelta_pY < n_Delta_pY_bins; iDelta_pY++)
	{
		integrated_bf_Dely[iDelta_pY] += Delta_pphi_binwidth * integrated_bf[iDelta_pphi * n_Delta_pY_bins + iDelta_pY];
	}

	for (int iDelta_pphi = 0; iDelta_pphi < n_Delta_pphi_bins; iDelta_pphi++)
	for (int iDelta_pY = 0; iDelta_pY < n_Delta_pY_bins; iDelta_pY++)
	{
		const double int_BF_val = integrated_bf[iDelta_pphi * n_Delta_pY_bins + iDelta_pY];
		integrated_bf_Dely[iDelta_pY] 		+= Delta_pphi_binwidth * int_BF_val;
		integrated_bf_Delphi[iDelta_pphi] 	+= Delta_pY_binwidth * int_BF_val;
	}


	return;
}





// End of file
