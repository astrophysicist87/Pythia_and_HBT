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



void BalanceFunction::Compute_spectra()
{
	Compute_dN_pTdpTdpphidpY();

	/*Compute_dN_2pipTdpTdpY();
	Compute_dN_pTdpTdpphi();
	Compute_dN_dpphidpY();

	Compute_dN_2pipTdpT();
	Compute_dN_dpphi();
	Compute_dN_2pidpY();*/

	Compute_dN_dpphidpY();

	Compute_N();

	return;
}



void BalanceFunction::Compute_dN_pTdpTdpphidpY()
{
	double Na = 0.0, Nb = 0.0;

	// Sum over all events
	int NEvents = allEvents.size();
	for (int iEvent = 0; iEvent < NEvents; ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];

			int ipT 	= floor( ( p.pT - pT_min ) / pT_bin_width );
			int ipphi 	= floor( ( p.pphi - pphi_min ) / pphi_bin_width );
			int ipY 	= floor( ( p.pY - pY_min ) / pY_bin_width );

			if ( ipT < 0 or ipT >= n_pT_bins
					or ipphi < 0 or ipphi >= n_pphi_bins
					or ipY < 0 or ipY >= n_pY_bins
				) 	{
						//err << "Warning: discard particle with (pT, pphi, pY) = " << p.pT << "   " << p.pphi << "   " << p.pY << endl;
						continue;
					}

			dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)] += 1.0;
			Na += 1.0;
		}
	}

	// Average over events
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);

		dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)]
				/= ( pT_bin_center * norm
						* pT_bin_width * pphi_bin_width * pY_bin_width );
	}

	return;
}



void BalanceFunction::Compute_dN_dpphidpY()
{
	//vector<double> dN_dpphidpY(n_pT_bins*n_pY_bins);

	int idx2D = 0;
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		dN_dpphidpY[idx2D] = 0.0;

		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		{
			double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
			dN_dpphidpY[idx2D] += pT_bin_center
									* dN_pTdpTdpphidpY[indexer(ipT, ipphi, ipY)];
		}

		dN_dpphidpY[idx2D++] *= pT_bin_width;
	}

	return;
}




void BalanceFunction::Compute_N()
{
	double N = 0.0;
	int NEvents = allEvents.size();

	int idx3D = 0;
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
		N += pT_bin_center * dN_pTdpTdpphidpY[idx3D++]
				* pT_bin_width * pphi_bin_width * pY_bin_width;
	}

	out << "Event-averaged multiplicity <N> = " << N << endl;
	out << "Total multiplicity N = " << NEvents*N << endl;

	return;
}



