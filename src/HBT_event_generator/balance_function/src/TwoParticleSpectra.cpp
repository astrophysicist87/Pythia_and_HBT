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



void BalanceFunction::Compute_2p_spectra(int aRefMCID, int aAssocMCID)
{
	Compute_dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y(aRefMCID, aAssocMCID);

	Compute_dN2_dp1phidp1Y_dp2phidp2Y();

	Compute_N2();

	return;
}



void BalanceFunction::Compute_dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y(int aRefMCID, int aAssocMCID)
{
	// Sum over all events
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			if ( iParticle == jParticle )
				continue;

			ParticleRecord p1 = event.particles[iParticle];
			ParticleRecord p2 = event.particles[jParticle];

			// if this is neither the reference particle nor its anti-partner, continue
			if ( not p1.aMCID == aRefMCID )
				continue;

			// if this is neither the associate particle nor its anti-partner, continue
			if ( not p2.aMCID == aAssocMCID )
				continue;

			int ip1T 	= floor( ( p1.pT - pT_min ) / pT_bin_width );
			int ip1phi 	= floor( ( p1.pphi - pphi_min ) / pphi_bin_width );
			int ip1Y 	= floor( ( p1.pY - pY_min ) / pY_bin_width );
			int ip2T 	= floor( ( p2.pT - pT_min ) / pT_bin_width );
			int ip2phi 	= floor( ( p2.pphi - pphi_min ) / pphi_bin_width );
			int ip2Y 	= floor( ( p2.pY - pY_min ) / pY_bin_width );

			if ( 	   ip1T   < 0 or ip1T   >= n_pT_bins
					or ip1phi < 0 or ip1phi >= n_pphi_bins
					or ip1Y   < 0 or ip1Y   >= n_pY_bins
					or ip2T   < 0 or ip2T   >= n_pT_bins
					or ip2phi < 0 or ip2phi >= n_pphi_bins
					or ip2Y   < 0 or ip2Y   >= n_pY_bins
				) 	{
						// skip particle pairs out of range
						continue;
					}

			int refID 	= static_cast<int>( p1.MCID < 0 );
			int assocID = static_cast<int>( p2.MCID < 0 );

			dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[refID][assocID][indexer(ip1T, ip1phi, ip1Y, ip2T, ip2phi, ip2Y)] += 1.0;
			N2[refID][assocID] += 1.0;
		}
	}

	// Normalize bin counts
	for (int iRefID = 0; iRefID < 2; iRefID++)
	for (int iAssocID = 0; iAssocID < 2; iAssocID++)
	for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
	for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
	for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
	{
		double p1T_bin_center = 0.5*(pT_pts[ip1T]+pT_pts[ip1T+1]);
		double p2T_bin_center = 0.5*(pT_pts[ip2T]+pT_pts[ip2T+1]);

		dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[iRefID][iAssocID][indexer(ip1T, ip1phi, ip1Y, ip2T, ip2phi, ip2Y)]
				/= ( p1T_bin_center * p2T_bin_center
						* pT_bin_width * pphi_bin_width * pY_bin_width
						* pT_bin_width * pphi_bin_width * pY_bin_width );
	}

	return;
}



void BalanceFunction::Compute_dN2_dp1phidp1Y_dp2phidp2Y()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	for (int iAssocID = 0; iAssocID < 2; iAssocID++)
	for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
	for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
	for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
	for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
	{
		for (int ip1T = 0; ip1T < n_pT_bins; ++ip1T)
		for (int ip2T = 0; ip2T < n_pT_bins; ++ip2T)
		{
			double p1T_bin_center = 0.5*(pT_pts[ip1T]+pT_pts[ip1T+1]);
			double p2T_bin_center = 0.5*(pT_pts[ip2T]+pT_pts[ip2T+1]);
			dN2_dp1phidp1Y_dp2phidp2Y[iRefID][iAssocID][indexer(ip1phi, ip1Y, ip2phi, ip2Y)]
				+= p1T_bin_center * p2T_bin_center
					* dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[iRefID][iAssocID][indexer(ip1T, ip1phi, ip1Y, ip2T, ip2phi, ip2Y)];
		}

		dN2_dp1phidp1Y_dp2phidp2Y[iRefID][iAssocID][indexer(ip1phi, ip1Y, ip2phi, ip2Y)] *= pT_bin_width * pT_bin_width;
	}

	return;
}




void BalanceFunction::Compute_N2()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	for (int iAssocID = 0; iAssocID < 2; iAssocID++)
	{
		double N2loc = 0.0;

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
			N2loc += p1T_bin_center * p2T_bin_center
					* dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[iRefID][iAssocID][idx6D++]
					* pT_bin_width * pphi_bin_width * pY_bin_width
					* pT_bin_width * pphi_bin_width * pY_bin_width;
		}

		out << "Total pairs N2 = " << N2loc << endl;
		out << "Check: N2 = " << N2[iRefID][iAssocID] << endl;
	}

	return;
}



