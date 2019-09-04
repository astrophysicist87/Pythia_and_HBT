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



void BalanceFunction::Compute_1p_spectra(int aRefMCID)
{
	Compute_dN_pTdpTdpphidpY(aRefMCID);

	Compute_dN_dpphidpY();

	Compute_N();

	return;
}



void BalanceFunction::Compute_dN_pTdpTdpphidpY(int aRefMCID)
{
	cout << "Starting here with allEvents.size() = " << allEvents.size() << " and particle_index = " << particle_index << endl;

	// Sum over all events
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];
		
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		{
			ParticleRecord p = event.particles[iParticle];

			// if this is neither the reference particle nor its anti-partner, continue
			if ( not p.aMCID == aRefMCID )
				continue;

			int ipT 	= floor( ( p.pT   - pT_min   ) / pT_bin_width   );
			int ipphi 	= floor( ( p.pphi - pphi_min ) / pphi_bin_width );
			int ipY 	= floor( ( p.pY   - pY_min   ) / pY_bin_width   );

			if ( 	   ipT   < 0 or ipT   >= n_pT_bins
					or ipphi < 0 or ipphi >= n_pphi_bins
					or ipY   < 0 or ipY   >= n_pY_bins
				) 	{
						//err << "Warning: discard particle with (pT, pphi, pY) = "
						//	<< p.pT << "   " << p.pphi << "   " << p.pY << endl;
						continue;
					}

			int refID = static_cast<int>( p.MCID < 0 );

			dN_pTdpTdpphidpY[refID][indexer(ipT, ipphi, ipY)] += 1.0;
			N[refID] += 1.0;
		}
	}

	// Normalize bin counts
	for (int ipT = 0; ipT < n_pT_bins; ++ipT)
	for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
	for (int ipY = 0; ipY < n_pY_bins; ++ipY)
	{
		double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);

		dN_pTdpTdpphidpY[refID][indexer(ipT, ipphi, ipY)]
				/= ( pT_bin_center
						* pT_bin_width * pphi_bin_width * pY_bin_width );
	}

	return;
}



void BalanceFunction::Compute_dN_dpphidpY()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	{
		int idx2D = 0;
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			for (int ipT = 0; ipT < n_pT_bins; ++ipT)
			{
				double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
				dN_dpphidpY[iRefID][idx2D] += pT_bin_center
										* dN_pTdpTdpphidpY[iRefID][indexer(ipT, ipphi, ipY)];
			}

			dN_dpphidpY[iRefID][idx2D++] *= pT_bin_width;
		}
	}

	return;
}




void BalanceFunction::Compute_N()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	{
		double N0 = 0.0;

		int idx3D = 0;
		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			double pT_bin_center = 0.5*(pT_pts[ipT]+pT_pts[ipT+1]);
			N0 += pT_bin_center * dN_pTdpTdpphidpY[iRefID][idx3D++]
					* pT_bin_width * pphi_bin_width * pY_bin_width;
		}

		out << "Total multiplicity N0 = " << N0 << endl;
		out << "Check: N[" << iRefID <<"] = " << N[iRefID] << endl;
	}

	return;
}



