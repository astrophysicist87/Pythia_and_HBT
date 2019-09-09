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



void BalanceFunction::Compute_rho1()
{
	Compute_rho1_pT_pphi_pY();

	Compute_rho1_pphi_pY();

	return;
}



void BalanceFunction::Compute_rho1_pT_pphi_pY()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	{
		// Normalize
		int idx3D = 0;
		for (int ipT = 0; ipT < n_pT_bins; ++ipT)
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			rho1_pT_pphi_pY[iRefID][idx3D]
				= dN_pTdpTdpphidpY[iRefID][idx3D] / N[iRefID];
			idx3D++;
		}
	}

	return;
}



void BalanceFunction::Compute_rho1_pphi_pY()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	{
		// Normalize
		int idx2D = 0;
		for (int ipphi = 0; ipphi < n_pphi_bins; ++ipphi)
		for (int ipY = 0; ipY < n_pY_bins; ++ipY)
		{
			rho1_pphi_pY[iRefID][idx2D]
				= dN_dpphidpY[iRefID][idx2D] / N[iRefID];
			idx2D++;
		}
	}

	return;
}


