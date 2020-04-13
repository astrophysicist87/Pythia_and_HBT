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



void BalanceFunction::Compute_rho2()
{
	Compute_rho2_p1Tp1phip1Y_p2Tp2phip2Y();

	Compute_rho2_p1phip1Y_p2phip2Y();

	return;
}



void BalanceFunction::Compute_rho2_p1Tp1phip1Y_p2Tp2phip2Y()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	for (int iAssocID = 0; iAssocID < 2; iAssocID++)
	{
		// Normalize
		int idx6D = 0;
		for (int ip1T   = 0; ip1T   < n_pT_bins;   ++ip1T)
		for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
		for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
		for (int ip2T   = 0; ip2T   < n_pT_bins;   ++ip2T)
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			rho2_p1Tp1phip1Y_p2Tp2phip2Y[iRefID][iAssocID][idx6D]
				= dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y[iRefID][iAssocID][idx6D] / N2[iRefID][iAssocID];
			idx6D++;
		}
	}

	return;
}



void BalanceFunction::Compute_rho2_p1phip1Y_p2phip2Y()
{
	for (int iRefID = 0; iRefID < 2; iRefID++)
	for (int iAssocID = 0; iAssocID < 2; iAssocID++)
	{
		// Normalize
		int idx4D = 0;
		for (int ip1phi = 0; ip1phi < n_pphi_bins; ++ip1phi)
		for (int ip1Y   = 0; ip1Y   < n_pY_bins;   ++ip1Y)
		for (int ip2phi = 0; ip2phi < n_pphi_bins; ++ip2phi)
		for (int ip2Y   = 0; ip2Y   < n_pY_bins;   ++ip2Y)
		{
			rho2_p1phip1Y_p2phip2Y[iRefID][iAssocID][idx4D]
				= dN2_dp1phidp1Y_dp2phidp2Y[iRefID][iAssocID][idx4D] / N2[iRefID][iAssocID];
			idx4D++;
		}
	}
	return;
}


