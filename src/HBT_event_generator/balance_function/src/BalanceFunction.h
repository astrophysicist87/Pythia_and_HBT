#ifndef HBTEG_H
#define HBTEG_H

#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <unordered_map>

#include "ParameterReader.h"
#include "Arsenal.h"
#include "EventRecord.h"
#include "ParticleRecord.h"
#include "gauss_quadrature.h"

using namespace std;

constexpr complex<double> i(0.0, 1.0);
constexpr double hbarC = 0.19733;	//GeV*fm

class BalanceFunction
{
	private:
		ParameterReader * paraRdr;

		int total_N_events, number_of_completed_events;

		std::unordered_map<int, int> particle_species_map;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		int n_pT_pts, n_pphi_pts, n_pY_pts;
		int n_pT_bins, n_pphi_bins, n_pY_bins;
		
		double pT_min, pT_max;
		double pphi_min, pphi_max;
		double pY_min, pY_max;

		double pT_bin_width, pphi_bin_width, pY_bin_width;

		vector<double> pT_pts, pphi_pts, pY_pts;

		vector<vector<double> > dN_pTdpTdpphidpY, dN_dpphidpY;
		vector<double> N;

		
		// miscellaneous
		string path;
		ostream & out;
		ostream & err;


	public:

		// Constructors, destructors, and initializers
		BalanceFunction( int particle_MCID1_in, int particle_MCID2_in,
								ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( particle_MCID1_in, particle_MCID2_in,
													paraRdr_in, allEvents_in ); };


		void initialize_all( int particle_MCID1_in, int particle_MCID2_in,
								ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in);

		~BalanceFunction();

		////////////////////
		// Library functions
		////////////////////
		inline int indexer(int ipT, int ipphi, int ipY)
		{
			return (
					( ipT * n_pphi_bins + ipphi )
										* n_pY_bins + ipY
					);
		}
		////////////////////

		////////////////////

		// Single-particle spectra
		void Compute_1p_spectra(int particle_index);
		void Compute_dN_pTdpTdpphidpY(int particle_index);
		void Compute_dN_dpphidpY(int particle_index);
		void Compute_N(int particle_index);

		
		// Single-particle distributions
		void Compute_rho1(int particle_index);
		void Compute_rho1_pT_pphi_pY(int particle_index);
		void Compute_rho1_pphi_pY(int particle_index);

		// Two-particle spectra
		void Compute_2p_spectra(int ip1, int ip2);
		void Compute_dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y(int ip1, int ip2);
		void Compute_dN2_dp1phidp1Y_dp2phidp2Y(int ip1, int ip2);
		void Compute_N2(int ip1, int ip2);

		// Two-particle distributions
		void Compute_rho2(int ip1, int ip2);
		void Compute_rho2_p1Tp1phip1Y_p2Tp2phip2Y(int ip1, int ip2);
		void Compute_rho2_p1phip1Y_p2phip2Y(int ip1, int ip2);
		






		// For subsequent chunks of events
		//void Update_distributions( const vector<EventRecord> & allEvents_in );

		// Correlation function itself
		//void Compute_balance_function();

		// Input/output
		//void Output_balance_function( string filename );
		void Output_1p_spectra( int particle_index, string filename );

};

#endif
