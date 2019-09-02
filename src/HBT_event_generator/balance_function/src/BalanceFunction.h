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

		int n_pT_bins, n_pphi_bins, n_pY_bins;
		
		double pT_min, pT_max;

		vecotr<double> pT_pts, pphi_pts, pY_pts;

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
		/*inline int indexer(int iKT, int iKphi, int iKL, int ijKT, int ijKphi, int ijKL, int iqo, int iqs, int iql)
		{
			return (
					( ( ( ( ( ( ( iKT * n_Kphi_bins + iKphi )
										* n_KL_bins + iKL )
										* n_KT_pts_per_bin + ijKT )
										* n_Kphi_pts_per_bin + ijKphi )
										* n_KL_pts_per_bin + ijKL )
										* n_qo_bins + iqo )
										* n_qs_bins + iqs )
										* n_ql_bins + iql
					);
		}*/
		////////////////////

		////////////////////

		// Single-particle spectra
		void Compute_spectra();
		void Compute_dN_pTdpTdpphidpY();
		/*void Compute_dN_2pipTdpTdpY();
		void Compute_dN_pTdpTdpphi();
		void Compute_dN_dpphidpY();
		void Compute_dN_2pipTdpT();
		void Compute_dN_dpphi();
		void Compute_dN_2pidpY();*/
		void Compute_N();

		// For subsequent chunks of events
		//void Update_distributions( const vector<EventRecord> & allEvents_in );

		// Correlation function itself
		//void Compute_balance_function();

		// Input/output
		//void Output_balance_function( string filename );

};

#endif
