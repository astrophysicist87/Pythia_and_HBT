#ifndef BALANCEFUNCTION_H
#define BALANCEFUNCTION_H

#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
//#include <unordered_map>

#include "ParameterReader.h"
#include "Arsenal.h"
#include "EventRecord.h"
#include "ParticleRecord.h"
#include "gauss_quadrature.h"
#include "EventMultiplicity.h"

using namespace std;

constexpr complex<double> i(0.0, 1.0);
constexpr double hbarC = 0.19733;								//GeV*fm

constexpr bool CONVERT_MM_TO_FM = true;							// needs to be true if running on Pythia output, false for Vishnu output
constexpr double MmPerFm = ( CONVERT_MM_TO_FM ) ? 1.e-12 : 1.0;	//mm-to-fm conversion


class BalanceFunction
{
	private:
		ParameterReader * paraRdr;

		int total_N_events, number_of_completed_events;

		int reference_MCID, associate_MCID;

		std::unordered_map<int, int> MCID_indices;
		vector<int> particle_MCIDs;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;
		vector<EventMultiplicity> ensemble;

		int n_pT_pts, n_pphi_pts, n_pY_pts;
		int n_pT_bins, n_pphi_bins, n_pY_bins;
		
		double pT_min, pT_max;
		double pphi_min, pphi_max;
		double pY_min, pY_max;

		double pT_bin_width, pphi_bin_width, pY_bin_width;

		vector<double> pT_pts, pphi_pts, pY_pts;

		vector<vector<double> >
			dN_pTdpTdpphidpY, dN_dpphidpY,
			rho1_pT_pphi_pY, rho1_pphi_pY;
		vector<double> N;

		vector<vector<vector<double> > >
			dN2_p1Tdp1Tdp1phidp1Y_p2Tdp2Tdp2phidp2Y,
			dN2_dp1phidp1Y_dp2phidp2Y,
			rho2_p1Tp1phip1Y_p2Tp2phip2Y,
			rho2_p1phip1Y_p2phip2Y;
		vector<vector<double> > N2;


		
		// miscellaneous
		string path;
		ostream & out;
		ostream & err;


	public:

		// Constructors, destructors, and initializers
		BalanceFunction( int reference_MCID_in, int associate_MCID_in,
								ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( reference_MCID_in, associate_MCID_in,
													paraRdr_in, allEvents_in ); };


		BalanceFunction( int reference_MCID_in, int associate_MCID_in,
								ParameterReader * paraRdr_in,
								const vector<string> & allEvents_filenames_in,
								const vector<EventMultiplicity> & ensemble_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( reference_MCID_in, associate_MCID_in,
													paraRdr_in, allEvents_filenames_in,
													ensemble_in ); };


		void initialize_all( int reference_MCID_in, int associate_MCID_in,
								ParameterReader * paraRdr_in,
								const vector<EventRecord> & allEvents_in );

		void initialize_all( int reference_MCID_in, int associate_MCID_in,
								ParameterReader * paraRdr_in,
								const vector<string> & allEvents_filenames_in,
								const vector<EventMultiplicity> & ensemble_in );

		~BalanceFunction();

		////////////////////
		// Library functions
		////////////////////
		inline int indexer(int ipphi, int ipY)
		{
			return (
					ipphi * n_pY_bins + ipY
					);
		}
		////////////////////
		inline int indexer(int ipT, int ipphi, int ipY)
		{
			return (
					( ipT * n_pphi_bins + ipphi )
										* n_pY_bins + ipY
					);
		}
		////////////////////
		inline int indexer(int ip1phi, int ip1Y, int ip2phi, int ip2Y)
		{
			return (
					( ( ip1phi * n_pY_bins + ip1Y )
								* n_pphi_bins + ip2phi )
								* n_pY_bins + ip2Y
					);
		}
		////////////////////
		inline int indexer(int ip1T, int ip1phi, int ip1Y, int ip2T, int ip2phi, int ip2Y)
		{
			return (
					( ( ( ( ip1T * n_pphi_bins + ip1phi )
								* n_pY_bins + ip1Y )
								* n_pT_bins + ip2T )
								* n_pphi_bins + ip2phi )
								* n_pY_bins + ip2Y
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
		
		void Check_normalizations(int particle_index, int ip1, int ip2);

		// For subsequent chunks of events
		void Get_spectra();

		// functions in FileReader.cpp
		void complete_particle(ParticleRecord & p);
		void read_in_file(string filename, vector<EventRecord> & eventsInFile, ParameterReader * paraRdr);
		void get_all_events(string file_name, vector<EventRecord> & allEvents, ParameterReader * paraRdr, bool verbose = false);

		// Correlation function itself
		void Get_balance_function();
		void Compute_event_averages(int particle_index, int ip1, int ip2);
		//void Compute_balance_function();

		// Input/output
		//void Output_balance_function( string filename );
		void Output_1p_spectra( int particle_index, string filename );
		void Output_2p_spectra( int ip1, int ip2, string filename );

};

#endif
