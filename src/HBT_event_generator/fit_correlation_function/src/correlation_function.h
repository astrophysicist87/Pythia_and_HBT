#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <omp.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_version.h>

#include "Arsenal.h"
#include "EventRecord.h"
#include "ParameterReader.h"
#include "ParticleRecord.h"

using namespace std;

const complex<double> i(0.0, 1.0);
const double hbarC = 0.197327053;	//GeV*fm

constexpr bool ignore_central_point = true;

const double fit_tolerance = 1e-6;
const int fit_max_iterations = 100;

struct correlationfunction_data
{
	size_t datalength;
	vector<double> qo, qs, ql;
	vector<double> A, B;	//numerator, denominator
};

struct Correlationfunction3D_data
{
	size_t data_length;
	vector<double> q_o;
	vector<double> q_s;
	vector<double> q_l;
	vector<double> y;
	vector<double> sigma;
};

#if GSL_MAJOR_VERSION < 2
	
	int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
	int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
	int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,	gsl_matrix *Jacobian_ptr);
	int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
	
	inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
	{
		return gsl_vector_get (solver_ptr->x, i);
	}
	
	inline double get_fit_err (int i, gsl_matrix * covariance_ptr)
	{
		return sqrt (gsl_matrix_get (covariance_ptr, i, i));
	}

#end if

class Correlation_function
{
	private:
		ParameterReader * paraRdr;

		//header info
		string particle_name;
		double particle_mass;

		int BE_mode, bin_mode, q_mode, fit_mode;
		bool include_cross_terms, use_slices_only, format_with_pairs;

		int n_Q_pts;
		int n_KT_pts, n_Kphi_pts, n_KL_pts;
		int n_qo_pts, n_qs_pts, n_ql_pts;
		int n_Kx_pts, n_Ky_pts, n_Kz_pts;
		int n_qx_pts, n_qy_pts, n_qz_pts;

		int n_Q_bins;
		int n_KT_bins, n_Kphi_bins, n_KL_bins;
		int n_qo_bins, n_qs_bins, n_ql_bins;
		int n_Kx_bins, n_Ky_bins, n_Kz_bins;
		int n_qx_bins, n_qy_bins, n_qz_bins;

		double KT_min, KT_max, Kphi_min, Kphi_max, KL_min, KL_max;
		double Kx_min, Kx_max, Ky_min, Ky_max, Kz_min, Kz_max;
		double qx_min, qx_max, qy_min, qy_max, qz_min, qz_max;

		double init_Q, delta_Q;
		double init_qo, init_qs, init_ql;
		double delta_qo, delta_qs, delta_ql;

		double KT_bin_width, Kphi_bin_width, KL_bin_width;
		double Kx_bin_width, Ky_bin_width, Kz_bin_width;
		double qx_bin_width, qy_bin_width, qz_bin_width;

		vector<string> all_file_names;
		vector<EventRecord> allEvents;

		vector<double> Q_pts;
		vector<double> KT_pts, Kphi_pts, KL_pts;
		vector<double> qo_pts, qs_pts, ql_pts;
		vector<double> Kx_pts, Ky_pts, Kz_pts;
		vector<double> qx_pts, qy_pts, qz_pts;
		
		//miscellaneous
		string path;
		ostream & out;
		ostream & err;

		vector<double> lambda_Correl, R2, R2_out, R2_side, R2_long,
						R2_outside, R2_outlong, R2_sidelong;
		vector<double> lambda_Correl_err, R2_err, R2_out_err, R2_side_err, R2_long_err,
						R2_outside_err, R2_outlong_err, R2_sidelong_err;
		vector<double> lambda_Correl_FRerr, R2_FRerr, R2_out_FRerr, R2_side_FRerr, R2_long_FRerr,
						R2_outside_FRerr, R2_outlong_FRerr, R2_sidelong_FRerr;

		vector<double> denominator, correlation_function, correlation_function_error;
		vector<bool> denominator_cell_was_filled;
		vector<complex<double> > numerator;
		vector<double> numCount, denCount;


	public:

		// Constructors, destructors, and initializers
		Correlation_function( ParameterReader * paraRdr_in,
								const string filename_in,
								ostream & out_stream = std::cout,
								ostream & err_stream = std::cerr )
								:
								out(out_stream),
								err(err_stream)
								{ initialize_all( paraRdr_in, filename_in ); };


		void initialize_all(ParameterReader * paraRdr_in,
								const string filename_in);

		~Correlation_function();

		////////////////////
		// Library functions
		inline int indexerK(int iKT, int iKphi, int iKL)
		{
			return ( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL );
		}
		////////////////////
		inline int indexer_Q_K(int iKT, int iKphi, int iKL, int iQ)
		{
			return (
					( ( iKT * n_Kphi_bins + iKphi )
							* n_KL_bins + iKL )
							* n_Q_bins + iQ
					);
		}
		////////////////////
		inline int indexer(int iKT, int iKphi, int iKL, int iqo, int iqs, int iql)
		{
			return (
					( ( ( ( iKT * n_Kphi_bins + iKphi )
								* n_KL_bins + iKL )
								* n_qo_bins + iqo )
								* n_qs_bins + iqs )
								* n_ql_bins + iql
					);
		}
		////////////////////


		// Input
		void Load_correlation_function( string filepath );

		// Output fit results
		void Output_HBTradii( string filepath );

		// Fit routines
		void Fit_correlation_function();
		void Fit_correlation_function_Q();
		void find_minimum_chisq_correlationfunction_full( int iKT, int iKphi, int iKL );
		void find_minimum_chisq_correlationfunction_full_FR( int iKT, int iKphi, int iKL, double qmin, double qmax, int nstep );
		void find_minimum_chisq_CFerr_full_FR( int iKT, int iKphi, int iKL, double qmax );
		void find_minimum_chisq_correlationfunction_Q( int iKT, int iKphi, int iKL );
		void find_minimum_chisq_correlationfunction_Q_FR( int iKT, int iKphi, int iKL, double Qmin, double Qmax, int nstep );
		void find_minimum_chisq_CFerr_Q_FR( int iKT, int iKphi, int iKL, double Qmax );

		// For the minimum log-likelihood approach
		void Fit_correlation_function_min_logL();
		//double LogL_PML_f(const gsl_vector *v, void *params);
		//void LogL_PML_df (const gsl_vector *v, void *params, gsl_vector *df);
		//void LogL_PML_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df);
		void set_CFdata(correlationfunction_data & CFdata, int iKT, int iKphi, int iKL);
		void fit_correlationfunction_minimum_log_likelihood(int iKT, int iKphi, int iKL);

		if ( GSL_MAJOR_VERSION < 2 )
		{
			void Fit_correlation_function_GF_LSQ();
			void fit_correlationfunction_GF_lsq( int iKT, int iKphi, int iKL );
		}
};

#endif
