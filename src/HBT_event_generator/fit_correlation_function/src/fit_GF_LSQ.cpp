#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "correlation_function.h"
#include "Stopwatch.h"
#include "Arsenal.h"

using namespace std;

void Correlation_function::Fit_correlation_function_GF_LSQ()
{

	out << "--> Getting HBT radii by Gaussian fit (least-squares) method" << endl;

	for (int iKT   = 0; iKT   < n_KT_bins;   ++iKT)
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
	for (int iKL   = 0; iKL   < n_KL_bins;   ++iKL)
	{
		out << "  --> Fitting in pair-momentum bin: "
			<< "KT=["   << KT_pts[iKT]     << ", " << KT_pts[iKT+1]     << "] GeV, "
			<< "Kphi=[" << Kphi_pts[iKphi] << ", " << Kphi_pts[iKphi+1] << "], "
			<< "KL=["   << KL_pts[iKL]     << ", " << KL_pts[iKL+1]     << "] GeV"
			<< endl;
		fit_correlationfunction_GF_lsq( iKT, iKphi, iKL );
	}

	out << "--> Finished getting HBT radii by Gaussian fit (least-squares) method" << endl;

	return;
}


int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		// digits in doubles

	int width = 15;		// setw width for output
	cout << scientific
		<< "iteration " << iteration << ": "
		<< "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 1)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 2)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 3)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 4)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 5)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 6)
		<< "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
		<< endl << endl;

	return 0;
}


int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	vector<double> q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	vector<double> q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	vector<double> q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	vector<double> y = ((struct Correlationfunction3D_data *) params_ptr)->y;
	vector<double> sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);
	double R2_ol = gsl_vector_get (xvec_ptr, 5);
	double R2_sl = gsl_vector_get (xvec_ptr, 6);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double Yi = 1.0 + lambda*exp( -q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o
                                      - 2.*q_o[i]*q_s[i]*R2_os - 2.*q_o[i]*q_l[i]*R2_ol - 2.*q_s[i]*q_l[i]*R2_sl);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
	}

	return GSL_SUCCESS;
}


int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	vector<double> q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	vector<double> q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	vector<double> q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	vector<double> sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);
	double R2_ol = gsl_vector_get (xvec_ptr, 5);
	double R2_sl = gsl_vector_get (xvec_ptr, 6);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double sig = sigma[i];

		//derivatives
		double common_elemt = exp( -q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o
                                   - 2.*q_o[i]*q_s[i]*R2_os - 2.*q_o[i]*q_l[i]*R2_ol - 2.*q_s[i]*q_l[i]*R2_sl);
      
		gsl_matrix_set (Jacobian_ptr, i, 0, common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q_o[i]*q_o[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 2, - lambda*q_s[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 3, - lambda*q_l[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 4, - 2.*lambda*q_o[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 5, - 2.*lambda*q_o[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 6, - 2.*lambda*q_s[i]*q_l[i]*common_elemt/sig);
	}

	return GSL_SUCCESS;
}


int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_correlfun3D_f_withlambda(xvec_ptr, params_ptr, f_ptr);
	Fittarget_correlfun3D_df_withlambda(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}


void Correlation_function::fit_correlationfunction_GF_lsq( int iKT, int iKphi, int iKL )
{
	const int VERBOSE = 10;
	const size_t n_para = 7;  // # of parameters

	// allocate space for a covariance matrix of size p by p
	gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

	// allocate and setup for generating gaussian distibuted random numbers
	gsl_rng_env_setup ();
	const gsl_rng_type *type = gsl_rng_default;
	gsl_rng *rng_ptr = gsl_rng_alloc (type);

	//set up test data
	struct Correlationfunction3D_data Correlfun3D_data;

	int n_usable_bins = 0;
	for (int i = 0; i < n_qo_bins; i++)
	for (int j = 0; j < n_qs_bins; j++)
	for (int k = 0; k < n_ql_bins; k++)
    {
		bool use_this_bin = true;
		if ( use_slices_only )
		{
			int i_in_center = int(i==(n_qo_bins-1)/2);
			int j_in_center = int(j==(n_qs_bins-1)/2);
			int k_in_center = int(k==(n_ql_bins-1)/2);
			use_this_bin = bool( i_in_center + j_in_center + k_in_center > 1 );
		}
		if ( not use_this_bin ) continue;

		int idx = indexer(iKT, iKphi, iKL, i, j, k);

        double q_out_local = 0.5*(qo_pts[i]+qo_pts[i+1]);
        double q_side_local = 0.5*(qs_pts[j]+qs_pts[j+1]);
        double q_long_local = 0.5*(ql_pts[k]+ql_pts[k+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if (correl_local < 1e-15) continue;

		if ( 	ignore_central_point
				and i==(n_qo_bins-1)/2
				and j==(n_qs_bins-1)/2
				and k==(n_ql_bins-1)/2)
			continue;

        Correlfun3D_data.q_o.push_back( q_out_local );
		Correlfun3D_data.q_s.push_back( q_side_local );
		Correlfun3D_data.q_l.push_back( q_long_local );
		Correlfun3D_data.y.push_back( correl_local );
		Correlfun3D_data.sigma.push_back( correl_err_local );

		// count this bin
		n_usable_bins++;
    }

	Correlfun3D_data.data_length = n_usable_bins;
	const size_t data_length     = Correlfun3D_data.data_length;

	// Make sure correlation function input is consistent (no reason it shouldn't be)
	if (    n_para                        >= data_length
         or Correlfun3D_data.q_o.size()   != data_length
         or Correlfun3D_data.q_s.size()   != data_length
         or Correlfun3D_data.q_l.size()   != data_length
         or Correlfun3D_data.y.size()     != data_length
         or Correlfun3D_data.sigma.size() != data_length )
	{
		cerr << "Error in fitting!  Exiting..." << endl;
		exit(8);
	}


	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f_withlambda;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df_withlambda;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf_withlambda;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &Correlfun3D_data;  // structure with the data and error bars

	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

	size_t iteration = 0;         // initialize iteration counter
	if (VERBOSE > 2) print_fit_state_3D_withlambda (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		if (VERBOSE > 2) cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		if (VERBOSE > 2) print_fit_state_3D_withlambda (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	if (VERBOSE > 2) cout << endl << "Covariance matrix: " << endl;
	if (VERBOSE > 2) gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	int width = 7;		// setw width for output

	const int iKT_iKphi_idx 			= indexerK(iKT, iKphi, iKL);
	lambda_Correl[iKT_iKphi_idx] 		= get_fit_results(0, solver_ptr);
	lambda_Correl_err[iKT_iKphi_idx] 	= c*get_fit_err(0, covariance_ptr);
	R2_out[iKT_iKphi_idx] 				= fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_side[iKT_iKphi_idx] 				= fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_long[iKT_iKphi_idx] 				= fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
	R2_outside[iKT_iKphi_idx] 			= get_fit_results(4, solver_ptr)*hbarC*hbarC;
	R2_outlong[iKT_iKphi_idx] 			= get_fit_results(5, solver_ptr)*hbarC*hbarC;
	R2_sidelong[iKT_iKphi_idx] 			= get_fit_results(6, solver_ptr)*hbarC*hbarC;
	R2_out_err[iKT_iKphi_idx] 			= c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_side_err[iKT_iKphi_idx] 			= c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_long_err[iKT_iKphi_idx] 			= c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[iKT_iKphi_idx] 		= c*get_fit_err(4, covariance_ptr)*hbarC*hbarC;
	R2_outlong_err[iKT_iKphi_idx] 		= c*get_fit_err(5, covariance_ptr)*hbarC*hbarC;
	R2_sidelong_err[iKT_iKphi_idx] 		= c*get_fit_err(6, covariance_ptr)*hbarC*hbarC;



	//clean up
	gsl_matrix_free (covariance_ptr);
	gsl_rng_free (rng_ptr);

	gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

	return;
}



//End of file
