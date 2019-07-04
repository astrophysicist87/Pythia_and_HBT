#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#include "correlation_function.h"
#include "Stopwatch.h"
#include "Arsenal.h"

using namespace std;


void Correlation_function::Fit_correlation_function_Q()
{
	if ( fit_mode == 0 )
	{

		out << "--> Getting HBT radii by Gaussian fit method" << endl;

		/*if ( use_fit_range_errors )
		{

			const double Qmax = -init_Q;	// choose maximum Q value to define maximum fit-range radius
			const double Qmin = 0.5 * Qmax;	// say
			const int nstep = 6;			// steps of 10% of Qmax

			for (int iKT = 0; iKT < n_KT_bins; ++iKT)
			for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
			for (int iKL = 0; iKL < n_KL_bins; ++iKL)
			{
				out << "  --> Doing fit-range analysis: "
					<< "KT=[" << KT_pts[iKT] << ", " << KT_pts[iKT+1] << "] GeV, "
					<< "Kphi=[" << Kphi_pts[iKphi] << ", " << Kphi_pts[iKphi+1] << "], "
					<< "KL=[" << KL_pts[iKL] << ", " << KL_pts[iKL+1] << "] GeV"
					<< endl;
				find_minimum_chisq_correlationfunction_Q_FR( iKT, iKphi, iKL, Qmin, Qmax, nstep );
			}
		}
		else*/
		{
			for (int iKT = 0; iKT < n_KT_bins; ++iKT)
			for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
			for (int iKL = 0; iKL < n_KL_bins; ++iKL)
			{
				out << "  --> Fitting in pair-momentum bin: "
					<< "KT=[" << KT_pts[iKT] << ", " << KT_pts[iKT+1] << "] GeV, "
					<< "Kphi=[" << Kphi_pts[iKphi] << ", " << Kphi_pts[iKphi+1] << "], "
					<< "KL=[" << KL_pts[iKL] << ", " << KL_pts[iKL+1] << "] GeV"
					<< endl;
				find_minimum_chisq_correlationfunction_Q( iKT, iKphi, iKL );
			}
		}

		out << "--> Finished getting HBT radii by Gaussian fit method" << endl;

	}
	else
	{
		err << "Fit_correlation_function_Q(): fit_mode = "
			<< fit_mode << " not currently supported!" << endl;
		exit(8);
	}

	return;
}


void Correlation_function::find_minimum_chisq_correlationfunction_Q( int iKT, int iKphi, int iKL )
{

	const size_t data_length = n_Q_bins;  // # of points

    double lambda, R;
    int dim = 2;	//# of fit parameters defined above
    int s_gsl;

    double * V = new double [dim];
    double * qweight = new double [dim];
    double ** T = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        V[i] = 0.0;
        T[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T[i][j] = 0.0;
    }

	gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
	gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
	gsl_permutation * perm = gsl_permutation_alloc (dim);

	int n_usable_bins = 0;
	for (int i = 0; i < n_Q_bins; i++)
    {
//cout << "Made it to i = " << i << endl;
		int idx = indexer_Q_K(iKT, iKphi, iKL, i);

        double Q_local = 0.5*(Q_pts[i]+Q_pts[i+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

//cout << "Check: " << correl_local << endl;

		if (correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_Q_bins-1)/2)
			continue;
//			correl_err_local = 1.0e10;	//ignore central point
        double sigma_k_prime = correl_err_local/correl_local;

        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
        double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;
//cout << "Check: " << Q_local << "   " << correl_local << endl;
//cout << "Check: " << inv_sigma_k_prime_sq << endl;
//cout << "Check: " << log_correl_over_sigma_sq << endl;

		qweight[0] = - 1.0;
		qweight[1] = Q_local*Q_local;

        for(int ij = 0; ij < dim; ij++)
        {
            V[ij] += qweight[ij]*log_correl_over_sigma_sq;
            T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
        }

        for(int ij = 1; ij < dim; ij++)
            T[ij][0] = T[0][ij];
            
        for(int ij = 1; ij < dim; ij++)
        {
            for(int lm = 1; lm < dim; lm++)
                T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
        }

		// if we made it here, go ahead and fit
		++n_usable_bins;
    }

	// if we don't have enough usable bins, don't fit!
	if ( n_usable_bins <= dim )
		return;

    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            gsl_matrix_set(T_gsl, i, j, T[i][j]);

    // Make LU decomposition of matrix T_gsl
    gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
    // Invert the matrix m
    gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

    double ** T_inverse = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        T_inverse[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
    }
    double * results = new double [dim];
    for(int i = 0; i < dim; i++)
    {
        results[i] = 0.0;
        for(int j = 0; j < dim; j++)
            results[i] += T_inverse[i][j]*V[j];
    }

	// compute results
	lambda_Correl[indexerK(iKT, iKphi, iKL)] = exp(results[0]);
	R2[indexerK(iKT, iKphi, iKL)] = results[1]*hbarC*hbarC;
out << "results: " << exp(results[0]) << "   " << results[1]*hbarC*hbarC << endl;

	// compute chi^2
    double chi_sq = 0.0;
	for (int i = 0; i < n_Q_bins; i++)
    {
		int idx = indexer_Q_K(iKT, iKphi, iKL, i);

        double Q_local = 0.5*(Q_pts[i]+Q_pts[i+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if(correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_Q_bins-1)/2)
			continue;
//			correl_err_local = 1.0e10;	//ignore central point
        double sigma_k_prime = correl_err_local/correl_local;

        chi_sq += pow( ( log(correl_local) - results[0] 
						+ results[1]*Q_local*Q_local )
						, 2 )
                  /sigma_k_prime/sigma_k_prime;
    }

    double chi_sq_per_dof = chi_sq/(data_length - dim);
    out << "chi_sq = " << chi_sq << endl;
    out << "Number d.o.f. = " << data_length - dim << endl;
    out << "chi_sq/d.o.f. = " << chi_sq_per_dof << endl;
	out << "Goodness-of-fit parameter Q = "
		<< gsl_sf_gamma_inc_Q (0.5*(data_length - dim), 0.5*chi_sq) << endl;

	//============================================
	// compute curvature and covariance matrices
    double ** curvature_mat = new double * [dim];
    double ** covariance_mat = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        curvature_mat[i] = new double [dim];
        covariance_mat[i] = new double [dim];
        for(int j = 0; j < dim; j++)
		{
            curvature_mat[i][j] = 0.0;
            covariance_mat[i][j] = 0.0;
		}
    }

	//=============================
	for (int i = 0; i < n_Q_bins; i++)
    {
		int idx = indexer_Q_K(iKT, iKphi, iKL, i);

        double Q_local = 0.5*(Q_pts[i]+Q_pts[i+1]);

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

		if (correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_Q_bins-1)/2)
			continue;
		//	correl_err_local = 1.0e10;	//ignore central point

        double sigma_k_prime = correl_err_local/correl_local;
        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);

		qweight[0] = - 1.0;
		qweight[1] = Q_local*Q_local;

		for(int ij = 0; ij < dim; ij++)
		for(int lm = 0; lm < dim; lm++)
            curvature_mat[ij][lm] += qweight[ij]*qweight[lm]
										* inv_sigma_k_prime_sq;
    }

	//=============================
	// covariance matrix is inverse of curvaure matrix
	gsl_matrix_invert( curvature_mat, covariance_mat, dim );

	//=============================
	// determine errors on fit parameters
	// (i.e., diagonal elements of covariance matrix)
	lambda_Correl_err[indexerK(iKT, iKphi, iKL)]
		= lambda_Correl[indexerK(iKT, iKphi, iKL)]*(exp(sqrt(covariance_mat[0][0]))-1);
	R2_err[indexerK(iKT, iKphi, iKL)] = sqrt(covariance_mat[1][1])*hbarC*hbarC;

out << "results errors: " << lambda_Correl_err[indexerK(iKT, iKphi, iKL)] << "   " << R2_err[indexerK(iKT, iKphi, iKL)] << endl;

	//=============================
    // clean up
    gsl_matrix_free (T_gsl);
    gsl_matrix_free (T_inverse_gsl);
    gsl_permutation_free (perm);

    delete [] qweight;
    delete [] V;
    for(int i = 0; i < dim; i++)
    {
        delete [] T[i];
        delete [] T_inverse[i];
		delete [] curvature_mat[i];
		delete [] covariance_mat[i];
    }
    delete [] T;
    delete [] T_inverse;
	delete [] curvature_mat;
	delete [] covariance_mat;
    delete [] results;

	return;
}


/*
void Correlation_function::find_minimum_chisq_correlationfunction_Q_FR( int iKT, int iKphi, int iKL, double Qmin, double Qmax, int nstep )
{
	// first set the radii themselves
	// this function also sets the ordinary fit errors (perhaps incorrectly)
	find_minimum_chisq_correlationfunction_full( iKT, iKphi, iKL );

	// check how much the fit radii change for different fir ranges and keep maximum deviation from radii from full grid
	for (int step = 0; step < nstep; step++)
		find_minimum_chisq_CFerr_Q_FR( iKT, iKphi, iKL, Qmin + step*(Qmax-Qmin)/static_case<double>(nstep-1.0) );

	return;
}


void Correlation_function::find_minimum_chisq_CFerr_Q_FR( int iKT, int iKphi, int iKL, double Qmax )
{

	const size_t data_length = n_Q_bins;  // # of points

    int dim = 2;	//# of fit parameters defined above
    int s_gsl;

    double * V = new double [dim];
    double * qweight = new double [dim];
    double ** T = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        V[i] = 0.0;
        T[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T[i][j] = 0.0;
    }

	gsl_matrix * T_gsl = gsl_matrix_alloc (dim, dim);
	gsl_matrix * T_inverse_gsl = gsl_matrix_alloc (dim, dim);
	gsl_permutation * perm = gsl_permutation_alloc (dim);

	int n_usable_bins = 0;
	for (int i = 0; i < n_Q_bins; i++)
    {
//cout << "Made it to i = " << i << endl;
		int idx = indexer_Q_K(iKT, iKphi, iKL, i);

        double Q_local = 0.5*(Q_pts[i]+Q_pts[i+1]);

		if ( Q_local*Q_local > Qmax*Qmax ) continue;

		double correl_local = correlation_function[idx]-1.0;
		double correl_err_local = correlation_function_error[idx];

//cout << "Check: " << correl_local << endl;

		if (correl_local < 1.0e-15) continue;

		bool ignore_central_point = true;
		if ( 	ignore_central_point
				and i==(n_Q_bins-1)/2)
			continue;
//			correl_err_local = 1.0e10;	//ignore central point
        double sigma_k_prime = correl_err_local/correl_local;

        double inv_sigma_k_prime_sq = 1./(sigma_k_prime*sigma_k_prime);
        double log_correl_over_sigma_sq = log(correl_local)*inv_sigma_k_prime_sq;
//cout << "Check: " << Q_local << "   " << correl_local << endl;
//cout << "Check: " << inv_sigma_k_prime_sq << endl;
//cout << "Check: " << log_correl_over_sigma_sq << endl;

		qweight[0] = - 1.0;
		qweight[1] = Q_local*Q_local;

        for(int ij = 0; ij < dim; ij++)
        {
            V[ij] += qweight[ij]*log_correl_over_sigma_sq;
            T[0][ij] += qweight[ij]*inv_sigma_k_prime_sq;
        }

        for(int ij = 1; ij < dim; ij++)
            T[ij][0] = T[0][ij];
            
        for(int ij = 1; ij < dim; ij++)
        {
            for(int lm = 1; lm < dim; lm++)
                T[ij][lm] += -qweight[ij]*qweight[lm]*inv_sigma_k_prime_sq;
        }

		// if we made it here, go ahead and fit
		++n_usable_bins;
    }

	// if we don't have enough usable bins, don't fit!
	if ( n_usable_bins <= dim )
		return;

    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            gsl_matrix_set(T_gsl, i, j, T[i][j]);

    // Make LU decomposition of matrix T_gsl
    gsl_linalg_LU_decomp (T_gsl, perm, &s_gsl);
    // Invert the matrix m
    gsl_linalg_LU_invert (T_gsl, perm, T_inverse_gsl);

    double ** T_inverse = new double * [dim];
    for(int i = 0; i < dim; i++)
    {
        T_inverse[i] = new double [dim];
        for(int j = 0; j < dim; j++)
            T_inverse[i][j] = gsl_matrix_get(T_inverse_gsl, i, j);
    }
    double * results = new double [dim];
    for(int i = 0; i < dim; i++)
    {
        results[i] = 0.0;
        for(int j = 0; j < dim; j++)
            results[i] += T_inverse[i][j]*V[j];
    }


	//=============================
	// determine errors on fit parameters
    double lambda 	= exp(results[0]);
	double R21D 	= results[1]*hbarC*hbarC;
	const int iK3D = indexerK(iKT, iKphi, iKL);

	//make error the value of the largest variation
	lambda_Correl_FRerr[iK3D] 	= max( abs( lambda - lambda_Correl[iK3D] ), lambda_Correl_FRerr[iK3D] );
	R2_FRerr[iK3D] 				= max( abs( R21D - R2[iK3D] ), R2_FRerr[iK3D] );



	//=============================
    // clean up
    gsl_matrix_free (T_gsl);
    gsl_matrix_free (T_inverse_gsl);
    gsl_permutation_free (perm);

    delete [] qweight;
    delete [] V;
    for(int i = 0; i < dim; i++)
    {
        delete [] T[i];
        delete [] T_inverse[i];
    }
    delete [] T;
    delete [] T_inverse;
    delete [] results;

	return;
}
*/

//End of file
