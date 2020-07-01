#include<iostream>
#include<sstream>
#include<string>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "correlation_function.h"
//#include "Stopwatch.h"
//#include "Arsenal.h"

using namespace std;

//https://github.com/astrophysicist87/iEBE-Plumberg/blob/fbec616f897940d49e6774b04c0b16252edaf42b/EBE-Node/HoTCoffeeh/cfwr/src/cfwr_GFroutines.cpp

//https://www.gnu.org/software/gsl/doc/html/multimin.html#examples

//https://arxiv.org/pdf/nucl-ex/0505014.pdf

//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.66.054906

//https://indico.cern.ch/category/6015/attachments/192/631/Statistics_Fitting_II.pdf

// Minimize log-likelihood function
double Correlation_function::LogL_PML_f(const gsl_vector *v, void *params)
{
	double norm, lambda, R2o, R2s, R2l, R2os, R2ol, R2sl;
	correlationfunction_data * p = (correlationfunction_data *)params;
	
	norm   = gsl_vector_get(v, 0);
	lambda = gsl_vector_get(v, 1);
	R2o    = gsl_vector_get(v, 2);
	R2s    = gsl_vector_get(v, 3);
	R2l    = gsl_vector_get(v, 4);
	R2os   = gsl_vector_get(v, 5);
	R2ol   = gsl_vector_get(v, 6);
	R2sl   = gsl_vector_get(v, 7);

	double chi2PML    = 0.0;
	size_t datalength = p->datalength;
	for (int k = 0; k < datalength; k++)
	{
		double qo = p->qo[k];
		double qs = p->qs[k];
		double ql = p->ql[k];
		double Ak = p->A[k];
		double Bk = p->B[k];

		double qi_qj_R2ij = qo*qo*R2o + qs*qs*R2s + ql*ql*R2l
                            +2.0*qo*qs*R2os+2.0*qo*ql*R2ol+2.0*qs*ql*R2sl;
		double exp_arg    = exp( -qi_qj_R2ij / (hbarC*hbarC) );
		double Ck         = norm * ( 1.0 + lambda * exp_arg );
		double arg1       = Ck*(Ak+Bk)/(Ak*(Ck+1.0));
		double arg2       = (Ak+Bk)/(Bk*(Ck+1.0));

		chi2PML          += Ak*log(arg1) + Bk*log(arg2);
	}

	return (-2.0 * chi2PML);
}



/* The gradient of f, df = (df/dx, df/dy). */
void Correlation_function::LogL_PML_df (const gsl_vector *v, void *params, gsl_vector *df)
{
	double norm, lambda, R2o, R2s, R2l, R2os, R2ol, R2sl;
	correlationfunction_data * p = (correlationfunction_data *)params;
	
	norm   = gsl_vector_get(v, 0);
	lambda = gsl_vector_get(v, 1);
	R2o    = gsl_vector_get(v, 2);
	R2s    = gsl_vector_get(v, 3);
	R2l    = gsl_vector_get(v, 4);
	R2os   = gsl_vector_get(v, 5);
	R2ol   = gsl_vector_get(v, 6);
	R2sl   = gsl_vector_get(v, 7);

	double dchi2PML_dnorm   = 0.0;
	double dchi2PML_dlambda = 0.0;
	double dchi2PML_dR2o    = 0.0;
	double dchi2PML_dR2s    = 0.0;
	double dchi2PML_dR2l    = 0.0;
	double dchi2PML_dR2os   = 0.0;
	double dchi2PML_dR2ol   = 0.0;
	double dchi2PML_dR2sl   = 0.0;

	size_t datalength = p->datalength;
	for (int k = 0; k < datalength; k++)
	{
		double qo = p->qo[k];
		double qs = p->qs[k];
		double ql = p->ql[k];
		double Ak = p->A[k];
		double Bk = p->B[k];

		double qi_qj_R2ij = qo*qo*R2o + qs*qs*R2s + ql*ql*R2l
                            +2.0*qo*qs*R2os+2.0*qo*ql*R2ol+2.0*qs*ql*R2sl;
		double exp_arg    = exp( -qi_qj_R2ij / (hbarC*hbarC) );
		double Ck         = norm * ( 1.0 + lambda * exp_arg );

		double dCdnorm   = 1.0+lambda*exp_arg;
		double dCdlambda = norm*exp_arg;
		double dCdR2o    = -norm*lambda*exp_arg*qo*qo/(hbarC*hbarC);
		double dCdR2s    = -norm*lambda*exp_arg*qs*qs/(hbarC*hbarC);
		double dCdR2l    = -norm*lambda*exp_arg*ql*ql/(hbarC*hbarC);
		double dCdR2os   = -2.0*norm*lambda*exp_arg*qo*qs/(hbarC*hbarC);
		double dCdR2ol   = -2.0*norm*lambda*exp_arg*qo*ql/(hbarC*hbarC);
		double dCdR2sl   = -2.0*norm*lambda*exp_arg*qs*ql/(hbarC*hbarC);

		double prefactor = -2.0*(Ak - Bk*Ck)/(Ck*(Ck+1.0));

		dchi2PML_dnorm   += prefactor*dCdnorm;
		dchi2PML_dlambda += prefactor*dCdlambda;
		dchi2PML_dR2o    += prefactor*dCdR2o;
		dchi2PML_dR2s    += prefactor*dCdR2s;
		dchi2PML_dR2l    += prefactor*dCdR2l;
		dchi2PML_dR2os   += prefactor*dCdR2os;
		dchi2PML_dR2ol   += prefactor*dCdR2ol;
		dchi2PML_dR2sl   += prefactor*dCdR2sl;

	}

	// Store results for gradient
	gsl_vector_set(df, 0, dchi2PML_dnorm);
	gsl_vector_set(df, 1, dchi2PML_dlambda);
	gsl_vector_set(df, 2, dchi2PML_dR2o);
	gsl_vector_set(df, 3, dchi2PML_dR2s);
	gsl_vector_set(df, 4, dchi2PML_dR2l);
	gsl_vector_set(df, 5, dchi2PML_dR2os);
	gsl_vector_set(df, 6, dchi2PML_dR2ol);
	gsl_vector_set(df, 7, dchi2PML_dR2sl);

	return;
}

/* Compute both f and df together. */
void Correlation_function::LogL_PML_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = LogL_PML_f(x, params);
	LogL_PML_df(x, params, df);
	return;
}

void Correlation_function::set_CFdata(correlationfunction_data & CFdata, int iKT, int iKphi, int iKL)
{
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
			use_this_bin    = bool( i_in_center + j_in_center + k_in_center > 1 );
		}
		if ( not use_this_bin ) continue;

		if ( 	ignore_central_point
				and i==(n_qo_bins-1)/2
				and j==(n_qs_bins-1)/2
				and k==(n_ql_bins-1)/2)
			continue;

		int idx             = indexer(iKT, iKphi, iKL, i, j, k);

        double q_out_local  = 0.5*(qo_pts[i]+qo_pts[i+1]);
        double q_side_local = 0.5*(qs_pts[j]+qs_pts[j+1]);
        double q_long_local = 0.5*(ql_pts[k]+ql_pts[k+1]);

		CFdata.qo.push_back( q_out_local );
		CFdata.qs.push_back( q_side_local );
		CFdata.ql.push_back( q_long_local );
		CFdata.A.push_back(  numCount[idx] );
		CFdata.B.push_back(  denCount[idx] );

		// count this bin
		n_usable_bins++;
    }

	CFdata.datalength = n_usable_bins;
	return;
}

void Correlation_function::fit_correlationfunction_minimum_log_likelihood(int iKT, int iKphi, int iKL)
{
	size_t iter = 0;
	int status;
	int dim = 8;
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	// Set data points here
	correlationfunction_data CFdata;
	set_CFdata(CFdata, iKT, iKphi, iKL);

	gsl_vector *x;
	gsl_multimin_function_fdf LogL_PM_func;
	
	LogL_PM_func.n      = dim;
	LogL_PM_func.f      = &Correlation_function::LogL_PML_f;
	LogL_PM_func.df     = &Correlation_function::LogL_PML_df;
	LogL_PM_func.fdf    = &Correlation_function::LogL_PML_fdf;
	LogL_PM_func.params = &CFdata;

	/* Starting point */
	x = gsl_vector_alloc (dim);
	gsl_vector_set (x, 0, 1.0);	// overall normalization
	gsl_vector_set (x, 1, 1.0);	// lambda
	gsl_vector_set (x, 2, 1.0);	// R2o  (fm)
	gsl_vector_set (x, 3, 1.0);	// R2s  (fm)
	gsl_vector_set (x, 4, 1.0);	// R2l  (fm)
	gsl_vector_set (x, 5, 0.0);	// R2os (fm)
	gsl_vector_set (x, 6, 0.0);	// R2ol (fm)
	gsl_vector_set (x, 7, 0.0);	// R2sl (fm)
	
	
	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc (T, dim);
	
	gsl_multimin_fdfminimizer_set (s, &LogL_PM_func, x, 0.01, 1e-4);
	
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		
		if (status)
			break;
		
		status = gsl_multimin_test_gradient (s->gradient, 1e-3);
		
		if (status == GSL_SUCCESS)
		printf ("Minimum found at:\n");
		
		printf ("%5d %.5f %.5f %10.5f\n", iter,
				gsl_vector_get (s->x, 0),
				gsl_vector_get (s->x, 1),
				gsl_vector_get (s->x, 2),
				gsl_vector_get (s->x, 3),
				gsl_vector_get (s->x, 4),
				gsl_vector_get (s->x, 5),
				gsl_vector_get (s->x, 6),
				gsl_vector_get (s->x, 7),
				s->f);
		
	}
	while (status == GSL_CONTINUE && iter < 100);
	
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
	
	return;
}

void Correlation_function::Fit_correlation_function_min_logL()
{
	out << "--> Getting HBT radii by minimum log-likelihood method" << endl;

	for (int iKT = 0; iKT < n_KT_bins; ++iKT)
	for (int iKphi = 0; iKphi < n_Kphi_bins; ++iKphi)
	for (int iKL = 0; iKL < n_KL_bins; ++iKL)
	{
		out << "  --> Fitting in pair-momentum bin: "
			<< "KT=[" << KT_pts[iKT] << ", " << KT_pts[iKT+1] << "] GeV, "
			<< "Kphi=[" << Kphi_pts[iKphi] << ", " << Kphi_pts[iKphi+1] << "], "
			<< "KL=[" << KL_pts[iKL] << ", " << KL_pts[iKL+1] << "] GeV"
			<< endl;
		fit_correlationfunction_minimum_log_likelihood( iKT, iKphi, iKL );
	}

	out << "--> Finished getting HBT radii by minimum log-likelihood method" << endl;

	return;
}


