#ifndef FIT_V1_H
#define FIT_V1_H

#if GSL_MAJOR_VERSION < 2
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_blas.h>
	#include <gsl/gsl_multifit_nlin.h>
	
	#define N 40
	
	struct data {
	   size_t n;
	   double * y;
	   double * sigma;
	 };
	 
	 int
	 expb_f (const gsl_vector * x, void *data, 
	         gsl_vector * f);
	 int
	 expb_df (const gsl_vector * x, void *data, 
	          gsl_matrix * J);
	 int
	 expb_fdf (const gsl_vector * x, void *data,
	           gsl_vector * f, gsl_matrix * J);
	
	void print_state (size_t iter, gsl_multifit_fdfsolver * s);
	
	void fit_driver();
#endif

#endif
