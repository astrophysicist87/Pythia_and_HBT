#ifndef FIT_V2_H
#define FIT_V2_H

#if GSL_MAJOR_VERSION < 2
	
	#include <stdlib.h>
	#include <stdio.h>
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
	
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_blas.h>
	#include <gsl/gsl_multifit_nlinear.h>
	
	#define N      100    /* number of data points to fit */
	#define TMAX   (40.0) /* time variable in [0,TMAX] */
	
	struct data {
	  size_t n;
	  double * t;
	  double * y;
	};
	
	int
	expb_f (const gsl_vector * x, void *data,
	        gsl_vector * f);
	
	int
	expb_df (const gsl_vector * x, void *data,
	         gsl_matrix * J);
	
	void
	callback(const size_t iter, void *params,
	         const gsl_multifit_nlinear_workspace *w);
	void fit_driver();

#endif

#endif
