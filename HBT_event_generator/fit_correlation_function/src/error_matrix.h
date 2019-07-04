#ifndef ERROR_MATRIX_H
#define ERROR_MATRIX_H

#include<iostream>
#include<string>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

namespace error_matrix
{

	int n, p;

	double * Y, * ThetaHat;
	double ** Phi, ** W;

	void get_errors(vector<double> & Y_in, vector<double> & W_in, vector<double> & ThetaHat_in, vector<double> & Phi_in)
	{
		// set Y, Phi, W, and ThetaHat in terms of input quantities
		n = Y_in.size();
		p = ThetaHat_in.size();

		Y = new double [n];
		ThetaHat = new double [p];
		Phi = new double * [n];
		W = new double * [n];

		for (int in = 0; in < n; in++)
		{
			Phi[in] = new double [p];
			W[in] = new double [n];
		}

		// evaluate Y - Phi * ThetaHat (result is vector)


		// evaluate ( Y - Phi * ThetaHat ) * W * ( Y - Phi * ThetaHat )^T


		// evaluate Phi^T * W * Phi


		// invert previous result


		// multiply two factors together



		return;

	}





}

#endif
