#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <random>
#include <algorithm>

#include "HBT_event_generator.h"
#include "Arsenal.h"
#include "Stopwatch.h"



///*
void HBT_event_generator::Compute_numerator_and_denominator_methodMode2_q_mode_3D()
{
	// ignore q-dependence in denominator
	//bool use_smoothness_approximation = false;
	const double SAfact = ( use_smoothness_approximation ) ? 0.0 : 1.0;

	constexpr bool perform_random_rotation = false;
	constexpr bool perform_random_shuffle  = false;
	//constexpr bool oneDim_slices           = true;
	//constexpr bool use_LCMS                = true;

	//int number_of_completed_events = 0;
	//out << "In file: " << __FILE__ << " and function: " << __FUNCTION__ << " at line: " << __LINE__ << endl;
	//out << "  * Computing numerator and denominator of correlation function with errors; qmode = 3D using bin-averaging" << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	const double KYmin = -0.5, KYmax = 0.5;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_num(numerator.size(), 0.0);
		vector<double> private_num2(numerator2.size(), 0.0);
		vector<double> private_den(denominator.size(), 0.0);
		vector<double> private_den2(denominator2.size(), 0.0);
		vector<double> private_num_den(numerator.size(), 0.0);

		//===================================
		//======== Doing numerator ==========
		//===================================

		double num_pairs_this_event = static_cast<double>(event.particles.size() * (event.particles.size() - 1));

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		//for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		for (int jParticle = iParticle + 1; jParticle < event.particles.size(); ++jParticle)
		{
			//if (iParticle == jParticle)
			//	continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// check which cell we're in
			double cx = pix - pjx;
			double cy = piy - pjy;
			double cz = piz - pjz;

			bool num_bin_false = 	   cx*cx > px_bin_width*px_bin_width
									or cy*cy > py_bin_width*py_bin_width
									or cz*cz > pz_bin_width*pz_bin_width;

			if ( num_bin_false )
				continue;

			double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
			double tj = pj.t, xj = pj.x, yj = pj.y, zj = pj.z;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);

			// boost pair if necessary
			if ( use_LCMS )
			{
				double betaL   = Kz / K0;
				double betaL2  = betaL*betaL;
				if (betaL2 >= 1.0) continue;
				double gamma   = 1.0/sqrt(1.0-betaL2);
				double piz_new = gamma*(piz - betaL*Ei);
				double Ei_new  = gamma*(Ei - betaL*piz);
				double pjz_new = gamma*(pjz - betaL*Ej);
				double Ej_new  = gamma*(Ej - betaL*pjz);

				piz            = piz_new;
				Ei             = Ei_new;
				pjz            = pjz_new;
				Ej             = Ej_new;

				K0 = 0.5*(Ei+Ej);
				Kz = 0.5*(piz+pjz);
				if ( abs(Kz) > 1e-4 )
				{
					#pragma omp critical
						err << "Something went wrong!!! Kz = " << Kz << endl;
					exit(8);
				}
			}

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double KL = Kz;

			// Get indices
			int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
			int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	= floor((KL - KL_min)/KL_bin_width);

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
			{
				//err << "KT; shouldn't have made it here!" << endl;
				continue;
			}

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
			{
				//err << "Kphi; shouldn't have made it here!" << endl;
				continue;
			}

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
			{
				//err << "KL; shouldn't have made it here!" << endl;
				continue;
			}

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				// allows to do only slices
				if ( oneDim_slices )
				{
					bool iqo_not_center = ( iqo != (n_qo_bins-1)/2 );
					bool iqs_not_center = ( iqs != (n_qs_bins-1)/2 );
					bool iql_not_center = ( iql != (n_ql_bins-1)/2 );
	
					// if we're not on an axis slice, skip this q-bin
					if (
							( iqo_not_center and iqs_not_center )
							or ( iqo_not_center and iql_not_center )
							or ( iqs_not_center and iql_not_center )
						)
						continue;
				}

				/*bool in_center = ( iqo == (n_qo_bins-1)/2
						and iqs == (n_qs_bins-1)/2
						and iql == (n_ql_bins-1)/2 );*/
				double overall_factor = 1.0;
				/*if ( not in_center and
						not use_smoothness_approximation )
							overall_factor = 2.0;*/

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);

				double arg =  q0 * (ti - tj)
							- qx * (xi - xj)
							- qy * (yi - yj)
							- qz * (zi - zj);

				double num_term = 2.0*cos(arg/hbarC);	// factor of 2 allows to collapse sum in half
				//						/ ( 1.0
				//							* px_bin_width
				//							* py_bin_width
				//							* pz_bin_width );

				private_num[index6D] += ( include_energy_factors ) ?
										overall_factor * num_term / ( Ei * Ej ) :
										overall_factor * num_term;
				//						/ ( K0*K0 - 0.25*q0*q0 );	// this factor is 1/(E1 E2)
				//						/ num_pairs_this_event;

			}
		}

		//=====================================
		//========= Doing denominator =========
		//=====================================

		// Sum over pairs of particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		//for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		for (int jParticle = iParticle + 1; jParticle < event.particles.size(); ++jParticle)
		{
			//if ( iParticle == jParticle )
			//	continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			// New method of binning
			double K0 = 0.5*(Ei+Ej), Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy), Kz = 0.5*(piz+pjz);

			// boost pair if necessary
			if ( use_LCMS )
			{
				double betaL   = Kz / K0;
				double betaL2  = betaL*betaL;
				if (betaL2 >= 1.0) continue;
				double gamma   = 1.0 / sqrt(1.0-betaL2);
				double piz_new = gamma*(piz - betaL*Ei);
				double Ei_new  = gamma*(Ei - betaL*piz);
				double pjz_new = gamma*(pjz - betaL*Ej);
				double Ej_new  = gamma*(Ej - betaL*pjz);

				piz            = piz_new;
				Ei             = Ei_new;
				pjz            = pjz_new;
				Ej             = Ej_new;

				K0 = 0.5*(Ei+Ej);
				Kz = 0.5*(piz+pjz);
				if ( abs(Kz) > 1e-4 )
				{
					#pragma omp critical
						err << "Something went wrong!!! Kz = " << Kz << endl;
					exit(8);
				}
			}

			double KT = sqrt(Kx*Kx+Ky*Ky);
			double Kphi = atan2(Ky, Kx);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double KL = Kz;

			// If pair survived cuts, get indices
			int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
			int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
				continue;

			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				continue;

			if ( KL_idx < 0 or KL_idx >= n_KL_bins )
				continue;

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);

			for (int iqo = 0; iqo < n_qo_bins; iqo++)
			for (int iqs = 0; iqs < n_qs_bins; iqs++)
			for (int iql = 0; iql < n_ql_bins; iql++)
			{

				// allows to do only slices
				if ( oneDim_slices )
				{
					bool iqo_not_center = ( iqo != (n_qo_bins-1)/2 );
					bool iqs_not_center = ( iqs != (n_qs_bins-1)/2 );
					bool iql_not_center = ( iql != (n_ql_bins-1)/2 );
	
					// if we're not on an axis slice, skip this q-bin
					if (
							( iqo_not_center and iqs_not_center )
							or ( iqo_not_center and iql_not_center )
							or ( iqs_not_center and iql_not_center )
						)
						continue;
				}

				//bool in_center = ( iqo == (n_qo_bins-1)/2 and iqs == (n_qs_bins-1)/2 and iql == (n_ql_bins-1)/2 );
				double overall_factor = 1.0;
				//if ( not use_smoothness_approximation )
				//{
				//	overall_factor = 0.5;
				//	if ( in_center )
				//		overall_factor = 1.0;
				//}

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);

				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				// check which cell we're in
				double cmx = pix - pjx - SAfact*qx;
				double cmy = piy - pjy - SAfact*qy;
				double cmz = piz - pjz - SAfact*qz;
				/*double cpx = pix - pjx + SAfact*qx;
				double cpy = piy - pjy + SAfact*qy;
				double cpz = piz - pjz + SAfact*qz;*/

				// modified binning condition
				bool this_pair_den_bin_false
						= 	   cmx*cmx > px_bin_width*px_bin_width
							or cmy*cmy > py_bin_width*py_bin_width
							or cmz*cmz > pz_bin_width*pz_bin_width;
				/*bool rev_pair_den_bin_false
						= 	   cpx*cpx > px_bin_width*px_bin_width
							or cpy*cpy > py_bin_width*py_bin_width
							or cpz*cpz > pz_bin_width*pz_bin_width;*/

				if ( this_pair_den_bin_false /*and rev_pair_den_bin_false*/ )
					continue;

				int index6D     = indexer(KT_idx, Kphi_idx, KL_idx, iqo, iqs, iql);
				int index6D_rev = indexer(KT_idx, Kphi_idx, KL_idx,
											(n_qo_bins-1)-iqo,
											(n_qs_bins-1)-iqs,
											(n_ql_bins-1)-iql );

				//double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

				//const double energy_factors = 1.0 / ( Ei * Ej );
				//const double energy_factors = 1.0 / ( K0*K0 - 0.25*q0*q0 );

				private_den[index6D]++;
				private_den[index6D_rev]++;
				//private_den[index6D] += ( include_energy_factors ) ?
				//						overall_factor * energy_factors :
				//						overall_factor;

			}
		}

		// Need this to avoid race conditions
		#pragma omp critical
		{
			int idx3D = 0, idx6D = 0;
			for (int iKT = 0; iKT < n_KT_bins; iKT++)
			for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
			for (int iKL = 0; iKL < n_KL_bins; iKL++)
			{
				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double numerator_this_event = private_num[idx6D];
					double denominator_this_event = private_den[idx6D];

					// first moments
					numerator[idx6D]
								+= numerator_this_event;
					denominator[idx6D]
								+= denominator_this_event;

					// second moments
					numerator2[idx6D]
								+= numerator_this_event
									* numerator_this_event;
					denominator2[idx6D]
								+= denominator_this_event
									* denominator_this_event;
					numerator_denominator[idx6D]
								+= numerator_this_event
									* denominator_this_event;

					++idx6D;
				}

				++idx3D;
			}

			++number_of_completed_events;

			out << "\t - finished "
					<< number_of_completed_events + total_N_events - allEvents.size()
					<< " of " << number_of_expected_events << endl;
		}

	}

	if ( number_of_completed_events == number_of_expected_events )
		out << "  * Finished " << total_N_events << " events so far!" << endl;

	return;
}
//*/



