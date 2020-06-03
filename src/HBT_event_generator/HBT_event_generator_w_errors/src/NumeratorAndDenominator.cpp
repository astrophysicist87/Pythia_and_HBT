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
#include <iterator>
#include <chrono>

#include "HBT_event_generator.h"
#include "Arsenal.h"
#include "Stopwatch.h"


void HBT_event_generator::Compute_numerator_and_denominator_methodMode0_q_mode_3D()
{
	constexpr bool average_over_Kphi = true;

	constexpr bool oneDim_slices = false;

	constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );

	out << "HBT_event_generator:  * Computing correlation function:" << endl
		<< "HBT_event_generator:  * - method_mode = " << method_mode << endl
		<< "HBT_event_generator:  * - q_mode = 3D" << endl;

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<complex<double> > sum1(numerator.size());
		vector<double> sum2(numerator.size());
		vector<double> sum3(denominator.size());
		vector<double> sum4(denominator.size());
		vector<double> sum5(denominator.size());

		vector<int> private_ABC(numerator_bin_count.size());
		vector<int> private_DBC(denominator_bin_count.size());

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KTmin = KT_pts[iKT];
			double KTmax = KT_pts[iKT+1];
			double KT = 0.5*(KTmin+KTmax);
			double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
			double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double Kx = KT * cKphi;
			double Ky = KT * sKphi;
			double Kz = KL;

			// Sum over particles
			for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
			{
				ParticleRecord pi = event.particles[iParticle];

				double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
				//err << "made it here: "
				//	<< Kx << "   " << pi.px << "   " << 0.5*px_bin_width << "   "
				//	<< Ky << "   " << pi.py << "   " << 0.5*py_bin_width << "   "
				//	<< Kz << "   " << pi.pz << "   " << 0.5*pz_bin_width << endl;

				bool num_bin_true = 	abs( Kx - pi.px ) <= 0.5*px_bin_width
									and abs( Ky - pi.py ) <= 0.5*py_bin_width
									and abs( Kz - pi.pz ) <= 0.5*pz_bin_width;

				if ( num_bin_true )
				{
					double num_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;
					//num_bin_factor = 1.0;

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

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
						double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
						double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
						double qx = qo * cKphi - qs * sKphi;
						double qy = qs * cKphi + qo * sKphi;
						double qz = ql;

						double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
						double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;
						double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
						double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);

						// rapidity cuts
						if ( impose_pair_rapidity_cuts
								and ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min )
								or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
							)
							continue;

						double q0 = get_q0(particle_mass, qo, qs, ql, KT, KL);

						double arg =  q0 * ti - qx * xi - qy * yi - qz * zi;

						complex<double> complex_num_term = exp(i*arg/hbarC) / num_bin_factor;
						//err << "CHECKCOMPLEX " << arg/hbarC << "   "
						//	<< i*arg/hbarC << "   " << exp(i*arg/hbarC) << endl;

						sum1[index6D] += complex_num_term;
						sum2[index6D] += 1.0 / (num_bin_factor*num_bin_factor);
						private_ABC[index6D]++;

					}	// end of q loops

				}		// end of if block

			}			// end of particle loop

		}				// end of K loops




		//=====================================
		//======== Doing denominator ==========
		//=====================================

		// Loop over q and K bins first
		for (int iKT 	= 0; iKT 	< n_KT_bins; 	iKT++)
		for (int iKphi 	= 0; iKphi 	< n_Kphi_bins; 	iKphi++)
		for (int iKL 	= 0; iKL 	< n_KL_bins; 	iKL++)
		{
			double KT = 0.5*(KT_pts[iKT]+KT_pts[iKT+1]);
			double Kphi = 0.5*(Kphi_pts[iKphi]+Kphi_pts[iKphi+1]);
			double KL = 0.5*(KL_pts[iKL]+KL_pts[iKL+1]);
			double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double Kx = KT * cKphi;
			double Ky = KT * sKphi;
			double Kz = KL;

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

				double qo = 0.5*(qo_pts[iqo]+qo_pts[iqo+1]);
				double qs = 0.5*(qs_pts[iqs]+qs_pts[iqs+1]);
				double ql = 0.5*(ql_pts[iql]+ql_pts[iql+1]);
				double qx = qo * cKphi - qs * sKphi;
				double qy = qs * cKphi + qo * sKphi;
				double qz = ql;

				double pax = Kx + 0.5 * qx, pay = Ky + 0.5 * qy, paz = Kz + 0.5 * qz;
				double pbx = Kx - 0.5 * qx, pby = Ky - 0.5 * qy, pbz = Kz - 0.5 * qz;
				double Ea = sqrt(particle_mass*particle_mass+pax*pax+pay*pay+paz*paz);
				double Eb = sqrt(particle_mass*particle_mass+pbx*pbx+pby*pby+pbz*pbz);
				//double paT = sqrt(pax*pax+pay*pay), pbT = sqrt(pbx*pbx+pby*pby);
				//double pa_phi = atan2(pay,pax), pb_phi = atan2(pby,pbx);

				// rapidity cuts
				if ( impose_pair_rapidity_cuts
						and ( ( 2.0*Kz/(Ea+Eb) < Kz_over_K0_min )
						or ( 2.0*Kz/(Ea+Eb) > Kz_over_K0_max ) )
					)
					continue;

				// Sum over particles for pa bin
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pax) <= 0.5*px_bin_width
										and abs(pi.py-pay) <= 0.5*py_bin_width
										and abs(pi.pz-paz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;
						//den_bin_factor = 1.0;

						denominator_cell_was_filled[index6D] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum3[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum3 particle loop

				// Sum over particles for pb bin
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pbx) <= 0.5*px_bin_width
										and abs(pi.py-pby) <= 0.5*py_bin_width
										and abs(pi.pz-pbz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width;
						//den_bin_factor = 1.0;

						denominator_cell_was_filled[index6D] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum4[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum4 particle loop

				// Sum over particles for pa and pb bins simultaneously
				for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
				{
					ParticleRecord pi = event.particles[iParticle];

					bool den_bin_true = 	abs(pi.px-pax) <= 0.5*px_bin_width
										and abs(pi.py-pay) <= 0.5*py_bin_width
										and abs(pi.pz-paz) <= 0.5*pz_bin_width
										and abs(pi.px-pbx) <= 0.5*px_bin_width
										and abs(pi.py-pby) <= 0.5*py_bin_width
										and abs(pi.pz-pbz) <= 0.5*pz_bin_width;

					if ( den_bin_true )
					{

						int index6D = indexer(iKT, iKphi, iKL, iqo, iqs, iql);

						// N.B. - use (bin volume)^2, here
						double den_bin_factor =
							px_bin_width*py_bin_width*pz_bin_width
							*px_bin_width*py_bin_width*pz_bin_width;
						//den_bin_factor = 1.0;

						denominator_cell_was_filled[index6D] = true;

						double den_term = 1.0 / den_bin_factor;	//no phase factor in denominator

						sum5[index6D] += den_term;
						private_DBC[index6D]++;

					}	// end of if block

				}		// end of sum5 particle loop

			}			// end of q loops

		}				// end of K loops

		// Need this to avoid race conditions
		#pragma omp critical
		{
			const int iqoC = (n_qo_pts - 1) / 2;
			const int iqsC = (n_qs_pts - 1) / 2;
			const int iqlC = (n_ql_pts - 1) / 2;

			int idx = 0;
			if ( average_over_Kphi )
			{
				for (int iKT = 0; iKT < n_KT_bins; iKT++)
				for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
				for (int iKL = 0; iKL < n_KL_bins; iKL++)
				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double abs_sum1 = abs(sum1[idx]);
					double numerator_contribution_from_this_event
							= abs_sum1*abs_sum1 - sum2[idx];
					double denominator_contribution_from_this_event
							= sum3[idx]*sum4[idx] - sum5[idx];

					int num_bin_count_from_this_event
							= private_ABC[idx];

					int den_bin_count_from_this_event
							= private_DBC[idx];

					for (int iKphiAVG = 0; iKphiAVG < n_Kphi_bins; iKphiAVG++)
					{
						int idxAVG = indexer(iKT, iKphiAVG, iKL, iqo, iqs, iql);

						// first moments
						numerator[idxAVG]
							+= numerator_contribution_from_this_event;
						denominator[idxAVG]
							+= denominator_contribution_from_this_event;

						// second moments
						numerator2[idxAVG]
							+= numerator_contribution_from_this_event
								* numerator_contribution_from_this_event;
						denominator2[idxAVG]
							+= denominator_contribution_from_this_event
								* denominator_contribution_from_this_event;
						numerator_denominator[idxAVG]
							+= numerator_contribution_from_this_event
								* denominator_contribution_from_this_event;

						// track total bin counts
						numerator_bin_count[idxAVG] += num_bin_count_from_this_event;
						denominator_bin_count[idxAVG] += den_bin_count_from_this_event;
					}

					++idx;
				}
			}
			else
			{
				for (int iKT = 0; iKT < n_KT_bins; iKT++)
				for (int iKphi = 0; iKphi < n_Kphi_bins; iKphi++)
				for (int iKL = 0; iKL < n_KL_bins; iKL++)
				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double abs_sum1 = abs(sum1[idx]);
					double numerator_contribution_from_this_event
							= abs_sum1*abs_sum1 - sum2[idx];
					double denominator_contribution_from_this_event
							= sum3[idx]*sum4[idx] - sum5[idx];
	//out << "den check: " << denominator_contribution_from_this_event
	//		 << "   " << sum3[idx] << "   " << sum4[idx] << "   " << sum5[idx] << endl;

					// first moments
					numerator[idx]
						+= numerator_contribution_from_this_event;
					denominator[idx]
						+= denominator_contribution_from_this_event;

					// second moments
					numerator2[idx]
						+= numerator_contribution_from_this_event
							* numerator_contribution_from_this_event;
					denominator2[idx]
						+= denominator_contribution_from_this_event
							* denominator_contribution_from_this_event;
					numerator_denominator[idx]
						+= numerator_contribution_from_this_event
							* denominator_contribution_from_this_event;

					// track number of events where bin count was non-vanishing
					//numerator_bin_count[idx] += int(private_ABC[idx] > 0);
					//denominator_bin_count[idx] += int(private_DBC[idx] > 0);
					// track total bin counts
					numerator_bin_count[idx] += private_ABC[idx];
					denominator_bin_count[idx] += private_DBC[idx];

					++idx;
				}
			}

			out << "HBT_event_generator:\t - finished "
				<< ++number_of_completed_events << " of "
				<< number_of_expected_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}

	}

	out << "HBT_event_generator:  * Finished!" << endl;

	return;
}














//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
//=====================================================================================
void HBT_event_generator::Compute_numerator_and_denominator_momentum_space_only_q_mode_3D()
{
	constexpr bool perform_random_rotation = false;
	constexpr bool perform_random_shuffle  = false;
	constexpr bool oneDim_slices           = true;
	constexpr bool use_LCMS                = true;

	out << "HBT_event_generator:  * Computing correlation function:" << endl
		<< "HBT_event_generator:  * - momentum-space only" << endl
		<< "HBT_event_generator:  * - 3D" << endl;

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	out << "HBT_event_generator:  * Starting at " << std::ctime(&start_time) << endl;

	double average_Npair_numerator = 0.0;
	double average_Nmixed_denominator = 0.0;

	/*constexpr bool impose_pair_rapidity_cuts = false;
	const double KYmin = -0.1, KYmax = 0.1;
	const double Kz_over_K0_min = tanh( KYmin );
	const double Kz_over_K0_max = tanh( KYmax );*/

	// Sum over all events
	#pragma omp parallel for schedule(static)
	for (int iEvent = 0; iEvent < allEvents.size(); ++iEvent)
	{
		EventRecord event = allEvents[iEvent];

		vector<double> private_num(numerator.size(), 0.0);
		vector<double> private_num2(numerator2.size(), 0.0);
		vector<double> private_numPair(numPair.size(), 0.0);
		vector<double> private_numPair2(numPair2.size(), 0.0);
		vector<double> private_den(denominator.size(), 0.0);
		vector<double> private_den2(denominator2.size(), 0.0);
		vector<double> private_denPair(denPair.size(), 0.0);
		vector<double> private_denPair2(denPair2.size(), 0.0);

		//===================================
		//======== Doing numerator ==========
		//===================================

		// Sum over pairs of distinct particles
		for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
		for (int jParticle = 0; jParticle < event.particles.size(); ++jParticle)
		{
			if (iParticle == jParticle)
				continue;

			ParticleRecord pi = event.particles[iParticle];
			ParticleRecord pj = event.particles[jParticle];

			//double ti = pi.t, xi = pi.x, yi = pi.y, zi = pi.z;
			//double tj = pj.t, xj = pj.x, yj = pj.y, zj = pj.z;
			double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;
			double Ej = pj.E, pjx = pj.px, pjy = pj.py, pjz = pj.pz;

			double K0 = 0.5*(Ei+Ej), Kz = 0.5*(piz+pjz);

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
				/*#pragma omp critical
				{
					out << "Check LCMS: "
                         << piz << "   " << Ei << "   "
                         << pjz << "   " << Ej << "   "
                         << piz_new << "   " << Ei_new << "   "
                         << pjz_new << "   " << Ej_new << "   "
                         << Kz << "   " << K0 << "   "
                         << 0.5*(piz+pjz) << "   " << 0.5*(Ei+Ej) << endl;
				}*/
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

			// New method of binning
			double Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy);
			double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

			double KT = sqrt(Kx*Kx+Ky*Ky);
			//double Kphi = atan2(Ky, Kx);
			//double cKphi = cos(Kphi), sKphi = sin(Kphi);
			double cKphi = Kx/(KT+1e-100), sKphi = Ky/(KT+1e-100);

			double qo = qx * cKphi + qy * sKphi;
			double qs = qy * cKphi - qx * sKphi;
			double ql = qz;

			// Get indices
			int KT_idx 	 = floor((KT - KT_min)/KT_bin_width);
			//int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
			int KL_idx 	 = floor((Kz - KL_min)/KL_bin_width);

			// Avoid problems if q goes out of range
			if (   qo < qo_min or qs < qs_min or ql < ql_min
				or qo > qo_max or qs > qs_max or ql > ql_max )
				continue;

			int qo_idx 	 = floor((qo - qo_min) / delta_qo);
			int qs_idx 	 = floor((qs - qs_min) / delta_qs);
			int ql_idx 	 = floor((ql - ql_min) / delta_ql);

			// allows to do only slices
			if ( oneDim_slices )
			{
				bool iqo_not_center = ( qo_idx != (n_qo_bins-1)/2 );
				bool iqs_not_center = ( qs_idx != (n_qs_bins-1)/2 );
				bool iql_not_center = ( ql_idx != (n_ql_bins-1)/2 );

				// if we're not on an axis slice, skip this q-bin
				if (
						( iqo_not_center and iqs_not_center )
						or ( iqo_not_center and iql_not_center )
						or ( iqs_not_center and iql_not_center )
					)
					continue;
			}

			// Momentum-space cuts
			if ( KT_idx < 0 or KT_idx >= n_KT_bins )
				continue;

			//if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
			//	continue;

			if ( !use_LCMS and ( KL_idx < 0 or KL_idx >= n_KL_bins ) )
				continue;

			if ( qo_idx < 0 or qo_idx >= n_qo_bins )
				continue;

			if ( qs_idx < 0 or qs_idx >= n_qs_bins )
				continue;

			if ( ql_idx < 0 or ql_idx >= n_ql_bins )
				continue;

			int Kphi_idx = floor((atan2(Ky, Kx) - Kphi_min)/Kphi_bin_width);
			if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				continue;

			int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
			int index6D = indexer(KT_idx, Kphi_idx, KL_idx, qo_idx, qs_idx, ql_idx);

			//double arg =  q0 * (ti - tj)
			//			- qx * (xi - xj)
			//			- qy * (yi - yj)
			//			- qz * (zi - zj);

			//double num_term = cos(arg/hbarC);
			double num_term = 1.0;

			private_num[index6D] += num_term;
			private_numPair[index3D]++;

		}


		//=====================================
		//========= Doing denominator =========
		//=====================================

		//------------------------------------
		// Randomly sample events to mix with.
		const unsigned int n_mixing_events = min( (int)allEvents.size()-1, n_mix_minimum );
		//const unsigned int n_mixing_events = allEvents.size()-1;
		//const unsigned int n_mixing_events = 100;

		vector<int> indices(allEvents.size());
		iota(indices.begin(), indices.end(), 0);

		/*if ( perform_random_shuffle
			or n_mixing_events < (int)allEvents.size()-1 )
			random_shuffle(indices.begin(), indices.end());*/
		vector<int> mixedEvents;
		if ( n_mixing_events < (int)allEvents.size()-1 )
			std::sample(indices.begin(), indices.end(), std::back_inserter(mixedEvents),
                n_mixing_events + 1, std::mt19937{std::random_device{}()});
		/*for (int mix_idx = 0; mix_idx <= n_mixing_events; ++mix_idx)
			if ( indices[mix_idx] != iEvent
					and mixedEvents.size() < n_mixing_events )
			{
				mixedEvents.push_back( indices[mix_idx] );
			}*/

		//--------------------------------
		// Randomly rotate sampled events.
		vector<double> random_angles(n_mixing_events, 0.0), cos_rand_angles, sin_rand_angles;
		if ( perform_random_rotation )
			get_random_angles(n_mixing_events, random_angles);
		for (int mix_idx = 0; mix_idx < random_angles.size(); ++mix_idx)
		{
			cos_rand_angles.push_back( cos( random_angles[mix_idx] ) );
			sin_rand_angles.push_back( sin( random_angles[mix_idx] ) );
		}

		//------------------------
		// Loop over mixed events.
		for (int jEvent = 0; jEvent < mixedEvents.size(); ++jEvent)
		{
			// Don't correlate this event with itself
			if ( iEvent == mixedEvents[jEvent] ) continue;
			
			#pragma omp critical
				cout << __FUNCTION__ << ":" << __LINE__ << ": " << iEvent << " mixing with " << mixedEvents[jEvent] << endl;

			EventRecord mixedEvent = allEvents[mixedEvents[jEvent]];

			double c_rand_phi = cos_rand_angles[jEvent],
					s_rand_phi = sin_rand_angles[jEvent];

			//-----------------------------
			// Sum over pairs of particles.
			for (int iParticle = 0; iParticle < event.particles.size(); ++iParticle)
			for (int jParticle = 0; jParticle < mixedEvent.particles.size(); ++jParticle)
			{

				ParticleRecord pi = event.particles[iParticle];
				ParticleRecord pj = mixedEvent.particles[jParticle];

				double Ei = pi.E, pix = pi.px, piy = pi.py, piz = pi.pz;

				//-----------------------------------------------------
				// Randomly rotate the mixed event in transverse plane.
				double Ej = pj.E,
						pjx = pj.px*c_rand_phi - pj.py*s_rand_phi,
						pjy = pj.px*s_rand_phi + pj.py*c_rand_phi,
						pjz = pj.pz;

	
				double K0 = 0.5*(Ei+Ej), Kz = 0.5*(piz+pjz);
	
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
					/*#pragma omp critical
					{
						out << "Check LCMS: "
                             << piz << "   " << Ei << "   "
                             << pjz << "   " << Ej << "   "
                             << piz_new << "   " << Ei_new << "   "
                             << pjz_new << "   " << Ej_new << "   "
                             << Kz << "   " << K0 << "   "
                             << 0.5*(piz+pjz) << "   " << 0.5*(Ei+Ej) << endl;
					}*/
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
	
				// New method of binning
				double Kx = 0.5*(pix+pjx), Ky = 0.5*(piy+pjy);
				double q0 = Ei-Ej, qx = pix-pjx, qy = piy-pjy, qz = piz-pjz;

				double KT = sqrt(Kx*Kx+Ky*Ky);
				//double Kphi = atan2(Ky, Kx);
				//double cKphi = cos(Kphi), sKphi = sin(Kphi);
				double cKphi = Kx/(KT+1e-100), sKphi = Ky/(KT+1e-100);

				double qo = qx * cKphi + qy * sKphi;
				double qs = qy * cKphi - qx * sKphi;
				double ql = qz;
	
				// If pair survived cuts, get indices
				int KT_idx 	= floor((KT - KT_min)/KT_bin_width);
				//int Kphi_idx = floor((Kphi - Kphi_min)/Kphi_bin_width);
				int KL_idx 	= floor((Kz - KL_min)/KL_bin_width);

				// Avoid problems if q goes out of range
				if (   qo < qo_min or qs < qs_min or ql < ql_min
					or qo > qo_max or qs > qs_max or ql > ql_max )
					continue;

				int qo_idx 	= floor((qo - qo_min) / delta_qo);
				int qs_idx 	= floor((qs - qs_min) / delta_qs);
				int ql_idx 	= floor((ql - ql_min) / delta_ql);

				// allows to do only slices
				if ( oneDim_slices )
				{
					bool iqo_not_center = ( qo_idx != (n_qo_bins-1)/2 );
					bool iqs_not_center = ( qs_idx != (n_qs_bins-1)/2 );
					bool iql_not_center = ( ql_idx != (n_ql_bins-1)/2 );
	
					// if we're not on an axis slice, skip this q-bin
					if (
							( iqo_not_center and iqs_not_center )
							or ( iqo_not_center and iql_not_center )
							or ( iqs_not_center and iql_not_center )
						)
						continue;
				}

				// Momentum-space cuts
				if ( KT_idx < 0 or KT_idx >= n_KT_bins )
					continue;

				//if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
				//	continue;

				if ( !use_LCMS and ( KL_idx < 0 or KL_idx >= n_KL_bins ) )
					continue;

				if ( qo_idx < 0 or qo_idx >= n_qo_bins )
					continue;

				if ( qs_idx < 0 or qs_idx >= n_qs_bins )
					continue;

				if ( ql_idx < 0 or ql_idx >= n_ql_bins )
					continue;

				int Kphi_idx = floor((atan2(Ky, Kx) - Kphi_min)/Kphi_bin_width);
				if ( Kphi_idx < 0 or Kphi_idx >= n_Kphi_bins )
					continue;

				int index3D = indexerK(KT_idx, Kphi_idx, KL_idx);
				int index6D = indexer(KT_idx, Kphi_idx, KL_idx, qo_idx, qs_idx, ql_idx);

				private_den[index6D]++;
				private_denPair[index3D]++;

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

				double private_numPair_val 		= private_numPair[idx3D];
				double private_denPair_val 		= private_denPair[idx3D];

				numPair[idx3D] 					+= private_numPair_val;
				numPair2[idx3D] 				+= private_numPair_val*private_numPair_val;
				denPair[idx3D] 					+= private_denPair_val;
				denPair2[idx3D] 				+= private_denPair_val*private_denPair_val;

				for (int iqo = 0; iqo < n_qo_bins; iqo++)
				for (int iqs = 0; iqs < n_qs_bins; iqs++)
				for (int iql = 0; iql < n_ql_bins; iql++)
				{
					double private_num_val 		= private_num[idx6D];
					double private_den_val 		= private_den[idx6D];

					numerator[idx6D] 			+= private_num_val;
					denominator[idx6D] 			+= private_den_val;

					numerator2[idx6D] 			+= private_num_val*private_num_val;
					numerator_numPair[idx6D] 	+= private_num_val*private_numPair_val;
					denominator2[idx6D] 		+= private_den_val*private_den_val;
					denominator_denPair[idx6D] 	+= private_den_val*private_denPair_val;

					++idx6D;
				}

				++idx3D;
			}

			out << "HBT_event_generator:\t - finished "
				<< ++number_of_completed_events << " of "
				<< number_of_expected_events << endl;
			//print_progressbar( static_cast<double>(++number_of_completed_events)
			//						/ static_cast<double>(total_N_events), err );
		}
		

	}

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    out << "HBT_event_generator:  * Ending at " << std::ctime(&end_time) << endl;

	out << "HBT_event_generator:  * Finished!" << endl;

	return;
}




