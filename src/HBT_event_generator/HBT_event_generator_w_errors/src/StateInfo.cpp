#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "HBT_event_generator.h"
#include "Arsenal.h"
#include "Stopwatch.h"


void HBT_event_generator::Dump_state(string path)
{

	const int q_space_size = ( q_mode == 0 ) ?
								n_qo_bins*n_qs_bins*n_ql_bins :
								n_Q_bins;
	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;

	string outPairsFilename = path + "/state/pairs.dat";
	ofstream outPairs( outPairsFilename.c_str() );
	for ( int iK = 0; iK < K_space_size; iK++ )
		outPairs
			<< numPair[iK]  << "   "
			<< numPair2[iK] << "   "
			<< denPair[iK]  << "   "
			<< denPair2[iK] << endl;
	outPairs.close();

	string outDistributionsFilename = path + "/state/distributions.dat";
	ofstream outDistributions( outDistributionsFilename.c_str() );
	for ( int iqK = 0; iqK < K_space_size*q_space_size; iqK++ )
		outDistributions
			<< numerator[iqK]                   << "   "
			<< denominator[iqK]                 << "   "
			<< numerator2[iqK]                  << "   "
			<< denominator2[iqK]                << "   "
			<< numerator_denominator[iqK]       << "   "
			<< denominator_cell_was_filled[iqK] << "   "
			<< numerator_bin_count[iqK]         << "   "
			<< denominator_bin_count[iqK]       << "   "
			<< numerator_numPair[iqK]           << "   "
			<< denominator_denPair[iqK]         << endl;
	outDistributions.close();

	return;
}

//
//void HBT_event_generator::Load_state(string path)
//{
//	const int q_space_size = ( q_mode == 0 ) ?
//								n_qo_bins*n_qs_bins*n_ql_bins :
//								n_Q_bins;
//	const int K_space_size = n_KT_bins*n_Kphi_bins*n_KL_bins;
//	double dummy = 0.0;
//
//	ifstream inPairs( path + "/state/pairs.dat" );
//	for ( int iK = 0; iK < K_space_size; iK++ )
//	{
//		/*inPairs
//			>> numPair[iK]
//			>> numPair2[iK]
//			>> denPair[iK]
//			>> denPair2[iK];*/
//		inPairs >> dummy; numPair[iK]  += static_cast<int>( dummy );
//		inPairs >> dummy; numPair2[iK] += static_cast<int>( dummy );
//		inPairs >> dummy; denPair[iK]  += static_cast<int>( dummy );
//		inPairs >> dummy; denPair2[iK] += static_cast<int>( dummy );
//	}
//	inPairs.close();
//
//	ifstream inDistributions( path + "/state/distributions.dat" );
//	for ( int iqK = 0; iqK < K_space_size*q_space_size; iqK++ )
//	{
//		/*inDistributions
//			>> numerator[iqK]
//			>> denominator[iqK]
//			>> numerator2[iqK]
//			>> denominator2[iqK]
//			>> numerator_denominator[iqK]
//			//>> denominator_cell_was_filled[iqK]
//			>> dummy
//			>> numerator_bin_count[iqK]
//			>> denominator_bin_count[iqK]
//			>> numerator_numPair[iqK]
//			>> denominator_denPair[iqK];
//		denominator_cell_was_filled[iqK] = static_cast<bool>( dummy );*/
//		inDistributions >> dummy; numerator[iqK]                   += dummy;
//		inDistributions >> dummy; denominator[iqK]                 += dummy;
//		inDistributions >> dummy; numerator2[iqK]                  += dummy;
//		inDistributions >> dummy; denominator2[iqK]                += dummy;
//		inDistributions >> dummy; numerator_denominator[iqK]       += dummy;
//		inDistributions >> dummy; denominator_cell_was_filled[iqK] += static_cast<bool>( dummy );
//		inDistributions >> dummy; numerator_bin_count[iqK]         += static_cast<int>( dummy );
//		inDistributions >> dummy; denominator_bin_count[iqK]       += static_cast<int>( dummy );
//		inDistributions >> dummy; numerator_numPair[iqK]           += dummy;
//		inDistributions >> dummy; denominator_denPair[iqK]         += dummy;
//	}
//	inDistributions.close();
//
//	return;
//}

// End of file
