// BoseEinstein.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BoseEinsten class.

#include "Pythia8/BoseEinstein.h"

namespace Pythia8 {

//==========================================================================

// The BoseEinstein class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

///===CJP(begin)===

const int npts = 7;
const double x_pts_7[7] = { -0.94910791234275852,
							-0.741531185599394440,
							-0.40584515137739717,
							0.0,
							0.40584515137739717,
							0.74153118559939444,
							0.94910791234275852 };

const double x_wts_7[7] = { 0.1294849661688696933,
							0.2797053914892766679,
							0.3818300505051189449,
							0.4179591836734693878,
							0.3818300505051189449,
							0.2797053914892766679,
							0.1294849661688696933 };

/*const double x_pts_51[51] = {-0.9989099908489035, -0.9942612604367526, -0.9859159917359030, 
							-0.9739033680193239, -0.9582678486139082, -0.9390675440029624, 
							-0.9163738623097802, -0.8902712180295273, -0.8608567111822924, 
							-0.8282397638230648, -0.7925417120993812, -0.7538953544853755, 
							-0.7124444575770366, -0.6683432211753701, -0.6217557046007233, 
							-0.5728552163513038, -0.5218236693661858, -0.4688509042860411, 
							-0.4141339832263039, -0.3578764566884095, -0.3002876063353319, 
							-0.2415816664477987, -0.1819770269570775, -0.1216954210188888, 
							-0.0609611001505787, 0.0, 0.0609611001505787, 0.1216954210188888, 
							0.1819770269570775, 0.2415816664477987, 0.3002876063353319, 
							0.3578764566884095, 0.4141339832263039, 0.4688509042860411, 
							0.521823669366186, 0.572855216351304, 0.621755704600723, 
							0.668343221175370, 0.712444457577037, 0.753895354485376, 
							0.792541712099381, 0.828239763823065, 0.860856711182292, 
							0.890271218029527, 0.916373862309780, 0.939067544002962, 
							0.958267848613908, 0.973903368019324, 0.985915991735903, 
							0.994261260436753, 0.998909990848903};
const double x_wts_51[51] = {0.002796807171089895576, 0.006500337783252600292, 
							0.010185191297821729939, 0.01383263400647782230, 
							0.01742871472340105226, 0.02095998840170321058, 
							0.02441330057378143427, 0.02777579859416247720, 
							0.03103497129016000845, 0.03417869320418833624, 
							0.03719526892326029284, 0.04007347628549645319, 
							0.04280260799788008665, 0.04537251140765006875, 
							0.04777362624062310200, 0.04999702015005740978, 
							0.05203442193669708756, 0.05387825231304556143, 
							0.05552165209573869302, 0.05695850772025866210, 
							0.05818347398259214060, 0.05919199392296154378, 
							0.05998031577750325209, 0.06054550693473779514, 
							0.06088546484485634388, 0.0609989248412058802, 
							0.06088546484485634388, 0.06054550693473779514, 
							0.05998031577750325209, 0.05919199392296154378, 
							0.05818347398259214060, 0.05695850772025866210, 
							0.05552165209573869302, 0.05387825231304556143, 
							0.05203442193669708756, 0.04999702015005740978, 
							0.04777362624062310200, 0.04537251140765006875, 
							0.04280260799788008665, 0.04007347628549645319, 
							0.03719526892326029284, 0.03417869320418833624, 
							0.03103497129016000845, 0.02777579859416247720, 
							0.02441330057378143427, 0.02095998840170321058, 
							0.01742871472340105226, 0.01383263400647782230, 
							0.010185191297821729939, 0.006500337783252600292, 
							0.002796807171089895576};*/

// Enumeration of id codes and table for particle species considered.
const int    BoseEinstein::IDHADRON[9] = { 211, -211, 111, 321, -321,
                                           130,  310, 221, 331 };
const int    BoseEinstein::ITABLE[9] = { 0, 0, 0, 1, 1, 1, 1, 2, 3 };

// Distance between table entries, normalized to min( 2*mass, QRef).
const double BoseEinstein::STEPSIZE  = 0.05;
//const double BoseEinstein::STEPSIZE  = 0.01;

// Skip shift for two extremely close particles, to avoid instabilities.
const double BoseEinstein::Q2MIN     = 1e-8;

// Parameters of energy compensation procedure: maximally allowed
// relative energy error, iterative stepsize, and number of iterations.
const double BoseEinstein::COMPRELERR = 1e-10;
const double BoseEinstein::COMPFACMAX = 1000.;
const int    BoseEinstein::NCOMPSTEP  = 10;

//--------------------------------------------------------------------------

// Find settings. Precalculate table used to find momentum shifts.

bool BoseEinstein::init(Info* infoPtrIn, Settings& settings,
  ParticleData& particleData) {

  // Save pointer.
  infoPtr         = infoPtrIn;

  // Main flags.
  doPion   = settings.flag("BoseEinstein:Pion");
  doKaon   = settings.flag("BoseEinstein:Kaon");
  doEta    = settings.flag("BoseEinstein:Eta");

  // Shape of Bose-Einstein enhancement/suppression.
  lambda   = settings.parm("BoseEinstein:lambda");
  QRef     = settings.parm("BoseEinstein:QRef");
  enhanceMode
           = settings.parm("BoseEinstein:enhanceMode");

  // Masses of particles with Bose-Einstein implemented.
  for (int iSpecies = 0; iSpecies < 9; ++iSpecies)
    mHadron[iSpecies] = particleData.m0( IDHADRON[iSpecies] );

  // Pair pi, K, eta and eta' masses for use in tables.
  mPair[0] = 2. * mHadron[0];
  mPair[1] = 2. * mHadron[3];
  mPair[2] = 2. * mHadron[7];
  mPair[3] = 2. * mHadron[8];

  // Loop over the four required tables. Local variables.
  for (int iTab = 0; iTab < 4; ++iTab)
    m2Pair[iTab]      = mPair[iTab] * mPair[iTab];

  number_of_pairs = 0;
  number_of_shifted_pairs = 0;
  number_of_too_close_pairs = 0;
  number_of_too_separated_pairs = 0;

  // Done.
  return true;

}

//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

bool BoseEinstein::shiftEvent( Event& event) {
  // Reset list of identical particles.
  hadronBE.resize(0);

    //for (int i = 0; i < event.size(); ++i)
	//	cout << "CHECK EVENT: " << event[i].id() << "   " << event[i].p();

  // Loop over all hadron species with BE effects.
  nStored[0] = 0;
  for (int iSpecies = 0; iSpecies < 9; ++iSpecies) {
    nStored[iSpecies + 1] = nStored[iSpecies];
    if (!doPion && iSpecies <= 2) continue;
    if (!doKaon && iSpecies >= 3 && iSpecies <= 6) continue;
    if (!doEta  && iSpecies >= 7) continue;

    // Properties of current hadron species.
    int idNow = IDHADRON[ iSpecies ];
    int iTab  = ITABLE[ iSpecies ];

    // Loop through event record to store copies of current species.
    for (int i = 0; i < event.size(); ++i)
      if ( event[i].id() == idNow && event[i].isFinal() )
		{
        hadronBE.push_back(
          BoseEinsteinHadron( idNow, i, event[i].p(), event[i].m(), event[i].vProd() ) );
//if (idNow==211) cout << "CHECK pi+'s: " << i << " of " << event.size() << ":   " << event[i].m() << "   " << event[i].p();
		}
    nStored[iSpecies + 1] = hadronBE.size();

/*cout << "CHECK sizes (pid=" << idNow << "):" << endl;
cout << "nStored[" << iSpecies << "] = " << nStored[iSpecies] << endl;
cout << "nStored[" << iSpecies+1 << "]-1 = " << nStored[iSpecies+1] - 1 << endl;
cout << "nStored[" << iSpecies+1 << "] = " << nStored[iSpecies+1] << endl;*/

/*if (idNow==211)
for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
	cout << i1 << "   " << hadronBE[i1].p;*/

	// ======================================================
	// This block to estimate pair shift using pair density
	vector< pair< double, pair <int,int> > > sortedPairs;
	vector<double> pairShifts;
	if ( enhanceMode == 1 )
	{
		// then need to loop over sorted vector of PAIRS
		// --> construct this here
		//vector< pair< double, double > > den_pdf, num_pdf;

		// get all values in vector first
		for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
		for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			sortedPairs.push_back(
				std::make_pair(
					m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab],
					std::make_pair(i1, i2)
				)
			);

		//cout << "sortedPairs.size() = " << sortedPairs.size() << endl << endl;

		if (sortedPairs.size() < 2)
		{
			cout << "Not enough sorted pairs of species with pid=" << idNow << endl;
			continue;
		}


		// THEN sort them (sorts on first column in ascending order automatically)
		sort( sortedPairs.begin(), sortedPairs.end() );

		// add fake first "pair"
		sortedPairs.insert(sortedPairs.begin(), std::make_pair( 0.0, std::make_pair(-1, -1) ) );

		// add fake last "pair" (new QVal is 10% larger than last one, just for instance)
		sortedPairs.push_back( std::make_pair( 1.1*sortedPairs.back().first, std::make_pair(-1, -1) ) );

/*cout << "Check sortedPairs: " << endl;
for (int iPair = 0; iPair < (int)sortedPairs.size(); ++iPair)
	cout << sortedPairs[iPair].first << "   "
			<< sortedPairs[iPair].second.first << "   "
			<< sortedPairs[iPair].second.second << endl;*/


		// Then construct PDF from CDF for CF numerator and denominator
		// as functions of Q, using these pairs
		// Set CDFs first
		vector< pair< double, double > > den_cdf, num_cdf;
		// Do phase space integral analytically
		vector<double> antiderivative;

		// Initialize at Q=0
		// format: Q value, CDF at that value
		//den_cdf.push_back( std::make_pair(0.0, 0.0) );
		//num_cdf.push_back( std::make_pair(0.0, 0.0) );
		double running_den_sum = 0.0;
		double running_num_sum = 0.0;

		// Loop over all pairs (except final fake "pairs")
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			const double QVal = sortedPairs[iPair].first;
			const double sqrt_Q2_plus_4m2 = sqrt(QVal*QVal + m2Pair[iTab]);
			antiderivative.push_back( 0.5*QVal*sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( QVal + sqrt_Q2_plus_4m2 ) );
			// boost particles to pair rest frame
			// and extract mag. of spatial separation
			const int iPair1 = sortedPairs[iPair].second.first;
			const int iPair2 = sortedPairs[iPair].second.second;
			Vec4 xDiff = ( hadronBE[iPair1].x - hadronBE[iPair2].x ) * MM2FM / HBARC;
			//xDiff.bstback( 0.5*(hadronBE[iPair1].p + hadronBE[iPair2].p) );

			den_cdf.push_back( std::make_pair(QVal, running_den_sum) );
			num_cdf.push_back( std::make_pair(QVal, running_num_sum) );

			// Note coordinate separation in GeV^{-1}.
			running_den_sum += 1.0;
			//running_num_sum += 1.0 + sphericalbesselj0( QVal * xDiff.pAbs() );
			const double arg = ( hadronBE[iPair1].p - hadronBE[iPair2].p ) * xDiff;
			running_num_sum += 1.0 + /*static_cast<double>(arg <= 16.0 * M_PI) * */sphericalbesselj0( arg );
cout << "Check QVal, etc.: " << QVal << "   " << arg << endl;
		//<< "dp = " << hadronBE[iPair1].p - hadronBE[iPair2].p
		//<< "dx = " << xDiff;
		}


		// Add the final fake "pair" to CDFs ensure last real pair is handled correctly
		// (QVal shifted by some positive number, running sums are unchanged)
		den_cdf.push_back( std::make_pair( den_cdf.back().first+1.0, den_cdf.back().second ) );
		num_cdf.push_back( std::make_pair( num_cdf.back().first+1.0, num_cdf.back().second ) );

		// Estimate PDFs using linear interpolation of CDFs
		// and evaluate running integral
		running_den_sum = 0.0;
		running_num_sum = 0.0;
		vector<double> running_den_integral(1, running_den_sum);
		vector<double> running_num_integral(1, running_num_sum);

		// Skip final fake "pair" again
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			double leftQVal = sortedPairs[iPair].first;
			double rightQVal = sortedPairs[iPair+1].first;
			// PDF ~ slope of CDF
			const double den_pdf = (den_cdf[iPair+1].second - den_cdf[iPair].second) / (rightQVal - leftQVal);
			const double num_pdf = (num_cdf[iPair+1].second - num_cdf[iPair].second) / (rightQVal - leftQVal);

			/*cout << "CDFs and PDFs: "
					<< iPair << "   "
					<< sortedPairs[iPair].first << "   "
					<< den_cdf[iPair].second << "   "
					<< num_cdf[iPair].second << "   "
					<< den_pdf << "   "
					<< num_pdf << "   "
					<< running_den_sum << "   "
					<< running_num_sum << endl;*/

			// Update running integrals
			running_den_sum += den_pdf * ( antiderivative[iPair+1] - antiderivative[iPair] );
			running_num_sum += num_pdf * ( antiderivative[iPair+1] - antiderivative[iPair] );

			// Store 
			running_den_integral.push_back( running_den_sum );
			running_num_integral.push_back( running_num_sum );
		}

		// Store one more entry for fake final "pair"
		running_den_integral.push_back( running_den_sum );
		running_num_integral.push_back( running_num_sum );

/*		cout << "CDFs and PDFs: "
				<< (int)sortedPairs.size()-1 << "   "
				<< sortedPairs[(int)sortedPairs.size()-1].first << "   "
				<< den_cdf[(int)sortedPairs.size()-1].second << "   "
				<< num_cdf[(int)sortedPairs.size()-1].second << "   "
				<< 0.0 << "   "
				<< 0.0 << "   "
				<< running_den_sum << "   "
				<< running_num_sum << endl;

cout << "Some size checks: "
		<< sortedPairs.size() << "   "
		<< den_cdf.size() << "   "
		<< num_cdf.size() << "   "
		<< antiderivative.size() << "   "
		<< running_den_integral.size() << "   "
		<< running_num_integral.size() << endl;

cout << "Element checks: " << endl;
for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
	cout << sortedPairs[iPair].first << "   "
			<< sortedPairs[iPair].second.first << "   "
			<< sortedPairs[iPair].second.second << "   "
			<< den_cdf[iPair].first << "   "
			<< den_cdf[iPair].second << "   "
			<< num_cdf[iPair].first << "   "
			<< num_cdf[iPair].second << "   "
			<< antiderivative[iPair] << "   "
			<< running_den_integral[iPair] << "   "
			<< running_num_integral[iPair] << endl;*/

//if (true) exit(8);


		// Finally, get shifts for each pair (except for first and last ones)
		// First "pair" is just the Q=0 point
		// Last "pair" is just a repeat to ensure PDF = 0 for Q > Qmax
		for (int iPair = 1; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			//int thisPair = iPair+1;
			int thisPair = iPair;
			const double thisPairDen = running_den_integral[thisPair];		// > 0
			double leftPairNum = running_num_integral[thisPair];
			double rightPairNum = running_num_integral[thisPair+1];
			//double leftPairNum0 = leftPairNum;
			//double rightPairNum0 = rightPairNum;

			// shift left and right numerator entries until they contain
			// denominator value; then invert using linear interpolation
			//if ( leftPairNum > thisPairDen )								// most likely
				while ( leftPairNum > thisPairDen )
				{
					// shift to the left
					leftPairNum = running_num_integral[--thisPair];
					rightPairNum = running_num_integral[thisPair+1];
				}
			//else if ( rightPairNum < thisPairDen )							// also possible
				while ( rightPairNum < thisPairDen )
				{
					// shift to the right
					leftPairNum = running_num_integral[++thisPair];
					rightPairNum = running_num_integral[thisPair+1];
				}

//====================================
//====================================
//Place restrictions on range of Q allowed?  etc.?
//double check the alignment of arrays, etc.
//====================================
//====================================

			if ( thisPair < 0 /*or thisPair >= (int)running_den_integral.size()-1*/ )
				break;

			// should now have leftPairNum < thisPairDen < rightPairNum
			// --> get shift Q from linear interpolation (Qnew === Q + deltaQ)
			const double Q0 = sortedPairs[iPair].first;	// original Q
			const int iPair1 = sortedPairs[iPair].second.first;		// first hadron
			const int iPair2 = sortedPairs[iPair].second.second;	// second hadron
			const double leftQ = sortedPairs[thisPair].first;
			const double rightQ = sortedPairs[thisPair+1].first;
//cout << "CHECK: " << iPair << "   " << thisPair << " (of " << running_den_integral.size() << ")   " << leftPairNum0 << "   " << rightPairNum0 << ";   "
//		<< leftPairNum << " < " << thisPairDen << " < " << rightPairNum << endl;
			double Qnew = leftQ + (thisPairDen-leftPairNum)*(rightQ-leftQ)/(rightPairNum-leftPairNum);
			//if (idNow==211) cout << "deltaQ = " << Qnew - Q0 << "; " << Q0 << "   " << Qnew << endl;
			//const double deltaQ = Qnew - Q0;				// shift for this pair

			// linear interpolation not good enough; need to solve numerically with Newton-Raphson
			// use Qnew as initial guess
			// need num_pdf in this interval
			const double num_pdf = (num_cdf[thisPair+1].second - num_cdf[thisPair].second) / (rightQ - leftQ);
			const double Numi = num_pdf * antiderivative[thisPair];
			// Solve antiderivative(Qnew) == c0
			const double c0 = (thisPairDen - leftPairNum + Numi) / num_pdf;	// DOUBLE CHECK THIS AND PREVIOUS STEPS

			//cout << "num_pdf = " << num_pdf << endl;
			//cout << "Numi = " << Numi << endl;

			constexpr double ACCURACY = 1.e-6;
			constexpr int MAXTRIES = 10;
			double Nnew = rightPairNum;
			int ntries = 0;
			Qnew = rightQ;
			while ( abs( Nnew - c0 ) > ACCURACY and ntries < MAXTRIES )
			{
				const double sqrt_Q2_plus_4m2 = sqrt(Qnew*Qnew + m2Pair[iTab]);
				const double deriv = Qnew*Qnew / sqrt_Q2_plus_4m2;
				Nnew = 0.5*Qnew*sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( Qnew + sqrt_Q2_plus_4m2 );			
				//cout << "NR method: " << Qnew << "   " << Nnew << "   " << deriv << "   " << c0 << endl;
				Qnew -= (Nnew - c0) / deriv;

				ntries++;
			}

			/*if (idNow==211)*/ cout << "(idNow=" << idNow << ") deltaQ = " << Qnew - Q0 << "; " << Q0 << "   " << Qnew << " (ntries=" << ntries << ")" << endl;

			const double arg = ( hadronBE[iPair1].p - hadronBE[iPair2].p )
								* ( hadronBE[iPair1].x - hadronBE[iPair2].x )
								* MM2FM / HBARC;
			//if ( abs(arg) > 16.0 * M_PI or Q0 > 2.5 )
			if ( arg > 16.0 * M_PI or Q0 > 2.5 )
				Qnew = Q0;

			pairShifts.push_back( Qnew - Q0 );

//if (true) exit(8);

		}

		/*cout << "<<<=======================================================>>>" << endl;
		cout << "CHECKS:" << endl;
		for (int iPair = 0; iPair < (int)sortedPairs.size(); ++iPair)
		{
			cout << "thisPair = " << iPair << " of " << sortedPairs.size() << endl;
			cout << sortedPairs[iPair].first << "   ("
					<< sortedPairs[iPair].second.first << "   "
					<< sortedPairs[iPair].second.second << ")" << endl;
			const int i1 = sortedPairs[iPair].second.first;
			const int i2 = sortedPairs[iPair].second.second;
			cout << "\t hadronBE[i1].p = " << hadronBE[i1].p;
			cout << "\t hadronBE[i1].x = " << hadronBE[i1].x;
			cout << "\t hadronBE[i2].p = " << hadronBE[i2].p;
			cout << "\t hadronBE[i2].x = " << hadronBE[i2].x;

			cout << "deltaQ = " << pairShifts[iPair] << endl;

			cout << endl;
		}
		cout << "<<<=======================================================>>>" << endl;*/

				
	}

//if (true) exit(8);

    // Loop through pairs of identical particles and find shifts.
	switch ( enhanceMode )
	{
		case 0: // the original and the default
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_fixedQRef( i1, i2, iTab);
			break;
		case 1:	// use a new Gaussian enhancement based on space-time interval between production points
			//for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			//for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			//  shiftPair_STint_GaussBE( i1, i2, iTab);
			
			// Finally, discard fake pairs at beginning and end of sortedPairs
			sortedPairs.erase ( sortedPairs.begin() );
			sortedPairs.erase ( sortedPairs.end() );
//cout << "Check sizes again: " << sortedPairs.size() << "   " << pairShifts.size() << endl;
			shiftPairs(sortedPairs, pairShifts);
			break;
		/*case 2:	// use a new cos(q*x) enhancement based on space-time interval between production points
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_STint_SphBesselBE( i1, i2, iTab);
			break;*/
		default:
			// Do nothing.
			break;
	}
  }

  // Must have at least two pairs to carry out compensation.
  if (nStored[9] < 2) return true;

  // Shift momenta and recalculate energies.
  double eSumOriginal = 0.;
  double eSumShifted  = 0.;
  double eDiffByComp  = 0.;
  for (int i = 0; i < nStored[9]; ++i) {
//cout << "(particle#=" << i << "): p = " << hadronBE[i].p;
//cout << "(particle#=" << i << "): pShift = " << hadronBE[i].pShift;
    eSumOriginal  += hadronBE[i].p.e();
    hadronBE[i].p += hadronBE[i].pShift;
/*cout 	<< setprecision(8)
		<< "Original shift: " << i << "   " << hadronBE[i].pShift;
    hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );*/
    eSumShifted   += hadronBE[i].p.e();
    eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                     / hadronBE[i].p.e();
cout 	<< setprecision(8)
		<< "Getting ready to balance energy budget (particle#=" << i << "): "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl;

  }


cout 	<< setprecision(8)
		<< "Balancing energy budget: "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;

  // Iterate compensation shift until convergence.
  int iStep = 0;
  while ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
    && abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
    && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
cout 	<< setprecision(8)
		<< "Balancing energy budget: "
		<< compFac << "   " << eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;
    eSumShifted      = 0.;
    eDiffByComp      = 0.;
    for (int i = 0; i < nStored[9]; ++i) {
      hadronBE[i].p += compFac * hadronBE[i].pComp;
      hadronBE[i].pShift += compFac * hadronBE[i].pComp;
/*cout 	<< setprecision(8)
		<< "Net shift at this point: " << i << "   " << hadronBE[i].pShift;*/
      hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
      eSumShifted   += hadronBE[i].p.e();
      eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                       / hadronBE[i].p.e();
    }
  }


  // Error if no convergence, and then return without doing BE shift.
  // However, not grave enough to kill event, so return true.
  if ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal ) {
    infoPtr->errorMsg("Warning in BoseEinstein::shiftEvent: "
      "no consistent BE shift topology found, so skip BE");
if (true) exit(8);
    return true;
  }


  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE[i].iPos, 99);
    event[ iNew ].p( hadronBE[i].p );
  }

  // Done.
  return true;

}

//---------------------------------------------------------------------------
// Calculate shift and (unnormalized) compensation for pair using fixed QRef.

void BoseEinstein::shiftPair_fixedQRef( int i1, int i2, int iTab) {

	//======================================
	// Start of initializations

  // Set relevant scales.
	// Multiples and inverses (= "radii") of distance parameters in Q-space.
	QRef2    = 2. * QRef;
	QRef3    = 3. * QRef;
	R2Ref    = 1. / (QRef * QRef);
	R2Ref2   = 1. / (QRef2 * QRef2);
	R2Ref3   = 1. / (QRef3 * QRef3);

  // Set various tables on a per-pair basis.
  double Qnow, Q2now, centerCorr;
    // Step size and number of steps in normal table.
    deltaQ[iTab]      = STEPSIZE * min(mPair[iTab], QRef);
    nStep[iTab]       = min( 199, 1 + int(3. * QRef / deltaQ[iTab]) );
    maxQ[iTab]        = (nStep[iTab] - 0.1) * deltaQ[iTab];
    centerCorr        = deltaQ[iTab] * deltaQ[iTab] / 12.;

    // Construct normal table recursively in Q space.
    shift[iTab][0]    = 0.;
    for (int i = 1; i <= nStep[iTab]; ++i) {
      Qnow            = deltaQ[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift[iTab][i]  = shift[iTab][i - 1] + exp(-Q2now * R2Ref)
        * deltaQ[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }

    // Step size and number of steps in compensation table.
    deltaQ3[iTab]     = STEPSIZE * min(mPair[iTab], QRef3);
    nStep3[iTab]      = min( 199, 1 + int(9. * QRef / deltaQ3[iTab]) );
    maxQ3[iTab]       = (nStep3[iTab] - 0.1) * deltaQ3[iTab];
    centerCorr        = deltaQ3[iTab] * deltaQ3[iTab] / 12.;

    // Construct compensation table recursively in Q space.
    shift3[iTab][0]   = 0.;
    for (int i = 1; i <= nStep3[iTab]; ++i) {
      Qnow            = deltaQ3[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift3[iTab][i] = shift3[iTab][i - 1] + exp(-Q2now * R2Ref3)
        * deltaQ3[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }
	// End of initializations
	//======================================


  // Calculate old relative momentum.
  double Q2old = m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab];
  if (Q2old < Q2MIN) return;
  double Qold  = sqrt(Q2old);
  double psFac = sqrt(Q2old + m2Pair[iTab]) / Q2old;

  // Calculate new relative momentum for normal shift.
  double Qmove = 0.;
  if (Qold < deltaQ[iTab]) Qmove = Qold / 3.;
  else if (Qold < maxQ[iTab]) {
    double realQbin = Qold / deltaQ[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove = ( shift[iTab][intQbin] + inter * (shift[iTab][intQbin + 1]
      - shift[iTab][intQbin]) ) * psFac;
  }
  else Qmove = shift[iTab][nStep[iTab]] * psFac;
  double Q2new = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove), 2. / 3.);

	/*cout << setprecision(6)
			<< "CHECK MOMENTUM SHIFT: "
			<< Qold << "   " << sqrt(Q2new) << "   "
			<< sqrt(m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab]) << "   "
			<< hadronBE[i1].p << "   " << hadronBE[i2].p << "   "
			<< sqrt(R2Ref) << endl;*/

  // Calculate corresponding three-momentum shift.
  double Q2Diff    = Q2new - Q2old;
  double p2DiffAbs = (hadronBE[i1].p - hadronBE[i2].p).pAbs2();
  double p2AbsDiff = hadronBE[i1].p.pAbs2() - hadronBE[i2].p.pAbs2();
  double eSum      = hadronBE[i1].p.e() + hadronBE[i2].p.e();
  double eDiff     = hadronBE[i1].p.e() - hadronBE[i2].p.e();
  double sumQ2E    = Q2Diff + eSum * eSum;
  double rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  double rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  double factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Add shifts to sum. (Energy component dummy.)
  Vec4   pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pShift += pDiff;
  hadronBE[i2].pShift -= pDiff;

//cout << "(idNow=" << hadronBE[i1].id << ") deltaQ = " << setprecision(8) << sqrt(Q2new)-sqrt(Q2old) << endl;

  // Calculate new relative momentum for compensation shift.
  double Qmove3 = 0.;
  if (Qold < deltaQ3[iTab]) Qmove3 = Qold / 3.;
  else if (Qold < maxQ3[iTab]) {
    double realQbin = Qold / deltaQ3[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove3 = ( shift3[iTab][intQbin] + inter * (shift3[iTab][intQbin + 1]
      - shift3[iTab][intQbin]) ) * psFac;
  }
  else Qmove3 = shift3[iTab][nStep3[iTab]] *psFac;
  double Q2new3 = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove3), 2. / 3.);

  // Calculate corresponding three-momentum shift.
  Q2Diff    = Q2new3 - Q2old;
  sumQ2E    = Q2Diff + eSum * eSum;
  rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Extra dampening factor to go from BE_3 to BE_32.
  factor   *= 1. - exp(-Q2old * R2Ref2);

  // Add shifts to sum. (Energy component dummy.)
  pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pComp += pDiff;
  hadronBE[i2].pComp -= pDiff;

}

//---------------------------------------------------------------------------
// Calculate shift and (unnormalized) compensation for pair using space-time
// interval and Gaussian form.
void BoseEinstein::shiftPairs( vector< pair< double, pair <int,int> > > & sortedPairs,
								vector<double> & pairShifts)
{

	/*
	//======================================
	// Start of modified initializations

  // Set relevant scales.
	// Set width of BE enhancement using pair's coordinate separation.
	Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;

	// Coordinate separation in GeV^{-1}.
	xDiff *= MM2FM / HBARC;

	number_of_pairs++;

	// Check that QRef will not be too large or too small: 0.05 <= QRef <= 1.0
	if ( abs( xDiff.mCalc() ) < 1.0 or abs( xDiff.mCalc() ) > 20.0 )
	{
		if ( abs( xDiff.mCalc() ) < 1.0 )
			number_of_too_close_pairs++;
		else
			number_of_too_separated_pairs++;
		return;
	}

	number_of_shifted_pairs++;

	R2Ref   = abs( xDiff * xDiff );
	QRef     = 1 / sqrt(R2Ref);
	QRef2    = 2. * QRef;
	QRef3    = 3. * QRef;
	R2Ref2   = R2Ref / 4.0;
	R2Ref3   = R2Ref / 9.0;

  // Set various tables on a per-pair basis.
  double Qnow, Q2now, centerCorr;
    // Step size and number of steps in normal table.
    deltaQ[iTab]      = STEPSIZE * min(mPair[iTab], QRef);
    nStep[iTab]       = min( 199, 1 + int(3. * QRef / deltaQ[iTab]) );
    maxQ[iTab]        = (nStep[iTab] - 0.1) * deltaQ[iTab];
    centerCorr        = deltaQ[iTab] * deltaQ[iTab] / 12.;

    // Construct normal table recursively in Q space.
    shift[iTab][0]    = 0.;
    for (int i = 1; i <= nStep[iTab]; ++i) {
      Qnow            = deltaQ[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift[iTab][i]  = shift[iTab][i - 1] + exp(-Q2now * R2Ref)
        * deltaQ[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }

    // Step size and number of steps in compensation table.
    deltaQ3[iTab]     = STEPSIZE * min(mPair[iTab], QRef3);
    nStep3[iTab]      = min( 199, 1 + int(9. * QRef / deltaQ3[iTab]) );
    maxQ3[iTab]       = (nStep3[iTab] - 0.1) * deltaQ3[iTab];
    centerCorr        = deltaQ3[iTab] * deltaQ3[iTab] / 12.;

    // Construct compensation table recursively in Q space.
    shift3[iTab][0]   = 0.;
    for (int i = 1; i <= nStep3[iTab]; ++i) {
      Qnow            = deltaQ3[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift3[iTab][i] = shift3[iTab][i - 1] + exp(-Q2now * R2Ref3)
        * deltaQ3[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }
	// End of modified initializations
	//======================================


  // Calculate old relative momentum.
  double Q2old = m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab];
  if (Q2old < Q2MIN) return;
  double Qold  = sqrt(Q2old);
  double psFac = sqrt(Q2old + m2Pair[iTab]) / Q2old;

  // Calculate new relative momentum for normal shift.
  double Qmove = 0.;
  if (Qold < deltaQ[iTab]) Qmove = Qold / 3.;
  else if (Qold < maxQ[iTab]) {
    double realQbin = Qold / deltaQ[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove = ( shift[iTab][intQbin] + inter * (shift[iTab][intQbin + 1]
      - shift[iTab][intQbin]) ) * psFac;
  }
  else Qmove = shift[iTab][nStep[iTab]] * psFac;
  double Q2new = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove), 2. / 3.);*/

	int pairIndex = 0;
	for (const auto & iPair : sortedPairs)
	{
		const int i1 	   = iPair.second.first;
		const int i2 	   = iPair.second.second;
		const double Qold  = iPair.first;
		const double Qnew  = Qold + pairShifts[pairIndex++];
		const double Q2old = Qold*Qold;
		const double Q2new = Qnew*Qnew;

		// Calculate corresponding three-momentum shift.
		double Q2Diff    = Q2new - Q2old;
		double p2DiffAbs = (hadronBE[i1].p - hadronBE[i2].p).pAbs2();
		double p2AbsDiff = hadronBE[i1].p.pAbs2() - hadronBE[i2].p.pAbs2();
		double eSum      = hadronBE[i1].p.e() + hadronBE[i2].p.e();
		double eDiff     = hadronBE[i1].p.e() - hadronBE[i2].p.e();
		double sumQ2E    = Q2Diff + eSum * eSum;
		double rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
		double rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
		double factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
							+ Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

		// Add shifts to sum. (Energy component dummy.)
		Vec4   pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
		hadronBE[i1].pShift += pDiff;
		hadronBE[i2].pShift -= pDiff;

		hadronBE[i1].pComp += pDiff;
		hadronBE[i2].pComp -= pDiff;

	}

	// turn off this compensation shift for now
	/*
  // Calculate new relative momentum for compensation shift.
  double Qmove3 = 0.;
  if (Qold < deltaQ3[iTab]) Qmove3 = Qold / 3.;
  else if (Qold < maxQ3[iTab]) {
    double realQbin = Qold / deltaQ3[iTab];
    int    intQbin  = int( realQbin );
    double inter    = (pow3(realQbin) - pow3(intQbin))
      / (3 * intQbin * (intQbin + 1) + 1);
    Qmove3 = ( shift3[iTab][intQbin] + inter * (shift3[iTab][intQbin + 1]
      - shift3[iTab][intQbin]) ) * psFac;
  }
  else Qmove3 = shift3[iTab][nStep3[iTab]] *psFac;
  double Q2new3 = Q2old * pow( Qold / (Qold + 3. * lambda * Qmove3), 2. / 3.);

  // Calculate corresponding three-momentum shift.
  Q2Diff    = Q2new3 - Q2old;
  sumQ2E    = Q2Diff + eSum * eSum;
  rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Extra dampening factor to go from BE_3 to BE_32.
  factor   *= 1. - exp(-Q2old * R2Ref2);

  // Add shifts to sum. (Energy component dummy.)
  pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pComp += pDiff;
  hadronBE[i2].pComp -= pDiff;
	*/

}

//==========================================================================

} // end namespace Pythia8
