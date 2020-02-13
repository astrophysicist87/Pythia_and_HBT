// BoseEinstein.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BoseEinsten class.

#include "Pythia8/BoseEinstein.h"

#include <ostream>
#include <chrono>

namespace Pythia8 {

//==========================================================================

// The BoseEinstein class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

///===CJP(begin)===

/*const int npts = 7;
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
							0.1294849661688696933 };*/

const int npts = 15;
const double x_pts[15] = { -0.98799251802048543,
							  -0.93727339240070590,
							  -0.84820658341042722,
							  -0.72441773136017005,
							  -0.57097217260853885,
							  -0.39415134707756337,
							  -0.20119409399743452,
							  0.0,
							  0.20119409399743452,
							  0.39415134707756337,
							  0.57097217260853885,
							  0.7244177313601700,
							  0.8482065834104272,
							  0.9372733924007059,
							  0.9879925180204854};

const double x_wts[15] = { 0.03075324199611726835,
							  0.07036604748810812471,
							  0.10715922046717193501,
							  0.1395706779261543144,
							  0.1662692058169939336,
							  0.1861610000155622110,
							  0.1984314853271115765,
							  0.2025782419255612729,
							  0.1984314853271115765,
							  0.1861610000155622110, 
							  0.1662692058169939336,
							  0.1395706779261543144,
							  0.10715922046717193501,
							  0.07036604748810812471,
							  0.03075324199611726835};

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
  infoPtr                = infoPtrIn;

  // Main flags.
  doPion 				 = settings.flag("BoseEinstein:Pion");
  doKaon 				 = settings.flag("BoseEinstein:Kaon");
  doEta 				 = settings.flag("BoseEinstein:Eta");
  useInvariantSize 		 = settings.flag("BoseEinstein:useInvariantSourceSize");
  useDistribution 		 = settings.flag("BoseEinstein:useDistribution");
  useRelativeDistance 	 = settings.flag("BoseEinstein:useRelativeDistance");
  useRestFrame 			 = settings.flag("BoseEinstein:useRestFrame");
  include_phase_space	 = true;	// for right now
  linear_interpolate_CDF = false;	// for right now
  include_posDelQ_in_compensation   = true;

  // Shape of Bose-Einstein enhancement/suppression.
  lambda 				= settings.parm("BoseEinstein:lambda");
  QRef 					= settings.parm("BoseEinstein:QRef");
  sourceDimension 		= settings.parm("BoseEinstein:sourceDimension");
  enhanceMode 			= settings.parm("BoseEinstein:enhanceMode");

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

void BoseEinstein::set_QRef(int iSpecies)
{
	if ( useDistribution )
	{
		cout << "Attempting to fix QRef from source distribution...";
		switch ( sourceDimension )
		{
			case 1:
				QRef = get_1D_source_size(iSpecies);
				break;
			default:
				infoPtr->errorMsg("Warning in BoseEinstein::set_QRef: "
      "choice of QRef makes no sense, defaulting to 1D source size");
				QRef = get_1D_source_size(iSpecies);
				break;
		}
	}

	// If out of range, just use default value (type this correctly later)
	if ( QRef < 0.05 )
	{
		// Reverting to minimum
		QRef = 0.05;
		cout << "failed!  Reverting to minimum QRef = " << QRef << "!" << endl;
	}
	else if ( QRef > 1.0 )
	{
		// Reverting to maximum
		QRef = 1.0;
		cout << "failed!  Reverting to maximum QRef = " << QRef << "!" << endl;
	}
	else
		cout << "success!  Using QRef = " << QRef << "!" << endl;
}

double BoseEinstein::get_1D_source_size(int iSpecies)
{
	double result = 0.0;
	double count = 0.0;

	if ( useInvariantSize )
	{
		if ( useRelativeDistance )
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			{
				Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
				result += xDiff*xDiff;
				count += 1.0;
			}
		else
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
			{
				Vec4 x = hadronBE[i1].x;
				result += x*x;
				count += 1.0;
			}
	}
	else
	{
		if ( useRelativeDistance )
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			{
				Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
				if ( useRestFrame ) xDiff.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
				result += xDiff.pAbs2();
				count += 1.0;
			}
		else
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
			{
				Vec4 x = hadronBE[i1].x;
				if ( useRestFrame ) x.bstback( hadronBE[i1].p );
				result += x.pAbs2();
				count += 1.0;
			}
	}

	double RMSsize = sqrt(result / (count+1.e-100)) * MM2FM / HBARC;

	// return RMS source size
	return ( 1.0 / RMSsize );
}


//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

bool BoseEinstein::shiftEvent( Event& event) {
  // Reset list of identical particles.
  hadronBE.resize(0);

	auto start = std::chrono::system_clock::now();

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
		}
    nStored[iSpecies + 1] = hadronBE.size();

	// ======================================================
	// Define these to estimate pair shift using pair density
	vector< pair< double, pair <int,int> > > sortedPairs;
	vector<double> pairShifts, pairCompensationShifts;
	// ======================================================

    // Loop through pairs of identical particles and find shifts.
	switch ( enhanceMode )
	{
		case 0: // the original and the default
			set_QRef(iSpecies);
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_fixedQRef( i1, i2, iTab );
			break;
		case 1:	// use a new cos(q*x) enhancement based on space-time interval between production points (mode 2)
			{
				bool enoughPairsToProceed = getSortedPairs( sortedPairs, iSpecies );
				if ( enoughPairsToProceed )
				{
					//cout << "BoseEinsteinCheck: NPair = " << sortedPairs.size()-2 << endl;
					shiftPairs_mode1( sortedPairs, pairShifts, pairCompensationShifts, iTab );
				}
			}
			break;
		default:
			// Do nothing.
			break;
	}
  }

//cout << "Made it here (1)" << endl;

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
/*cout    << setprecision(12) << setw(16)
                << "Original p and pShift: " << i
                << " (id=" << hadronBE[i].id << ")" << endl
                << "p=" << hadronBE[i].p
                << "pShift=" << hadronBE[i].pShift;*/
    hadronBE[i].p += hadronBE[i].pShift;
    hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
    eSumShifted   += hadronBE[i].p.e();
    eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                     / hadronBE[i].p.e();
/*cout << "eDiffByComp at this step = " << eDiffByComp << endl;
cout << "pComp at this step = " << hadronBE[i].pComp;
cout << "p at this step = " << hadronBE[i].p;
cout << "e at this step = " << hadronBE[i].p.e() << endl;
cout 	<< setprecision(8)
		<< "Getting ready to balance energy budget (particle#=" << i << "): "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl;*/

  }

//cout << "Made it here (2)" << endl;


/*cout 	<< setprecision(8)
		<< "Balancing energy budget: "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;*/

  // Iterate compensation shift until convergence.
  int iStep = 0;
  while ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
    && abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
    && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
/*cout 	<< setprecision(8)
		<< "Balancing energy budget: "
		<< compFac << "   " << eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;*/
    eSumShifted      = 0.;
    eDiffByComp      = 0.;
    for (int i = 0; i < nStored[9]; ++i) {
      hadronBE[i].p += compFac * hadronBE[i].pComp;
/*      hadronBE[i].pShift += compFac * hadronBE[i].pComp;
cout 	<< setprecision(8)
		<< "Net shift at this point: " << i << "   " << hadronBE[i].pShift;*/
      hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
      eSumShifted   += hadronBE[i].p.e();
      eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                       / hadronBE[i].p.e();
    }
  }

//cout << "Made it here (3)" << endl;

  // Error if no convergence, and then return without doing BE shift.
  // However, not grave enough to kill event, so return true.
  if ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal ) {
    infoPtr->errorMsg("Warning in BoseEinstein::shiftEvent: "
      "no consistent BE shift topology found, so skip BE");
	cout << setprecision(16) << "BoseEinsteinCheck: This event did not pass! Check: "
			<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	infoPtr->setBECShifts( false );
    return true;
  }
  else
  {
	cout << setprecision(16) << "BoseEinsteinCheck: This event passes! Check: "
			<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	infoPtr->setBECShifts( true );
  }


  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE[i].iPos, 99);
    event[ iNew ].p( hadronBE[i].p );
  }

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
std::cout << "BoseEinsteinCheck: elapsed time: " << elapsed_seconds.count() << " s" << "\n";

  // Done.
  return true;

}


//--------------------------------
// Construct list of sorted pairs.

bool BoseEinstein::getSortedPairs(
		vector< pair< double, pair <int,int> > > & sortedPairs,
		int iSpecies )
{
	// Reset.
    int iTab  = ITABLE[ iSpecies ];
	sortedPairs.clear();
	//sortedPairs_xDiffPRF.clear();

	// get all values in vector first
	for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
	for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
		sortedPairs.push_back(
			std::make_pair(
				sqrt( m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab] ),
				std::make_pair(i1, i2)
			)
		);

	// check if there are enough pairs of this species to do shift
	if (sortedPairs.size() < 2)
	{
		//cout << "Not enough sorted pairs of species with pid=" << idNow << endl;
		return false;
	}

	// THEN sort them (sorts on first column in ascending order automatically)
	sort( sortedPairs.begin(), sortedPairs.end() );

	// add fake first "pair"
	sortedPairs.insert(sortedPairs.begin(), std::make_pair( 0.0, std::make_pair(-1, -1) ) );

/*
cout << "Check sortedPairs: " << endl;
for (const auto & iPair : sortedPairs)
{
	const int i1 = iPair.second.first;
	const int i2 = iPair.second.second;
	Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
	xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );

	double thisQ = iPair.first;
	double nextQ = (&iPair == &sortedPairs.back() ) ? 1.1*iPair.first : (*(&iPair+1)).first;

	cout << setprecision(12) << thisQ << "   " << nextQ - thisQ << "   " << xDiffPRF.pAbs() * HBARC << "   " << xDiffPRF.pAbs() << "   " << m2Pair[iTab] << endl;
}

if (1) exit(8);
*/

	// add fake last "pair" (new QVal is 10% larger than last one, just for definiteness)
	sortedPairs.push_back( std::make_pair( 1.1*sortedPairs.back().first, std::make_pair(-1, -1) ) );

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


///*
double BoseEinstein::compute_integral_with_phasespace(double a_in, double b_in, vector<double> & cvec_in, double d_in)
{
	// For each c in cvec, computes integral given by
	// Integrate[Q^2 j0[c Q] / Sqrt[Q^2 + d], {Q, a, b}]
	double a = a_in, b = b_in, d = d_in;
	double overallSign = 1.0;
	if (b<a)
	{
		a=b_in;
		b=a_in;
		overallSign = -1.0;
	}
	double result = 0.0;

	const double cen = 0.5 * ( a + b );
	//const double sq_ad = sqrt(a*a+d), sq_bd = sqrt(b*b+d);
	const double sqrt_cen2_plus_d = sqrt(cen*cen + d);

	const double cen2 = cen*cen;
	const double cen3 = cen2*cen;
	const double cen5 = cen2*cen3;
	const double a2 = a*a, b2 = b*b, d2 = d*d;
	const double sqrt_cen2_plus_d_5 = sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d;
	const double approx_lhs
					= ( 2.0*cen5 - cen3*d + 8.0*cen2*d*a
						+ 2.0*d2*a - 3.0*cen*d*a2
						) / ( 2.0*sqrt_cen2_plus_d_5 );
	const double approx_rhs
					= ( 2.0*cen5 - cen3*d + 8.0*cen2*d*b
						+ 2.0*d2*b - 3.0*cen*d*b2
						) / ( 2.0*sqrt_cen2_plus_d_5 );
	const double exact_lhs = a / sqrt(a2 + d), exact_rhs = b / sqrt(b2 + d);
	const double eps_check1 = 0.001;

	if ( 		abs( 1.0-approx_lhs/exact_lhs ) < eps_check1
			and abs( 1.0-approx_rhs/exact_rhs ) < eps_check1 )
	{
		cout << "doing check 1" << endl;
		int localCount = 0;
		for (const auto & c : cvec_in)
		{
			// Enforce 'cluster' radius of maximum size.
			//if ( a*c > M_PI ) break;
			if ( b*c > 0.5*M_PI ) break;

			const double ac = a*c, bc = b*c;
			/*const double ac2 = ac*ac, bc2 = bc*bc;
			const double ac3 = ac*ac2, bc3 = bc*bc2;
			const double ac4 = ac*ac3, bc4 = bc*bc3;
			const double cos_ac = (15120 - 6900*ac2 + 313*ac4)/(15120 + 660*ac2 + 13*ac4);
			const double sin_ac = (5880*ac - 620*ac3)/(5880 + 360*ac2 + 11*ac4);
			const double cos_bc = (15120 - 6900*bc2 + 313*bc4)/(15120 + 660*bc2 + 13*bc4);
			const double sin_bc = (5880*bc - 620*bc3)/(5880 + 360*bc2 + 11*bc4);*/

			const double cos_ac = 1.00271 + ac*(-0.0428685 + ac*(-0.3822 + ac*(-0.12786 + ac*(0.110319 + ac*(-0.0191261 + 0.00101467*ac)))));
			const double sin_ac = 0.00635213 + ac*(0.901504 + ac*(0.243202 + ac*(-0.383693 + ac*(0.0854396 - 0.00543925*ac))));
			const double cos_bc = 1.00271 + bc*(-0.0428685 + bc*(-0.3822 + bc*(-0.12786 + bc*(0.110319 + bc*(-0.0191261 + 0.00101467*bc)))));
			const double sin_bc = 0.00635213 + bc*(0.901504 + bc*(0.243202 + bc*(-0.383693 + bc*(0.0854396 - 0.00543925*bc))));

			const double c2 = c*c;
			const double c4 = c2*c2;
			result +=	 ( 1.0/(2.0*c4*sqrt_cen2_plus_d_5) )
						*(    ( 6.0*cen*d + c2*( 2.0*cen5 - cen*(3.0*a2 - 8.0*a*cen + cen2)*d + 2.0*a*d2) )*cos_ac
							+ (-6.0*cen*d + c2*(-2.0*cen5 + cen*(3.0*b2 - 8.0*b*cen + cen2)*d - 2.0*b*d2) )*cos_bc
							+ 2.0*c*d*( - (-3.0*a*cen + 4.0*cen2 + d)*sin_ac
										+ (-3.0*b*cen + 4.0*cen2 + d)*sin_bc ) );
			//cout << "CHECK result: " << c << "   " << result << endl;
			localCount++;
		}
		cout << "Broke on localCount = " << localCount << endl;
	}
	else if ( a > 1000.0*d )	// if d is negligible compared to Q
	{
		cout << "doing check 2" << endl;
		for (const auto & c : cvec_in)
		{
			// Enforce 'cluster' radius of maximum size.
			//if ( a*c > 2.0*M_PI ) break;

			result += ( cos(a*c) - cos(b*c) ) / (c*c);
			//cout << "CHECK result: " << c << "   " << result << endl;
		}
	}
	else
	{
		cout << "doing full version" << endl;
		for (const auto & c : cvec_in)
		{
			// Enforce 'cluster' radius of maximum size.
			//if ( a*c > 2.0*M_PI ) break;

			const double period = 2.0*M_PI/c;
			double tmp_result = 0.0;

			// While limits contain more than one full period,
			// do one period at a time
			while ( b - a > period )
			{
				double a0 = a, a1 = a0 + period;
				double hw = 0.5*period, center = 0.5*(a0+a1);
				for (int i = 0; i < npts; ++i)
				{
					double Qloc = center + hw * x_pts[i];
					tmp_result += hw * x_wts[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
				}
				//cout << setprecision(16) << "CHECK result: " << c << "   " << tmp_result << endl;
				a += period;
			}
	
			// Finally, do remaining (fraction of a) period
			{
				double a0 = a, a1 = b;
				double hw = 0.5*(a1-a0), center = 0.5*(a0+a1);
				for (int i = 0; i < npts; ++i)
				{
					double Qloc = center + hw * x_pts[i];
					tmp_result += hw * x_wts[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
				}
				//cout << setprecision(16) << "CHECK result: " << c << "   " << tmp_result << endl;
			}

			//cout << setprecision(16) << "CHECK result: " << a << "   " << b << "   " << c << "   " << tmp_result << endl;
			result += tmp_result;
		}
	}
//if (1) exit(8);

	return ( overallSign*result );
}






double BoseEinstein::compute_integral_with_phasespace(double a_in, double b_in, double c_in, double d_in)
{
	// Computes integral given by Integrate[Q^2 j0[c Q] / Sqrt[Q^2 + d], {Q, a, b}]
	double a = a_in, b = b_in, c = c_in, d = d_in;
	double overallSign = 1.0;
	if (b<a)
	{
		a=b_in;
		b=a_in;
		overallSign = -1.0;
	}
	double result = 0.0;

	// Enforce 'cluster' radius of maximum size.
	//if ( a*c > 2.0*M_PI ) return ( 0.0 );

	///*
	//---------
	// Check 0.
	// Taylor expand part of phase space to second order
	// and do resulting integral analytically.
	{
		const double cen = 0.5 * ( a + b );
		const double sqrt_cen2_plus_d = sqrt(cen*cen + d);
		const double cen2 = cen*cen;
		const double cen3 = cen2*cen;
		const double cen5 = cen2*cen3;
		const double a2 = a*a, b2 = b*b, d2 = d*d;
		const double sqrt_cen2_plus_d_5 = sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d;
		const double approx_lhs
						= ( 2.0*cen5 - cen3*d + 8.0*cen2*d*a
							+ 2.0*d2*a - 3.0*cen*d*a2
							) / ( 2.0*sqrt_cen2_plus_d_5 );
		const double approx_rhs
						= ( 2.0*cen5 - cen3*d + 8.0*cen2*d*b
							+ 2.0*d2*b - 3.0*cen*d*b2
							) / ( 2.0*sqrt_cen2_plus_d_5 );
		const double exact_lhs = a / sqrt(a2 + d), exact_rhs = b / sqrt(b2 + d);
		const double eps_check0 = 0.001;
		if ( 		abs( 1.0-approx_lhs/exact_lhs ) < eps_check0
				and abs( 1.0-approx_rhs/exact_rhs ) < eps_check0 )
		{
//cout << "check 0" << endl;
			const double c2 = c_in*c_in;
			const double c4 = c2*c2;
/*cout 
	//<< scientific
	<< setprecision(12)
	<< a << "   "
	<< b << "   "
	<< c << "   "
	<< d << "   "
	<< cen << "   "
	<< c4 << "   "
	<< sqrt_cen2_plus_d_5 << "   "
	<< cos(a*c) << "   "
	<< cos(b*c) << "   "
	<< sin(a*c) << "   "
	<< sin(b*c) << "   "
	<< ( 6.0*cen*d + c2*( 2.0*cen5 - cen*(3.0*a2 - 8.0*a*cen + cen2)*d + 2.0*a*d2) )*cos(a*c) << "   "
	<< (-6.0*cen*d + c2*(-2.0*cen5 + cen*(3.0*b2 - 8.0*b*cen + cen2)*d - 2.0*b*d2) )*cos(b*c) << "   "
	<< (-3.0*a*cen + 4.0*cen2 + d)*sin(a*c) << "   "
	<< (-3.0*b*cen + 4.0*cen2 + d)*sin(b*c) << "   "
	<< overallSign *
		 ( 1.0/(2.0*c4*sqrt_cen2_plus_d_5) )
		*(    ( 6.0*cen*d + c2*( 2.0*cen5 - cen*(3.0*a2 - 8.0*a*cen + cen2)*d + 2.0*a*d2) )*cos(a*c)
			+ (-6.0*cen*d + c2*(-2.0*cen5 + cen*(3.0*b2 - 8.0*b*cen + cen2)*d - 2.0*b*d2) )*cos(b*c)
			+ 2.0*c*d*( - (-3.0*a*cen + 4.0*cen2 + d)*sin(a*c)
						+ (-3.0*b*cen + 4.0*cen2 + d)*sin(b*c) ) ) << endl;

if (1) exit(8);*/
			return ( overallSign *
								 ( 1.0/(2.0*c4*sqrt_cen2_plus_d_5) )
								*(    ( 6.0*cen*d + c2*( 2.0*cen5 - cen*(3.0*a2 - 8.0*a*cen + cen2)*d + 2.0*a*d2) )*cos(a*c)
									+ (-6.0*cen*d + c2*(-2.0*cen5 + cen*(3.0*b2 - 8.0*b*cen + cen2)*d - 2.0*b*d2) )*cos(b*c)
									+ 2.0*c*d*( - (-3.0*a*cen + 4.0*cen2 + d)*sin(a*c)
												+ (-3.0*b*cen + 4.0*cen2 + d)*sin(b*c) ) )
					);
		}
	}

	//---------
	// Check 1.
	// Lowest order Riemann 'sum' if b - a is small enough
	{
//cout << "check 1" << endl;
		const double cen = 0.5 * ( a + b );
		const double delta = b - a;

		const double cot_cenc = 1.0/tan(c*cen);
		const double cen2_plus_d = cen*cen + d;
		const double eps = delta*delta
							* ( 2.0*c*d*cen2_plus_d*cot_cenc
								- cen*(3.0*d+c*c*cen2_plus_d*cen2_plus_d) )
							/ ( 24.0*cen*cen2_plus_d*cen2_plus_d );
/*cout 
	//<< scientific
	<< setprecision(12)
	<< a << "   "
	<< b << "   "
	<< c << "   "
	<< d << "   "
	<< cen << "   "
	<< eps << "   "
	<< cen2_plus_d << "   "
	<< delta << "   "
	<< sphericalbesselj0(c*cen) << "   "
	<< sqrt(cen * cen + d) << "   "
	<< overallSign * delta * cen * cen * sphericalbesselj0(c*cen) / sqrt(cen * cen + d) << endl;*/

//if (1) exit(8);

		if ( abs(eps) < 0.00001 )
			return ( overallSign * delta * cen * cen * sphericalbesselj0(c*cen) / sqrt(cen * cen + d) );
	}

	//---------
	// Check 2.
	if ( a > 1000.0*d )	// do the integral approximately here
	{
//cout << "check 2" << endl;
		result = ( cos(a*c) - cos(b*c) ) / (c*c);
		return ( overallSign*result );
	}
	
	//---------
	// Check 3.
	const double sin_ac = sin(a*c), sin_bc = sin(b*c),
				 cos_ac = cos(a*c), cos_bc = cos(b*c);
	const double sq_ad = sqrt(a*a+d), sq_bd = sqrt(b*b+d);
	const double comp1 =   a*a*sq_bd*sq_bd*sq_bd*sin_ac
						 - b*b*sq_ad*sq_ad*sq_ad*sin_bc;
	const double comp2 = sq_ad*sq_ad*sq_bd*sq_bd*c
						 *(a*sq_bd*cos_ac - b*sq_ad*cos_bc);
	if ( abs(comp1) < 0.0001*abs(comp2) )		// do it approximately here too
	{
//cout << "check 3" << endl;
		const double term1 =  ( a * c * cos_ac - sin_ac ) / sq_ad;
		const double term2 = -( b * c * cos_bc - sin_bc ) / sq_bd;
		result = ( term1 + term2 ) / (c*c*c);
		return ( overallSign*result );
	}
	else
	{
//cout << "exact" << endl;
		const double period = 2.0*M_PI/c;

		// While limits contain more than one full period,
		// do one period at a time
		while ( b - a > period )
		{
			double a0 = a, a1 = a0 + period;
			double hw = 0.5*period, cen = 0.5*(a0+a1);
			for (int i = 0; i < npts; ++i)
			{
				double Qloc = cen + hw * x_pts[i];
				result += hw * x_wts[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
			}
			a += period;
		}
	
		// Finally, do remaining (fraction of a) period
		{
			double a0 = a, a1 = b;
			double hw = 0.5*(a1-a0), cen = 0.5*(a0+a1);
			for (int i = 0; i < npts; ++i)
			{
				double Qloc = cen + hw * x_pts[i];
				result += hw * x_wts[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
			}
		}
	}

/*
// Limit in which to do integral approximately
	//---------------------------------------
	// Stuff for check #1
	const double cen = 0.5 * ( a_in + b_in );
	//const double delta = abs(b_in - a_in);

	//const double cot_cenc = 1.0/tan(c_in*cen);
	const double cen2_plus_d = cen*cen + d_in;
	//const double eps = delta*delta
	//					* ( 2.0*c_in*d_in*cen2_plus_d*cot_cenc
	//						- cen*(3.0*d_in+c_in*c_in*cen2_plus_d*cen2_plus_d) )
	//					/ ( 24.0*cen*cen2_plus_d*cen2_plus_d );

	//---------------------------------------
	// Stuff for check #0
	const double sqrt_cen2_plus_d = sqrt(cen2_plus_d);
	const double cen2 = cen*cen;
	const double cen3 = cen2*cen;
	const double cen5 = cen2*cen3;
	const double sqrt_cen2_plus_d_5 = sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d*sqrt_cen2_plus_d;
	const double approx_lhs
					= ( 2.0*cen5 - cen3*d_in + 8.0*cen2*d_in*a_in
						+ 2.0*d_in*d_in*a_in - 3.0*cen*d_in*a_in*a_in
						) / ( 2.0*sqrt_cen2_plus_d_5 );
	const double approx_rhs
					= ( 2.0*cen5 - cen3*d_in + 8.0*cen2*d_in*b_in
						+ 2.0*d_in*d_in*b_in - 3.0*cen*d_in*b_in*b_in
						) / ( 2.0*sqrt_cen2_plus_d_5 );
	const double exact_lhs = a_in / sqrt(a_in*a_in + d_in), exact_rhs = b_in / sqrt(b_in*b_in + d_in);
	const double eps_check0 = 0.001;
if ( abs( 1.0-approx_lhs/exact_lhs ) < eps_check0 and abs( 1.0-approx_rhs/exact_rhs ) < eps_check0 )
{
	const double c2 = c_in*c_in;
	const double d2 = d_in*d_in;
	const double c4 = c2*c2;
	const double approximate_result
					= //overallSign *
						 ( 1.0/(2.0*c4*sqrt_cen2_plus_d_5) )
						*(    ( 6.0*cen*d + c2*( 2.0*cen5 - cen*(3.0*a_in*a_in - 8.0*a_in*cen + cen2)*d + 2.0*a_in*d2) )*cos(a_in*c)
							+ (-6.0*cen*d + c2*(-2.0*cen5 + cen*(3.0*b_in*b_in - 8.0*b_in*cen + cen2)*d - 2.0*b_in*d2) )*cos(b_in*c)
							+ 2.0*c*d*( - (-3.0*a_in*cen + 4.0*cen2 + d)*sin(a_in*c)
										+ (-3.0*b_in*cen + 4.0*cen2 + d)*sin(b_in*c) ) );
	const double error_estimate = 100.0*(approximate_result-result)/result;
	//if ( abs( error_estimate ) > 1.0 )
	{
		cout << "Test integration: " << setprecision(12) << setw(16) << a_in << "   " << b_in << "   " << c_in << "   " << d_in;
		cout << "  COMPARE (approx 0): " << result << "   " << approximate_result << "   "
				<< abs( 1.0-approx_lhs/exact_lhs ) << "   " << abs( 1.0-approx_rhs/exact_rhs )
				<< " (error of " << error_estimate << "%)" << endl;
	}

}
else if ( abs(eps) < 0.00001 )
{
	const double approximate_result = overallSign * delta * cen * cen * sphericalbesselj0(c_in*cen) / sqrt(cen * cen + d_in);
	const double error_estimate = 100.0*(approximate_result-result)/result;
	if ( abs( error_estimate ) > 1.0 )
	{
		cout << "Test integration: " << setprecision(12) << setw(16) << a_in << "   " << b_in << "   " << c_in << "   " << d_in;
		cout << "  COMPARE (approx 1): " << result << "   " << approximate_result  << " (error of " << error_estimate << "%)" << endl;
	}

}
else if ( a_in > 1000.0*d_in )
{
	const double approximate_result = ( cos(a_in * c_in) - cos(b_in * c_in) ) / (c_in*c_in);
	const double error_estimate = 100.0*(approximate_result-result)/result;
	if ( abs( error_estimate ) > 1.0 )
	{
		cout << "Test integration: " << setprecision(12) << setw(16) << a_in << "   " << b_in << "   " << c_in << "   " << d_in;
		cout << "  COMPARE (approx 2): " << result << "   " << approximate_result  << " (error of " << error_estimate << "%)" << endl;
	}
}
else
{
	const double sin_ac = sin(a_in*c_in), sin_bc = sin(b_in*c_in),
				 cos_ac = cos(a_in*c_in), cos_bc = cos(b_in*c_in);
	const double sq_ad = sqrt(a_in*a_in+d_in), sq_bd = sqrt(b_in*b_in+d_in);
	const double comp1 =   a_in*a_in*sq_bd*sq_bd*sq_bd*sin_ac
						 - b_in*b_in*sq_ad*sq_ad*sq_ad*sin_bc;
	const double comp2 = sq_ad*sq_ad*sq_bd*sq_bd*c_in
						 *(a_in*sq_bd*cos_ac - b_in*sq_ad*cos_bc);
	if ( abs(comp1) < 0.0001*abs(comp2) )		// do it approximately here too
	{
		const double term1 =  ( a_in * c_in * cos(a_in * c_in) - sin(a_in * c_in) ) / sq_ad;
		const double term2 = -( b_in * c_in * cos(b_in * c_in) - sin(b_in * c_in) ) / sq_bd;
		const double approximate_result = ( term1 + term2 ) / (c_in*c_in*c_in);
		const double error_estimate = 100.0*(approximate_result-result)/result;
		if ( abs( error_estimate ) > 1.0 )
		{
			cout << "Check comps: " << comp1 << "   " << comp2 << endl;
			cout << "Test integration: " << setprecision(12) << setw(16) << a_in << "   " << b_in << "   " << c_in << "   " << d_in;
			cout << "  COMPARE (approx 3): " << result << "   " << approximate_result  << " (error of " << error_estimate << "%)" << endl;
		}
	}
}
*/
//if (1) exit(8);

	return ( overallSign*result );
}


double BoseEinstein::compute_integral_without_phasespace(double a_in, double b_in, double c_in)
{
	// Computes integral given by Integrate[j0[c Q], {Q, a, b}]
	double a = a_in, b = b_in, c = c_in;
	double result = 0.0;
	const double period = 2.0*M_PI/c;

	// While limits contain more than one full period,
	// do one period at a time
	while ( b - a > period )
	{
		double a0 = a, a1 = a0 + period;
		double hw = 0.5*period, cen = 0.5*(a0+a1);
		for (int i = 0; i < npts; ++i)
		{
			double Qloc = cen + hw * x_pts[i];
			result += hw * x_wts[i] * sphericalbesselj0(c*Qloc);
		}
		a += period;
	}
	
	// Finally, do remaining (fraction of a) period
	{
		double a0 = a, a1 = b;
		double hw = 0.5*(a1-a0), cen = 0.5*(a0+a1);
		for (int i = 0; i < npts; ++i)
		{
			double Qloc = cen + hw * x_pts[i];
			result += hw * x_wts[i] * sphericalbesselj0(c*Qloc);
		}
	}

	return ( result );
}


//---------------------------------------------------------------------------
// Compute the derivative of the linearly interpolated pair CDF vs. Q.
void BoseEinstein::set_pair_density(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector<double> & denBar
			)
{
	// Reset.
	denBar.clear();

	if ( linear_interpolate_CDF )
	{
		// Take the derivative of the linear interpolant.
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
			denBar.push_back( 1.0 / ( sortedPairs.at(iPair+1).first - sortedPairs.at(iPair).first ) );
	}
	else
	{
		// Do not estimate pair density.
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
			denBar.push_back( 1.0 );
	}
	
	// Density must go to zero eventually.
	//denBar.push_back( 0.0 );
	denBar.back() = 0.0;

	return;
}


//-------------------------------------------------
// Set lefthand side of shift relation.
void BoseEinstein::set_sorted_xDiffs(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector<double> & sorted_xDiffs
			)
{
	int xDiff_pairCount = 0;
	for (const auto & eachPair : sortedPairs)
	{
		    // Skip unphysical dummy pairs.
		    if (   &eachPair == &sortedPairs.front()
		            or &eachPair == &sortedPairs.back() )
		            continue;

		    const int i1 = eachPair.second.first;
		    const int i2 = eachPair.second.second;
		    Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
		    xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );

		sorted_xDiffs[xDiff_pairCount++] = xDiffPRF.pAbs();
	}
	sort( sorted_xDiffs.begin(), sorted_xDiffs.end() );
	return;
}


//-------------------------------------------------
// Set lefthand side of shift relation.
void BoseEinstein::set_LHS(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector< pair< double, double > > & LHS,
			vector<double> & denBar, int iTab
			)
{
	int pairCount = 0;
	double result = 0.0;
	if ( include_phase_space )
		for (const auto & thisPair : sortedPairs)
		{
			auto nextPair = ( pairCount == (int)sortedPairs.size()-1 ) ? thisPair : *(&thisPair+1);
			const double thisQ = thisPair.first;
			const double nextQ = nextPair.first;

			LHS.push_back( std::make_pair( thisQ, result ) );

			if ( pairCount == (int)denBar.size() )
				continue;

			const double this_sqrt_Q2_plus_4m2 = sqrt(thisQ*thisQ + m2Pair[iTab]);
			const double next_sqrt_Q2_plus_4m2 = sqrt(nextQ*nextQ + m2Pair[iTab]);
			const double this_antideriv = 0.5*thisQ*this_sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( thisQ + this_sqrt_Q2_plus_4m2 );
			const double next_antideriv = 0.5*nextQ*next_sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( nextQ + next_sqrt_Q2_plus_4m2 );

			result += denBar.at(pairCount) * ( next_antideriv - this_antideriv );

			pairCount++;

		}
	else
		for (const auto & thisPair : sortedPairs)
		{
			auto nextPair = ( pairCount == (int)sortedPairs.size()-1 ) ? thisPair : *(&thisPair+1);
			const double thisQ = thisPair.first;
			const double nextQ = nextPair.first;

			LHS.push_back( std::make_pair( thisPair.first, result ) );

			if ( pairCount == (int)denBar.size() )
				continue;

			result += denBar.at(pairCount) * ( nextQ - thisQ );

			pairCount++;
		}


// estimate scale for Q sampling
double size_estimate_sum = 0.0, size_estimate_sqsum = 0.0, size_estimate_cubsum = 0.0;
double mean = 0.0, rms = 0.0, avg_log = 0.0;

const double npairs = sortedPairs.size()-2;

vector<double> effSource(100001);
const double dQ = 0.001;
//for (const auto & eachPair : sortedPairs)
for (int iPair = 1; iPair < (int)sortedPairs.size()-1; iPair++)
{
	//cerr << "iPair = " << iPair << " of " << (int)sortedPairs.size() << "\n";
	auto eachPair = sortedPairs[iPair];
	// Skip unphysical dummy pairs.
	//if (   &eachPair == &sortedPairs.front()
	//	or &eachPair == &sortedPairs.back() )
	//	continue;

	const int i1 = eachPair.second.first;
	const int i2 = eachPair.second.second;
	Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
	xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );

	const double xDiffval = xDiffPRF.pAbs();

	size_estimate_sum += 1.0 / xDiffval;
	size_estimate_sqsum += 1.0 / (xDiffval*xDiffval);
	size_estimate_cubsum += 1.0 / (xDiffval*xDiffval*xDiffval);
	mean += xDiffval;
	rms += xDiffval*xDiffval;
	avg_log += log(1.0 / xDiffval);

	///*
	double currentQ = 0.0;
	for (int iQ = 0; iQ < (int)effSource.size(); iQ++)
	{
		effSource[iQ] += sphericalbesselj0(currentQ * xDiffval);
		currentQ += dQ;
	}
	//*/
}

//cout << setprecision(12) << "size estimate (inverse linear) = " << npairs / size_estimate_sum << endl;
//cout << setprecision(12) << "size estimate (inverse quadratic) = " << sqrt(npairs / size_estimate_sqsum) << endl;
//cout << setprecision(12) << "size estimate (inverse cubic) = " << sqrt(npairs / size_estimate_cubsum) << endl;
//cout << setprecision(12) << "size estimate (mean) = " << mean / npairs << endl;
//cout << setprecision(12) << "size estimate (rms) = " << sqrt(rms / npairs) << endl;
cout << "npairs = " << (int)npairs << endl;
cout << setprecision(12) << "size estimate (linear halfwidth) = " << 2.0 * size_estimate_sqsum / (M_PI * size_estimate_sum) << endl;
//cout << setprecision(12) << "size estimate (gaussian approximation) = " << sqrt(3.0 * size_estimate_cubsum / size_estimate_sum) << endl;
cout << setprecision(12) << "size estimate (gaussian approximation / 2) = " << 0.5 * sqrt(3.0 * size_estimate_cubsum / size_estimate_sum) << endl;
cout << setprecision(12) << "size estimate (geometric mean) = " << exp(avg_log / npairs) << endl;

//if (1) exit (8);
cout << "<<<============================================>>>" << endl;
//*/
double currentQ = 0.0;
for (int iQ = 0; iQ < (int)effSource.size(); iQ++)
{
	cout << "effSource: " << setprecision(12) << currentQ << "   " << effSource[iQ] / (double)(sortedPairs.size()-2) << endl;
	currentQ += dQ;
}
//if (1) return;
if (1) exit (8);
	return;
}


//-------------------------------------------------
// Set righthand side of shift relation.
void BoseEinstein::set_RHS(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector<double> & sorted_xDiffs,
			vector< pair< double, double > > & LHS,
			vector< pair< double, double > > & RHS,
			vector<double> & denBar, int iTab
			)
{
//cout << "RHS:" << endl;
	int pairCount = 0;
	double result = 0.0;
//cout << "Sizes: " << sortedPairs.size() << "   " << LHS.size() << "   " << RHS.size() << "   "
//		<< denBar.size() << "   " << sorted_xDiffs.size() << endl;
	if ( include_phase_space )
		for (const auto & thisPair : sortedPairs)
		{
			auto nextPair = ( pairCount == (int)sortedPairs.size()-1 ) ? thisPair : *(&thisPair+1);
			const double thisQ = thisPair.first;
			const double nextQ = nextPair.first;

//cout << "Check Qs: " << setprecision(12) << setw(16) << thisQ << "   " << nextQ << endl;

			RHS.push_back( std::make_pair( thisQ, result ) );
//cout << setprecision(12) << thisQ << "   " << result << endl;

			if ( pairCount == (int)sorted_xDiffs.size() )
				continue;

//cout << "pairCount = " << pairCount << " of " << sortedPairs.size() - 2 << endl;

//if (pairCount > 10000) exit(8);

			const double one_by_N = 1.0 / static_cast<double>(sortedPairs.size() - 2);
			int eachPairIndex = 0;
			for (const auto & eachPair : sortedPairs)
			{
//if (not pairCount == eachPairIndex) continue;
				// Skip unphysical dummy pairs.
				if (   &eachPair == &sortedPairs.front()
					or &eachPair == &sortedPairs.back() )
					continue;
				
				//const int i1 = eachPair.second.first;
				//const int i2 = eachPair.second.second;
				//Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
				//xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
				//const double xDiffVal = xDiffPRF.pAbs();
				const double xDiffPRFVal = sorted_xDiffs.at(eachPairIndex++);

//cout << "Check Qs(again): " << setprecision(12) << setw(16) << thisQ << "   " << nextQ << endl;
				// Add in physical results
//cout << "need to access " << pairCount << "th element of " << denBar.size() << "-length denBar at line=" << __LINE__ << "...";
//cout << "Line = " << __LINE__ << endl;
				result += one_by_N * denBar.at(pairCount)
							* compute_integral_with_phasespace(
								thisQ, nextQ, xDiffPRFVal, m2Pair[iTab]);
/*cout << "check result here: " << pairCount << "   " << one_by_N << "   " << denBar.at(pairCount) << "   " << thisQ << "   " << nextQ << "   " << xDiffPRFVal << "   " << m2Pair[iTab] << "   " << result << "   " << compute_integral_with_phasespace(
								thisQ, nextQ, xDiffPRFVal, m2Pair[iTab]) << endl;*/
//if (1) exit(8);
//cout << "Line = " << __LINE__ << endl;
//cout << "success!" << endl;
			}

			// whole vector inside at once
			//result += one_by_N * denBar.at(pairCount)
            //                      * compute_integral_with_phasespace( thisQ, nextQ, sorted_xDiffs, m2Pair[iTab]);
//cout << "Line = " << __LINE__ << endl;
			// DOUBLE CHECK INDICES HERE!!!
/*cout << "pairCount = " << pairCount << endl;
cout << "Line = " << __LINE__ << endl;
cout << "denBar.at(pairCount) = " << denBar.at(pairCount) << endl;
cout << "Line = " << __LINE__ << endl;
cout << "sorted_xDiffs.at(pairCount) = " << sorted_xDiffs.at(pairCount) << endl;
cout << "Line = " << __LINE__ << endl;*/
/*			result += one_by_N * denBar.at(pairCount)
						* compute_integral_with_phasespace(
							thisQ, nextQ, sorted_xDiffs.at(pairCount), m2Pair[iTab]);
cout << "check result here: "
		<< pairCount << "   " << one_by_N << "   " << denBar.at(pairCount) << "   "
		<< thisQ << "   " << nextQ << "   " << sorted_xDiffs.at(pairCount) << "   " << m2Pair[iTab] << "   "
		<< result << "   " << compute_integral_with_phasespace(thisQ, nextQ, sorted_xDiffs.at(pairCount), m2Pair[iTab]) << endl;*/

//cout << "Line = " << __LINE__ << endl;
			pairCount++;
		}
	else
		for (const auto & thisPair : sortedPairs)
		{
			auto nextPair = ( pairCount == (int)sortedPairs.size()-1 ) ? thisPair : *(&thisPair+1);
			const double thisQ = thisPair.first;
			const double nextQ = nextPair.first;

			RHS.push_back( std::make_pair( thisQ, result ) );

			if ( pairCount == (int)denBar.size() )
				continue;

			const double one_by_N = 1.0 / static_cast<double>(sortedPairs.size() - 2);
			for (const auto & eachPair : sortedPairs)
			{
				// Skip unphysical dummy pairs.
				if (   &eachPair == &sortedPairs.front()
					or &eachPair == &sortedPairs.back() )
					continue;

				const int i1 = eachPair.second.first;
				const int i2 = eachPair.second.second;
				Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
				xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
				
				// Add in physical results
				result += one_by_N * denBar.at(pairCount)
							* compute_integral_without_phasespace(
								thisQ, nextQ, xDiffPRF.pAbs());
			}

			pairCount++;
		}

	// Add in constant piece to RHS result.
	for (int iPair = 0; iPair < (int)LHS.size(); ++iPair)
		RHS.at(iPair).second += LHS.at(iPair).second;
//cout << "Line = " << __LINE__ << endl;
//cout << "Sizes: " << sortedPairs.size() << "   " << LHS.size() << "   " << RHS.size() << "   " << sorted_xDiffs.size() << endl;
//for (int iPair = 0; iPair < (int)LHS.size(); ++iPair)
//	cout << "CHECK: "  << setprecision(12) << sortedPairs.at(iPair).first << "   " << LHS.at(iPair).second << "   " << RHS.at(iPair).second << endl;
//if (1) exit(8);

	return;
}


//-------------------------------------------------
// Compute the unshifted pair integrals at each Qi.
void BoseEinstein::evaluate_shift_relation_at_Qi(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector< pair< double, double > > & LHS,
			vector< pair< double, double > > & RHS,
			vector<double> & denBar, int iTab
			)
{
	// Do some timing operations.
	auto start = std::chrono::system_clock::now();

	//-------
	// Reset.
	LHS.clear();
	RHS.clear();

//cout << "This is just a test: "
//		<< compute_integral_with_phasespace( 1.0, 2.0, 3.0, 4.0 )
//		<< endl;

/*cout << "First check of sorted pairs:" << endl;
for (const auto & thisPair : sortedPairs)
	cout << thisPair.first << "   "
		<< thisPair.second.first << "   "
		<< thisPair.second.second << endl;
cout << endl;*/

	vector<double> sorted_xDiffs(sortedPairs.size()-2);
	set_sorted_xDiffs( sortedPairs, sorted_xDiffs );

	//------------------
	// Set LHS integral.
	set_LHS(sortedPairs, LHS, denBar, iTab);

	//------------------
	// Set RHS integral.
	set_RHS(sortedPairs, sorted_xDiffs, LHS, RHS, denBar, iTab);

//cout << "Done setting RHS" << endl;
//if (1) exit(8);

   auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "BoseEinsteinCheck: spent " << elapsed_seconds.count() << "s in " << __FUNCTION__ << endl;

	return;
}


//---------------------------------------------------------------------------
// Calculate shift and (unnormalized) compensation for pair using space-time
// interval and and mode 1 of evaluating BE enhancement.
void BoseEinstein::shiftPairs_mode1(
					vector< pair< double, pair <int,int> > > & sortedPairs,
					vector<double> & pairShifts,
					vector<double> & pairCompensationShifts, int iTab)
{
//	auto start = std::chrono::system_clock::now();

	//--------------------------------
	// Construct and set pair density.
	vector<double> denBar;
	set_pair_density(sortedPairs, denBar);


	//--------------------------------------------------
	// Set LHS and RHS of shift relation at each pair Q.
	vector< pair< double, double > > LHS, RHS;
	evaluate_shift_relation_at_Qi( sortedPairs, LHS, RHS, denBar, iTab );
/*
cout << "<<<============================================>>>" << endl;
cout << "Check sizes: " << setprecision(12) << setw(16) << sortedPairs.size() << "   " << LHS.size() << "   " << RHS.size() << "   " << denBar.size() << endl;
cout << "Check sortedPairs: " << endl;
for (const auto & iPair : sortedPairs) cout << iPair.first << "   " << iPair.second.first << "   " << iPair.second.second << endl;
//cout << "Check LHS: " << endl;
//for (const auto & iLHS : LHS) cout << iLHS.first << "   " << iLHS.second << endl;
//cout << "Check RHS: " << endl;
//for (const auto & iRHS : RHS) cout << iRHS.first << "   " << iRHS.second << endl;
cout << "Check denBar: " << endl;
for (const auto & iPair : denBar) cout << iPair << endl;
cout << "Check xDiffs: " << endl;
vector<double> effSource(1001);
for (const auto & eachPair : sortedPairs)
{
	// Skip unphysical dummy pairs.
	if (   &eachPair == &sortedPairs.front()
		or &eachPair == &sortedPairs.back() )
		continue;

	const int i1 = eachPair.second.first;
	const int i2 = eachPair.second.second;
	Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
	xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );

	//cout << "xDiffval: " << xDiffPRF.pAbs() << endl;
	double currentQ = 0.0;
	for (int iQ = 0; iQ < (int)effSource.size(); iQ++)
	{
		effSource[iQ] += sphericalbesselj0(currentQ * xDiffPRF.pAbs());
		currentQ += 0.01;
	}
}
cout << "<<<============================================>>>" << endl;
*/
/*double currentQ = 0.0;
for (int iQ = 0; iQ < (int)effSource.size(); iQ++)
{
	cout << "effSource: " << currentQ << "   " << effSource[iQ] / (double)(sortedPairs.size()-2) << endl;
	currentQ += 0.01;
}
//if (1) return;
if (1) exit (8);*/

	// keep track of which pairs are shifted
	int n_skipped_pairs = 0;
	vector<bool> this_pair_shifted(sortedPairs.size()-2, false);

	// -------------------------------------------
	// Finally, get shifts for each physical pair.
	for (int iPair = 1; iPair < (int)sortedPairs.size()-1; ++iPair)
	{
		int thisPair = iPair;
		const double Q0          = sortedPairs.at(iPair).first;				// original Q
		const int iPair1         = sortedPairs.at(iPair).second.first;		// first hadron
		const int iPair2         = sortedPairs.at(iPair).second.second;	// second hadron
		const double thisPairLHS = LHS.at(thisPair).second;		// > 0
		double RHS_lower         = RHS.at(thisPair).second;
		double RHS_upper         = RHS.at(thisPair+1).second;

//cout << "Test: " << Q0 << "   " << iPair << "   " << iPair1 << "   " << iPair2 << "   " << RHS_lower << "   " << thisPairLHS << "   " << RHS_upper << endl;

		//--------------------------------------------
		// Find RHS interval which contains LHS value.
		while ( RHS_lower > thisPairLHS
			and thisPair  >= 0 )
		{
			// shift to the left
			RHS_lower = RHS.at(--thisPair).second;
			RHS_upper = RHS.at(thisPair+1).second;
		}
		while ( RHS_upper < thisPairLHS
			and thisPair  < (int)RHS.size()-2 )
		{
			// shift to the right
			RHS_lower = RHS.at(++thisPair).second;
			RHS_upper = RHS.at(thisPair+1).second;
		}

//cout << "Test again: " << Q0 << "   " << iPair << "   " << iPair1 << "   " << iPair2 << "   " << RHS_lower << "   " << thisPairLHS << "   " << RHS_upper << endl;
//if (1) exit(8);

		//-------------------------------------------------
		// If search ended up out of range, skip this pair.
		if ( thisPair < 0
				or thisPair  > (int)sortedPairs.size() - 2
				or RHS_lower > thisPairLHS
				or RHS_upper < thisPairLHS )
		{
			n_skipped_pairs++;
			pairShifts.push_back( 0.0 );
			pairCompensationShifts.push_back( 0.0 );
			continue;
		}

//cout << "Test thisPair = " << thisPair << "   " << sortedPairs.size() << endl;

		const double leftQ        = sortedPairs.at(thisPair).first;
		const double rightQ       = sortedPairs.at(thisPair+1).first;
		const double thisdenBar   = ( linear_interpolate_CDF ) ? 1.0 / ( rightQ - leftQ ) : 1.0;
		const double one_by_N     = 1.0 / static_cast<double>(sortedPairs.size() - 2);

		//------------------------------------------
		// Estimate pair shift using Newton-Raphson.
		constexpr double ACCURACY = 1.e-6;
		constexpr int MAXTRIES    = 10;
		const double c0           = thisPairLHS;	// The target value
		double Qlower             = leftQ;
		double Qupper             = rightQ;

		//--------------------------------
		// Initialize to closest endpoint.
		double Qnew = Qlower;
		double RHSnew = RHS_lower;
		if ( abs(c0-RHS_lower) > abs(c0-RHS_upper) )
		{
			Qnew   = Qupper;
			RHSnew = RHS_upper;
		}

		//-----------------
		// Start iteration.
		int ntries = 0;
		if ( include_phase_space )
			while ( abs( RHSnew - c0 ) > ACCURACY and ntries < MAXTRIES )
			{
				//cout << "with phase space: ntries = " << ntries << endl;

				const double sqrt_Q2_plus_4m2 = sqrt(Qnew*Qnew + m2Pair[iTab]);
				const double PSfactor         = Qnew*Qnew / sqrt_Q2_plus_4m2;
				double deriv                  = thisdenBar * PSfactor;

				if (ntries > 0)
				{
					//-----------------------------
					// Set RHSnew to left endpoint.
					RHSnew = RHS_lower;

					//---------------------------------------
					// Add in the constant phase space piece.
					const double sqrt_Qlower2_plus_4m2 = sqrt(Qlower*Qlower + m2Pair[iTab]);
					RHSnew += thisdenBar * (
								0.5*Qnew*sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( Qnew + sqrt_Q2_plus_4m2 )
								-0.5*Qlower*sqrt_Qlower2_plus_4m2 + 0.5*m2Pair[iTab]*log( Qlower + sqrt_Qlower2_plus_4m2 )
									);

					//---------------------------------
					// Add in the BE enhancement piece.
					/*for (const auto & eachPair : sortedPairs)
					{
						// Skip unphysical dummy pairs.
						if (   &eachPair == &sortedPairs.front()
							or &eachPair == &sortedPairs.back() )
							continue;

						const int i1  = eachPair.second.first;
						const int i2  = eachPair.second.second;
						Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
						xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
		
						// Add in physical results
						const double xDiffMag = xDiffPRF.pAbs();
						RHSnew += one_by_N * thisdenBar
									* compute_integral_with_phasespace( Qlower, Qnew, xDiffMag, m2Pair[iTab] );
						deriv  += one_by_N * thisdenBar
									* PSfactor * sphericalbesselj0( Qnew * xDiffMag );
					}*/
// DOUBLE CHECK INDICES HERE!!!
{
	const int i1  = sortedPairs.at(iPair).second.first;
	const int i2  = sortedPairs.at(iPair).second.second;
	Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
	xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
	const double xDiffMag = xDiffPRF.pAbs();
	RHSnew += one_by_N * thisdenBar
				* compute_integral_with_phasespace(
					Qlower, Qnew, xDiffMag, m2Pair[iTab]);
	deriv  += one_by_N * thisdenBar
				* PSfactor * sphericalbesselj0( Qnew * xDiffMag );
}

				}

				//if ( ntries>0 and ( RHSnew > RHS_upper or RHSnew < RHS_lower ) )
				//	cout << "WARNING: RHSnew = " << RHSnew << " not in range (" << RHS_lower << ", " << RHS_upper << ")!" << endl;
				if ( Qnew > Qupper or Qnew < Qlower )
					cout << "WARNING (ntries=" << ntries << "): Qnew = " << Qnew << " not in range (" << Qlower << ", " << Qupper << ")!" << endl;
				//cout << "NR method (ntries=" << ntries << "): " << Qnew << "   " << RHSnew << "   " << deriv << "   " << c0 << endl;	

				//----------------------------
				// Get next estimate for Qnew.
				Qnew -= (RHSnew - c0) / deriv;

				ntries++;
			}
		else
			while ( abs( RHSnew - c0 ) > ACCURACY and ntries < MAXTRIES )
			{
				//cout << "without phase space: ntries = " << ntries << endl;

				double deriv = thisdenBar;

				if (ntries > 0)
				{
					//-----------------------------
					// Set RHSnew to left endpoint.
					RHSnew = RHS_lower;

					//---------------------------------------
					// Add in the constant phase space piece.
					RHSnew += thisdenBar * ( Qnew - Qlower );

					//---------------------------------
					// Add in the BE enhancement piece.
					for (const auto & eachPair : sortedPairs)
					{
						// Skip unphysical dummy pairs.
						if (   &eachPair == &sortedPairs.front()
							or &eachPair == &sortedPairs.back() )
							continue;

						const int i1 = eachPair.second.first;
						const int i2 = eachPair.second.second;
						Vec4 xDiffPRF = ( hadronBE[i1].x - hadronBE[i2].x ) * MM2FM / HBARC;
						xDiffPRF.bstback( 0.5*(hadronBE[i1].p + hadronBE[i2].p) );
			
						// Add in physical results
						const double xDiff = xDiffPRF.pAbs();
						RHSnew += one_by_N * thisdenBar
									* compute_integral_without_phasespace( Qlower, Qnew, xDiff );
						deriv += one_by_N * thisdenBar * sphericalbesselj0( Qnew * xDiff );

						//pairCount++;
					}
				}

				//if ( ntries>0 and ( RHSnew > RHS_upper or RHSnew < RHS_lower ) )
				//	cout << "WARNING: RHSnew = " << RHSnew << " not in range (" << RHS_lower << ", " << RHS_upper << ")!" << endl;
				if ( Qnew > Qupper or Qnew < Qlower )
					cout << "WARNING (ntries=" << ntries << "): Qnew = " << Qnew << " not in range (" << Qlower << ", " << Qupper << ")!" << endl;
				//cout << "NR method (ntries=" << ntries << "): " << Qnew << "   " << RHSnew << "   " << deriv << "   " << c0 << endl;	

				//----------------------------
				// Get next estimate for Qnew.
				Qnew -= (RHSnew - c0) / deriv;

				ntries++;
			}


		//-------------------------------------------
		// Store results based on size of BE effects.
		Vec4 p1 = hadronBE[iPair1].p;
		Vec4 p2 = hadronBE[iPair2].p;
		Vec4 x1 = hadronBE[iPair1].x;
		Vec4 x2 = hadronBE[iPair2].x;
		Vec4 xDiffPRF = ( x1 - x2 ) * MM2FM / HBARC;
		xDiffPRF.bstback( 0.5*(p1 + p2) );

		/*Vec4 qLocPRF = p1 - p2;
		qLocPRF.bstback( 0.5*(p1 + p2) );
		cout << setprecision(12) << setw(16) << " --> Checks Q: " << Q0 << "   " << sqrt(abs(( p1 - p2 ) * ( p1 - p2 ))) << endl;
		cout << " --> Check qLocPRF: " << qLocPRF;
		cout << " --> Check xDiffPRF: " << xDiffPRF;
		cout << " --> Check arg: " << arg << "   " << qLocPRF*xDiffPRF << endl;*/

		// Do this for shifts and compensation shifts together
		//const double arg = ( p1 - p2 ) * ( x1 - x2 ) * MM2FM / HBARC;
		/*if ( abs(arg) > 16.0 * M_PI )
		{
			pairShifts.push_back( 0.0 );
			pairCompensationShifts.push_back( Qnew - Q0 );
		}
		else*/
		//{
			pairShifts.push_back( Qnew - Q0 );
			pairCompensationShifts.push_back( 0.0 );
		//}
		const double Qnew_LININT = Qlower + (Qupper-Qlower)*(thisPairLHS-RHS_lower)/(RHS_upper-RHS_lower);
		cout << "BoseEinsteinCheck: CHECK shift is "
				<< Qnew - Q0 << " vs. " << Qnew_LININT - Q0 << " of " << Q0 << endl;

//cout << "Line = " << __LINE__ << endl;
		// Shift pairs closer together; use pairs
		// shifted further away to compensate.
		if ( not include_posDelQ_in_compensation or Qnew <= Q0 )
			this_pair_shifted.at(iPair-1) = true;
//cout << "Line = " << __LINE__ << endl;

	}

//	cout << "BoseEinsteinCheck: n_skipped_pairs = " << n_skipped_pairs << endl;


//if (1) exit(8);



//cout << "MAde it ere..." << endl;



	// Finally, discard fake pairs at beginning and end of sortedPairs
	sortedPairs.erase ( sortedPairs.begin() );
	sortedPairs.erase ( sortedPairs.end() );

	int pairIndex = 0;
	for (const auto & iPair : sortedPairs)
	{

		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		if ( this_pair_shifted.at(pairIndex) )
		{
			const double Qold  = iPair.first;
			double Qnew  = Qold + pairShifts.at(pairIndex);
			const double Q2old = Qold*Qold;
			double Q2new = Qnew*Qnew;

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
		}
		else
		{
			// Now get compensation shifts
			///*
			//cout << "Line = " << __LINE__ << endl;
			const double Qold  = iPair.first;
			double Qnew  = Qold + pairShifts.at(pairIndex);
			const double Q2old = Qold*Qold;
			double Q2new = Qnew*Qnew;
			//cout << "Line = " << __LINE__ << endl;

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
			Vec4 pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
			//*/
			//Vec4 pDiff = hadronBE[i1].p - hadronBE[i2].p;
			hadronBE[i1].pComp += pDiff;
			hadronBE[i2].pComp -= pDiff;
		}


		pairIndex++;

	}

//cout << "Finished here" << endl;
//	auto end = std::chrono::system_clock::now();

//    std::chrono::duration<double> elapsed_seconds = end-start;
//	std::cout << "BoseEinsteinCheck: spent " << elapsed_seconds.count() << "s in " << __FUNCTION__ << endl;


	return;
}
//*/


//==========================================================================

} // end namespace Pythia8
