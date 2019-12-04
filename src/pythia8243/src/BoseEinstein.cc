// BoseEinstein.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BoseEinsten class.

#include "Pythia8/BoseEinstein.h"

#include <ostream>

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
  //useInv   = settings.flag("BoseEinstein:useInvariantSourceSize");
  //useDist  = settings.flag("BoseEinstein:useDistribution");

  // Shape of Bose-Einstein enhancement/suppression.
  lambda   = settings.parm("BoseEinstein:lambda");
  QRef     = settings.parm("BoseEinstein:QRef");
  //dim      = settings.parm("BoseEinstein:sourceDimension");
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

void BoseEinstein::set_QRef(int iSpecies)
{
	if ( useDist )
	{
		switch (dim)
		{
			case 1:
				QRef = get_1D_source_size(iSpecies, useInv);
				break;
			/*case 3:
				QRef = get_3D_source_size(iSpecies);
				break;
			case 4:
				QRef = get_4D_source_size(iSpecies);
				break;*/
			default:
				infoPtr->errorMsg("Warning in BoseEinstein::set_QRef: "
      "choice of QRef makes no sense, defaulting to 1D source size");
				QRef = get_1D_source_size(iSpecies, useInv);
				break;
		}
	}
	// if not, do nothing
}

double BoseEinstein::get_1D_source_size(int iSpecies, bool useInvariant)
{
	double result = 0.0;
	double count = 0.0;

	if (useInvariant)
		for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
		//for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
		{
			//Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
			Vec4 x = hadronBE[i1].x;
			result += x*x;
			count += 1.0;
		}
	else
		for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
		//for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
		{
			//Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
			Vec4 x = hadronBE[i1].x;
			result += x.pAbs2();
			count += 1.0;
		}

	// return RMS source size
	return (sqrt(result / (count+1.e-100)));
}

void BoseEinstein::get_3D_source_size(int iSpecies, double & xSize, double & ySize, double & zSize)
{
	xSize = 0.0, ySize = 0.0, zSize = 0.0;
	double count = 0.0;

	for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
	//for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
	{
		//Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
		Vec4 x = hadronBE[i1].x;
		xSize += x.px()*x.px();
		ySize += x.py()*x.py();
		zSize += x.pz()*x.pz();
		count += 1.0;
	}

	xSize = sqrt(xSize / (count+1.e-100));
	ySize = sqrt(ySize / (count+1.e-100));
	zSize = sqrt(zSize / (count+1.e-100));

	return;
}

void BoseEinstein::get_4D_source_size(int iSpecies, double & xSize, double & ySize, double & zSize, double & tSize)
{
	xSize = 0.0, ySize = 0.0, zSize = 0.0, tSize = 0.0;
	double count = 0.0;

	for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
	//for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
	{
		//Vec4 xDiff = hadronBE[i1].x - hadronBE[i2].x;
		Vec4 x = hadronBE[i1].x;
		xSize += x.px()*x.px();
		ySize += x.py()*x.py();
		zSize += x.pz()*x.pz();
		tSize += x.e()*x.e();
		count += 1.0;
	}

	xSize = sqrt(xSize / (count+1.e-100));
	ySize = sqrt(ySize / (count+1.e-100));
	zSize = sqrt(zSize / (count+1.e-100));
	tSize = sqrt(tSize / (count+1.e-100));

	return;
}

//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

bool BoseEinstein::shiftEvent( Event& event) {
  // Reset list of identical particles.
  hadronBE.resize(0);

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
		case 1:	// use a new cos(q*x) enhancement based on space-time interval between production points (mode 1)
			{
				bool enoughPairsToProceed = getSortedPairs( sortedPairs, iSpecies );
				if ( enoughPairsToProceed )
					shiftPairs_mode1( sortedPairs, pairShifts, pairCompensationShifts, iTab );
			}
			break;
		/*case 2:	// use a new cos(q*x) enhancement based on space-time interval between production points (mode 2)
			{
				bool enoughPairsToProceed = getSortedPairs( sortedPairs, iSpecies );
				if ( enoughPairsToProceed )
					shiftPairs_mode2( sortedPairs, pairShifts, pairCompensationShifts );
			}
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
cout 	<< setprecision(8)
		<< "Original p and pShift: " << i << endl << "p=" << hadronBE[i].p << "pShift=" << hadronBE[i].pShift;
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
      hadronBE[i].pShift += compFac * hadronBE[i].pComp;
cout 	<< setprecision(8)
		<< "Net shift at this point: " << i << "   " << hadronBE[i].pShift;
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
cout << setprecision(16) << "This event did not pass! Check: " << abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << endl;
    return true;
  }
else cout << setprecision(16) << "This event passes! Check: " << abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << endl;


  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE[i].iPos, 99);
    event[ iNew ].p( hadronBE[i].p );
  }

  // Done.
  return true;

}


//--------------------------------
// Construct list of sorted pairs.

bool BoseEinstein::getSortedPairs( vector< pair< double, pair <int,int> > > & sortedPairs, int iSpecies )
{
    int iTab  = ITABLE[ iSpecies ];
	sortedPairs.clear();	//reset

	// get all values in vector first
	for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
	for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
		sortedPairs.push_back(
			std::make_pair(
				m2(hadronBE[i1].p, hadronBE[i2].p) - m2Pair[iTab],
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

	// add fake last "pair" (new QVal is 10% larger than last one, just for definiteness)
	sortedPairs.push_back( std::make_pair( 1.1*sortedPairs.back().first, std::make_pair(-1, -1) ) );

cout << "<<<==========================================================>>>" << endl;
cout << "CHECK sortedPairs (size = " << sortedPairs.size() << "): " << endl;
	for (const auto & iPair : sortedPairs)
		cout << iPair.first << "   " << iPair.second.first << "   " << iPair.second.second << endl;
cout << "<<<==========================================================>>>" << endl;

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

//---------------------------------------------------------------------------
// Calculate shift and (unnormalized) compensation for pair using space-time
// interval and and mode 1 of evaluating BE enhancement.
void BoseEinstein::shiftPairs_mode1(  vector< pair< double, pair <int,int> > > & sortedPairs,
								vector<double> & pairShifts,
								vector<double> & pairCompensationShifts, int iTab)
{

		// Construct PDF from CDF for CF numerator and denominator
		// as functions of Q, using these sortedPairs
		// Set CDFs first
		vector< pair< double, double > > den_cdf, num_cdf;
		// Do phase space integral analytically
		vector<double> antiderivative;

		// Initialize at Q=0
		double running_den_sum = 0.0;
		double running_num_sum = 0.0;

		// Loop over all pairs (except final fake "pair")
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			const double QVal = sortedPairs[iPair].first;
			const double sqrt_Q2_plus_4m2 = sqrt(QVal*QVal + m2Pair[iTab]);
			antiderivative.push_back( 0.5*QVal*sqrt_Q2_plus_4m2 - 0.5*m2Pair[iTab]*log( QVal + sqrt_Q2_plus_4m2 ) );

			den_cdf.push_back( std::make_pair(QVal, running_den_sum) );
			num_cdf.push_back( std::make_pair(QVal, running_num_sum) );

cout << "CHECK: " << running_den_sum << "   " << running_num_sum << endl;

			// Note coordinate separation in GeV^{-1}.
			running_den_sum += 1.0;
			if ( iPair == 0 )
				running_num_sum += 2.0;
			else
			{
				const int iPair1 = sortedPairs[iPair].second.first;
				const int iPair2 = sortedPairs[iPair].second.second;

				Vec4 xDiff = ( hadronBE[iPair1].x - hadronBE[iPair2].x ) * MM2FM / HBARC;
				const double arg = ( hadronBE[iPair1].p - hadronBE[iPair2].p ) * xDiff;
				running_num_sum += 1.0 + sphericalbesselj0( arg );
			}
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
			<< running_num_integral[iPair] << endl;

		// Finally, get shifts for each pair (except for first and last ones)
		// First "pair" is just the Q=0 point
		// Last "pair" is just a repeat to ensure PDF = 0 for Q > Qmax
		for (int iPair = 1; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			int thisPair = iPair;
			const double thisPairDen = running_den_integral.at(thisPair);		// > 0
			double leftPairNum = running_num_integral.at(thisPair);
			double rightPairNum = running_num_integral.at(thisPair+1);

			if ( thisPairDen > running_num_integral.back() )
			{
				cout << "WARNING: no solution for this pair (iPair = " << iPair << ")!  Skipping..." << endl;
				break;
			}

//cout << "Do shifting if necessary: " << thisPairDen << "   " << leftPairNum << "   " << rightPairNum << endl;
			// shift left and right numerator entries until they contain
			// denominator value; then invert using linear interpolation
			while ( leftPairNum > thisPairDen )					// most likely
			{
				// shift to the left
				leftPairNum = running_num_integral.at(--thisPair);
				rightPairNum = running_num_integral.at(thisPair+1);
			}
			while ( rightPairNum < thisPairDen )				// also possible
			{
				// shift to the right
				leftPairNum = running_num_integral.at(++thisPair);
				rightPairNum = running_num_integral.at(thisPair+1);
			}
//cout << "Finished shifting" << endl;

			// check if needed pair interval is included (skip shift of this pair if not)
			if ( thisPair < 0 /*or thisPair >= (int)running_den_integral.size()-1*/ )
			{
				cout << "WARNING: thisPair = " << thisPair << " was out of range!" << endl;
				break;
			}

			// should now have leftPairNum < thisPairDen < rightPairNum
			// --> get shift Q from linear interpolation (Qnew === Q + deltaQ)
			const double Q0 = sortedPairs[iPair].first;	// original Q
			const int iPair1 = sortedPairs[iPair].second.first;		// first hadron
			const int iPair2 = sortedPairs[iPair].second.second;	// second hadron
/*cout << "Q0=" << Q0 << "; (thisPair, iPair) = (" << thisPair << ", " << iPair << ") of "
		<< "(" << sortedPairs.size() << ", " << num_cdf.size() << ", "
		<< den_cdf.size() << ", " << antiderivative.size() << ")" << endl;//std::flush;
cout << "\t" << thisPairDen << "   " << leftPairNum << "   " << rightPairNum << endl;*/
			const double leftQ = sortedPairs.at(thisPair).first;
			const double rightQ = sortedPairs.at(thisPair+1).first;
//cout << "   " << leftQ << "   " << rightQ << endl;
			double Qnew = leftQ + (thisPairDen-leftPairNum)*(rightQ-leftQ)/(rightPairNum-leftPairNum);

			// linear interpolation not good enough; need to solve numerically with Newton-Raphson
			// use Qnew as initial guess
			// need num_pdf in this interval
			const double num_pdf = (num_cdf[thisPair+1].second - num_cdf[thisPair].second) / (rightQ - leftQ);
			const double Numi = num_pdf * antiderivative[thisPair];
			// Solve antiderivative(Qnew) == c0
			const double c0 = (thisPairDen - leftPairNum + Numi) / num_pdf;	// DOUBLE CHECK THIS AND PREVIOUS STEPS

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

			const double arg = ( hadronBE[iPair1].p - hadronBE[iPair2].p )
								* ( hadronBE[iPair1].x - hadronBE[iPair2].x )
								* MM2FM / HBARC;
			// Do this for shifts and compensation shifts together
			/*if ( abs(arg) > 2.0 * M_PI )
			{
				pairShifts.push_back( 0.0 );
				//pairCompensationShifts.push_back( Qnew - Q0 );
				pairCompensationShifts.push_back( 0.0 );
			}
			else if ( abs(arg) > 1.0 * M_PI )
			{
				pairShifts.push_back( 0.0 );
				pairCompensationShifts.push_back( Qnew - Q0 );
			}
			else
			{
				pairShifts.push_back( Qnew - Q0 );
				pairCompensationShifts.push_back( 0.0 );
			}*/
			if ( abs(arg) > 16.0 * M_PI )
			{
				pairShifts.push_back( 0.0 );
				pairCompensationShifts.push_back( Qnew - Q0 );
			}
			else
			{
				pairShifts.push_back( Qnew - Q0 );
				pairCompensationShifts.push_back( 0.0 );
			}

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


	// Finally, discard fake pairs at beginning and end of sortedPairs
	sortedPairs.erase ( sortedPairs.begin() );
	sortedPairs.erase ( sortedPairs.end() );

	int pairIndex = 0;
	for (const auto & iPair : sortedPairs)
	{
		const int i1 	   = iPair.second.first;
		const int i2 	   = iPair.second.second;
		const double Qold  = iPair.first;
		double Qnew  = Qold + pairShifts[pairIndex];
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

		// Now get compensation shifts
		Qnew  = Qold + pairCompensationShifts[pairIndex];
		Q2new = Qnew*Qnew;

		// Calculate corresponding three-momentum shift.
		Q2Diff    = Q2new - Q2old;
		p2DiffAbs = (hadronBE[i1].p - hadronBE[i2].p).pAbs2();
		p2AbsDiff = hadronBE[i1].p.pAbs2() - hadronBE[i2].p.pAbs2();
		eSum      = hadronBE[i1].p.e() + hadronBE[i2].p.e();
		eDiff     = hadronBE[i1].p.e() - hadronBE[i2].p.e();
		sumQ2E    = Q2Diff + eSum * eSum;
		rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
		rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
		factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
						+ Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

		// Add shifts to sum. (Energy component dummy.)
		pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
		hadronBE[i1].pComp += pDiff;
		hadronBE[i2].pComp -= pDiff;

		pairIndex++;

	}

	return;
}

double BoseEinstein::compute_integral(double a_in, double b_in, double c, double d)
{
	// Computes integral given by Integrate[Q^2 ( 0 + j0[c Q] ) / Sqrt[Q^2 + d], {Q, a, b}]
	double result = 0.0;
	const double period = 2.0*M_PI/c;
	double a = a_in, b = b_in;

	while ( b - a > period )
	{
		double a0 = a, a1 = a0 + period;
		double hw = 0.5*period, cen = 0.5*(a0+a1);
		for (int i = 0; i < npts; ++i)
		{
			double Qloc = cen + hw * x_pts_7[i];
			result += hw * x_wts_7[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
		}
		a += period;
	}
	
	// finally, do this
	{
		double a0 = a, a1 = b;
		double hw = 0.5*period, cen = 0.5*(a0+a1);
		for (int i = 0; i < npts; ++i)
		{
			double Qloc = cen + hw * x_pts_7[i];
			result += hw * x_wts_7[i] * Qloc * sin(c*Qloc) / ( c*sqrt(Qloc*Qloc + d) );
		}
	}

	// don't recalculate this after denominator!
	/*
	// from the constant part
	const double sq_a = sqrt(d + a*a);
	const double sq_b = sqrt(d + b*b);
	const double constantTerm = 0.5*( b*sq_b - a*sq_a
									+ d*log( (a + sq_a) / (b + sq_b) ) );

	result += constantTerm;*/

	return ( result );
}


/*
//---------------------------------------------------------------------------
// Calculate shift and (unnormalized) compensation for pair using space-time
// interval and and mode 2 of evaluating BE enhancement.
void BoseEinstein::shiftPairs_mode2(  vector< pair< double, pair <int,int> > > & sortedPairs,
								vector<double> & pairShifts,
								vector<double> & pairCompensationShifts, int iTab)
{
		// ----------------------------------------------
		// Construct PDF from CDF for CF denominator only
		// as function of Q, using sortedPairs
		// Set CDF first
		vector< pair< double, double > > den_cdf;
		// Do phase space integral analytically
		vector<double> antiderivative;

		// -----------------------------------------------
		// Loop over all pairs (except final fake "pairs")
		// Initialize at Q=0
		double running_den_sum = 0.0;
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			const double QVal = sortedPairs[iPair].first;
			den_cdf.push_back( std::make_pair(QVal, running_den_sum) );
			running_den_sum += 1.0;
		}

		// ---------------------------------------------------------------------------
		// Add the final fake "pair" to CDF ensure last real pair is handled correctly
		// (QVal shifted by some positive number, running sums are unchanged)
		den_cdf.push_back( std::make_pair( den_cdf.back().first+1.0, den_cdf.back().second ) );


		// -----------------------------------------------------------
		// Estimate denominator PDF using linear interpolation
		// of denominator CDF and evaluate running integral
		// with extra factor for numerator
		running_den_sum = 0.0;
		running_num_sum = 0.0;
		vector<double> running_den_integral(1, running_den_sum);
		vector<double> running_num_integral(1, running_num_sum);

		// ----------------------------
		// Skip final fake "pair" again
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			const int iPair1 = sortedPairs[iPair].second.first;		// first hadron
			const int iPair2 = sortedPairs[iPair].second.second;	// second hadron
			double leftQVal = sortedPairs[iPair].first;
			double rightQVal = sortedPairs[iPair+1].first;
			// PDF ~ slope of CDF
			const double den_pdf = (den_cdf[iPair+1].second - den_cdf[iPair].second) / (rightQVal - leftQVal);

			Vec4 xDiffPRF = ( hadronBE[iPair1].x - hadronBE[iPair2].x ) * MM2FM / HBARC;
			xDiffPRF.bstback( 0.5*(hadronBE[iPair1].p + hadronBE[iPair2].p) );

			// Update running integrals
			const double increment = den_pdf * ( antiderivative[iPair+1] - antiderivative[iPair] );
			running_den_sum += increment;
			running_num_sum += increment + den_pdf * compute_integral(leftQVal, rightQVal, xDiffPRF.pAbs(), m2Pair[iTab]);

			// Store 
			running_den_integral.push_back( running_den_sum );
			running_num_integral.push_back( running_num_sum );
		}

		// ------------------------------------------
		// Store one more entry for fake final "pair"
		running_den_integral.push_back( running_den_sum );
		running_num_integral.push_back( running_num_sum );

		// ------------------------------------------------------------------
		// Finally, get shifts for each pair (except for first and last ones)
		// First "pair" is just the Q=0 point
		// Last "pair" is just a repeat to ensure PDF = 0 for Q > Qmax
		for (int iPair = 1; iPair < (int)sortedPairs.size()-1; ++iPair)
		{
			int thisPair = iPair;
			const double thisPairDen = running_den_integral[thisPair];		// > 0
			double leftPairNum = running_num_integral[thisPair];
			double rightPairNum = running_num_integral[thisPair+1];

			// shift left and right numerator entries until they contain
			// denominator value; then invert using linear interpolation
			while ( leftPairNum > thisPairDen )					// most likely
			{
				// shift to the left
				leftPairNum = running_num_integral[--thisPair];
				rightPairNum = running_num_integral[thisPair+1];
			}
			while ( rightPairNum < thisPairDen )				// also possible
			{
				// shift to the right
				leftPairNum = running_num_integral[++thisPair];
				rightPairNum = running_num_integral[thisPair+1];
			}

			if ( thisPair < 0 )
				break;

			// should now have leftPairNum < thisPairDen < rightPairNum
			// --> get shift Q from linear interpolation (Qnew === Q + deltaQ)
			const double Q0 = sortedPairs[iPair].first;	// original Q
			const int iPair1 = sortedPairs[iPair].second.first;		// first hadron
			const int iPair2 = sortedPairs[iPair].second.second;	// second hadron
			const double leftQ = sortedPairs[thisPair].first;
			const double rightQ = sortedPairs[thisPair+1].first;
			double Qnew = leftQ + (thisPairDen-leftPairNum)*(rightQ-leftQ)/(rightPairNum-leftPairNum);

			// linear interpolation not good enough; need to solve numerically with Newton-Raphson
			// use Qnew as initial guess
			// need num_pdf in this interval
			const double num_pdf = (num_cdf[thisPair+1].second - num_cdf[thisPair].second) / (rightQ - leftQ);
			const double Numi = num_pdf * antiderivative[thisPair];
			// Solve antiderivative(Qnew) == c0
			const double c0 = (thisPairDen - leftPairNum + Numi) / num_pdf;	// DOUBLE CHECK THIS AND PREVIOUS STEPS

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
				Qnew -= (Nnew - c0) / deriv;

				ntries++;
			}

			const double arg = ( hadronBE[iPair1].p - hadronBE[iPair2].p )
								* ( hadronBE[iPair1].x - hadronBE[iPair2].x )
								* MM2FM / HBARC;
			// Do this for shifts and compensation shifts together
			if ( abs(arg) > 16.0 * M_PI )
			{
				pairShifts.push_back( 0.0 );
				pairCompensationShifts.push_back( Qnew - Q0 );
			}
			else
			{
				pairShifts.push_back( Qnew - Q0 );
				pairCompensationShifts.push_back( 0.0 );
			}

		}

	// Finally, discard fake pairs at beginning and end of sortedPairs
	sortedPairs.erase ( sortedPairs.begin() );
	sortedPairs.erase ( sortedPairs.end() );

	int pairIndex = 0;
	for (const auto & iPair : sortedPairs)
	{
		const int i1 	   = iPair.second.first;
		const int i2 	   = iPair.second.second;
		const double Qold  = iPair.first;
		double Qnew  = Qold + pairShifts[pairIndex];
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

		// Now get compensation shifts
		Qnew  = Qold + pairCompensationShifts[pairIndex];
		Q2new = Qnew*Qnew;

		// Calculate corresponding three-momentum shift.
		Q2Diff    = Q2new - Q2old;
		p2DiffAbs = (hadronBE[i1].p - hadronBE[i2].p).pAbs2();
		p2AbsDiff = hadronBE[i1].p.pAbs2() - hadronBE[i2].p.pAbs2();
		eSum      = hadronBE[i1].p.e() + hadronBE[i2].p.e();
		eDiff     = hadronBE[i1].p.e() - hadronBE[i2].p.e();
		sumQ2E    = Q2Diff + eSum * eSum;
		rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
		rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
		factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
						+ Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

		// Add shifts to sum. (Energy component dummy.)
		pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
		hadronBE[i1].pComp += pDiff;
		hadronBE[i2].pComp -= pDiff;

		pairIndex++;

	}

	return;
}*/


//==========================================================================

} // end namespace Pythia8
