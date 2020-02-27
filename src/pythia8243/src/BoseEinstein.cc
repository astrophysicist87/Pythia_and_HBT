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

const double BoseEinstein::dQ = 1e-4;
const double BoseEinstein::Qmaximum = 1.0;

const bool BoseEinstein::BE_VERBOSE = false;

//--------------------------------------------------------------------------

// Find settings. Precalculate table used to find momentum shifts.

bool BoseEinstein::init(Info* infoPtrIn, Settings& settings,
  ParticleData& particleData) {

  // Save pointer.
  infoPtr                = infoPtrIn;

  // Main flags.
  doPion 				          = settings.flag("BoseEinstein:Pion");
  doKaon 				          = settings.flag("BoseEinstein:Kaon");
  doEta 				          = settings.flag("BoseEinstein:Eta");
  useInvariantSize 		          = settings.flag("BoseEinstein:useInvariantSourceSize");
  useDistribution 		          = settings.flag("BoseEinstein:useDistribution");
  useRelativeDistance 	          = settings.flag("BoseEinstein:useRelativeDistance");
  useRestFrame 			          = settings.flag("BoseEinstein:useRestFrame");
  include_phase_space	          = settings.flag("BoseEinstein:includePhaseSpace");
  linear_interpolate_CDF          = settings.flag("BoseEinstein:linearInterpolateCDF");
  include_posDelQ_in_compensation = settings.flag("BoseEinstein:usePositiveShiftsForCompensation");
  compute_BE_enhancement_exactly  = settings.flag("BoseEinstein:computeBEEnhancementExactly");

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
				Vec4 xDiff = hadronBE.at(i1).x - hadronBE.at(i2).x;
				result += xDiff*xDiff;
				count += 1.0;
			}
		else
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
			{
				Vec4 x = hadronBE.at(i1).x;
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
				Vec4 xDiff = hadronBE.at(i1).x - hadronBE.at(i2).x;
				if ( useRestFrame ) xDiff.bstback( 0.5*(hadronBE.at(i1).p + hadronBE.at(i2).p) );
				result += xDiff.pAbs2();
				count += 1.0;
			}
		else
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1]; ++i1)
			{
				Vec4 x = hadronBE.at(i1).x;
				if ( useRestFrame ) x.bstback( hadronBE.at(i1).p );
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






//===========Added by Chris Plumberg================
	// if using debugging version, reset pion momentum to random value
	///*
	constexpr bool debugging = true;
	if ( debugging )
	{
		//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		//std::default_random_engine generator(seed);
		std::default_random_engine generator;

		const double scale = 5.0;
		std::normal_distribution<double> distribution(0.0, scale);
		// loop over particles in this event
		int lastPion = -1;
		for (int i = 0; i < event.size(); ++i)
		{
			Particle & p = event[i];
			if ( p.isHadron() and p.id() == 211 )
			{
				// record this particle index
				lastPion = i;

				//================================
				// random number generation here

				// reset positions
				p.tProd( 0.0 );
				p.xProd( 0.0*FM2MM * distribution( generator ) );
				p.yProd( 0.0*FM2MM * distribution( generator ) );
				p.zProd( FM2MM * distribution( generator ) );
				//const double tmp_rand = FM2MM * distribution( generator );
				//p.tProd( tmp_rand );
				//p.zProd( tmp_rand );

				const double pm = p.m();
				const double px = 0.0 * distribution( generator );
				const double py = 0.0 * distribution( generator );
				const double pz = 1.0 * distribution( generator );
				//const double Epi = sqrt(0.13957*0.13957 + px*px + py*py + pz*pz);
				const double Epi = sqrt(pm*pm + px*px + py*py + pz*pz);

				// reset momenta
				p.e( Epi );
				p.px( px );
				p.py( py );
				p.pz( pz );

				p.statusPos();

			}
		}

		if ( lastPion < 0 ) lastPion = event.size()-1;

		// Then add a bunch more pions...
		if (true)
		for (int i = 0; i < 10; i++)
		{
			// same "mother" for all
			int jNew = event.copy( lastPion, 99 );
			//event[ iNew ].p( hadronBE.at(i).p );

			// reset positions
			event[ jNew ].tProd( 0.0 );
			event[ jNew ].xProd( 0.0*FM2MM * distribution( generator ) );
			event[ jNew ].yProd( 0.0*FM2MM * distribution( generator ) );
			event[ jNew ].zProd( FM2MM * distribution( generator ) );
			//const double tmp_rand = FM2MM * distribution( generator );
			//event[ jNew ].tProd( tmp_rand );
			//event[ jNew ].zProd( tmp_rand );

			const double pm = event[ jNew ].m();
			const double px = 0.0 * distribution( generator );
			const double py = 0.0 * distribution( generator );
			const double pz = 1.0 * distribution( generator );
			//const double Epi = sqrt(0.13957*0.13957 + px*px + py*py + pz*pz);
			const double Epi = sqrt(pm*pm + px*px + py*py + pz*pz);

			// reset momenta
			event[ jNew ].e( Epi );
			event[ jNew ].px( px );
			event[ jNew ].py( py );
			event[ jNew ].pz( pz );

			event[ jNew ].statusPos();

			//cout << "Added particle in position = " << iNew << endl;
		}
	}

	//*/
//===========End of Chris Plumberg's addition================






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
      if ( event[i].id() == idNow && (event[i].isFinal() or debugging) )
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
					cout << "BoseEinsteinCheck: NPair = " << sortedPairs.size()-2 << " in iSpecies = " << iSpecies << endl;
//if (1) exit(0);
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
//cout << "(particle#=" << i << "): p = " << hadronBE.at(i).p;
//cout << "(particle#=" << i << "): pShift = " << hadronBE.at(i).pShift;
    eSumOriginal  += hadronBE.at(i).p.e();
/*cout    << setprecision(12) << setw(16)
                << "Original p and pShift: " << i
                << " (id=" << hadronBE.at(i).id << ")" << endl
                << "p=" << hadronBE.at(i).p
                << "pShift=" << hadronBE.at(i).pShift;*/
    hadronBE.at(i).p += hadronBE.at(i).pShift;
    hadronBE.at(i).p.e( sqrt( hadronBE.at(i).p.pAbs2() + hadronBE.at(i).m2 ) );
    eSumShifted   += hadronBE.at(i).p.e();
    eDiffByComp   += dot3( hadronBE.at(i).pComp, hadronBE.at(i).p)
                     / hadronBE.at(i).p.e();
/*
cout << "eDiffByComp at this step = " << eDiffByComp << endl;
cout << "pComp at this step = " << hadronBE.at(i).pComp;
cout << "p at this step = " << hadronBE.at(i).p;
cout << "e at this step = " << hadronBE.at(i).p.e() << endl;
cout 	<< setprecision(8)
		<< "Getting ready to balance energy budget (particle#=" << i << "): "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl;
*/

  }

//cout << "Made it here (2)" << endl;


///*
cout 	<< setprecision(8)
		<< "Balancing energy budget: "
		<< eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;
//*/

constexpr bool perform_compensation = false;

  // Iterate compensation shift until convergence.
  int iStep = 0;
  while ( perform_compensation
    && abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
    && abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
    && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
///*
cout 	<< setprecision(8)
		<< "Balancing energy budget(2): "
		<< compFac << "   " << eSumOriginal << "   " << eSumShifted << "   " << eDiffByComp << "   "
		<< COMPRELERR * eSumOriginal << "   " << COMPFACMAX * abs(eDiffByComp) << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " > " << COMPRELERR * eSumOriginal << endl
		<< "TEST: " << abs(eSumShifted - eSumOriginal) << " < " << COMPFACMAX * abs(eDiffByComp) << endl;
//*/
    eSumShifted      = 0.;
    eDiffByComp      = 0.;
    for (int i = 0; i < nStored[9]; ++i) {
      hadronBE.at(i).p += compFac * hadronBE.at(i).pComp;
///*
      hadronBE.at(i).pShift += compFac * hadronBE.at(i).pComp;
cout 	<< setprecision(8)
		<< "Net shift at this point: " << i << "   " << hadronBE.at(i).pShift;
//*/
      hadronBE.at(i).p.e( sqrt( hadronBE.at(i).p.pAbs2() + hadronBE.at(i).m2 ) );
      eSumShifted   += hadronBE.at(i).p.e();
      eDiffByComp   += dot3( hadronBE.at(i).pComp, hadronBE.at(i).p)
                       / hadronBE.at(i).p.e();
    }
  }

//cout << "Made it here (3)" << endl;

constexpr bool check_for_bad_events = false;

  // Error if no convergence, and then return without doing BE shift.
  // However, not grave enough to kill event, so return true.
  if ( perform_compensation
	and check_for_bad_events
	and abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal )
  {
    infoPtr->errorMsg("Warning in BoseEinstein::shiftEvent: "
      "no consistent BE shift topology found, so skip BE");
	cout << setprecision(16) << "BoseEinsteinCheck: This event did not pass! Check: "
			<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	infoPtr->setBECShifts( false );
//if (1) exit(8);
    return true;
  }
  else
  {
	cout << setprecision(16) << "BoseEinsteinCheck: This event passes! Check: "
			<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	infoPtr->setBECShifts( true );
//if (1) exit(8);
  }


  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE.at(i).iPos, 99);
    event[ iNew ].p( hadronBE.at(i).p );
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

	// get all values in vector first
	for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
	for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
	{
		Vec4 qPRF = hadronBE.at(i1).p - hadronBE.at(i2).p;
		qPRF.bstback( 0.5*( hadronBE.at(i1).p + hadronBE.at(i2).p ) );

		sortedPairs.push_back(
			std::make_pair(
				sqrt( m2(hadronBE.at(i1).p, hadronBE.at(i2).p) - m2Pair[iTab] ),
				std::make_pair(i1, i2)
			)
		);
		pairs_sorted_by_qxPRF.push_back( std::make_pair( qPRF.px(), std::make_pair(i1, i2) ) );
		pairs_sorted_by_qyPRF.push_back( std::make_pair( qPRF.py(), std::make_pair(i1, i2) ) );
		pairs_sorted_by_qzPRF.push_back( std::make_pair( qPRF.pz(), std::make_pair(i1, i2) ) );
	}

	// check if there are enough pairs of this species to do shift
	if (sortedPairs.size() < 2)
	{
		//cout << "Not enough sorted pairs of species with pid=" << idNow << endl;
		return false;
	}

	// THEN sort them (sorts on first column in ascending order automatically)
	sort( sortedPairs.begin(), sortedPairs.end() );
	sort( pairs_sorted_by_qxPRF.begin(), pairs_sorted_by_qxPRF.end() );
	sort( pairs_sorted_by_qyPRF.begin(), pairs_sorted_by_qyPRF.end() );
	sort( pairs_sorted_by_qzPRF.begin(), pairs_sorted_by_qzPRF.end() );

cout << "CHECK pairs:" << endl;
for (const auto & eachPair : sortedPairs)
{
	const int i1 = eachPair.second.first;
	const int i2 = eachPair.second.second;
	Vec4 q = hadronBE.at(i1).p - hadronBE.at(i2).p;
	Vec4 qPRF = q;
	qPRF.bstback( 0.5*( hadronBE.at(i1).p + hadronBE.at(i2).p ) );
	cout << i1 << "   " << i2 << "   " << setprecision(12) << eachPair.first << "   " << qPRF.pAbs() << endl;
	cout << i1 << "   " << i2 << "   " << q;
	cout << i1 << "   " << i2 << "   " << qPRF;
}

cout << endl << endl << endl << "CHECK qzPRF pairs:" << endl;
for (const auto & eachPair : pairs_sorted_by_qzPRF)
{
	const int i1 = eachPair.second.first;
	const int i2 = eachPair.second.second;
	cout << i1 << "   " << i2 << "   " << setprecision(12) << eachPair.first << endl;
}


if (1) exit(8);

	// add fake first "pair"
	sortedPairs.insert(sortedPairs.begin(), std::make_pair( 0.0, std::make_pair(-1, -1) ) );

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
  double Q2old = m2(hadronBE.at(i1).p, hadronBE.at(i2).p) - m2Pair[iTab];
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
  double p2DiffAbs = (hadronBE.at(i1).p - hadronBE.at(i2).p).pAbs2();
  double p2AbsDiff = hadronBE.at(i1).p.pAbs2() - hadronBE.at(i2).p.pAbs2();
  double eSum      = hadronBE.at(i1).p.e() + hadronBE.at(i2).p.e();
  double eDiff     = hadronBE.at(i1).p.e() - hadronBE.at(i2).p.e();
  double sumQ2E    = Q2Diff + eSum * eSum;
  double rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  double rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  double factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

  // Add shifts to sum. (Energy component dummy.)
  Vec4   pDiff     = factor * (hadronBE.at(i1).p - hadronBE.at(i2).p);
  hadronBE.at(i1).pShift += pDiff;
  hadronBE.at(i2).pShift -= pDiff;

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
  pDiff     = factor * (hadronBE.at(i1).p - hadronBE.at(i2).p);
  hadronBE.at(i1).pComp += pDiff;
  hadronBE.at(i2).pComp -= pDiff;

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
		if ( BE_VERBOSE )
			cout << "linearly interpolating!" << endl;
		// Take the derivative of the linear interpolant.

		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
			denBar.push_back( 1.0 / ( sortedPairs.at(iPair+1).first - sortedPairs.at(iPair).first ) );

	}
	else
	{
		if ( BE_VERBOSE )
			cout << "NOT linearly interpolating!" << endl;
		// Do not estimate pair density.
		for (int iPair = 0; iPair < (int)sortedPairs.size()-1; ++iPair)
			denBar.push_back( 1.0 );
	}
	
	// Density must go to zero eventually.
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
		    Vec4 xDiffPRF = ( hadronBE.at(i1).x - hadronBE.at(i2).x ) * MM2FM / HBARC;
		    xDiffPRF.bstback( 0.5*(hadronBE.at(i1).p + hadronBE.at(i2).p) );

		sorted_xDiffs[xDiff_pairCount++] = xDiffPRF.pAbs();
	}
	sort( sorted_xDiffs.begin(), sorted_xDiffs.end() );
	return;
}



//-------------------------------------------------
// Set lefthand side of shift relation.
void BoseEinstein::set_effective_source(
			vector< pair< double, pair <int,int> > > & sortedPairs, int iTab
			)
{
	//dummy statement
	if (iTab){;}

	//-------------------------------
	// Estimate scale for Q sampling.
	//double size_estimate_sum = 0.0;
	//double size_estimate_sqsum = 0.0
	//double size_estimate_cubsum = 0.0;
	//double mean = 0.0
	//double rms = 0.0
	//double avg_log = 0.0;

	//-----------
	// Set grids.
	const int Qgridsize = static_cast<int>(Qmaximum / dQ) + 1;
	Qgrid.assign( Qgridsize, 0.0 );
	phase_space.assign( Qgridsize, 0.0 );
	effSource.assign( Qgridsize, 0.0 );
	integrated_effective_source.assign( Qgridsize, 0.0 );

	//const double npairs = sortedPairs.size()-2;
	//cout << "npairs = " << (int)npairs << endl;

	for (int iPair = 1; iPair < (int)sortedPairs.size()-1; iPair++)
	{

		auto eachPair = sortedPairs.at(iPair);

		const int i1 = eachPair.second.first;
		const int i2 = eachPair.second.second;
		Vec4 xDiffPRF = ( hadronBE.at(i1).x - hadronBE.at(i2).x ) * MM2FM / HBARC;
		xDiffPRF.bstback( 0.5*(hadronBE.at(i1).p + hadronBE.at(i2).p) );

		const double xDiffval = xDiffPRF.pAbs();

		/*size_estimate_sum += 1.0 / xDiffval;
		size_estimate_sqsum += 1.0 / (xDiffval*xDiffval);
		size_estimate_cubsum += 1.0 / (xDiffval*xDiffval*xDiffval);
		mean += xDiffval;
		rms += xDiffval*xDiffval;
		avg_log += log(1.0 / xDiffval);*/

		double currentQ = 0.0;
		for (int iQ = 0; iQ < (int)effSource.size(); iQ++)
		{
			//effSource[iQ] += sphericalbesselj0(currentQ * xDiffval);
			effSource[iQ] += cos(currentQ * xDiffval);
			currentQ += dQ;
		}
	}

	/*
	cout << setprecision(12) << "size estimate (linear halfwidth) = "
			<< 2.0 * size_estimate_sqsum / (M_PI * size_estimate_sum) << endl;
	cout << setprecision(12) << "size estimate (gaussian approximation / 2) = "
			<< 0.5 * sqrt(3.0 * size_estimate_cubsum / size_estimate_sum) << endl;
	cout << setprecision(12) << "size estimate (geometric mean) = "
			<< exp(avg_log / npairs) << endl;
	*/

	double currentQ = 0.0;
	/*if (include_phase_space)
		for (int iQ = 0; iQ < (int)Qgrid.size(); iQ++)
		{
			phase_space[iQ] = currentQ*currentQ / sqrt(currentQ*currentQ + m2Pair[iTab]);
			currentQ += dQ;
		}
	else*/
		for (int iQ = 0; iQ < (int)Qgrid.size(); iQ++)
			phase_space[iQ] = 1.0;

	//reset
	double previousIntegrand = 0.0, currentIntegrand = 0.0;
	currentQ = 0.0;
	double runningIntegral = 0.0;
	for (int iQ = 0; iQ < (int)Qgrid.size(); iQ++)
	{
		Qgrid[iQ] = currentQ;
		effSource[iQ] /= (double)(sortedPairs.size()-2);
		previousIntegrand = currentIntegrand;
		currentIntegrand = phase_space[iQ] * effSource[iQ];
		runningIntegral += 0.5*dQ*(currentIntegrand + previousIntegrand);
		integrated_effective_source[iQ] = runningIntegral;
		cout << "effSource: " << setprecision(12)
				<< currentQ << "   "
				<< phase_space[iQ] << "   "
				<< effSource[iQ] << "   "
				<< integrated_effective_source[iQ] << endl;
		currentQ += dQ;
	}
	//if (1) return;
	if (1) exit (8);

}


//-------------------------------------------------
// Set lefthand side of shift relation.
void BoseEinstein::set_LHS(
			vector< pair< double, pair <int,int> > > & sortedPairs,
			vector< pair< double, double > > & LHS,
			vector<double> & denBar, int iTab
			)
{
	//dummy statement
	if (iTab){;}

	int pairCount = 0;
	double result = 0.0;
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


	// Allows me to compute effective 
	if (not compute_BE_enhancement_exactly)
		set_effective_source( sortedPairs, iTab );

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
	if (iTab){;}

	int pairCount = 0;
	double result = 0.0;

	// Compute two-body enhancement here
	for (const auto & thisPair : sortedPairs)
	{
		//auto nextPair = ( pairCount == (int)sortedPairs.size()-1 ) ? thisPair : *(&thisPair+1);
		const double thisQ = thisPair.first;
		//const double nextQ = nextPair.first;

		RHS.push_back( std::make_pair( thisQ, result ) );

		if ( pairCount == (int)denBar.size() )
			continue;

		if (sorted_xDiffs.size() > 0){;}

		//const double one_by_N = 1.0 / static_cast<double>(sortedPairs.size() - 2);
		//int eachPairIndex = 0;
		for (const auto & eachPair : sortedPairs)
		{
			// Skip unphysical dummy pairs.
			if (   &eachPair == &sortedPairs.front()
				or &eachPair == &sortedPairs.back() )
				continue;

			//const double xDiffPRFVal = sorted_xDiffs.at(eachPairIndex++);
			
			// Add in results from this pair
			/*result += one_by_N * denBar.at(pairCount)
						* compute_integral_without_phasespace(
							thisQ, nextQ, xDiffPRFVal);*/

		}

		pairCount++;
	}

	// Add in constant piece to RHS result.
	for (int iPair = 0; iPair < (int)LHS.size(); ++iPair)
		RHS.at(iPair).second += LHS.at(iPair).second;


	// check LHS and RHS at this point
	//cout << "Line = " << __LINE__ << endl;
	/*cout << "Sizes: " << sortedPairs.size() << "   "
			<< LHS.size() << "   " << RHS.size() << "   "
			<< sorted_xDiffs.size() << endl;

	for (int iPair = 0; iPair < (int)LHS.size(); ++iPair)
		cout << "CHECK: "  << setprecision(12)
				<< sortedPairs.at(iPair).first << "   "
				<< LHS.at(iPair).second << "   "
				<< RHS.at(iPair).second << endl;*/

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

	vector<double> sorted_xDiffs(sortedPairs.size()-2);
	set_sorted_xDiffs( sortedPairs, sorted_xDiffs );

	//------------------
	// Set LHS integral.
	set_LHS(sortedPairs, LHS, denBar, iTab);

	//------------------
	// Set RHS integral.
	set_RHS(sortedPairs, sorted_xDiffs, LHS, RHS, denBar, iTab);

   auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
	if ( BE_VERBOSE )
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


	// vector to track sorted pair separations
	vector<double> sorted_xDiffs(sortedPairs.size()-2);
	set_sorted_xDiffs( sortedPairs, sorted_xDiffs );


	// keep track of which pairs are shifted
	int n_skipped_pairs = 0;
	vector<bool> this_pair_shifted(sortedPairs.size()-2, false);


	// -------------------------------------------
	// Finally, get shifts for each physical pair.
	for (int iPair = 1; iPair < (int)sortedPairs.size()-1; ++iPair)
	{
		int thisPair = iPair;

		const double Q0          = sortedPairs.at(iPair).first;          // original Q
		const int iPair1         = sortedPairs.at(iPair).second.first;   // first hadron
		const int iPair2         = sortedPairs.at(iPair).second.second;	 // second hadron
		const double thisPairLHS = LHS.at(thisPair).second;              // > 0
		double RHS_lower         = RHS.at(thisPair).second;
		double RHS_upper         = RHS.at(thisPair+1).second;

		//-------------------------------------------------
		// If thisPairLHS < 0, skip this pair.
		if ( thisPairLHS < 0.0 )
		{
			n_skipped_pairs++;
			pairShifts.push_back( 0.0 );
			pairCompensationShifts.push_back( 0.0 );
			continue;
		}



		//cout << "Test: " << Q0 << "   " << iPair << "   " << thisPair << "   "
		//		<< RHS.size() << "   " << iPair1 << "   " << iPair2 << "   "
		//		<< RHS_lower << "   " << thisPairLHS << "   " << RHS_upper << endl;


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


		//cout << "Test again: " << Q0 << "   " << iPair << "   " << thisPair << "   "
		//		<< RHS.size() << "   " << iPair1 << "   " << iPair2 << "   "
		//		<< RHS_lower << "   " << thisPairLHS << "   " << RHS_upper << endl;
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


		// If Qlower is too large, skip this pair.
		if ( Qlower >= Qmaximum )
		{
			n_skipped_pairs++;
			pairShifts.push_back( 0.0 );
			pairCompensationShifts.push_back( 0.0 );
			continue;
		}


		//--------------------------------
		// Initialize to closest endpoint.
		double Qnew = Qlower;
		double RHSnew = RHS_lower;
		if ( abs(c0-RHS_lower) > abs(c0-RHS_upper) or Qnew < 1.e-10 )
		{
			Qnew   = Qupper;
			RHSnew = RHS_upper;
		}


		//-----------------
		// Start iteration.
		int ntries = 0;
		bool success = true;	// assume true unless something fails
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
				int eachPairIndex = 0;
				for (const auto & eachPair : sortedPairs)
				{
					// Skip unphysical dummy pairs.
					if (   &eachPair == &sortedPairs.front()
						or &eachPair == &sortedPairs.back() )
						continue;

					const double xDiffMag = sorted_xDiffs.at(eachPairIndex++);
		
					// Add in results for this pair
					//RHSnew += one_by_N * thisdenBar * compute_integral_without_phasespace( Qlower, Qnew, xDiffMag );
					deriv += one_by_N * thisdenBar * cos( Qnew * xDiffMag );

					//pairCount++;
				}
			}

			//if ( ntries>0 and ( RHSnew > RHS_upper or RHSnew < RHS_lower ) )
			//	cout << "WARNING: RHSnew = " << RHSnew << " not in range (" << RHS_lower << ", " << RHS_upper << ")!" << endl;
			if ( BE_VERBOSE and (Qnew > Qupper or Qnew < Qlower) )
				cout << "WARNING (ntries=" << ntries << "): Qnew = " << Qnew << " not in range (" << Qlower << ", " << Qupper << ")!" << endl;
			//cout << "NR method (ntries=" << ntries << "): " << Qnew << "   " << RHSnew << "   " << deriv << "   " << c0 << endl;	

			//----------------------------
			// Get next estimate for Qnew.
			Qnew -= (RHSnew - c0) / deriv;

			ntries++;
		}


		//-------------------------------------------------
		// If shifting this pair was unsuccessful, move on.
		if ( not success )
		{
			n_skipped_pairs++;
			pairShifts.push_back( 0.0 );
			pairCompensationShifts.push_back( 0.0 );
			continue;
		}


		//-------------------------------------------
		// Store results based on size of BE effects.
		Vec4 p1 = hadronBE.at(iPair1).p;
		Vec4 p2 = hadronBE.at(iPair2).p;
		Vec4 x1 = hadronBE.at(iPair1).x;
		Vec4 x2 = hadronBE.at(iPair2).x;
		Vec4 xDiffPRF = ( x1 - x2 ) * MM2FM / HBARC;
		xDiffPRF.bstback( 0.5*(p1 + p2) );


		pairShifts.push_back( Qnew - Q0 );
		pairCompensationShifts.push_back( 0.0 );

		if ( BE_VERBOSE )
			cout << "BoseEinsteinCheck: CHECK shift is "
					<< Qnew - Q0 /*<< " vs. " << Qnew_LININT - Q0*/ << " of " << Q0 << endl;

		// Shift pairs closer together; use pairs
		// shifted further away to compensate.

		// Use this pair in compensation instead if, either
		// (1) - the shift was positive and we want to use these for compensation, or
		// (2) - the magnitude of the shift was larger than some % of the original Q0
		bool shift_failure = ( include_posDelQ_in_compensation and Qnew > Q0 )
								//or abs(Qnew - Q0) > 0.05 * Q0
								;

		bool shift_success = not shift_failure;

		//if ( not include_posDelQ_in_compensation or Qnew <= Q0 )
		if ( shift_success )
			this_pair_shifted.at(iPair-1) = true;

	}



	// Finally, discard fake pairs at beginning and end of sortedPairs
	sortedPairs.erase ( sortedPairs.begin() );
	sortedPairs.erase ( sortedPairs.end() );

	int pairIndex = 0;
	int number_of_pairs_shifted = 0, number_of_pairs_not_shifted = 0;
	for (const auto & iPair : sortedPairs)
	{

		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		constexpr bool rescale_pair_momenta = true;

		// Get the computed shifts
		const double Qold  = iPair.first;
		double Qnew  = Qold + pairShifts.at(pairIndex);
		const double Q2old = Qold*Qold;
		double Q2new = Qnew*Qnew;

		// Calculate corresponding three-momentum shift.
		double Q2Diff    = Q2new - Q2old;
		double p2DiffAbs = (hadronBE.at(i1).p - hadronBE.at(i2).p).pAbs2();
		double p2AbsDiff = hadronBE.at(i1).p.pAbs2() - hadronBE.at(i2).p.pAbs2();
		double eSum      = hadronBE.at(i1).p.e() + hadronBE.at(i2).p.e();
		double eDiff     = hadronBE.at(i1).p.e() - hadronBE.at(i2).p.e();
		double sumQ2E    = Q2Diff + eSum * eSum;
		double rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
		double rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
		double factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
							+ Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

		// Add shifts to sum. (Energy component dummy.)
		Vec4   pDiff     = factor * (hadronBE.at(i1).p - hadronBE.at(i2).p);


		if ( rescale_pair_momenta or this_pair_shifted.at(pairIndex) )
		{
			// Compute appropriate shift for pair
			number_of_pairs_shifted++;

			hadronBE.at(i1).pShift += pDiff;
			hadronBE.at(i2).pShift -= pDiff;

			if ( rescale_pair_momenta )
			{
				// add symmetrically to both momenta
				pDiff = factor * (hadronBE.at(i1).p + hadronBE.at(i2).p);
				hadronBE.at(i1).pComp += pDiff;
				hadronBE.at(i2).pComp += pDiff;
			}
		}
		else
		{
			// Use computed shift for compensation
			number_of_pairs_not_shifted++;

			//*/
			//Vec4 pDiff = hadronBE.at(i1).p - hadronBE.at(i2).p;
			hadronBE.at(i1).pComp += pDiff;
			hadronBE.at(i2).pComp -= pDiff;
		}

		pairIndex++;

	}

	if ( BE_VERBOSE or true )
	{
		cout << "number_of_pairs_shifted = " << number_of_pairs_shifted << endl;
		cout << "number_of_pairs_not_shifted = " << number_of_pairs_not_shifted << endl;
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
