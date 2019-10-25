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
//int checkpoint = 0;
  // Reset list of identical particles.
  hadronBE.resize(0);

	//===========Added by Chris Plumberg================
	// if using debugging version, reset pion momentum to random value
	///*
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	const double scale = 5.0;
	std::normal_distribution<double> distribution(0.0, scale);
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;
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
			p.xProd( FM2MM * distribution( generator ) );
			p.yProd( FM2MM * distribution( generator ) );
			p.zProd( FM2MM * distribution( generator ) );

			const double pm = p.m();
			const double px = 0.1 * distribution( generator );
			const double py = 0.1 * distribution( generator );
			const double pz = 0.1 * distribution( generator );
			//const double Epi = sqrt(0.13957*0.13957 + px*px + py*py + pz*pz);
			const double Epi = sqrt(pm*pm + px*px + py*py + pz*pz);

			// reset momenta
			p.e( Epi );
			p.px( px );
			p.py( py );
			p.pz( pz );

		}
	}
	
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

if ( lastPion < 0 ) lastPion = event.size()-1;

	// Then add a bunch more pions...
	if (true)
	for (int i = 0; i < 1000; i++)
	{
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;
		// same "mother" for all
		int jNew = event.copy( lastPion, 99 );
		//event[ iNew ].p( hadronBE[i].p );

//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;
//cout << lastPion << "   " << jNew << "   " << event.size() << endl;
		// reset positions
		event[ jNew ].tProd( 0.0 );
		event[ jNew ].xProd( FM2MM * distribution( generator ) );
		event[ jNew ].yProd( FM2MM * distribution( generator ) );
		event[ jNew ].zProd( FM2MM * distribution( generator ) );
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

		const double pm = event[ jNew ].m();
		const double px = 0.1 * distribution( generator );
		const double py = 0.1 * distribution( generator );
		const double pz = 0.1 * distribution( generator );
		//const double Epi = sqrt(0.13957*0.13957 + px*px + py*py + pz*pz);
		const double Epi = sqrt(pm*pm + px*px + py*py + pz*pz);

//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

		// reset momenta
		event[ jNew ].e( Epi );
		event[ jNew ].px( px );
		event[ jNew ].py( py );
		event[ jNew ].pz( pz );

//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

		//cout << "Added particle in position = " << iNew << endl;
	}

//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;
	//*/
	//===========End of Chris Plumberg's addition================


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
        hadronBE.push_back(
          BoseEinsteinHadron( idNow, i, event[i].p(), event[i].m(), event[i].vProd() ) );
    nStored[iSpecies + 1] = hadronBE.size();

//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;
    // Loop through pairs of identical particles and find shifts.
	switch ( enhanceMode )
	{
		case 0: // the original and the default
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_fixedQRef( i1, i2, iTab);
			break;
		case 1:	// use a new Gaussian enhancement based on space-time interval between production points
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_STint_GaussBE( i1, i2, iTab);
			break;
		case 2:	// use a new cos(q*x) enhancement based on space-time interval between production points
			for (int i1 = nStored[iSpecies]; i1 < nStored[iSpecies+1] - 1; ++i1)
			for (int i2 = i1 + 1; i2 < nStored[iSpecies+1]; ++i2)
			  shiftPair_STint_SphBesselBE( i1, i2, iTab);
			break;
		default:
			// Do nothing.
			break;
	}
  }
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

  // Must have at least two pairs to carry out compensation.
  if (nStored[9] < 2) return true;

  // Shift momenta and recalculate energies.
  double eSumOriginal = 0.;
  double eSumShifted  = 0.;
  double eDiffByComp  = 0.;
  for (int i = 0; i < nStored[9]; ++i) {
    eSumOriginal  += hadronBE[i].p.e();
    hadronBE[i].p += hadronBE[i].pShift;
    hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
    eSumShifted   += hadronBE[i].p.e();
    eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                     / hadronBE[i].p.e();
  }
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

  // Iterate compensation shift until convergence.
  int iStep = 0;
  while ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
    && abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
    && iStep < NCOMPSTEP ) {
    ++iStep;
    double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
    eSumShifted      = 0.;
    eDiffByComp      = 0.;
    for (int i = 0; i < nStored[9]; ++i) {
      hadronBE[i].p += compFac * hadronBE[i].pComp;
      hadronBE[i].p.e( sqrt( hadronBE[i].p.pAbs2() + hadronBE[i].m2 ) );
      eSumShifted   += hadronBE[i].p.e();
      eDiffByComp   += dot3( hadronBE[i].pComp, hadronBE[i].p)
                       / hadronBE[i].p.e();
    }
  }

/*cout << "CHECK ERROR: "
		<< setprecision(16)
		<< abs(eSumShifted - eSumOriginal) << "   "
		<< COMPRELERR * eSumOriginal << "   "
		<< COMPFACMAX * abs(eDiffByComp) << endl;*/
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

  // Error if no convergence, and then return without doing BE shift.
  // However, not grave enough to kill event, so return true.
  if ( abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal ) {
    infoPtr->errorMsg("Warning in BoseEinstein::shiftEvent: "
      "no consistent BE shift topology found, so skip BE");
    return true;
  }

  // Store new particle copies with shifted momenta.
  for (int i = 0; i < nStored[9]; ++i) {
    int iNew = event.copy( hadronBE[i].iPos, 99);
    event[ iNew ].p( hadronBE[i].p );
  }
//cout << "checkpoint #" << checkpoint++ << " at " << __LINE__ << endl;

  /*cout << "BoseEinstein(): pair counts = "
		<< number_of_pairs << "   "
		<< number_of_shifted_pairs << "   "
		<< number_of_too_close_pairs << "   "
		<< number_of_too_separated_pairs << endl;*/

	//===========Added by Chris Plumberg================
	// just for fun
	/*for (int i = 0; i < 100; i++)
	{
		//int iNew = event.append( Particle( 211 ) );
		int iNew = event.copy( event.size()-1, 101 );
		cout << "Added particle in position = " << iNew << endl;
	}*/
	/*for (int i = 0; i < event.size(); ++i)
	{
		Particle & p = event[i];
		if ( p.isHadron() )
		{
	 		if ( p.id() == 211 )	// i.e., is pi^+
			{
				cout << "this is another check: "
						<< setprecision(16)
						<< "   " << p.e()
						<< "   " << p.px()
						<< "   " << p.py()
						<< "   " << p.pz()
						<< "   " << (p.vProd())[0]
						<< "   " << (p.vProd())[1]
						<< "   " << (p.vProd())[2]
						<< "   " << (p.vProd())[3]
						<< endl;
			}
		}
	}*/
	//===========End of Chris Plumberg's addition================


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
void BoseEinstein::shiftPair_STint_GaussBE( int i1, int i2, int iTab) {

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
// interval and spherical Bessel form.
void BoseEinstein::shiftPair_STint_SphBesselBE( int i1, int i2, int iTab) {

	//======================================
	// Start of modified initializations

	// phase space info for this pair
	Vec4 x1 = hadronBE[i1].x;
	Vec4 x2 = hadronBE[i2].x;
	Vec4 p1 = hadronBE[i1].p;
	Vec4 p2 = hadronBE[i2].p;

	// define pair and relative momenta and coordinate separation
	Vec4 xDiff = x1 - x2;
	Vec4 pairMomentumK = 0.5 * ( p1 + p2 );
	Vec4 relMomentumq  = p1 - p2;

	// evaluate momentum shift in pair rest frame (PRF)
	xDiff.bstback( pairMomentumK );
	relMomentumq.bstback( pairMomentumK );

	// Coordinate separation in GeV^{-1}.
	xDiff *= MM2FM / HBARC;

	// set relevant scales
	RRef = xDiff.pAbs();
	//RRef = 5.0;	// just for debugging
	QRef = 1.0 / RRef;

	// --------------------------------------------------------------
	// tally how many pairs we're looking at of the total possible
	number_of_pairs++;

	// Check that QRef will not be too large or too small: 0.05 <= QRef <= 1.0
	/*if ( RRef < 1.0 or RRef > 20.0 )
	{
		if ( RRef < 1.0 )
			number_of_too_close_pairs++;
		else
			number_of_too_separated_pairs++;
		return;
	}*/

	number_of_shifted_pairs++;
	// --------------------------------------------------------------


	// if we passed cuts on QRef, then proceed...

	const int NSTEP = 199;
	//const int NSTEP = 999;

	// Set various tables on a per-pair basis.
	double Qnow, Q2now, centerCorr;
    // Step size and number of steps in normal table.
    deltaQ[iTab]      = STEPSIZE * min(mPair[iTab], QRef);
    nStep[iTab]       = min( NSTEP, 1 + int(3. * QRef / deltaQ[iTab]) );
    maxQ[iTab]        = (nStep[iTab] - 0.1) * deltaQ[iTab];
    centerCorr        = deltaQ[iTab] * deltaQ[iTab] / 12.;

    // Construct normal table recursively in Q space.
    shift[iTab][0]    = 0.;
    for (int i = 1; i <= nStep[iTab]; ++i) {
      Qnow            = deltaQ[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift[iTab][i]  = shift[iTab][i - 1]
						+ sphericalbesselj0(Qnow * RRef)
						* deltaQ[iTab] * (Q2now + centerCorr)
						/ sqrt(Q2now + m2Pair[iTab]);
						// recall: Q is Lorentz scalar and can
						//         be evaluated in any frame
    }

    // Step size and number of steps in compensation table.
    //deltaQ3[iTab]     = STEPSIZE * min(mPair[iTab], QRef3);
    deltaQ3[iTab]     = STEPSIZE * mPair[iTab];
    nStep3[iTab]      = min( NSTEP, 1 + int(9. * QRef / deltaQ3[iTab]) );
    maxQ3[iTab]       = (nStep3[iTab] - 0.1) * deltaQ3[iTab];
    centerCorr        = deltaQ3[iTab] * deltaQ3[iTab] / 12.;

    // Construct compensation table recursively in Q space.
    shift3[iTab][0]   = 0.;
    for (int i = 1; i <= nStep3[iTab]; ++i) {
      Qnow            = deltaQ3[iTab] * (i - 0.5);
      Q2now           = Qnow * Qnow;
      shift3[iTab][i] = shift3[iTab][i - 1] + 1.0//*exp(-Q2now * R2Ref3)
        * deltaQ3[iTab] * (Q2now + centerCorr) / sqrt(Q2now + m2Pair[iTab]);
    }
	// End of modified initializations
	//======================================

	// Calculate old relative momentum.
	//double Q2old = m2(p1, p2) - m2Pair[iTab]; // THIS IS THE WRONG Q2!!!!!
	double Q2old = relMomentumq.pAbs2();
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

	double my_alternate_result = 0.0;
	//double c0_glob = 0.0;

	bool try_alternate_shift_calculation = true;
	if ( try_alternate_shift_calculation )
	{
		const double Q0 = Qold;

		// compute constant integral
		double c0 = 0.0;
		/*for (int ipt = 0; ipt < 7; ipt++)
		{
			const double cen = 0.5 * Q0;
			const double hw = 0.5 * Q0;
			double qloc = cen + hw * x_pts_7[ipt];
			c0 += hw * x_wts_7[ipt] * pow2(qloc)
				* sphericalbesselj0(qloc*RRef)
				/ sqrt(pow2(qloc) + m2Pair[iTab]);
		}*/
		const double max = Q0 * RRef;
		const int nPeriods = static_cast<int>( 0.5*max/M_PI );

		// do any whole periods which will fit
		for ( int iPeriod = 0; iPeriod < nPeriods; iPeriod++ )
		{
			double tmp = 0.0;
			const double cen = (2.0 * iPeriod + 1.0) * M_PI;
			const double hw = M_PI;
			for ( int ipt = 0; ipt < 7; ipt++ )
			{
				double xloc = cen + hw * x_pts_7[ipt];
				tmp += hw * x_wts_7[ipt] * xloc * sin(xloc)
					/ sqrt(pow2(xloc) + m2Pair[iTab]*pow2(RRef));
			}
			//if (iTab == 2)
			//	cout << setprecision(8) << "CHECK: " << iPeriod << "(" << cen - hw << "," << cen + hw << "): "
			//			<< tmp << "   " << c0 << endl;
			c0 += tmp;
		}
		
		// now do remainder
		double tmp2 = 0.0;
		const double cen0 = 0.5 * ( max + 2.0 * nPeriods * M_PI );
		const double hw0 = 0.5 * ( max - 2.0 * nPeriods * M_PI );
		for ( int ipt = 0; ipt < 7; ipt++ )
		{
			double xloc = cen0 + hw0 * x_pts_7[ipt];
			tmp2 += hw0 * x_wts_7[ipt] * xloc * sin(xloc)
				/ sqrt(pow2(xloc) + m2Pair[iTab]*pow2(RRef));
		}
		c0 += tmp2;
		//if (iTab == 2)
		//	cout << setprecision(8) << "CHECK: " << nPeriods << "(" << cen0 - hw0 << "," << cen0 + hw0 << "): "
		//			<< tmp2 << "   " << c0 << "   " << c0/pow2(RRef) << endl;

		
		c0 /= pow2(RRef);

		//if (iTab == 2)
		//	cout << setprecision(8) << "CHECK MY VERSION: "
		//			<< Qold*RRef << "   " << mPair[iTab]*RRef << ";   "
		//			<< Qold << "   " << 0.5*mPair[iTab] << "   " << RRef << ":   "
		//			<< c0 << endl;

		// use the result to get initial estimate for delta Q
		//c0_glob = c0;
		double initial_estimate
			= -c0 * sqrt(pow2(Q0)+m2Pair[iTab])
				/ ( pow2(Q0) *(1.0+sphericalbesselj0(Q0*RRef)) );

		double current_estimate = initial_estimate;

		/*cout << setprecision(8) << "CHECK MY VERSION: "
				<< Qold << "   " << sqrt(Q2new) << "   "
				<< 0.0 << "   " << initial_estimate << "   "
				<< Qold + initial_estimate << endl;*/

		// sum must be 0 (solve with Newton's method)
		// iterate as needed
		const int NITER = 100;
		const double ACCURACYGOAL = 1.e-8;
		double ratio = 0.0;
		int ITER = 0;
		do
		{
			const double a = Q0;
			const double b = Q0 + current_estimate;
			const double cen = 0.5 * (b + a);
			const double hw = 0.5 * (b - a);

			// get integral evaluation at current point
			double f_val = c0;
			for (int ipt = 0; ipt < 7; ipt++)
			{
				double qloc = cen + hw * x_pts_7[ipt];
				f_val += hw * x_wts_7[ipt] * pow2(qloc)
					* (1.0 + sphericalbesselj0(qloc*RRef))
					/ sqrt(pow2(qloc) + m2Pair[iTab]);
			}

			// get derivative evaluation at current point
			double fp_val = pow2(b)
					* (1.0 + sphericalbesselj0(b * RRef) )
					/ sqrt(pow2(b) + m2Pair[iTab]);

			ratio = f_val / fp_val;
			current_estimate -= ratio;
			/*cout << setprecision(8) << "CHECK MY VERSION: "
					<< Qold << "   " << sqrt(Q2new) << "   "
					<< c0 << "   " << ITER << "   " << abs(ratio) << "   "
					<< f_val << "   " << current_estimate << "   "
					<< Qold + current_estimate << endl;*/
			ITER++;
		} while ( ITER < NITER and abs(ratio) > ACCURACYGOAL );

		my_alternate_result = current_estimate;

		// I don't think I need to do this
		//Q2new = Q2old * pow( Qold / (Qold + 3. * lambda * my_alternate_result), 2. / 3.);

		// just do this instead
		Q2new = pow2(Qold + my_alternate_result);
	}

	/*
	//if (iTab == 2)
		cout << setprecision(8)
				<< "CHECK MOMENTUM SHIFT (EnMode = 2): "
				<< Qold << "   "
				//<< relMomentumq.pAbs() << "   "
				//<< sqrt(m2(p1, p2) - m2Pair[iTab]) << "   "
				<< sqrt(Q2new) << "   "
				<< Qold + my_alternate_result << "   "
				//<< sqrt(Q2old * pow( Qold / (Qold - 3. * my_alternate_result), 2. / 3.)) << "   "
				//<< Qmove << "   "
				<< 0.5*mPair[iTab]
				//<< ";   ("
				//<< (p1 + p2).m2Calc() << " =?= "
				//<< m2Pair[iTab] - (p1 - p2).m2Calc()<< ")   " 
				//<< relMomentumq << "   " << p1 << "   " << p2 << "   "
				//<< c0_glob << "   "
				//<< RRef
				<< endl;
	*/

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

	//=================================
	// Using old method here, so keep the correction (stuff)^(2/3) part
	//=================================
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
  double Q2new3 = Q2old * pow( Qold / (Qold + 3. * Qmove3), 2. / 3.);

  // Calculate corresponding three-momentum shift.
  Q2Diff    = Q2new3 - Q2old;
  sumQ2E    = Q2Diff + eSum * eSum;
  rootA     = eSum * eDiff * p2AbsDiff - p2DiffAbs * sumQ2E;
  rootB     = p2DiffAbs * sumQ2E - p2AbsDiff * p2AbsDiff;
  factor    = 0.5 * ( rootA + sqrtpos(rootA * rootA
    + Q2Diff * (sumQ2E - eDiff * eDiff) * rootB) ) / rootB;

//cout << "CHECK THIS (EnMode = 2): " << Q2Diff << "   " << sumQ2E << "   " << rootA << "   " << rootB << "   " << factor << endl;

  // Extra dampening factor to go from BE_3 to BE_32.
  //factor   *= 1. - exp(-Q2old * R2Ref2);

  // Add shifts to sum. (Energy component dummy.)
  pDiff     = factor * (hadronBE[i1].p - hadronBE[i2].p);
  hadronBE[i1].pComp += pDiff;
  hadronBE[i2].pComp -= pDiff;

}

///===CJP(end)===

//==========================================================================

} // end namespace Pythia8
