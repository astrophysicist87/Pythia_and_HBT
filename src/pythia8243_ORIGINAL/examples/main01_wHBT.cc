// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main()
{
	// Generator. Process selection. LHC initialization. Histogram.
	Pythia pythia;
	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");
	pythia.readString("Beams:eCM = 7000.0");
	pythia.readString("Beams:frameType = 1");
	pythia.readString("SoftQCD:all = on");

	pythia.readString("Fragmentation:setVertices = on");
	pythia.readString("PartonVertex:setVertex = off");
	pythia.readString("HadronLevel:BoseEinstein = on");
	pythia.readString("BoseEinstein:QRef = 0.2");
	pythia.readString("ColourReconnection:reconnect = off");
	pythia.readString("Ropewalk:RopeHadronization = off");
	pythia.readString("Ropewalk:doShoving = off");
	pythia.readString("Ropewalk:doFlavour = off");

	pythia.init();

	// Begin event loop. Generate event. Skip if error. List first one.
	for (int iEvent = 0; iEvent < 10; ++iEvent)
	{
		if (!pythia.next()) continue;

		int event_multiplicity = 0;
		int charged_multiplicity = 0;
		int hadron_multiplicity = 0;
		int charged_hadron_multiplicity = 0;    // N_{ch}
		int Nch_absEta_lt_0_5 = 0;              // N_{ch}|_{|\eta|<=0.5}
		int Nch_absEta_lt_1_2 = 0;              // N_{ch}|_{|\eta|<=1.2}


		for (int i = 0; i < pythia.event.size(); ++i)
		{
			Particle & p = pythia.event[i];
			if ( p.isFinal() )
			{
				event_multiplicity++;
				if ( p.isCharged() )
					charged_multiplicity++;

				if ( p.isHadron() )
				{
					//count all final hadrons in multiplicity
					hadron_multiplicity++;
					if ( p.isCharged() )
					{
						charged_hadron_multiplicity++;
						const double abs_p_eta = abs(p.eta());
						if ( abs_p_eta < 1.2 )
						{
							Nch_absEta_lt_1_2++;
							if ( abs_p_eta < 0.5 )
								Nch_absEta_lt_0_5++;
						}
					}

					// i.e., is pion(+) or another HBT particle in the HBT_particle_IDs map
			 		if ( p.id() == 211 )
					{

						// save the final (shifted) particles, no matter what
						//particles_to_output.push_back( p );

						cout
							<< "results: "
							<< fixed
							<< setprecision(6)
							<< iEvent << "   " << i
							<< "   " << p.e()
							<< "   " << p.px()
							<< "   " << p.py()
							<< "   " << p.pz()
							<< scientific
							<< "   " << p.tProd()
							<< "   " << p.xProd()
							<< "   " << p.yProd()
							<< "   " << p.zProd()
							<< endl;

					}
				}
			}
		}
	}

	return 0;
}
