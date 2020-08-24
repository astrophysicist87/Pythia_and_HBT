#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

/*typedef struct
{
	int eventID;		//which event did this particle come from
	int pdgID;
	int thermal_or_decay;

	double t, x, y, z;
	double E, px, py, pz;
	
} ParticleRecord;*/

int main(int argc, char ** argv)
{
	bool include_thermal = true;
	bool include_decay = false;

	int mult_min = 1, mult_max = 11;
	double KTmin = 0.0, KTmax = 0.2;
	double Kphimin = 0.0, Kphimax = 2.0*M_PI;
	double KZmin = -1.0, KZmax = 1.0;

	//vector<ParticleRecord> particles;
	ifstream inMultiplicities( argv[1] );
	int nEvents = 30000000;
	int nCols = 9;
	int columnToUse = 5;	//zero-indexed
	vector<int> multiplicities (nEvents);
	int count = 0;
	//for (int iEvent = 0; iEvent < nEvents; iEvent++)
	while (inMultiplicities.good())
	{
		int dummy, tmp;
		for (int iCol = 0; iCol < columnToUse; iCol++)
			inMultiplicities >> dummy;
		inMultiplicities >> multiplicities[count++];
		for (int iCol = columnToUse + 1; iCol < nCols; iCol++)
			inMultiplicities >> dummy;
	}
	inMultiplicities.close();

	cerr << "Read in " << count << " events." << endl;

	// read file(s)
	for (int iFile = 2; iFile < argc; iFile++)
	{
		cerr << "Reading in " << argv[iFile] << "..." << endl;
		ifstream infile( argv[iFile] );
		count = 0;
		//int nEventsPerFile = 100000;
		//for (int iEvent = 0; iEvent < nEventsPerFile; iEvent++)
		while (infile.good())
		{
			int eventID, multiplicity;
			infile >> eventID >> multiplicity;
			double evdummy;
			if ( multiplicities[eventID] < mult_min or multiplicities[eventID] > mult_max )
			{
				for (int iParticle = 0; iParticle < multiplicity; iParticle++)
				for (int iCol = 0; iCol < 10; iCol++)
					infile >> evdummy;
				continue;
			}
			if (count==0)
			{
				cerr << " --> " << argv[iFile] << " was good" << endl;
				cerr << " --> " << eventID << "   " << multiplicity << "   "
						<< multiplicities[eventID] << endl;
			}
			count++;
			for (int iParticle = 0; iParticle < multiplicity; iParticle++)
			{
				/*ParticleRecord particle;
				particle.eventID = eventID;
				infile >> particle.pdgID >> particle.thermal_or_decay
						>> particle.E >> particle.px >> particle.py >> particle.pz
						>> particle.t >> particle.x >> particle.y >> particle.z;
				particle.t *= 1.0e12;
				particle.x *= 1.0e12;
				particle.y *= 1.0e12;
				particle.z *= 1.0e12;*/
				double dummy;
				int thermal_or_decay;
				double E, px, py, pz;
				double t, x, y, z;
				infile >> dummy >> thermal_or_decay
						>> E >> px >> py >> pz >> t >> x >> y >> z;

				if ( thermal_or_decay==0 and !include_thermal ) continue;
				else if ( thermal_or_decay==1 and !include_decay ) continue;
				if ( px*px+py*py > KTmax*KTmax or px*px+py*py < KTmin*KTmin ) continue;
				if ( pz > KZmax or pz < KZmin ) continue;
				double this_Kphi = atan2(py, px);
				if ( this_Kphi < Kphimin or this_Kphi > Kphimax ) continue;

				t *= 1.0e12;
				x *= 1.0e12;
				y *= 1.0e12;
				z *= 1.0e12;

				// if we made it here, output particle (space-time) info
				cout << t << "   " << x << "   " << y << "   " << z << endl;
			}
		}

		infile.close();
	}

	return 0;
}
