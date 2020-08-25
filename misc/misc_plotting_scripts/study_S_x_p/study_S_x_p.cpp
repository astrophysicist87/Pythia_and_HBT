#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

//using namespace std;

int main(int argc, char ** argv)
{
	bool include_thermal = true;
	bool include_decay = false;

	//int mult_min = 1, mult_max = 11;
	//double KTmin = 0.0, KTmax = 0.2;
	double Kphimin = -M_PI, Kphimax = M_PI;
	double KZmin = -1.0, KZmax = 1.0;

	int mult_min = std::stoi( argv[1] );
	int mult_max = std::stoi( argv[2] );
	double KTmin = std::stod( argv[3] );
	double KTmax = std::stod( argv[4] );

	std::ifstream inMultiplicities( argv[5] );
	int nEvents = 60000000;
	int nCols = 8;
	int columnToUse = 5;	//zero-indexed
	std::vector<int> multiplicities (nEvents);
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

	std::cerr << "Read in " << count << " events." << std::endl;

	// read file(s)
	for (int iFile = 6; iFile < argc; iFile++)
	{
		std::cerr << "Reading in " << argv[iFile] << "..." << std::endl;
		std::ifstream infile( argv[iFile] );
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
				std::cerr << " --> " << argv[iFile] << " was good" << std::endl;
				std::cerr << " --> " << eventID << "   " << multiplicity << "   "
						<< multiplicities[eventID] << std::endl;
			}
			count++;
			for (int iParticle = 0; iParticle < multiplicity; iParticle++)
			{
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
				//if ( this_Kphi < Kphimin or this_Kphi > Kphimax ) continue;
				if ( this_Kphi < 0.0 ) this_Kphi += 2.0*M_PI;

				t *= 1.0e12;
				x *= 1.0e12;
				y *= 1.0e12;
				z *= 1.0e12;

				// if we made it here, output particle (space-time) info
				//std::cout << t << "   " << x << "   " << y << "   " << z << std::endl;
				printf("%f   %f   %f   %f   %f\n", this_Kphi, t, x, y, z);
			}
		}

		infile.close();
	}

	return 0;
}
