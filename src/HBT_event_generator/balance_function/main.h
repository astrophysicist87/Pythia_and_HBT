#ifndef MAIN_H
#define MAIN_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>


using namespace std;

//this is just to give this file a reason to exist for the moment...
const double plumbergtest = 0.;

// function to read in catalogue of event files and return number of files to read
int read_file_catalogue(string catalogue_name, vector<string> & allLines)
{
	ifstream catalogue_in(catalogue_name.c_str());

	string line;
	while (getline(catalogue_in, line))
		allLines.push_back(line);

	return ( allLines.size() );
}


#endif
