#include <stdio.h>
#include <math.h>

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iterator>
#include <netcdf.h>
#include <algorithm>
#include <vector>

#define pi 3.14159265f


class param {
public:
	double LonMin, LonMax, dlon, dlat, LatMin, LatMax;
	std::string Trackfile, Outputncfile;

};

class TCparam {
public:
	double TClon, TClat, cP, eP, dP, rMax, vFm, thetaFm, beta, rho;
//double TClat, double TClon, double cP, double eP, double rMax, double vFm, double thetaFm, double beta, double rho
};



extern "C" void creatncfile(std::string ncfileout, int nx, int ny, float totaltime, float * xval, float * yval, float *R, float *V, float *Z);
extern "C" void writestep2nc(std::string ncfileout, int nx, int ny, float totaltime, float *R, float *V, float *Z);
std::string findparameter(std::string parameterstr, std::string line);
param readparamstr(std::string line, param grid);
std::vector<TCparam> readBSHfile(std::string BSHfilename);
TCparam readBSHline(std::string line);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);