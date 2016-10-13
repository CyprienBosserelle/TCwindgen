#include <stdio.h>
#include <math.h>

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cmath>
#include <fstream>
#include <netcdf.h>
#include <algorithm>

#define pi 3.14159265f


class param {
public:
	double LonMin, LonMax, dlon, dlat, LatMin, LatMax;
	std::string Trackfile, Outputncfile;

};


extern "C" void creatncfile(std::string ncfileout, int nx, int ny, float totaltime, float * xval, float * yval, float *R, float *V, float *Z);
extern "C" void writestep2nc(std::string ncfileout, int nx, int ny, float totaltime, float *R, float *V, float *Z);
std::string findparameter(std::string parameterstr, std::string line);
param readparamstr(std::string line, param grid);