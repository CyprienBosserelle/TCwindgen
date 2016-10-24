#include <stdio.h>
#include <math.h>


#include <iostream>
#include <fstream>
#include <iomanip>


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cmath>

#include <sstream>
#include <iterator>
#include <netcdf.h>
#include <algorithm>
#include <vector>
#include <ctime>

#define pi 3.14159265f


class param {
public:
	double LonMin, LonMax, dlon, dlat, LatMin, LatMax, dt;
	tm datestart, dateend;
	std::string Trackfile, Outputncfile, SWANout, Delft3Dout;

};

class TCparam {
public:
	double TClon, TClat, cP, eP, dP, rMax, vFm, thetaFm, beta, rho;
	tm datetime;
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
std::string trim(const std::string& str, const std::string& whitespace);
template <class T> const T& max(const T& a, const T& b);
template <class T> const T& min(const T& a, const T& b);
double interptime(double next, double prev, double timenext, double time);
std::vector<TCparam> checkTCtrack(std::vector<TCparam> TCtrack);
void createSWANwindfile(std::string SWANfileOut, int nx, int ny, float *U, float *V);
void createSWANwindfileFPF(std::string SWANfileOut, int nx, int ny, float *U, float *V);
void writeSWANWindstep(std::string SWANfileOut, int nx, int ny, float *U, float *V);
void writeSWANWindstepFPF(std::string SWANfileOut, int nx, int ny, float *U, float *V);
void createD3DAtmfile(std::string D3dfileOut, int nx, int ny, double lonmin, double latmin, double dlon, double dlat, tm  datestart, float *P, float *U, float *V);
void writeD3DAtmstep(std::string D3dfileOut, int nx, int ny, tm datestart, double totaltime, float *P, float *U, float *V);
void createD3DAtmfileFPF(std::string D3dfileOut, int nx, int ny, double lonmin, double latmin, double dlon, double dlat, tm  datestart, float *P, float *U, float *V);
void writeD3DAtmstepFPF(std::string D3dfileOut, int nx, int ny, tm datestart, double totaltime, float *P, float *U, float *V);

std::vector<TCparam> readtrackfile(std::string trackfilename);
TCparam readtrackline(std::string line);