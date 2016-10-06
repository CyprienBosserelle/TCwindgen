
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <fstream>
//#include <netcdf.h>
#include <algorithm>

#define pi 3.14159265f


__global__ void Rdist(int nx, int ny, float *Gridlon, float *Gridlat, double TClon, double TClat, float *R)
{
	//
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	float Rearth = 6372797.560856f;
	float dlat = (Gridlat[iy] - TClat)*pi / 180.0f;
	float lat1 = TClat * pi / 180.0f;
	float lat2 = Gridlat[iy] * pi / 180.0f;
	float dlon = (Gridlon[ix] - TClon)*pi / 180.0f;
	float a = sin(dlat / 2)*sin(dlat / 2) + cos(lat1)*cos(lat2)*sin(dlon / 2)*sin(dlon / 2);
	float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));
	R[i] = c*Rearth;
}


double Vmax_models(int model, double cP, double eP, double beta, double rho)
{
	//Calculate max wind speed in model
	// model is the type of model used :
	//0: Willoughby & Rahn(2004)
	//1: Holland (1980)
	//2: Atkinson and Holliday(1977)
	double dP = eP - cP;

	switch (model) {
	case 0: return 0.6252 * sqrt(dP);
			break;
			//Willoughby Wind Speed
			/*Willoughby & Rahn (2004), Parametric Representation of the Primary
			Hurricane Vortex. Part I: Observations and Evaluation of the
			Holland (1980) Model.  Mon. Wea. Rev., 132, 3033-3048*/
					
	case 1: return sqrt(beta * dP / (exp(1.0) * rho));
			break;
		
			//Holland model
			/*.. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

			Holland (1980), An Analytic Model of the Wind and Pressure Profiles
			in Hurricanes. Mon. Wea. Rev, 108, 1212-1218 Density of air is
			assumed to be 1.15 kg/m^3.  |beta| is assumed to be 1.3. Other values
			can be specified.  Gradient level wind (assumed maximum).*/
			//beta = 1.3;
			//rho = 1.15;
			
	case 2 : return 3.04 * pow(1010.0 - cP / 100.0, 0.644);
			 break;
		
			//Atkinson and Holliday(1977), *Tropical Cyclone Minimum SeaLevel Pressure Maximum Sustained Wind Relationship for
			//the Western North Pacific*.Mon.Wea.Rev., **105 * *, 421 - 427
			//Maximum 10m, 1 - minute wind speed.Uses ``pEnv`` as 1010 hPa.

			
		

	default : return sqrt(beta * dP / (exp(1.0) * rho));
				  break;
	}
}

int main()
{
	// Grid parameters
	double LonMin = 177.0;
	double LonMax = 180.0;

	double dlon = 0.001;
	double dlat = dlon;
	
	double LatMin = -19.0;
	double LatMax = -17.0;
	
	
	// Generate grid parameters;
	float * Gridlon; // Longitude vector is length of nx
	float * Gridlat; // Latitude vector isd length of ny
	float * R; // array of distance  is size of nx*ny;
	int nx, ny; //grid dimension


	nx = ceil((LonMin - LonMax) / dlon); // in case not an exact match then LonMax is extended
	ny = ceil((LatMin - LatMax) / dlat);

	// Allocate on the CPU
	Gridlon = (float *)malloc(nx*sizeof(float));
	Gridlat = (float *)malloc(ny*sizeof(float));

	for (int i = 0; i < nx; i++)
	{
		Gridlon[i] = LonMin + i*dlon;
	}

	for (int j = 0; j < ny; j++)
	{
		Gridlat[j] = LatMin + j*dlat;
	}

	// Results parameters 
	R = (float *)malloc(nx*ny*sizeof(float));
	


	//Cyclone parameters
	double TClat=-18.0;//Latitude of TC centre
	double TClon=178.0;//Longitude of TC centre

	double cP = 900; //central pressure hpa
	double eP = 1013; //Env pressure hpa
	double rMax = 40; // Radius of maximum wind (km)

	//Calculated parameters
	double dP; //Pressure difference
	double Vmax; // Max speed
	double beta=1.3;
	double rho=1.15;
	float Rearth = 6372797.560856f;
	float E, d2Vm,f;
	// convert from hPa to Pa
	cP = cP * 100;
	eP = eP * 100;

	
	dP = eP - cP;
		
	Vmax = Vmax_models(0, cP, eP, beta, rho);

	double TClatrad = TClat*pi / 180.0;
	double wearth = pi*(1.0 / 24.0) / 1800.0;
	f = 2.0f*wearth*sin(TClatrad);

	// Calculate R from the present cyclone position
	//calc distance between each i  and i and the TC using haversine formula
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			float dlat = (Gridlat[j] - TClat)*pi / 180.0f;
			float lat1 = TClat * pi / 180.0f;
			float lat2 = Gridlat[j] * pi / 180.0f;
			float dlon = (Gridlon[i] - TClon)*pi / 180.0f;
			float a = sin(dlat / 2)*sin(dlat / 2) + cos(lat1)*cos(lat2)*sin(dlon / 2)*sin(dlon / 2);
			float c = 2 * atan2f(sqrtf(a), sqrtf(1 - a));
			R[i + j*nx] = c*Rearth;
		}
	}

	//Calculate velocity and Vorticity

	E = exp(1);
	d2Vm = ((beta * dP * (-4 * pow(beta, 3) * dP / rho -
		(-2 + beta * beta) * E * pow(f * rMax, 2))) /
		(E * rho * sqrt((4 * beta * dP) / (E * rho)
		+ pow(f * rMax, 2)) * (4 * beta * dP * rMax *rMax / rho
		+ E * pow(f * rMax * rMax, 2))));

	double aa, bb, cc;
	aa = ((d2Vm / 2.0 - (-1.0*Vmax / rMax) / rMax) / rMax);
	bb = (d2Vm - 6 * aa * rMax) / 2.0;
	cc = -3 * aa * rMax * rMax - 2 * bb * rMax;


	//delta = pow(rMax / R, beta);
	//edelta = exp(-1.0*delta)


}
