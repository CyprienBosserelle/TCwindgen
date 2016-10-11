
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
#include <netcdf.h>
#include <algorithm>

#define pi 3.14159265f


extern "C" void creatncfile(char outfile[], int nx, int ny, float totaltime, float * xval, float * yval, float *R, float *V, float *Z)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim;
	size_t nxx, nyy, nnpart;
	int  var_dimids[3], var_dimzb[2];

	int R_id, time_id, xx_id, yy_id, V_id, Z_id;
	
	nxx = nx;
	nyy = ny;
	//nnpart=npart;

	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { 1, ny, nx };
	static size_t zbstart[] = { 0, 0 }; // start at first value 
	static size_t zbcount[] = { ny, nx };
	//static size_t pstart[] = {0, 0}; // start at first value 
	// 	static size_t pcount[] = {1, npart};
	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value 
	static size_t xcount[] = { nx };
	
	static size_t ystart[] = { 0 }; // start at first value 
	static size_t ycount[] = { ny };


	
	//create the netcdf dataset
	status = nc_create(outfile, NC_NOCLOBBER, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "x", nxx, &xx_dim);
	status = nc_def_dim(ncid, "y", nyy, &yy_dim);
	//status = nc_def_dim(ncid, "npart",nnpart,&p_dim);
	status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	//int pdim[2];
	//pdim[0]=time_dim;
	//pdim[1]=p_dim;
	//define variables: Name, Type,...
	var_dimids[0] = time_dim;
	var_dimids[1] = yy_dim;
	var_dimids[2] = xx_dim;
	var_dimzb[0] = yy_dim;
	var_dimzb[1] = xx_dim;

	status = nc_def_var(ncid, "time", NC_FLOAT, 1, tdim, &time_id);
	status = nc_def_var(ncid, "x", NC_FLOAT, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "y", NC_FLOAT, 1, ydim, &yy_id);



	status = nc_def_var(ncid, "R", NC_FLOAT, 3, var_dimids, &R_id);
	status = nc_def_var(ncid, "V", NC_FLOAT, 3, var_dimids, &V_id);
	status = nc_def_var(ncid, "Z", NC_FLOAT, 3, var_dimids, &Z_id);

	//put attriute: assign attibute values
	//nc_put_att

	//End definitions: leave define mode
	status = nc_enddef(ncid);

	//Provide values for variables
	status = nc_put_var1_float(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_float(ncid, xx_id, xstart, xcount, xval);
	status = nc_put_vara_float(ncid, yy_id, ystart, ycount, yval);

	
	status = nc_put_vara_float(ncid, R_id, start, count, R);
	status = nc_put_vara_float(ncid, V_id, start, count, V);
	status = nc_put_vara_float(ncid, Z_id, start, count, Z);
	// U, V, P are also needed

	//close and save new file
	status = nc_close(ncid);
}

__global__ void Rdist(int nx, int ny, float *Gridlon, float *Gridlat, double TClon, double TClat, float *R, float *lam)
{
	//
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	if (ix < nx && iy < ny)
	{
		float Rearth = 6372797.560856f;
		float dlat = (Gridlat[iy] - TClat)*pi / 180.0f;
		float lat1 = TClat * pi / 180.0f;
		float lat2 = Gridlat[iy] * pi / 180.0f;
		
		float dlon = (Gridlon[ix] - TClon)*pi / 180.0f;
		float a = sinf(dlat / 2.0f)*sinf(dlat / 2.0f) + cosf(lat1)*cosf(lat2)*sinf(dlon / 2.0f)*sinf(dlon / 2.0f);
		float c = 2.0f * atan2f(sqrtf(a), sqrtf(1.0f - a));
		R[i] = c*Rearth/100.0f;//convert to km

		float x = sinf(dlon)*cosf(lat2);
		float y = cosf(lat1)*sinf(lat2) - sinf(lat1)*cosf(lat2)*cosf(dlon);
		lam[i] = atan2f(y, x)*180.0f / pi;

	}
}

__global__ void lonGrid(int nx, float dx, float LonMin, float *Gridlon)
{
	//
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < nx)
	{
		Gridlon[i] = LonMin+dx*i;
	}
}

__global__ void JelesnianskiWindProfile(int nx, int ny, float f, float vMax, float rMax, float *R, float *V,float *Z)
{
	//
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	float Vi, Ri,Zi, sf;
	sf = (f / abs(f));

	if (ix < nx && iy < ny)
	{
		//
		Ri = R[i];
		Vi = 2.0f * vMax * rMax * Ri / (rMax *rMax + Ri * Ri) * sf;
		V[i] = Vi;

		Zi = (sf * 2.0f * vMax * rMax / (rMax *rMax + Ri * Ri) + sf * 2.0f * vMax * rMax * (rMax *rMax - Ri * Ri) /((rMax *rMax + Ri * Ri) * (rMax *rMax + Ri * Ri)));
		Z[i] = Zi;

	}
}
__global__ void HollandWindProfile(int nx, int ny, float f, float vMax, float rMax, float dP, float rho, float beta, float *R, float *V, float *Z)
{
	//
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	float Vi, Ri, Zi;
	float E, d2Vm,aa,bb,cc; 
	float delta, edelta;
	if (ix < nx && iy < ny)
	{
		//

		E = expf(1.0f);
		d2Vm = ((beta * dP * (-4.0f * beta *beta *beta * dP / rho - (-2.0f + beta *beta) * E * (f * rMax) *(f * rMax))) / (E * rho * sqrt((4.0f * beta * dP) / (E * rho) + (f * rMax) *(f * rMax)) * (4.0f * beta * dP * rMax *rMax / rho + E * (f * rMax *rMax) *(f * rMax *rMax))));

		aa = ((d2Vm / 2.0f - (-vMax / rMax) / rMax) / rMax);

		bb = (d2Vm - 6.0f * aa * rMax) / 2.0f;

		cc = -3.0f * aa * rMax *rMax - 2.0f * bb * rMax;

		delta = powf(rMax / Ri, beta);
		edelta = expf(-delta);
		if (Ri <= rMax)
		{
			Vi = (Ri * (Ri * (Ri * aa + bb) + cc));
		}
		else
		{
			Vi =  (sqrtf((dP * beta / rho)* delta * edelta + (Ri * f / 2.0f)*(Ri * f / 2.0f)) - Ri *abs(f) / 2.0f);
			Zi = ((sqrt((dP * beta / rho) * delta * edelta + (Ri * f / 2..0f) ** 2)) / R -
				np.abs(self.f) + edelta *
				(2 * (beta ** 2) * self.dP * (delta - 1) * delta +
				self.rho * edelta * (self.f * R) ** 2) /
				(2 * self.rho * R *
				np.sqrt(4 * (beta * self.dP / self.rho) * delta * edelta
				+ (self.f * R) ** 2)))
		}
		V[i] = Vi*f / abs(f);
	}

}

__global__ void HubbertWindField(int nx, int ny, float rMax, float vFm, float thetaFm, float *R, float *lam, float *V, float *Uw, float *Vw)
{
	//
	/*lam: Direction(geographic bearing, positive clockwise)
	from storm centre to the grid.
	: type  lam : : class :`numpy.ndarray`
	:param float vFm : Foward speed of the storm(m / s).
	: param float thetaFm : Forward direction of the storm(geographic
	bearing, positive clockwise).
	: param float thetaMax : Bearing of the location of the maximum
	wind speed, relative to the direction of
	motion.*/
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	float Km = 0.70;
	float inflow = 25.0f;
	float Ri,Vi;
	float lami;
	float thetaMax = 0.0f;

	float thetaMaxAbsolute, asym, Vsf, phi, Ux, Vy;

	if (ix < nx && iy < ny)
	{
		//V = self.velocity(R)
		Ri = R[i];
		lami = lam[i]*pi/180.0f;
		Vi = V[i];

		if (Ri < rMax)
		{
			inflow = 0;
		}

		inflow = inflow * pi / 180.0f;

		thetaMaxAbsolute = thetaFm + thetaMax;
		asym = vFm * cosf(thetaMaxAbsolute - lami + pi);
		Vsf = Km * Vi +asym;
		phi = inflow - lami;
		Uw[i] = Vsf *sinf(phi);
		Vw[i] = Vsf *cosf(phi);


	}
}

void CUDA_CHECK(cudaError CUDerr)
{


	if( cudaSuccess != CUDerr) {

		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \

		__FILE__, __LINE__, cudaGetErrorString( CUDerr) );

		exit(EXIT_FAILURE);

	}
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

int main(int argc, char **argv)
{
	// Grid parameters
	double LonMin = 177.0;
	double LonMax = 180.0;

	double dlon = 0.001;
	double dlat = dlon;
	
	double LatMin = -19.0;
	double LatMax = -17.0;
	
	
	// Generate grid parameters;
	float * Gridlon, * Gridlon_g; // Longitude vector is length of nx
	float * Gridlat, * Gridlat_g; // Latitude vector isd length of ny
	float * R, * R_g; // array of distance  is size of nx*ny;
	float *lam, *lam_g; // array of bearing (forward azimuth) from TC center to each grid point
	float * V, *V_g; // array for Wind velocity
	float * Z, *Z_g; // array for TC vorticity
	float *Uw, *Vw, *Uw_g, *Vw_g; // Array of U and V wind from the cyclone
	int nx, ny; //grid dimension


	nx = ceil((LonMax - LonMin) / dlon); // in case not an exact match then LonMax is extended
	ny = ceil((LatMax - LatMin) / dlat);

	printf("nx=%i; ny=%i\n", nx, ny);
	// Allocate on the CPU
	Gridlon = (float *)malloc(nx*sizeof(float));
	Gridlat = (float *)malloc(ny*sizeof(float));

	//Allocate on the GPU
	CUDA_CHECK(cudaMalloc((void **)&Gridlon_g, nx*sizeof(float)));
	CUDA_CHECK(cudaMalloc((void **)&Gridlat_g, ny*sizeof(float)));

	//Following should be quicker on GPU too
	/*for (int i = 0; i < nx; i++)
	{
		Gridlon[i] = LonMin + i*dlon;
	}

	for (int j = 0; j < ny; j++)
	{
		Gridlat[j] = LatMin + j*dlat;
	}*/
	dim3 blockDim1D(32, 1, 1);// This means that the grid has to be a factor of 16 on both x and y
	dim3 gridDim1Dlon(ceil((nx*1.0f) / blockDim1D.x), 1, 1);
	dim3 gridDim1Dlat(ceil((ny*1.0f) / blockDim1D.x), 1, 1);
	lonGrid <<<gridDim1Dlon, blockDim1D, 0 >>>(nx,(float) dlon, (float) LonMin, Gridlon_g);
	lonGrid <<<gridDim1Dlat, blockDim1D, 0 >>>(ny, dlat, LatMin, Gridlat_g);

	CUDA_CHECK(cudaMemcpy(Gridlon, Gridlon_g, nx*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Gridlat, Gridlat_g, ny*sizeof(float), cudaMemcpyDeviceToHost));

	// Results parameters 
	R = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&R_g, nx*ny*sizeof(float)));
	lam = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&lam_g, nx*ny*sizeof(float)));
	V = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&V_g, nx*ny*sizeof(float)));
	Z = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&Z_g, nx*ny*sizeof(float)));
	Uw = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&Uw_g, nx*ny*sizeof(float)));
	Vw = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&Vw_g, nx*ny*sizeof(float)));


	printf("Gridlon[0]=%f\tGridlon[nx-1]=%f\n", Gridlon[0], Gridlon[1]);
	//printf("Gridlat[0]=%f\tGridlat[ny-1]=%f\n", Gridlat[0], Gridlat[1]);

	//Cyclone parameters
	double TClat=-18.0;//Latitude of TC centre
	double TClon=178.0;//Longitude of TC centre

	double cP = 900.0; //central pressure hpa
	double eP = 1013.0; //Env pressure hpa
	double rMax = 40.0; // Radius of maximum wind (km)
	double vFm = 15.0; //Foward speed of the storm(m / s)
	double thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);

	//Calculated parameters
	double dP; //Pressure difference
	double Vmax; // Max speed
	double beta=1.30;
	double rho=1.15;
	float Rearth = 6372797.560856f;
	float E, d2Vm,f;
	// convert from hPa to Pa
	cP = cP * 100.0;
	eP = eP * 100.0;

	
	dP = eP - cP;
		
	Vmax = Vmax_models(0, cP, eP, beta, rho);

	double TClatrad = TClat*pi / 180.0;
	double wearth = pi*(1.0 / 24.0) / 1800.0;
	f = 2.0f*wearth*sin(TClatrad);

	// Calculate R from the present cyclone position
	//calc distance between each i  and i and the TC using haversine formula
	/*for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			float dlat = (Gridlat[j] - TClat)*pi / 180.0f;
			float lat1 = TClat * pi / 180.0f;
			float lat2 = Gridlat[j] * pi / 180.0f;
			float dlon = (Gridlon[i] - TClon)*pi / 180.0f;
			float a = sinf(dlat / 2.0f)*sinf(dlat / 2.0f) + cosf(lat1)*cosf(lat2)*sinf(dlon / 2.0f)*sinf(dlon / 2.0f);
			float c = 2.0f * atan2f(sqrtf(a), sqrtf(1 - a));
			R[i + j*nx] = c*Rearth;
		}
	}*/

	//printf("R[0]=%f\tR[nx*ny-1]=%f\n", R[0], R[nx*ny - 1]);


	dim3 blockDim(16, 16, 1);// This means that the grid has to be a factor of 16 on both x and y
	dim3 gridDim(ceil((nx*1.0f) / blockDim.x), ceil((ny*1.0f) / blockDim.y), 1);


	Rdist <<<gridDim, blockDim, 0 >>>(nx, ny, Gridlon_g, Gridlat_g, TClon, TClat, R_g, lam_g);
	CUDA_CHECK(cudaThreadSynchronize());

	JelesnianskiWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, R_g, V_g, Z_g);
	CUDA_CHECK(cudaThreadSynchronize());

	//CUDA_CHECK(cudaMemcpy(R, R_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(V, V_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(lam, lam_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//printf("R[0]=%f\tR[nx*ny-1]=%f\n", V[0], V[nx*ny - 1]);


	HubbertWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, (float)vFm, (float)thetaFm, R_g, lam_g, V_g, Uw_g, Vw_g);
	CUDA_CHECK(cudaThreadSynchronize());


	CUDA_CHECK(cudaMemcpy(V, V_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Vw, Vw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Uw, Uw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));

	creatncfile("test.nc", nx, ny, 0.0f, Gridlon, Gridlat, V, Uw, Vw);
	//Calculate velocity and Vorticity
	/*
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
	*/

	//delta = pow(rMax / R, beta);
	//edelta = exp(-1.0*delta)


}
