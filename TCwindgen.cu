
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


// Global parameters;
float * Gridlon, *Gridlon_g; // Longitude vector is length of nx
float * Gridlat, *Gridlat_g; // Latitude vector isd length of ny
float * R, *R_g; // array of distance  is size of nx*ny;
float *lam, *lam_g; // array of bearing (forward azimuth) from TC center to each grid point
float * V, *V_g; // array for Wind velocity
float * Z, *Z_g; // array for TC vorticity
float *Uw, *Vw, *Uw_g, *Vw_g; // Array of U and V wind from the cyclone
int nx, ny; //grid dimension


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
		Ri = R[i];

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
			Zi = Ri * (Ri * 4.0f * aa + 3.0f * bb) + 2.0f * cc;
		}
		else
		{
			Vi =  (sqrtf((dP * beta / rho) * delta * edelta + (Ri * f / 2.0f)*(Ri * f / 2.0f)) - Ri *abs(f) / 2.0f);
			Zi = ((sqrtf((dP * beta / rho) * delta * edelta + (Ri * f / 2.0f)*(Ri * f / 2.0f))) / Ri - abs(f) + edelta * (2.0f * (beta * beta) * dP * (delta - 1.0f) * delta + rho * edelta * (f * Ri) *(f * Ri)) / (2.0f * rho * Ri * sqrtf(4.0f * (beta * dP / rho) * delta * edelta + (f * Ri) *(f * Ri))));
		}
		V[i] = Vi*f / abs(f);
		Z[i] = Zi*f / abs(f);
	}

}

__global__ void NewHollandWindProfile(int nx, int ny, float f, float rMax, float dP, float rho, float TClat, float *R, float *V, float *Z)
{
	//Holland et al. 2010.  In this version, the exponent is allowed to
	//vary linearly outside the radius of maximum wind.i.e.rather than
	//	take the sqare root, the exponent varies around 0.5.Currently
	//	this version does not have a corresponding vorticity profile set up
	//	in windVorticity, so it cannot be applied in some wind field modelling.
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	float Vi, Ri, Zi;
	float E, d2Vm, aa, bb, cc;
	float delta, edelta;

	float Bs, deltag, edeltag, rgterm, xn, xx;

	
	float rGale = 150.0; // Radius for gale force wind. This should be user defined

	if (ix < nx && iy < ny)
	{
		//
		Ri = R[i];
		Bs = (-0.000044f * powf(dP / 100.0f, 2.0f) + 0.01 * (dP / 100.0f) - 0.014f * abs(TClat) + 1.0);
		deltag = powf(rMax / rGale, Bs);
		edeltag = exp(-1.0f * deltag);
		rgterm = Bs * dP * deltag * edeltag / rho;
		xn = log(17.0f) / log(rgterm);
		xx = 0.5;

		if (Ri > rMax)
		{
			xx = (0.5 + (Ri - rMax) * (xn - 0.5) / (rGale - rMax));
		}

		delta = powf(rMax / Ri, Bs);
		edelta = exp(-delta);

		V[i] = (f / abs(f)) * pow((dP * Bs / rho) *delta * edelta, xx);
		Z[i] = 0.0f;// Warning dummy value
		
	}
}

__global__ void DoubleHollandWindProfile(int nx, int ny, float f, float vMax, float rMax, float dP, float cP, float rho, float beta, float *R, float *V, float *Z)
{
	//McConochie *et al*'s double Holland vortex model (based on Cardone *et
	//al*, 1994).This application is the Coral Sea adaptation of the
	//double vortex model(it can also be used for concentric eye - wall
	//configurations).
	//
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;
	float Vi, Ri, Zi;
	float E, d2Vm, aa, bb, cc;
	float delta, edelta;

	float beta1, beta2;

	float rMax2 = 150.0f;
	float rMax1 = rMax;
	float gradientV1, gradientV2;

	float chi, psi;

	float dp2,dp1,nu,mu,enu,emu;

	if (dP < 1500.0f)
	{
		dp2 = ((dP / 1500.0f) * (800.0f + (dP - 800.0f) / 2000.0f));
	}
	else
	{
		dp2 = 800.0f + (dP - 800.0f) / 2000.0f;
	}

	dp1 = dP - dp2;


	//Second derivative of the profile
	beta1 = beta;
	beta2 = 7.2f - cP / 16000.0f;

	E = exp(1.0f);
	
	nu = pow((rMax2 / rMax1), beta2);

	d2Vm = (-1.0f / (8.0f * pow(4.0f * beta1 * dp1 / (rho * E) + (4.0f * beta2 * dp2 / rho) * nu * exp(-nu) + powf(rMax1 * f, 2.0f), 1.5f))*
		(-(4.0f * (beta1 *beta1) * dp1 / (rho * rMax1 * E)) + (4.0f * (beta1 * beta1) * dp1 / (rho * rMax1 * E)) - (4 * (beta2 *beta2) * dp2 / rho) *
		(nu / rMax1) * exp(-nu) + (4.0f * (beta2 *beta2) * dp2 / rho) *((nu *nu) / rMax1) * exp(-nu) + 2.0f * rMax1 * f *f, 2.0f)
		+ 1.0f / (4.0f * sqrt((4 * beta1 * dp1 / (rho * E)) +
		(4.0f * beta2 * dp2 / rho) * nu * 2.0f +
		exp(-nu) + pow(rMax1 * f,2.0f)))
		* ((4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))
		+ (4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))
		- (12.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))
		- (4.0f * (beta1 *beta1) * dp1 / (rho * (rMax1 *rMax1) * E))
		+ (4.0f * (beta1 *beta1*beta1) * dp1 / (rho * (rMax1 *rMax1) * E))
		+ (4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu / (rMax1 *rMax1)) * exp(-nu)
		+ (4.0f * (beta2 *beta2) * dp2 / rho) *
		(nu / (rMax1 *rMax1)) * exp(-nu)
		- (12.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu *nu) / (rMax1 *rMax1) * exp(-nu)
		- (4.0f * (beta2 *beta2) * dp2 / rho) *
		(nu *nu) / (rMax1 *rMax1) * exp(-nu)
		+ (4.0f * (beta2 *beta2*beta2) * dp2 / rho) *
		(nu *nu*nu) / (rMax1 *rMax1) * exp(-nu)
		+ 2.0f * f *f));


	if (ix < nx && iy < ny)
	{
		//
		Ri = R[i];
		mu = powf(rMax / Ri, beta1);
		nu = powf(rMax2 / Ri, beta2);
		emu = exp(-mu);
		enu = exp(-nu);

		chi = beta1 * dp1 / rho;
		psi = beta2 * dp2 / rho;

		gradientV1 = (chi) * mu * emu;
		gradientV2 = (psi) * nu * enu;


		aa = (d2Vm / 2.0f - (-vMax / rMax) / rMax) / rMax;
		bb = (d2Vm - 6.0f * aa * rMax) / 2.0f;
		cc = -3.0f * aa * rMax * rMax - 2.0f * bb * rMax;

		Vi = (f / abs(f) * sqrt(gradientV1 + gradientV2 + (Ri *f / 2.0f) *(Ri *f / 2.0f)) - Ri * abs(f) / 2.0f);

		if (dP >= 1500.0f && Ri <= rMax)
		{
			Vi = (f / abs(f) * Ri * (Ri * (Ri * aa + bb) + cc));
		}

		V[i] = Vi;
		Z[i] = 0.0f;
			//(f / abs(f) * sqrtf(chi * delta * edelta + psi * nu * enu + (f * Ri / 2.0f) *(f * Ri / 2.0f)) / Ri -
			//abs(f) + (0.5f) *
			//(chi * ddelta * edelta * (1 - delta) +
			//psi * dgamma * egamma * (1 - gamma) +
			//R * self.f ** 2) /
			//np.sqrt(chi * delta * edelta + psi * gamma *
			//egamma + (self.f * R / 2) ** 2))


	}
}





__global__ void HubbertWindField(int nx, int ny, float rMax, float vFm, float thetaFm, float *R, float *lam, float *V, float *Uw, float *Vw)
{
	//
	//Hubbert, G.D., G.J.Holland, L.M.Leslie and M.J.Manton, 1991:
	//A Real - Time System for Forecasting Tropical Cyclone Storm Surges.
	//	*Weather and Forecasting*, **6 * *, 86 - 97

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

__global__ void McConochieWindField(int nx, int ny, float rMax, float vMax, float vFm, float thetaFm, float *R, float *lam, float *V, float *Uw, float *Vw)
{
	//
	//McConochie, J.D., T.A.Hardy and L.B.Mason, 2004:
	//Modelling tropical cyclone over - water wind and pressure fields.
	//	Ocean Engineering, 31, 1757 - 1782.
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	float Km = 0.70;
	float inflow = 25.0f;
	float Ri, Vi;
	float lami;
	float thetaMax = 0.0f;

	float thetaMaxAbsolute, asym, Vsf, phi, Ux, Vy,swrf;

	if (ix < nx && iy < ny)
	{
		//
		//V = self.velocity(R)
		Ri = R[i];
		lami = lam[i] * pi / 180.0f;
		Vi = V[i];
		if (Ri < 1.2f*rMax)
		{
			//
			inflow = 10.0f + 75.0f * (Ri / rMax - 1.0f);
		}
		if (Ri < rMax)
		{
			//
			inflow = 10.0f * Ri / rMax;
		}
		inflow = inflow*pi / 180.0f;

		thetaMaxAbsolute = thetaFm + thetaMax;
		phi = inflow - lami;

		asym = (0.5f * (1.0f + cosf(thetaMaxAbsolute - lami)) * vFm * (Vi / vMax));
		Vsf = Vi + asym;

		swrf = 0.81f;
		if (Vsf >= 6.0f);
		{
			swrf = 0.81f - (2.93f * (Vsf - 6.0f) / 1000.0f);
		}
		if (Vsf >= 19.5f)
		{
			swrf = 0.77f - (4.31f * (Vsf - 19.5f) / 1000.0f);
		}
		if (Vsf > 45.0f)
		{
			swrf = 0.66f;
		}
		Uw[i] = swrf * Vsf * sinf(phi);
		Vw[i] = swrf * Vsf * cosf(phi);
	}
}

__global__ void KepertWindField(int nx, int ny, float rMax, float vMax, float vFm, float thetaFm, float f, float *R, float *lam, float *V, float *Z, float *Uw, float *Vw)
{
	// Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
	//Tropical Cyclone Core.Part I : Linear Theory.J.Atmos.Sci., 58,
	//	2469 - 2484
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	float Ri, Vi, Zi;
	float lami;
	float K = 50.0f; //diffusivity
	float Cd = 0.002f; // Constant drag coeff

	float Vt,al,be,gam, albe;
	float chi, eta, psi;

	float A0r, A0i,u0s,v0s,Amr,Ami,ums,vms,Apr,Api,ups,vps;
	float us, vs, usf, vsf,phi;

	if (ix < nx && iy < ny)
	{
		//
		//V = self.velocity(R)
		Ri = R[i];
		lami = lam[i] * pi / 180.0f;
		Vi = V[i];
		Zi = Z[i];

		Vt = vFm;
		if (Ri>4.0f*rMax)
		{
			Vt = vFm * expf(-pow((Ri / rMax) - 4.0f, 2.0f));
		}

		al = ((2.0f * Vi / Ri) + f) / (2.0f * K);
		be = (f + Zi) / (2.0f * K);
		gam = (-1.0f * f/abs(f)) * Vi / (2.0f * K * Ri);

		albe = sqrtf(al / be);

		chi = (Cd / K) * Vi / sqrt(sqrt(al * be));
		eta = (Cd / K) * Vi / sqrt(sqrt(al * be) + abs(gam));
		psi = (Cd / K) * Vi / sqrt(abs(sqrt(al * be) - gam));

		// converted from complex number formula to this
		A0r = (-chi * Vi) / (2.0f * chi *chi + 3.0f* chi + 2.0f);
		A0i = (-chi * Vi * (1. + chi) / (2.0f * chi *chi + 3.0f* chi + 2.0f));

		//Symmetric surface wind component
		u0s = albe*A0r;
		v0s = A0i;

		// converted from complex number formula to this
		Amr = (-(psi * (1.0f + 2.0f * albe + (1.0f + albe) * eta)) * Vt / (albe * (2.0f * (1 + eta * psi) + 3.0f * psi + 3.0f * i * eta)));
		Ami = (-(psi * ( (1.0f + albe) * eta)) * Vt / (albe * (2.0f * (1.0f + eta * psi))));

		if (abs(gam) > sqrt(al*be))
		{
			Amr = (-(psi * (1.0f + 2.0f * albe +  (1.0f + albe) * eta) * Vt) / (albe * (2.0f + 3.0f * (eta + psi) + 2.0f * eta * psi)));
			Ami = (-(psi * ( (1.0f + albe) * eta) * Vt) / (albe * (- 2.0f + 2.0f  * eta * psi)));
		}

		//First asymmetric surface component
		ums = albe * Amr;
		vms = (Ami * exp(-(lami)));
		
		Apr = (-(eta * (1.0f - 2.0f * albe + (1.0f - albe) * psi)) * Vt / (albe * (2.0f * (1.0f + eta * psi) + 3.0f * eta )));
		Api = (-(eta * ((1.0f - albe) * psi)) * Vt / (albe * (2.0f * (1.0f + eta * psi) + 3.0f * psi)));
		
		if (abs(gam) > sqrt(al*be))
		{
			Apr = (-(eta * (1.0f - 2.0f * albe + (1.0f - albe) * psi) * Vt) / (albe * (2.0f + 3.0f * (eta + psi) + 2.0f * eta * psi)));
			Api = (-(eta * ((- 1.0f) * (1.0f - albe) * psi) * Vt) / (albe * ( 2.0f - 2.0f * eta * psi)));
		}

		//Second asymmetric surface component
		ups = albe * Apr;
		vps = Api * exp(lami);

		//Total surface wind in (moving coordinate system)
		us = u0s + ups + ums;
		vs = Vi + v0s + vps + vms;

		usf = us + Vt * cosf(lami - thetaFm);
		vsf = vs - Vt * sinf(lami - thetaFm);
		phi = atan2(usf, vsf);

		Uw[i] = (sqrtf(usf*usf + vsf *vsf) * sinf(phi - lami));
		Vw[i] = (sqrtf(usf*usf + vsf *vsf) * cosf(phi - lami));

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

int GenPUV(int Profile, int Field, int Vmaxmodel, double TClat, double TClon, double cP, double eP, double rMax, double vFm, double thetaFm, double beta, double rho)
{
	//Cyclone parameters
	//double TClat = -18.0;//Latitude of TC centre
	//double TClon = 178.0;//Longitude of TC centre

	//double cP = 900.0; //central pressure hpa
	//double eP = 1013.0; //Env pressure hpa
	//double rMax = 40.0; // Radius of maximum wind (km)
	//double vFm = 15.0; //Foward speed of the storm(m / s)
	//double thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);

	//Calculated parameters
	double dP; //Pressure difference
	double Vmax; // Max speed
	//double beta = 1.30;
	//double rho = 1.15;
	float Rearth = 6372797.560856f;
	float E, d2Vm, f;
	// convert from hPa to Pa
	cP = cP * 100.0;
	eP = eP * 100.0;


	dP = eP - cP;

	Vmax = Vmax_models(1, cP, eP, beta, rho);

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


	Rdist << <gridDim, blockDim, 0 >> >(nx, ny, Gridlon_g, Gridlat_g, TClon, TClat, R_g, lam_g);
	CUDA_CHECK(cudaThreadSynchronize());

	//JelesnianskiWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	//HollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, dP, rho, beta, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	//NewHollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, rMax, dP, rho, TClat, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	DoubleHollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, dP, cP, rho, beta, R_g, V_g, Z_g);
	CUDA_CHECK(cudaThreadSynchronize());
	//CUDA_CHECK(cudaMemcpy(R, R_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(V, V_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(lam, lam_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//printf("R[0]=%f\tR[nx*ny-1]=%f\n", V[0], V[nx*ny - 1]);


	//HubbertWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, (float)vFm, (float)thetaFm, R_g, lam_g, V_g, Uw_g, Vw_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	McConochieWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, Vmax, vFm, thetaFm, R_g, lam_g, V_g, Uw_g, Vw_g);
	CUDA_CHECK(cudaThreadSynchronize());

	//KepertWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, Vmax, vFm, thetaFm, f, R_g, lam_g, V_g, Z_g, Uw_g, Vw_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	return 0;

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


	//printf("Gridlon[0]=%f\tGridlon[nx-1]=%f\n", Gridlon[0], Gridlon[1]);
	//printf("Gridlat[0]=%f\tGridlat[ny-1]=%f\n", Gridlat[0], Gridlat[1]);

	// First iteration
	int iteration = 0;

	double dt = 600;// needs to be small enough to make a smooth model

	int Profilemodeltype, WindFieldmodeltype, Vmaxmodeltype; //not yet set. Default is Double Holland for the wind profile, McConochie for the wind field, and Holland for vmax

	double TClat = -18.0;//Latitude of TC centre
	double TClon = 178.0;//Longitude of TC centre

	double cP = 900.0; //central pressure hpa
	double eP = 1013.0; //Env pressure hpa
	double rMax = 40.0; // Radius of maximum wind (km)
	double vFm = 15.0; //Foward speed of the storm(m / s)
	double thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);
	double beta = 1.30;
	double rho = 1.15;


	// Below function modifies and call global parameters 
	int dummy = GenPUV(Profilemodeltype, WindFieldmodeltype, Vmaxmodeltype, TClat, TClon,cP, eP, rMax, vFm, thetaFm, beta, rho);

	CUDA_CHECK(cudaMemcpy(V, V_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Vw, Vw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Uw, Uw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));

	creatncfile("test.nc", nx, ny, 0.0f, Gridlon, Gridlat, V, Uw, Vw);

	 
	



}
