
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "TCwindgen.h"


// Global parameters;
float * Gridlon, *Gridlon_g; // Longitude vector is length of nx
float * Gridlat, *Gridlat_g; // Latitude vector isd length of ny
float * R, *R_g; // array of distance  is size of nx*ny;
float *lam, *lam_g; // array of bearing (forward azimuth) from TC center to each grid point
float * V, *V_g; // array for Wind velocity
float * Z, *Z_g; // array for TC vorticity
float *Uw, *Vw, *Uw_g, *Vw_g, *P, *P_g; // Array of U and V wind and P pressure from the cyclone
int nx, ny; //grid dimension





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
__global__ void HollandPressureProfile(int nx, int ny, float rMax, float dP, float cP, float beta, float *R, float *P)
{
	//Holland pressure profile
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	float Ri;

	if (ix < nx && iy < ny)
	{
		//
		Ri = R[i];

		P[i] = cP + dP*exp(-1.0f*pow(rMax / Ri, beta));
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


__global__ void DoubleHollandPressureProfile(int nx, int ny,  float rMax, float dP, float cP,  float beta, float *R, float *P)
{
	//Holland pressure profile
	unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int i = ix + iy*nx;

	float Ri;
	float dp1,dp2;
	float beta1, beta2;
	float nu, mu, enu, emu;

	float rMax2 = 150.0f;


	if (dP < 1500.0f)
	{
		dp2 = (dP / 1500.0f)*(800.0f + (dP - 800.0f) / 2000.0f);
	}
	else
	{
		dp2 = 800.0f + (dP - 800.0f) / 2000.0f;
	}

	dp1 = dP - dp2;

	beta1 = beta;
	//beta1 = 7.3f - cP / 16000.0f;
	beta2 = 7.2 - cP / 16000.0f;

	if (ix < nx && iy < ny)
	{
		//
		Ri = R[i];
		mu = powf(rMax / Ri,beta1);
		nu = powf(rMax2 / Ri,beta2);
		emu = exp(-mu);
		enu = exp(-nu);
		P[i] = cP + dp1*emu + dp2*enu;
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

int GenPUV(int Profile, int Field, int Vmaxmodel, TCparam TCpara)
{
	//Wrapper for TC generation of Pressure, and Wind feilds
	
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
	TCpara.cP = TCpara.cP * 100.0;
	TCpara.eP = TCpara.eP * 100.0;


	TCpara.dP = TCpara.eP - TCpara.cP;

	Vmax = Vmax_models(1, TCpara.cP, TCpara.eP, TCpara.beta, TCpara.rho);

	double TClatrad = TCpara.TClat*pi / 180.0;
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


	Rdist << <gridDim, blockDim, 0 >> >(nx, ny, Gridlon_g, Gridlat_g, TCpara.TClon, TCpara.TClat, R_g, lam_g);
	CUDA_CHECK(cudaThreadSynchronize());

	//JelesnianskiWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	//HollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, rMax, dP, rho, beta, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	//NewHollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, rMax, dP, rho, TClat, R_g, V_g, Z_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	DoubleHollandWindProfile << <gridDim, blockDim, 0 >> >(nx, ny, f, Vmax, TCpara.rMax, TCpara.dP, TCpara.cP, TCpara.rho, TCpara.beta, R_g, V_g, Z_g);
	CUDA_CHECK(cudaThreadSynchronize());

	DoubleHollandPressureProfile << <gridDim, blockDim, 0 >> >(nx, ny, TCpara.rMax, TCpara.dP, TCpara.cP, TCpara.beta, R_g, P_g);
	CUDA_CHECK(cudaThreadSynchronize());
	//CUDA_CHECK(cudaMemcpy(R, R_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(V, V_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//CUDA_CHECK(cudaMemcpy(lam, lam_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	//printf("R[0]=%f\tR[nx*ny-1]=%f\n", V[0], V[nx*ny - 1]);


	//HubbertWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, (float)vFm, (float)thetaFm, R_g, lam_g, V_g, Uw_g, Vw_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	McConochieWindField << <gridDim, blockDim, 0 >> >(nx, ny, TCpara.rMax, Vmax, TCpara.vFm, TCpara.thetaFm, R_g, lam_g, V_g, Uw_g, Vw_g);
	CUDA_CHECK(cudaThreadSynchronize());

	//KepertWindField << <gridDim, blockDim, 0 >> >(nx, ny, rMax, Vmax, vFm, thetaFm, f, R_g, lam_g, V_g, Z_g, Uw_g, Vw_g);
	//CUDA_CHECK(cudaThreadSynchronize());

	return 0;

}


int main(int argc, char **argv)
{
	param grid;
	// initialise parameters
	grid.LonMin = 177.0;
	grid.LonMax = 180.0;

	grid.dlon = 0.005;
	grid.dlat = grid.dlon;
	
	grid.LatMin = -19.0;
	grid.LatMax = -17.0;
	
	grid.Trackfile = "";
	grid.Outputncfile = "";


	
	std::ifstream fs("TC_param.txt");

	if (fs.fail()){
		std::cerr << "TC_param.txt file could not be opened" << std::endl;
		exit(1);
	}

	std::string line;
	while (std::getline(fs, line))
	{		
		//std::cout << line << std::endl;
		
		//Get param or skip empty lines
		if (!line.empty())
		{
			grid = readparamstr(line, grid);
			//std::cout << line << std::endl;
		}

	}
	fs.close();
	
	
	//std::cout.precision(7);
	//std::cout << "Class test: "<< std::fixed << grid.LonMin << std::endl;

	


	nx = ceil((grid.LonMax - grid.LonMin) / grid.dlon); // in case not an exact match then LonMax is extended
	ny = ceil((grid.LatMax - grid.LatMin) / grid.dlat);

	std::cout << "nx=" << nx << " ny="<< ny << std::endl;
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
	lonGrid <<<gridDim1Dlon, blockDim1D, 0 >>>(nx,(float) grid.dlon, (float) grid.LonMin, Gridlon_g);
	lonGrid <<<gridDim1Dlat, blockDim1D, 0 >>>(ny, grid.dlat, grid.LatMin, Gridlat_g);

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
	P = (float *)malloc(nx*ny*sizeof(float));
	CUDA_CHECK(cudaMalloc((void **)&P_g, nx*ny*sizeof(float)));



	//printf("Gridlon[0]=%f\tGridlon[nx-1]=%f\n", Gridlon[0], Gridlon[1]);
	//printf("Gridlat[0]=%f\tGridlat[ny-1]=%f\n", Gridlat[0], Gridlat[1]);

	// First iteration
	int iteration = 0;

	double dt = 600;// needs to be small enough to make a smooth model

	int Profilemodeltype = 1;//not yet implemented. Default is Double Holland 
	int WindFieldmodeltype = 1;//not yet  implemented. Default is McConochie for the wind field
	int Vmaxmodeltype = 1;; //not yet  implemented. Default is Holland 


	TCparam TCinit;


	std::vector<TCparam> TCparamlist = readBSHfile(grid.Trackfile);

	/*
	TCinit.TClat = -18.0;//Latitude of TC centre
	TCinit.TClon = 178.0;//Longitude of TC centre

	TCinit.cP = 900.0; //central pressure hpa
	TCinit.eP = 1013.0; //Env pressure hpa
	TCinit.rMax = 40.0; // Radius of maximum wind (km)
	TCinit.vFm = 15.0; //Foward speed of the storm(m / s)
	TCinit.thetaFm = 180.0; //Forward direction of the storm(geographic bearing, positive clockwise);
	TCinit.beta = 1.30;
	TCinit.rho = 1.15;
	*/

	TCinit = TCparamlist[0];

	std::cout << TCinit.TClat << std::endl;

	int dummy;

	// Below function modifies and call global parameters 
	dummy = GenPUV(Profilemodeltype, WindFieldmodeltype, Vmaxmodeltype, TCinit);

	CUDA_CHECK(cudaMemcpy(P, P_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Vw, Vw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(Uw, Uw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));

	if (!grid.Outputncfile.empty())// empty string means no netcdf output
	{
		creatncfile(grid.Outputncfile, nx, ny, 0.0f, Gridlon, Gridlat, P, Uw, Vw);
				
	}

	for (int i = 0; i < 20; i++)
	{
		//dummy main loop
		TCinit.TClat = TCinit.TClat + 0.01;
		dummy = GenPUV(Profilemodeltype, WindFieldmodeltype, Vmaxmodeltype, TCinit);

		CUDA_CHECK(cudaMemcpy(P, P_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(Vw, Vw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(Uw, Uw_g, nx*ny*sizeof(float), cudaMemcpyDeviceToHost));

		if (!grid.Outputncfile.empty())// empty string means no netcdf output
		{
			writestep2nc(grid.Outputncfile, nx, ny, (i + 1)*600.0f, P, Uw, Vw);
		}
	}
	
	
}
