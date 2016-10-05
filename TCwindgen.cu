
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>


int main()
{
	double cP = 900; //current pressure
	double eP = 1013; //Env pressure

	double beta,beta1,beta2;
	double rho;

	// convert from hPa to Pa
	cP = cP * 100;
	eP = eP * 100;

	double dP; //Pressure difference

	double Vmax; // Max speed
	dP = eP - cP;


	//Vmax = MaxVelModels()
	//Willoughby Wind Speed
	/*Willoughby & Rahn (2004), Parametric Representation of the Primary
    Hurricane Vortex. Part I: Observations and Evaluation of the
    Holland (1980) Model.  Mon. Wea. Rev., 132, 3033-3048*/
	Vmax = 0.6252 * sqrt(dP);

	printf("Vmax= %f\n",Vmax);
	//Holland model
	/*.. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    Holland (1980), An Analytic Model of the Wind and Pressure Profiles
    in Hurricanes. Mon. Wea. Rev, 108, 1212-1218 Density of air is
    assumed to be 1.15 kg/m^3.  |beta| is assumed to be 1.3. Other values
    can be specified.  Gradient level wind (assumed maximum).*/
	beta = 1.3;
	rho = 1.15;
	Vmax = sqrt(beta * dP / (exp(1.0) * rho));
	printf("Vmax= %f\n", Vmax);

}
