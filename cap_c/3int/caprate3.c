#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cubature.h"

#define VERBOSE 0

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

//GLOBAL CONSTANTS
#define ESCAPE_VEL sqrt((2 * 6.67408E-11 * 1.4 * 2E30) / (10E3))	//NS escape velocity in m/s. Need to make into function of radius
#define TEMP (1E3 * 1E-9) / (1.16E4)	//NS temp in GeV
#define NM 0.939	//neutron mass in GeV
#define FERMI_ENERGY 0.085 // in GeV (will need to make a variable)
#define FERMI_VEL (2 * FERMI_ENERGY) / NM
#define SOL 299792458.0 // speed of light in vacuum m/s
#define VELDISP 270e3 //DM velocity dispersion for MB dist
#define NSVEL 200e3 // NS velocity in galactic frame

//GENERAL FUNCTIONS
double mu(double dm){
	/* RATIO OF DARK MATTER TO NEUTRON MASS */

	return dm / NM;
}

double mu_plus(double dm){

	return ( 1 + mu(dm)) / 2.0;
}

//TRIG FUNCTIONS

double cosine(double s, double t, double v){
	/* v = escape_vel for initial */

	double num = v * v - s * s - t * t;
	double den = 2 * s * t;

	return num / den;
}

//STEP FUNCTIONS
double step(double f){
	/* a step function to constrain integration region */
	if ( fabs(f) > 1){
		return 0;
	}
	else{
		return 1;
	}
}

double heaviside_product(double s, double t, double v){
	/* product of relevant step functions */
	return step(cosine(s, t, ESCAPE_VEL)) * step(cosine(s, t, v));
}

//ENERGIES


double velsqr(double s, double t, double v, double dm){
	/* square of velocity: v = ESCAPE_VEL for initial particle */
	double a = 2 * mu(dm) * mu_plus(dm) * t * t + 2 * mu_plus(dm) * s*s- mu(dm) * v * v;
	return a / (SOL*SOL);
}

double energy(double s, double t, double vel, double dm){
	return 0.5 * NM * velsqr(s, t, vel, dm);
}


//FERMI-DIRAC DISTRIBUTION FUNCTIONS

double FD(double s, double t, double v, double dm){


 return 1 / (double)( 1 + exp( ( energy(s, t, v , dm) - FERMI_ENERGY ) / TEMP ) );

}


// SETTING UP INTEGRATION

struct int_params {double dm;};

int myintegrand(unsigned ndim, const double *x, void *p, unsigned fdim, double *fval){

    struct int_params *params = (struct int_params *)p;

    double dm = (params->dm);

    fval[0] = heaviside_product(x[1], x[2], x[0]) * FD(x[1], x[2], ESCAPE_VEL, dm) * (1 - FD(x[1], x[2], x[0], dm));

    return 0;
}

double tbound(double dm){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL * ESCAPE_VEL) / ( 2.0 * mu(dm) * mu_plus(dm) ) );
}

double sbound(double dm){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL * ESCAPE_VEL) / ( 2.0 * mu(dm) ) );
}



double doing_integral(double dm){

    double val, error;

    struct int_params params = {dm};

    double smax = sbound(dm);
    double tmax = tbound(dm);
    double vmax = ESCAPE_VEL;

    double xl[3] = {0, 0, 0}; // lower bounds for (v, s, t)
    double xu[3] = {vmax, smax, tmax}; // upper bounds for (v, s, t)

    hcubature(1, &myintegrand, &params, 3, xl, xu, 0, 1e-6, 1e-6, ERROR_INDIVIDUAL, &val, &error);

    return val;
}

//MAIN
int main( ){

	double a = doing_integral(1);
	printf("%0.8e\n", a);

	return 0;
}
