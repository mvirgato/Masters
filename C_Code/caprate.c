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
#define FERMI_ENERGY 0.085 // in GeV
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

	return ( 1 + mu(dm)) / (double)2;
}

//TRIG FUNCTIONS
double initial_cos(double s, double t){
	/* angle between incoming DM and neutron */
    double num = ESCAPE_VEL * ESCAPE_VEL - s * s - t * t;
	double den = 2 * s * t;

	return num / den;
}

double final_cos(double s, double t, double v){
	/* angle between outgoing DM and neutron */

	double num = v * v - s * s - t * t;
	double den = 2 * s * t;

	return num / den;
}

//STEP FUNCTIONS
double step(double f){
	/* a step function to constrain integration region */
	if ( f > 1 || f < -1){
		return 0;
	}
	else{
		return 1;
	}
}

int heaviside_product(double s, double t, double v){
	/* product of relevant step functions */
	return step(initial_cos(s, t)) * step(final_cos(s, t, v));
}

//ENERGIES


double velsqr(double s, double t, double v, double dm){
	/*centre of mass value of final neutron velocity */
	double a = 2 * mu(dm) * mu_plus(dm) * t * t + 2 * mu_plus(dm) * s*s- mu(dm) * v * v;
	return a / (SOL*SOL);
}

double energy(double s, double t, double v, double dm){
	return 0.5 * NM * velsqr(s, t, v, dm);
}


//FERMI-DIRAC DISTRIBUTION FUNCTIONS

double FD(double s, double t, double v, double dm){
	return 1 / (double)( 1 + exp( ( energy(s, t, v , dm) - FERMI_ENERGY ) / TEMP ) );
}

struct int_params {double dm;};

int myintegrand(unsigned ndim, const double *x, void *p, unsigned fdim, double *fval){

    struct int_params *params = (struct int_params *)p;
    
    double dm = (params->dm);
    
    fval[0] = heaviside_product(x[1], x[2], x[0]) * FD(x[1], x[2], ESCAPE_VEL, dm) * (1 - FD(x[1], x[2], x[0], dm));
    
    return 0;
}

double tbound(double dm){
    return 1.05 * sqrt( (SOL*SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL*ESCAPE_VEL) / (2 * mu(dm) * mu_plus(dm) ));
}

double sbound(double dm){
    return 1.05 * sqrt( (SOL*SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL*ESCAPE_VEL) / (2 * mu(dm)));
}



double doing_integral(double dm){
    
    double val, error;
    
    struct int_params params = {dm};
    
    double smax = sbound(dm);
    double tmax = tbound(dm);
    double vmax = ESCAPE_VEL;
    
    double xl[3] = {0, 0, 0};
    double xu[3] = {vmax, smax, tmax};
    
    hcubature(1, &myintegrand, &params, 3, xl, xu, 0, 1e-6, 1e-6, ERROR_INDIVIDUAL, &val, &error);
    
    return val;
    

}

//MAIN
int main(/* int argc, char *argv[] */){
	
	//double mass = atof(argv[1]);

	double a = doing_integral(1.0);
	printf("%f\n", a);
	return 0;
}
