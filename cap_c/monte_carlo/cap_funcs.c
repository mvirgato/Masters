#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "NSinterp.c"

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

double esc_vel(double radius, int npts){ //rad in km

	return sqrt((2 * 6.67408E-11 * mass_interp(radius, npts) * 2E30) / (radius*1E3));
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

//INTEGRAND
//for double int: s = x[0], t = x[1]
//for tripple int: v = x[0], s = x[1], t = x[2]

double tbound(double dm){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL * ESCAPE_VEL) / ( 2.0 * mu(dm) * mu_plus(dm) ) );
}

double sbound(double dm){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL + mu(dm) * ESCAPE_VEL * ESCAPE_VEL) / ( 2.0 * mu(dm) ) );
}

struct int_params {double dm_mass; double radius; int npoints;};

double myintegrand(double *x, size_t dim, void *p){

    struct int_params *params = (struct int_params *)p;

    double dm = (params->dm_mass);
		double radius = (params->radius);
		int npts = (params->npoints);


    return heaviside_product(x[1], x[2], x[0]) * FD(x[1], x[2], ESCAPE_VEL, dm) * (1 - FD(x[1], x[2], x[0], dm));

}


// LOGSPACE

double* logspace(double a, double b, int n, double u[])
{
    double c;
    int i;

    /* make sure number of points and array are valid */
    if(n < 2 || u == 0)
        return (void*)0;

    /* step size */
    c = (b - a)/(n - 1);

    /* fill vector */
    for(i = 0; i < n - 1; ++i)
        u[i] = pow(10., a + i*c);

    /* fix last entry to 10^b */
    u[n - 1] = pow(10., b);

    /* done */
    return u;
}
