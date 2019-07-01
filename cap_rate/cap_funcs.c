#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "cap_funcs.h"




//=========================================================

//GENERAL FUNCTIONS

//=========================================================

double mu(double dm){
	/* RATIO OF DARK MATTER TO NEUTRON MASS */

	return dm / NM;
}

//=========================================================

double mu_plus(double dm){

	return ( 1.0 + mu(dm)) / 2.0;
}


//=========================================================


double FERMI_VEL(double muF){

	return  (2.0 * muF ) / NM ;

}

//=========================================================

//TRIG FUNCTIONS

//=========================================================

double cosine(double s, double t, double v){
	/* v = escape_vel for initial */

	double num = v * v - s * s - t * t;
	double den = 2.0 * s * t;

	return num / den;
}

//=========================================================

//STEP FUNCTIONS

//=========================================================

double step(double f){
	/* a step function to constrain integration region */
	if ( f > 1){
		return 1;
	}
	else{
		return 0;
	}
}

//=========================================================

double heaviside_product(double s, double t, double vi, double vf){
	/* product of relevant step functions */
	return step( vi - fabs(s-t) ) * step( s + t - vi ) * step( vf - fabs(s-t) ) * step(s + t -vf) ;
}

//=========================================================

//ENERGIES

//=========================================================

double velsqr(double s, double t, double vel, double dm){
	/* square of velocity: v = ESCAPE_VEL for initial particle */
	double a = 2.0 * mu(dm) * mu_plus(dm) * t * t + 2.0 * mu_plus(dm) * s*s - mu(dm) * vel * vel;
	return a / (SOL*SOL);
}

//=========================================================

double energy(double s, double t, double vel, double dm){
	return 0.5 * NM * velsqr(s, t, vel, dm);
}

//=========================================================

//FERMI-DIRAC DISTRIBUTION FUNCTIONS

double FD(double s, double t, double vel, double chempot, double dm){


 return 1.0 /(  1.0 + exp( ( energy(s, t, vel , dm) -  chempot)  / TEMP ) );

}

//=========================================================

// INITIAL VELOCITY

double w_init(double escvel, double DMvel) {
    return sqrt(escvel*escvel+DMvel*DMvel);
}


//=========================================================

// DM halo velocity distribution

double fvel(double DMvel) {
   double prefactor = sqrt(3./2./M_PI)*DMvel/NSVEL/VELDISP;

   return prefactor*(exp(-3.*(DMvel-NSVEL)*(DMvel-NSVEL)/2./VELDISP/VELDISP) -  exp(-3.*(DMvel-NSVEL)*(DMvel+NSVEL)/2./VELDISP/VELDISP));

}

//=========================================================

//INTEGRAND

//=========================================================

/*
for double int: s = x[0], t = x[1]
for tripple int: v = x[0], s = x[1], t = x[2]
*/
double tbound( double dm, double escvel, double muF, double DMvel){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL(muF) + mu(dm) * w_init(escvel, DMvel) * w_init(escvel, DMvel)) / ( 2.0 * mu(dm) * mu_plus(dm) ) );
}

//=========================================================

double sbound( double dm, double escvel, double muF, double DMvel){
    return 1.05 * sqrt( (SOL * SOL * FERMI_VEL(muF) + mu(dm) * w_init(escvel, DMvel) * w_init(escvel, DMvel) )/ ( 2.0 * mu(dm) ) );
}

//=========================================================



//=========================================================

double myintegrand(double *x, size_t dim, void *p){

    struct int_params *params = (struct int_params *)p;

    double dm = (params->dm_mass);
    double chempot = (params->muF);
    double escvel = (params->escvel);
		double dmvel = (params->DMvel);

    return (16.0 / dm ) * mu_plus(dm) * mu_plus(dm) * mu_plus(dm) * mu_plus(dm) * (x[0] ) * x[2] * heaviside_product(x[1], x[2], w_init(escvel, dmvel), x[0]) *
		 FD(x[1], x[2], w_init(escvel, dmvel), chempot, dm) * (1.0 - FD(x[1], x[2], x[0], chempot, dm));

}

//=========================================================

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
