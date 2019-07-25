#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "cap_funcs.h"




//=========================================================

//GENERAL FUNCTIONS

//=========================================================

double mu(double dm){

	return dm / NM;
}

//=========================================================

double mu_plus(double dm){

	return ( 1.0 + dm / NM) / 2.0;
}


//=========================================================


double FERMI_VEL(double muF){
//square of fermi vel in natural units

	return  (2.0 * muF ) / NM ;

}


double prefactors(double dm){

	return 64*M_PI* mu_plus(dm) * mu_plus(dm) * mu_plus(dm) * mu_plus(dm) *
					erf(sqrt(3/2)*NSVEL/VELDISP)/NSVEL*(1e6/dm); // rho_chi/dm in GeV/m^2
}
//=========================================================

//TRIG FUNCTIONS

//=========================================================

double cosine(double s, double t, double v){
	/* v = escape_vel for initial */

	double num = v * v - s * s - t * t;
	double den = 2.0 * s * t;

	return num*num / den/den;
}

//=========================================================

double sineSqr(double s, double t, double v){

	return (1 - cosine(s, t, v)*cosine(s, t, v) );
}

//=========================================================

//STEP FUNCTIONS

//=========================================================

double step(double f){
	/* a step function to constrain integration region */
	if ( f >= 1){
		return 1;
	}
	else{
		return 0;
	}
}

//=========================================================

double heaviside_product(double s, double t, double vi, double vf){
	/* product of relevant step functions */
	return ( step( vi - fabs(s-t) ) * step( s + t - vi ) *
					step( vf - fabs(s-t) ) * step(s + t - vf) );
}

//=========================================================

//ENERGIES

//=========================================================

double velsqr(double s, double t, double vel, double dm){
	/* square of velocity: v = ESCAPE_VEL for initial particle */
	double a = (2.0 * mu(dm) * mu_plus(dm) * t * t + 2.0 * mu_plus(dm) *
							s*s - mu(dm) * vel * vel);
	return a;
}

//=========================================================

double energy(double s, double t, double vel, double dm){
	return (0.5 * NM * velsqr(s, t, vel, dm));
}

//=========================================================

//FERMI-DIRAC DISTRIBUTION FUNCTIONS

//=========================================================

double FD(double s, double t, double vel, double chempot, double dm){

	double MUP = mu_plus(dm);
	double MU  = mu(dm);


 return (1.0 /(  1.0 + exp( ( 0.5*NM*( 2.*MU*MUP*t*t + 2.*MUP*s*s - MU*vel*vel)/SOL/SOL -  chempot)  / TEMP ) ) );

}

//=========================================================

// INITIAL VELOCITY

//=========================================================


//=========================================================

// DM halo velocity distribution

//=========================================================

double fvel(double DMvel) {
   double prefactor = sqrt(3./2./M_PI)*DMvel/NSVEL/VELDISP;

   return (prefactor*(exp(-3.*(DMvel-NSVEL)*(DMvel-NSVEL)/2./VELDISP/VELDISP)
	 					-  exp(-3.*(DMvel+NSVEL)*(DMvel+NSVEL)/2./VELDISP/VELDISP)));

}

//=========================================================

//INTEGRANDS

//=========================================================

/*
for double int: s = x[0], t = x[1]
for tripple int: v = x[0], s = x[1], t = x[2]
*/
double tbound( double dm, double escvel, double muF){
    return 1.05 * sqrt( ( FERMI_VEL(muF)*SOL*SOL + mu(dm) * escvel * escvel) / ( 2.0 * mu(dm) * mu_plus(dm) ) );
}

//=========================================================

double sbound( double dm, double escvel, double muF){
    return 1.05 * sqrt( (FERMI_VEL(muF)*SOL*SOL + mu(dm) * escvel * escvel )/ ( 2.0 * mu_plus(dm) ) );
}

//=========================================================


double OmegaIntegrand(double *x, size_t dim, void *p){

    struct omega_params *params = (struct omega_params *)p;

    double dm      = (params->dm_mass);
    double chempot = (params->muF);
    double escvel  = (params->escvel);

		double sp = x[1];
		double tp = x[2];


		// double sp = x[1]/(1-x[1]);
		// double tp = x[2]/(1-x[2]);
		//
		// double jacobian = 1./((1. - x[1])*(1. - x[1])*(1. - x[2])*(1. - x[2]));



    return x[0] * tp * heaviside_product(sp, tp, escvel, x[0]) *
		 FD(sp, tp, escvel, chempot, dm) * (1.0 - FD(sp, tp, x[0], chempot, dm));

}

//=========================================================

// DIFFERENTIAL CROSS SECTIONS

//=========================================================

double constCS(){
	return 1e-45/2 * 1e-4; // m^2
}
//=========================================================

double mom4CS(double s, double t, double vinit, double vfin){

return 2*M_PI*(1 - 2*cosine(s, t, vfin)*cosine(s, t, vinit) + cosine(s, t, vfin)*cosine(s, t, vinit)*cosine(s, t, vfin)*cosine(s, t, vinit) + 0.5*sineSqr(s, t, vfin)*sineSqr(s, t, vinit) );
}

//=========================================================

// LOGSPACE

//=========================================================

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
