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

	return (0.5*( 1.0 + dm / NM) );
}

//=========================================================

double FermiVel(double chempot){
  return sqrt(0.5*chempot/NM);
}


//=========================================================


double prefactors(double dm){

	return 64*M_PI* mu_plus(dm) * mu_plus(dm) * mu_plus(dm) * mu_plus(dm) *erf(sqrt(3/2)*NSVEL/VELDISP)/NSVEL*(1e6/dm); // rho_chi/dm in GeV/m^2
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

//FERMI-DIRAC DISTRIBUTION FUNCTIONS

//=========================================================

double FD(double s, double t, double vel, double chempot, double dm){

	double MUP = mu_plus(dm);
	double MU  = mu(dm);


 return (1.0 /(  1.0 + exp( ( 0.5*NM*( 2.*MU*MUP*t*t + 2.*MUP*s*s - MU*vel*vel) -  chempot)  / TEMP ) ) );

}

//=========================================================

//STEP FUNCTION

double step(double x){
  if (x <0){
    return ((double) 0.);
  }
  else{
    return ((double) 1.);
  }
}

//=========================================================

// INITIAL VELOCITY

//=========================================================

double w_init(double escvel, double DMvel) {
    return sqrt(escvel*escvel+DMvel*DMvel);
}


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
//=========================================================

//=========================================================

// DIFFERENTIAL CROSS SECTIONS

//=========================================================

double constCS(){
	return 0.5e-49; // m^2
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
