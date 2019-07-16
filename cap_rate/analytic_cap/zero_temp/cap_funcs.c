#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "cap_funcs.h"


// =========================================================

double mu(double DMmass){
  return DMmass/NM;
}

double muPlus(double DMmass){
  return (0.5*((DMmass/NM) + 1));
}

double muMinus(double DMMass){
  return (0.5*((DMMass/NM) - 1));
}

//=========================================================



double prefactors(double dm){

	return (32/3)*M_PI* muPlus(dm) * muPlus(dm) * muPlus(dm) *erf(sqrt(3/2)*NSVEL/VELDISP)/NSVEL*(1e6/dm); // rho_chi/dm in GeV/m^2
}

//=========================================================

double constCS(){
	return 1e-45/2 * 1e-4; // m^2
}

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
