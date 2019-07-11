#include <stdio.h>
#include <stdlib.h>
#include<time.h>

//=========================================================

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

//=========================================================

#include "constants.h"
#include "NSinterp.h"
#include "cap_funcs.h"

//=========================================================

#define Nrpts 100

double radint[Nrpts];
double cap_full[Nrpts];
double dCdr[Nrpts];

//=========================================================

struct ItParams {double svel; double finalvel; double initvel; double chempot; double DMmass;};

//=========================================================

double ItIntegrand(double tvel, void *p){

  struct ItParams *params = (struct ItParams *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double s        =  (params->svel);
  double v        =  (params->finalvel);
  double w        =  (params->initvel);

  return tvel * FD(s, tvel, w, chempot, dm) * ( 1 - FD(s, tvel, v, chempot, dm) );
}

//=========================================================

double It1Integral(double svel, double finalvel, double initvel, double chempot, double DMmass){

  struct ItParams params = {svel, finalvel, initvel, chempot, DMmass};

  double tmin = initvel - svel;
  double tmax = finalvel + svel;

  double result1, error1;

  gsl_integration_workspace * wp1 = gsl_integration_workspace_alloc (1000);

  gsl_function F1;

  F1.function = &ItIntegrand;
  F1.params   = &params;

  gsl_integration_qag(&F1, tmin, tmax, 1.e-6, 1.e-6, 1000, 6, wp1, &result1, &error1);
  gsl_integration_workspace_free(wp1);

  return result1;
}

//=========================================================

double It2Integral(double svel, double finalvel, double initvel, double chempot, double DMmass){

  struct ItParams params2 = {svel, finalvel, initvel, chempot, DMmass};

  double tmin = initvel - svel;
  double tmax = svel - finalvel;

  double result2, error2;

  gsl_integration_workspace * wp2 = gsl_integration_workspace_alloc (1000);

  gsl_function F2;

  F2.function = &ItIntegrand;
  F2.params   = &params2;

  gsl_integration_qag(&F2, tmin, tmax, 1.e-6, 1.e-6, 1000, 6, wp2, &result2, &error2);
  gsl_integration_workspace_free(wp2);

  return result2;
}

//=========================================================

struct IParams {double finalvel; double initvel; double chempot; double DMmass;};

//=========================================================

double I1Integrand(double Svel, void *p){

  struct IParams *params = (struct IParams *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double v        =  (params->finalvel);
  double w        =  (params->initvel);

  return It1Integral(Svel, v, w, chempot, dm);
}


//=========================================================

double I2Integrand(double Svel, void *p){

  struct IParams *params = (struct IParams *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double v        =  (params->finalvel);
  double w        =  (params->initvel);

  return It2Integral(Svel, v, w, chempot, dm);
}

//=========================================================

double I1Integral(double finalvel, double initvel, double chempot, double DMmass){

  struct IParams params3 = {finalvel, initvel, chempot, DMmass};

  double smin = 0.5*(initvel - finalvel);
  double smax = 0.5*(initvel + finalvel);

  double res3, err3;

  gsl_integration_workspace * wp3 = gsl_integration_workspace_alloc (1000);

  gsl_function F3;

  F3.function = &I1Integrand;
  F3.params   = &params3;

  gsl_integration_qag(&F3, smin, smax, 1.e-6, 1.e-6, 1000, 6, wp3, &res3, &err3);
  gsl_integration_workspace_free(wp3);

  return res3;

}

//=========================================================

double I2Integral(double finalvel, double initvel, double chempot, double DMmass){

  struct IParams params4 = {finalvel, initvel, chempot, DMmass};

  double smin = 0.5*(initvel + finalvel);

  double res4, err4;

  gsl_integration_workspace * wp4 = gsl_integration_workspace_alloc (1000);

  gsl_function F4;

  F4.function = &I2Integrand;
  F4.params   = &params4;

  size_t limit;

  gsl_integration_qagiu(&F4, smin, 1.e-6, 1.e-6, limit, wp4, &res4, &err4);
  gsl_integration_workspace_free(wp4);

  return res4;

}

//=========================================================

struct OmegaParms {double initvel; double chempot; double DMmass;};

//=========================================================

double OmegaIntegrand(double vvel, void *p){

  struct OmegaParms *params = (struct OmegaParms *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double w        =  (params->initvel);

  return vvel*(I1Integral(vvel, w, chempot, dm) + I2Integral(vvel, w, chempot, dm));

}

//=========================================================

double OmegaIntegral(double initvel, double chempot, double DMmass){

  struct OmegaParms params5 = {initvel, chempot, DMmass};

  double res5, err5;

  gsl_integration_workspace * wp5 = gsl_integration_workspace_alloc (1000);

  gsl_function F5;

  F5.function = &OmegaIntegrand;
  F5.params   = &params5;

  gsl_integration_qag(&F5, 0, initvel, 1.e-6, 1.e-6, 1000, 6, wp5, &res5, &err5);
  gsl_integration_workspace_free(wp5);
  printf("result = %0.8E\n", res5);
  return res5;
}

//=========================================================

int main()
{

  int i,j;

  double test_mass   = 1.e0;

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  FILE *outfile = fopen("analytic_rad_cap.dat", "w");

  for (i = 0; i<Nrpts; i++){

    radint[i] = rmin + ((double) i)*(rmax-rmin)/(Nrpts-1);

    double initialvel  = esc_vel_full(radint[i], npts);
    double nd          = nd_interp(radint[i], npts)* 1.e45 ; // m^-3
    double chempot     = muFn_interp(radint[i], npts);
    double ndfree      = pow(2.*NM*chempot,1.5)/3./M_PI/M_PI/hbarc/hbarc/hbarc * 1.e45 ;


    dCdr[i] = prefactors(test_mass)* constCS() * OmegaIntegral(initialvel, chempot, test_mass)* nd*nd/ndfree;
    fprintf(outfile, "%0.10E\t%0.10E\t%0.10E\n", radint[i], dCdr[i], nd*nd/ndfree);

  }

  fclose(outfile);
}
