#include <stdio.h>
#include <stdlib.h>
#include<time.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

#include "constants.h"
#include "NSinterp.h"
#include "cap_funcs.h"

//=========================================================

struct ItParams {double svel; double finalvel; double initvel; double chempot; double DMmass;};

//=========================================================

double It1Integrand(double tvel, void *p){

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

  double result, error;

  gsl_integration_workspace * wp = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &It1Integrand;
  F.params   = &params;

  gsl_integration_qag(&F, tmin, tmax, 1.e-6, 1.e-6, 1000, 6, wp, &result, &error);
  gsl_integration_workspace_free(wp);

  return result;
}

//=========================================================

double It2Integral(double svel, double finalvel, double initvel, double chempot, double DMmass){

  struct ItParams params = {svel, finalvel, initvel, chempot, DMmass};

  double tmin = initvel - svel;
  double tmax = svel - finalvel;

  double result, error;

  gsl_integration_workspace * wp = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &It1Integrand;
  F.params   = &params;

  gsl_integration_qag(&F, tmin, tmax, 1.e-6, 1.e-6, 1000, 6, wp, &result, &error);
  gsl_integration_workspace_free(wp);

  return result;
}

//=========================================================

int main()
{
  int npts;
  npts = readdata("eos_24_lowmass.dat");

  double test_radius = 11.3;
  double test_mass   = 1.e0;

  double chempot     = muFn_interp(test_radius, npts);
  double initialvel  = esc_vel_full(test_radius, npts);

  double test_result = It2Integral(initialvel/2, initialvel/2, initialvel, chempot, test_mass);
  printf("result = %0.6E\n", test_result);
}
