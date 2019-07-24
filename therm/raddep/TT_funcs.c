#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// =========================================================

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

// =========================================================

#include "TT_funcs.h"
#include "NSinterp.h"
#include "constants.h"

// =========================================================

double constME(double DMmass){

return 16*NM*NM*DMmass*DMmass;

}

// =========================================================

double mom4supME(double q){

  return q*q*q*q;
}

// =========================================================

double prefac(double DMmass){

  return TEMP/64./M_PI/M_PI/M_PI/DMmass/DMmass;

}

// =========================================================

double cosIntegrand(double x, void *p){

  struct qParams *params1 = (struct qParams *)p;


  double ki         = (params1->initmom);
  double kf         = (params1->finmom);
  double dm         = (params1->DMmass);
  double muFn       = (params1->chempot);

  double q          = sqrt(ki*ki + kf*kf -2.*ki*kf*x);

  double q0         = 0.5*(ki*ki - kf*kf)/dm;

  double z          = q0/TEMP;

  double eMinus     = 0.2*NM*(q0 - 0.5*q*q/NM)/q/q;

  double zetaMinus  = log( (1 + exp((eMinus - muFn)/TEMP))/(1 + exp((eMinus+q0-muFn)/TEMP)) );

  double brace      = (z*(1 + zetaMinus/z))/(1 + exp(-z));

  return constME(dm)*brace*kf*kf/q;
}

// =========================================================

double cosIntegral(double initmom, double finmom, double DMmass, double chempot){

  struct qParams params1 = {initmom, finmom, DMmass, chempot};

  double qmin = -1.;
  double qmax = 1.;

  double result1, error1;

  gsl_integration_workspace * wp1 = gsl_integration_workspace_alloc (5000);

  gsl_function F1;

  F1.function = &cosIntegrand;
  F1.params   = &params1;

  gsl_integration_qag(&F1, qmin, qmax, 1.e-6, 1.e-6, 5000, 6, wp1, &result1, &error1);
  gsl_integration_workspace_free(wp1);

  return result1;

}

// =========================================================

double kfIntegrand(double kf, void *p){

  struct kfParams *params2 = (struct kfParams *)p;

  double ki         = (params2->initmom);
  double dm         = (params2->DMmass);
  double muFn       = (params2->chempot);

  return cosIntegral(ki, kf, dm, muFn);

}

// =========================================================

double kfIntegral(double initmom, double DMmass, double chempot){

  struct kfParams params2 = {initmom, DMmass, chempot};

  double kfmin = 0;
  double kfmax = initmom;

  double result2, error2;

  gsl_integration_workspace * wp2 = gsl_integration_workspace_alloc (5000);

  gsl_function F2;

  F2.function = &kfIntegrand;
  F2.params   = &params2;

  gsl_integration_qag(&F2, kfmin, kfmax, 1.e-6, 1.e-6, 5000, 6, wp2, &result2, &error2);
  gsl_integration_workspace_free(wp2);

  return result2;

}

// =========================================================

double Gamma(double initmom, double DMmass, double chempot){
  //Gamma at a specific r value
  return prefac(DMmass)*kfIntegral(initmom, DMmass, chempot);
}

// =========================================================

double volumeIntegrand(double rad, void *p){

  return rad*rad;

}

// =========================================================

double volumeIntegral(){

  double result3, error3;

  gsl_integration_workspace * wp3 = gsl_integration_workspace_alloc (5000);

  gsl_function F3;

  F3.function = &volumeIntegrand;
  F3.params   = 0;

  gsl_integration_qag(&F3, rmin, rmax, 1.e-6, 1.e-6, 5000, 6, wp3, &result3, &error3);
  gsl_integration_workspace_free(wp3);

  return result3;


}

// =========================================================

double volAvgRateIntegrand(double r, void *p){

  struct rGammaParams *params4 = (struct rGammaParams *)p;


  double dm   = (params4->DMmass);
  double np   = (params4->npts);

  double eval = esc_vel_full(dm, np)/SOL;
  double ki   = eval*dm;
  // double ki   = eval*dm/sqrt(1 - eval*eval);
  double muFn = muFn_interp(r, np)*1e9;

  return r*r*Gamma(ki, dm, muFn);

}

// =========================================================

double volAvgRateIntegral(double DMmass, int npts){

  struct rGammaParams params4 = {DMmass, npts};

  double result4, error4;

  gsl_integration_workspace * wp4 = gsl_integration_workspace_alloc (5000);

  gsl_function F4;

  F4.function = &volAvgRateIntegrand;
  F4.params   = &params4;


  gsl_integration_qag(&F4, rmin, rmax, 1.e-6, 1.e-6, 5000, 6, wp4, &result4, &error4);
  gsl_integration_workspace_free(wp4);

  return result4;
}

// =========================================================

// Average final energy

// =========================================================

double numeratorGammaIntegrand(double kf, void *p){

  struct kfParams *params5 = (struct kfParams *)p;

  double ki         = (params5->initmom);
  double dm         = (params5->DMmass);
  double muFn       = (params5->chempot);

  return cosIntegral(ki, kf, dm, muFn)*kf*kf*0.5/dm/dm;

}

// =========================================================

double numeratorGammaIntegral(double initmom, double DMmass, double chempot){

  struct kfParams params5 = {initmom, DMmass, chempot};

  double nummin = 0;
  double nummax = initmom;

  double result5, error5;

  gsl_integration_workspace * wp5 = gsl_integration_workspace_alloc (5000);

  gsl_function F5;

  F5.function = &numeratorGammaIntegrand;
  F5.params   = &params5;

  gsl_integration_qag(&F5, nummin, nummax, 1.e-6, 1.e-6, 5000, 6, wp5, &result5, &error5);
  gsl_integration_workspace_free(wp5);

  return result5*0.5/DMmass/DMmass;

}

// =========================================================

double volAvgEnergyIntegrand(double r, void *p){

  struct rGammaParams *params6 = (struct rGammaParams *)p;


  double dm   = (params6->DMmass);
  double np   = (params6->npts);

  double eval = esc_vel_full(dm, np)/SOL;
  double ki   = eval*dm;
  // double ki   = eval*dm/sqrt(1 - eval*eval);
  double muFn = muFn_interp(r, np)*1e9;

  return r*r*numeratorGammaIntegral(ki, dm, muFn);

}

// =========================================================

double volAvgEnergyIntegral(double DMmass, int npts){

  struct rGammaParams params6 = {DMmass, npts};

  double result6, error6;

  gsl_integration_workspace * wp6 = gsl_integration_workspace_alloc (5000);

  gsl_function F6;

  F6.function = &volAvgEnergyIntegrand;
  F6.params   = &params6;


  gsl_integration_qag(&F6, rmin, rmax, 1.e-6, 1.e-6, 5000, 6, wp6, &result6, &error6);
  gsl_integration_workspace_free(wp6);

  return result6;
}
