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

double prefac(double DMmass){

  return TEMP/64./M_PI/M_PI/M_PI/DMmass/DMmass;

}

// =========================================================

double constME(double DMmass){

return 16*NM*NM*DMmass*DMmass;

}

// =========================================================

double mom4supME(double q){

  return q*q*q*q;
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
double finalEnergyNumIntegrand(double kf, void *p){

  struct kfParams *params5 = (struct kfParams *)p;

  double ki         = (params5->initmom);
  double dm         = (params5->DMmass);
  double muFn       = (params5->chempot);

  return 0.5*cosIntegral(ki, kf, dm, muFn)*kf*kf/dm;

}

// =========================================================

double finalEnergyNumIntegral(double initmom, double DMmass, double chempot){

  struct kfParams params5 = {initmom, DMmass, chempot};

  double kfmin = 0;
  double kfmax = initmom;

  double result5, error5;

  gsl_integration_workspace * wp5 = gsl_integration_workspace_alloc (5000);

  gsl_function F5;

  F5.function = &finalEnergyNumIntegrand;
  F5.params   = &params5;

  gsl_integration_qag(&F5, kfmin, kfmax, 1.e-6, 1.e-6, 5000, 6, wp5, &result5, &error5);
  gsl_integration_workspace_free(wp5);

  return result5;

}

// =========================================================

double nextEnergy(double initmom, double DMmass, double chempot){

  double numerator   = finalEnergyNumIntegral(initmom, DMmass, chempot);
  double denominator = kfIntegral(initmom, DMmass, chempot);

  return numerator/denominator;


}

// =========================================================

double EtoMom(double energy, double DMmass){

  return sqrt(2*DMmass*energy);
}

// =========================================================

double momToEnergy(double mom, double DMmass){

  return (0.5*mom*mom/DMmass);
}

// =========================================================

double TTintegrand(double rad, void *p){

  struct rGammaParams *params6 = (struct rGammaParams *)p;

  double dm      = (params6->DMmass);
  int    npts    = (params6->npts);

  double muFn    = muFn_interp(rad, npts)*1e9;

  double v0      = esc_vel_full(rad, npts)/SOL;
  double k0      = dm*v0/sqrt( 1 - v0*v0);
  double E0      = sqrt(k0*k0 + dm*dm);

  double sum    = 0;
  double Gdummy = 0;

  Gdummy = (Gamma(k0, dm, muFn))/dm;
  sum   += 1./Gdummy;

  double initialE = E0;
  double finalE   = nextEnergy(k0, dm, muFn);

  while ( (initialE - finalE)>TEMP ) {


    double initialk = EtoMom(finalE, dm);


    Gdummy = Gamma(initialk, dm, muFn)/dm;
    sum   += 1./Gdummy;

    printf("%0.8e\n", sum);



    initialE = finalE;
    finalE   = nextEnergy(initialk, dm, muFn);
  }

  return rad*rad*sum/dm;

}

// =========================================================

double TTintegral(double DMmass, double npts){

  struct kfParams params6 = {DMmass, npts};

  double result6, error6;

  gsl_integration_workspace * wp6 = gsl_integration_workspace_alloc (5000);

  gsl_function F6;

  F6.function = &TTintegrand;
  F6.params   = &params6;


  gsl_integration_qag(&F6, rmin, rmax, 1.e-6, 1.e-6, 5000, 6, wp6, &result6, &error6);
  gsl_integration_workspace_free(wp6);

  return result6;


}
