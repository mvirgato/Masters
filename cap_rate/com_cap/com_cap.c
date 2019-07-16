#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

double FD(double nvel, double initvel, double finvel, double chempot, double dmmass){

  double num   = (double) 1.0;
  double denom = exp((0.5*NM*nvel*nvel - chempot - 0.5*dmmass*(initvel*initvel - finvel*finvel))/TEMP) + 1;

  return num/denom;
}

//=========================================================

struct xparams {double nvel; double initvel; double finvel; double chempot; double dmmass;};

//=========================================================

double xintegrand(double x, void *p){
  struct xparams *params1 = (struct xparams *)p;

  double u       = (params1->nvel);
  double w       = (params1->initvel);
  double v       = (params1->finvel);
  double chempot = (params1->chempot);
  double dm      = (params1->dmmass);

  return sqrt(w*w + u*u - 2.*w*u*x)*FD(u, 0., 0., chempot, dm)*FD(u, w, v, chempot, dm);

}

//=========================================================

double xintegral(double nvel, double initvel, double finvel, double chempot, double dmmass ){

  struct xparams params1 = {nvel, initvel, finvel, chempot, dmmass};

  double res1, err1;

  gsl_integration_workspace * wp1 = gsl_integration_workspace_alloc (5000);

  gsl_function F1;

  F1.function = &xintegrand;
  F1.params   = &params1;


  gsl_integration_qag(&F1, -1, 1, 1.e-6, 1.e-6, 5000, 6, wp1, &res1, &err1);
  gsl_integration_workspace_free(wp1);

  return res1;
}

//=========================================================

struct uparams { double initvel; double finvel; double chempot; double dmmass;};

//=========================================================

double uintegrand(double u, void *p){

  struct uparams *params2 = (struct uparams *)p;

  double w       = (params2->initvel);
  double v       = (params2->finvel);
  double chempot = (params2->chempot);
  double dm      = (params2->dmmass);

  return u*u*xintegral(u, w, v, chempot, dm);

}

//=========================================================



double uintegral(double initvel, double finvel, double chempot, double dmmass ){

  struct xparams params2 = {initvel, finvel, chempot, dmmass};

  double res2, err2;

  gsl_integration_workspace * wp2 = gsl_integration_workspace_alloc (5000);

  gsl_function F2;

  F2.function = &uintegrand;
  F2.params   = &params2;

  size_t limit;

  gsl_integration_qagiu(&F2, 0, 1.e-6, 1.e-6, limit, wp2, &res2, &err2);
  gsl_integration_workspace_free(wp2);

  return res2;
}

//=========================================================

struct vparams {double initvel; double chempot; double dmmass;};

//=========================================================

double vintegrand(double v, void *p){

  struct vparams *params3 = (struct vparams *)p;

  double w       = (params3->initvel);
  double chempot = (params3->chempot);
  double dm      = (params3->dmmass);

  return uintegral(w, v, chempot, dm);

}

//=========================================================

double vintegral(double initvel, double chempot, double dmmass ){

  struct xparams params3 = {initvel, chempot, dmmass};

  double res3, err3;

  gsl_integration_workspace * wp3 = gsl_integration_workspace_alloc (5000);

  gsl_function F3;

  F3.function = &vintegrand;
  F3.params   = &params3;


  gsl_integration_qag(&F3, 0., initvel, 1.e-6, 1.e-6, 5000, 6, wp3, &res3, &err3);
  gsl_integration_workspace_free(wp3);

  return res3;
}


//=========================================================

struct rparams {double dmmass; int npts;};

//=========================================================

double rintegrand(double r, void *p){

  struct rparams *params4 = (struct rparams *)p;

  double dm   = (params4->dmmass);
  int    npts = (params4->npts);

  double initialvel  = esc_vel_full(r, npts)/SOL;
  double nd          = nd_interp(r, npts)* 1.e45 ; // m^-3
  double chempot     = muFn_interp(r, npts);
  double ndfree      = pow(2.*NM*chempot,1.5)/3./M_PI/M_PI/hbarc/hbarc/hbarc * 1.e45 ;


  return initialvel*r*r*vintegral(initialvel, chempot, dm)*SOL*SOL * nd*nd/ndfree ;
}

//=========================================================

double rintegral(double dmmass, int npts){

  struct rparams params4 = {dmmass, npts};

  double res4, err4;

  gsl_integration_workspace * wp4 = gsl_integration_workspace_alloc (5000);

  gsl_function F4;

  F4.function = &rintegrand;
  F4.params   = &params4;


  gsl_integration_qag(&F4, rmin, rmax, 1.e-6, 1.e-6, 5000, 6, wp4, &res4, &err4);
  gsl_integration_workspace_free(wp4);

  return res4;

}
//=========================================================

int main() {

  double testmass = 1e-3;

  int i;

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  int range = 100;
  double mass_vals[range];
  logspace(-9, 4, range, mass_vals);

  FILE *outfile = fopen("com_cap.dat", "w");

  for (i = 0; i<Nrpts; i++){

        cap_full[i] = 4*M_PI*M_PI*M_PI*0.5e-49*SOL*SOL * rintegral(mass_vals[i], npts);
        fprintf(outfile, "%0.10E\t%0.10E\n", mass_vals[i], cap_full[i]);

      }

      fclose(outfile);


  // double test = vintegral(w, chempot, testmass);
  // printf("%0.8e\n", test);



  return 0;
}
