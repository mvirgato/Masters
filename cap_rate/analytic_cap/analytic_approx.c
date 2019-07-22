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

struct ItParams {double svel; double finalvel; double initvel; double chempot; double DMmass;};

//=========================================================

double ItIntegrand(double tvel, void *p){

  struct ItParams *params = (struct ItParams *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double s        =  (params->svel);
  double v        =  (params->finalvel);
  double w        =  (params->initvel);

  double MU       =  mu(dm);
  double FV       =  FermiVel(chempot);
  double MUP      =  mu_plus(dm);

  return tvel * FD(s, tvel, w, chempot, dm) * ( 1 - FD(s, tvel, v, chempot, dm) );
  // return tvel*step(FV*FV + w*w - 2.*MU*MUP*tvel*tvel - 2.*MUP*s*s)*step(2.*MU*MUP*tvel*tvel +2.*MUP*s*s - FV*FV - MU*v*v);
}

//=========================================================

double It1Integral(double svel, double finalvel, double initvel, double chempot, double DMmass){

  struct ItParams params = {svel, finalvel, initvel, chempot, DMmass};

  double tmin = initvel - svel;
  double tmax = finalvel + svel;

  double result1, error1;

  gsl_integration_workspace * wp1 = gsl_integration_workspace_alloc (5000);

  gsl_function F1;

  F1.function = &ItIntegrand;
  F1.params   = &params;

  gsl_integration_qag(&F1, tmin, tmax, 1.e-6, 1.e-6, 5000, 6, wp1, &result1, &error1);
  gsl_integration_workspace_free(wp1);

  return result1;
}

//=========================================================

double It2Integral(double svel, double finalvel, double initvel, double chempot, double DMmass){

  struct ItParams params2 = {svel, finalvel, initvel, chempot, DMmass};

  double tmin = initvel - svel;
  double tmax = svel - finalvel;

  double result2, error2;

  gsl_integration_workspace * wp2 = gsl_integration_workspace_alloc (5000);

  gsl_function F2;

  F2.function = &ItIntegrand;
  F2.params   = &params2;

  gsl_integration_qag(&F2, tmin, tmax, 1.e-6, 1.e-6, 5000, 6, wp2, &result2, &error2);
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

  gsl_integration_workspace * wp3 = gsl_integration_workspace_alloc (5000);

  gsl_function F3;

  F3.function = &I1Integrand;
  F3.params   = &params3;

  gsl_integration_qag(&F3, smin, smax, 1.e-6, 1.e-6, 5000, 6, wp3, &res3, &err3);
  gsl_integration_workspace_free(wp3);

  return res3;

}

//=========================================================

double I2Integral(double finalvel, double initvel, double chempot, double DMmass){

  struct IParams params4 = {finalvel, initvel, chempot, DMmass};

  double smin = 0.5*(initvel + finalvel);

  double res4, err4;

  gsl_integration_workspace * wp4 = gsl_integration_workspace_alloc (5000);

  gsl_function F4;

  F4.function = &I2Integrand;
  F4.params   = &params4;


  gsl_integration_qagiu(&F4, smin, 1.e-6, 1.e-6, 5000, wp4, &res4, &err4);
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

  return (vvel * (I1Integral(vvel, w, chempot, dm) + I2Integral(vvel, w, chempot, dm)));

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
  printf("result        = %0.8E\n", res5);
  printf("error         = %0.8E\n", err5);
  printf("==============================\n");
  return res5;
}

//=========================================================

double dCdr_interp(double r) {
  // printf("called nd\n");

   double dCdr_r;

   if (r >= radint[0] && r <= radint[Nrpts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Nrpts);

      gsl_spline_init (spline, radint, dCdr, Nrpts);

      dCdr_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return dCdr_r;
   }
   else
      return 0.;

}

//=========================================================

double dCdr_integrand(double r, void *p){

  return dCdr_interp(r);
}

//=========================================================


double capture_rate (double min, double max){

   double result6, error6;

   gsl_integration_workspace * wp = gsl_integration_workspace_alloc (1000);

   gsl_function F6;


   F6.function = &dCdr_integrand;
   F6.params = 0;

//   size_t limit; gsl_integration_qagiu (&F, 0, 1e-6, 1e-6, limit, w, &result, &error);

   gsl_integration_qag (&F6, min, max, 1.e-6, 1.e-6,1000,6, wp, &result6, &error6);
   gsl_integration_workspace_free (wp);

   return result6;
}

//=========================================================

int main()
{

  double total_time;
  clock_t start, end;
  start = clock();
  srand(time(NULL));

  int i,j;

  int range = 50;
  double mass_vals[range];
  logspace(-1, 4, range, mass_vals);

  // double test_mass   = 1.e2;

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  double test_rad   = 11.3;
  // double test_mass  = 1e-3;
  // double muFn       = muFn_interp(test_rad, npts);
  // double initialvel = esc_vel_full(test_rad, npts)/SOL;



  FILE *outfile = fopen("analytic_complete_cap.dat", "w");
  // FILE *outfile = fopen("analytic_rad_cap.dat", "w");
  //
  for (j =0; j<range;j++){
    // double test = 0;
    // test = 3.*(mass_vals[j]/NM + 1)*(I1Integral(initialvel, initialvel, muFn, mass_vals[j]) + I2Integral(initialvel, initialvel, muFn, mass_vals[j]));
    // fprintf(outfile, "%0.8e\t%0.8E\n",mass_vals[j], test);

    for (i = 0; i<Nrpts; i++){

      radint[i] = rmin + ((double) i)*(rmax-rmin)/(Nrpts-1);

      double initialvel  = esc_vel_full(radint[i], npts);
      double nd          = nd_interp(radint[i], npts); // m^-3
      double chempot     = muFn_interp(radint[i], npts);
      double ndfree      = pow(2.*NM*chempot,1.5)/3./M_PI/M_PI/hbarc/hbarc/hbarc;


      dCdr[i] = prefactors(mass_vals[j])*constCS()* OmegaIntegral(initialvel, chempot, mass_vals[j]) * nd*nd/ndfree*radint[i]*radint[i]*1.e54/SOL/SOL/SOL;
      // fprintf(outfile, "%0.10E\t%0.10E\t%0.10E\n", radint[i], dCdr[i], initialvel);

    }

 // double test = capture_rate(rmin, rmax);
 // printf("%0.8e\n", test);

    cap_full[j] = capture_rate(rmin, rmax);
    fprintf(outfile, "%0.10e\t%0.10e\n", mass_vals[j], cap_full[j]);

    for (i = 0; i<Nrpts; i++){
      dCdr[i] = 0;
    }
  }

  fclose(outfile);

  end = clock();
  total_time = ((double) (end - start)) / CLOCKS_PER_SEC;

  int hrs, min;
  float sec;

  if (total_time < 60){
    hrs =0; min =0; sec = total_time;
  } else if (total_time < 3600){
    hrs = 0; min = floor(total_time/60); sec = total_time - 60 * min;
  } else {
    hrs = floor(total_time/3600); min = floor(total_time/60 - 60 * hrs); sec = total_time - 60 * min;
  }

  printf("\nTime taken is: %d hrs : %d min : %f sec\n",hrs, min, sec);
}
