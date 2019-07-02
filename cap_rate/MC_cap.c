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


#define Nrpts 100

double radint[Nrpts];
double dCdr[Nrpts];

//=========================================================
void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result           = % .8e\n", result);
  printf ("sigma            = % .8e\n", error);
  printf ("percent error = % .3f\n", error/result * 100);
}


//=========================================================

double OmegaIntegral(double dm, double muFn, double vmax, double DMvel){
  double res, err;

  struct omega_params params = {dm, muFn, vmax, DMvel };

  double smax = sbound(dm, vmax, muFn, DMvel);
  double tmax = tbound(dm, vmax, muFn, DMvel);


  double xl[3] = {0, 0, 0}; // lower bounds for (v, s, t)
  double xu[3] = {vmax, smax, tmax}; // upper bounds for (v, s, t)

  const gsl_rng_type *T;
  gsl_rng *r;


  gsl_monte_function G = { &OmegaIntegrand, 3, &params }; // {function, dimension, params}

  size_t calls = 500000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // {
  //   gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
  //   gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s,
  //                              &res, &err);
  //   gsl_monte_plain_free (s);
  //
  //   // display_results ("plain", res, err);
  // }

  // {
  //   gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
  //   gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
  //                              &res, &err);
  //   gsl_monte_miser_free (s);
  //
  //   display_results ("miser", res, err);
  // }
  //
  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xl, xu, 3, 100000, r, s,
                               &res, &err);
    // display_results ("vegas warm-up", res, err);
    //
    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                   &res, &err);
        // printf ("result = % .6e sigma = % .6e "
        //         "chisq/dof = %.1e\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);

  return res;
}

//=========================================================

double DMvel_integrand(double DMvel, void *p){

  struct DMvelint_params *params2 = (struct DMvelint_params *)p;

  double dm = (params2->dm_mass);
  double muF = (params2->muF);
  double escvel = (params2->escvel);

  return fvel(DMvel)*OmegaIntegral(dm, muF, escvel, DMvel)/DMvel;
}

//=========================================================

double DMvel_integral (double dm, double muFn, double vmax){

   double result, error;

   struct DMvelint_params params2 = {dm, muFn, vmax};

   gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

   gsl_function F;


   F.function = &DMvel_integrand;
   F.params = &params2;

   size_t limit;
   gsl_integration_qagiu (&F, 0, 1e-6, 1e-6, limit, w, &result, &error);



   gsl_integration_workspace_free (w);

   return result;
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

      return exp(dCdr_r);
   }
   else
      return 0.;

}

//=========================================================

double dCdr_integrand(double r, void *p){

  return dCdr_interp(r);
}

//=========================================================


double capture_rate (double rmin, double rmax){

   double result, error;

   gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

   gsl_function F;


   F.function = &dCdr_integrand;
   F.params = 0;

//   size_t limit; gsl_integration_qagiu (&F, 0, 1e-6, 1e-6, limit, w, &result, &error);

   gsl_integration_qag (&F, rmin, rmax, 1.e-6, 1.e-6,1000,6, w, &result, &error);

   gsl_integration_workspace_free (w);

   return result;
}


//=========================================================

int main ()
  {

    double total_time;
  	clock_t start, end;
  	start = clock();
  	srand(time(NULL));

    double ev[Nrpts];

    int npts;
    npts = readdata("eos_24_lowmass.dat");

    // double test = potnl(rmin, rmax, npts);
    // printf("%0.10E\n", test);
   int i,j;
   //double mdm = 1.e0;

    double ev_out = (2.0 * Grav * mass_interp(rmax, npts) * 2E30) /rmax/1e3;

    // FILE *outfile = fopen("esc_vel.dat", "w");
    //
    // for (i = 0; i < Nrpts; i++){
    //   radint[i] = rmin + ((double) i)*(rmax-rmin)/(Nrpts-1);
    //   ev[i] = sqrt( potnl(radint[i], rmax, npts) + ev_out)/SOL;
    //   fprintf(outfile, "%0.8e\t%0.8e\n", radint[i], ev[i] );
    // }
    //
    // fclose(outfile);




    int range = 40;
    double mass_vals[range];


    logspace(-9, 3, range, mass_vals);

    //FILE *outfile = fopen("cap_rate.dat", "w");

    //for (j = 0; j < range; j++){

    FILE *outfile = fopen("cap_rate_rad.dat", "w");
    // fprintf(outfile, "%s\t%s\t%s\t%s\t%s\t%s\n", "radius(km)", "dCdr", "nd", "nresc", "mu", "vmax" );
    for (i = 0; i < Nrpts; i++){
       radint[i] = rmin + ((double) i)*(rmax-rmin)/(Nrpts-1);
       double nd = nd_interp(radint[i], npts);
       double muFn = muFn_interp(radint[i], npts);
       double vmax = sqrt( potnl(radint[i], rmax, npts) + ev_out);
       double nresc = 2.*NM*pow(muFn,1.5)/3./M_PI/M_PI/hbarc/hbarc/hbarc;
       //double Br = B_r(vmax);

       //dCdr[i] = OmegaIntegral( 1., muFn, vmax )*sqrt(1.-Br)/Br/Br;
       dCdr[i] = prefactors(1e-3) * DMvel_integral( 1e-3, muFn, vmax , 0) * nd*nd/nresc * (0.5e-41) /4./M_PI;

       fprintf(outfile,"%0.10E\t%.10E\t%.10E\t%0.10E\t%0.10E\n", radint[i], dCdr[i] , nd_interp(radint[i], npts)/ nresc, muFn, vmax/SOL);
       //dCdr[i] = log(radint[i] * radint[i] *1e6 * nd_interp(radint[i], npts) /nresc*OmegaIntegral( mass_vals[j], muFn, vmax ));
    }


 //   double test = capture_rate(rmin, rmax);
 //   printf("%0.10E\n", test);






//       fprintf(outfile,"%0.10E\t%.10E\n", mass_vals[j], capture_rate(rmin, rmax)*1e3);

//    }

    fclose(outfile);

    end = clock();
  	total_time = ((double) (end - start)) / CLOCKS_PER_SEC/60;
  	printf("\nTime taken is: %f min\n", total_time);

    return 0;
  }
