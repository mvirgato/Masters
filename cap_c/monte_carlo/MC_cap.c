#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "cap_funcs.c"

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

double doing_integral_1(double dm, double radius, int npts){

  double res, err;

  struct int_params params = {dm, npts, radius };

  double smax = sbound(dm, radius, npts);
  double tmax = tbound(dm, radius, npts);
  double vmax = esc_vel(radius, npts);


  double xl[3] = {0, 0, 0}; // lower bounds for (v, s, t)
  double xu[3] = {vmax, smax, tmax}; // upper bounds for (v, s, t)

  const gsl_rng_type *T;
  gsl_rng *r;


  gsl_monte_function G = { &myintegrand, 3, &params }; // {function, dimension, params}

  size_t calls = 500000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
    gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s,
                               &res, &err);
    gsl_monte_plain_free (s);

    display_results ("plain", res, err);
  }

  // {
  //   gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
  //   gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
  //                              &res, &err);
  //   gsl_monte_miser_free (s);
  //
  //   display_results ("miser", res, err);
  // }
  //
  // {
  //   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
  //
  //   gsl_monte_vegas_integrate (&G, xl, xu, 3, 100000, r, s,
  //                              &res, &err);
  //   display_results ("vegas warm-up", res, err);
  //
  //   printf ("converging...\n");
  //
  //   do
  //     {
  //       gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
  //                                  &res, &err);
  //       printf ("result = % .6e sigma = % .6f "
  //               "chisq/dof = %.1e\n", res, err, gsl_monte_vegas_chisq (s));
  //     }
  //   while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
  //
  //   display_results ("vegas final", res, err);
  //
  //   gsl_monte_vegas_free (s);
  // }

  gsl_rng_free (r);

  return res;
}

//=========================================================

struct int_params2 {double dm_mass2; double npoints2;};

//=========================================================

double r_integrand(double *x, size_t dim, void *p){

  struct int_params2 *params2 = (struct int_params2 *)p;

  double dm = (params2->dm_mass2);
  int npts = (params2->npoints2);

  return x[0] * x[0] * doing_integral_1(dm, x[0], npts);
}

//=========================================================

double all_integrals(double dm, int npoints){

    double res, err;

    struct int_params2 params2 = {dm, npoints};

    double xl[1] = {1};
    double xu[1] = {12}; // r in km

    const gsl_rng_type *T;
    gsl_rng *r;


    gsl_monte_function F = { &r_integrand, 1, &params2 }; // {function, dimension, params}

    size_t calls = 500000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    {
      gsl_monte_plain_state *s = gsl_monte_plain_alloc (1);
      gsl_monte_plain_integrate (&F, xl, xu, 1, calls, r, s,
                                 &res, &err);
      gsl_monte_plain_free (s);

      display_results ("plain", res, err);
    }

    // {
    //   gsl_monte_miser_state *s = gsl_monte_miser_alloc (1);
    //   gsl_monte_miser_integrate (&F, xl, xu, 1, calls, r, s,
    //                              &res, &err);
    //   gsl_monte_miser_free (s);
    //
    //   display_results ("miser", res, err);
    // }


  // {
  //   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
  //
  //   gsl_monte_vegas_integrate (&F, xl, xu, 1, 100000, r, s,
  //                              &res, &err);
  //   display_results ("vegas warm-up", res, err);
  //
  //   printf ("converging...\n");
  //
  //   do
  //     {
  //       gsl_monte_vegas_integrate (&F, xl, xu, 1, calls/5, r, s,
  //                                  &res, &err);
  //       printf ("result = % .6e sigma = % .6f "
  //               "chisq/dof = %.1e\n", res, err, gsl_monte_vegas_chisq (s));
  //     }
  //   while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
  //
  //   display_results ("vegas final", res, err);
  //
  //   gsl_monte_vegas_free (s);
  // }

  gsl_rng_free (r);

  return res;


}

//=========================================================

int main ()
  {
    int npts;
    npts = readdata("eos_24_lowmass.dat");

    double test = all_integrals(1, npts);
    printf("%0.10e\n", test);

    // int range = 300;
    // double mass_vals[range];
    //
    // logspace(-9, 1, range, mass_vals);
    //
    //
    // int i;
    //
    // for (i = 0; i < range; i++){
    //   printf("%0.10e\n", mass_vals[i]);
    // }
    // FILE *outfile = fopen("cap_rate_MC.dat", "w");
    //
    // for (i = 0; i < range; i++) {
    //    // linear interpolation
    //    fprintf(outfile,"%.10E\t%.10E\n",mass_vals[i], doing_integral( mass_vals[i] ) );
    // }
    //
    // fclose(outfile);

  return 0;
}
