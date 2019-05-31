#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "cap_funcs.c"


void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result           = % .8e\n", result);
  printf ("sigma            = % .8e\n", error);
  printf ("percent error = % .3f\n", error/result * 100);
}



double doing_integral(double dm){
  printf("mass = %0.3e GeV\n", dm);

  double res, err;

  struct int_params params = {dm};

  double smax = sbound(dm);
  double tmax = tbound(dm);
  double vmax = ESCAPE_VEL;

  double xl[3] = {0, 0, 0}; // lower bounds for (v, s, t)
  double xu[3] = {vmax, smax, tmax}; // upper bounds for (v, s, t)

  const gsl_rng_type *T;
  gsl_rng *r;


  gsl_monte_function G = { &myintegrand, 3, &params }; // {function, dimension, params}

  size_t calls = 700000;

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

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xl, xu, 3, 500000, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                   &res, &err);
        // printf ("result = % .6e sigma = % .6f "
        //         "chisq/dof = %.1e\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);

  return 0;
}

  int
  main ()
  {

    doing_integral(1e2);

  return 0;
}
