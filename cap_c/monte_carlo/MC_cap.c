#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "cap_funcs.c"

/* Computation of the integral,

      I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */

/* For simplicity we compute the integral over the region
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */

// double exact = 120.174;

// struct int_params { double dm_mass; double v_e; };

// double
// g (double *x, size_t dim, void *p)
// {
//   // (void)(dim); /* avoid unused parameter warnings */
//   // (void)(params);
//   struct int_params * fp = (struct int_params *)p;
//
//   double a = (fp->a);
//   double b = (fp->b);
//
//   return (x[0] * x[1] * x[2])*mu(a+b+c);
// }
struct int_params {double dm_mass; double vel;};

double myintegrand(double *x, size_t dim, void *p){

    struct int_params *params = (struct int_params *)p;

    double dm = (params->dm_mass);
		double vel = (params->vel);


    return heaviside_product(x[0], x[1], vel) * FD(x[0], x[1], ESCAPE_VEL, dm) * (1 - FD(x[0], x[1], vel, dm));

}


void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result           = % .6e\n", result);
  printf ("sigma            = % .6e\n", error);
  printf ("percentage error = % .6e\n", error/result );
  // printf ("exact  = % .6f\n", exact);
  // printf ("error  = % .6f = %.2g sigma\n", result - exact,
          // fabs (result - exact) / error);
}

int
main ()
{
  struct int_params fp = {10, ESCAPE_VEL/2};

  double dm = (int_params->dm_mass);


  double res, err;

  double xl[2] = { 0, 0};
  double xu[2] = { sbound(10), tbound(10) };

  const gsl_rng_type *T;
  gsl_rng *r;


  gsl_monte_function G = { &myintegrand, 2, &fp }; // {function, dimension, params}

  size_t calls = 500000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&G, xl, xu, 2, calls, r, s,
                               &res, &err);
    gsl_monte_plain_free (s);

    display_results ("plain", res, err);
  }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
    gsl_monte_miser_integrate (&G, xl, xu, 2, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,
                                   &res, &err);
        printf ("result = % .6e sigma = % .6f "
                "chisq/dof = %.1e\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);

  return 0;
}
