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

#include "NSinterp.h"
#include "constants.h"
#include "TT_funcs.h"

// =========================================================



// =========================================================

int main(){

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  double testmass = 1e15; //in eV
  double testrad  = 11.;

  double initvel     = esc_vel_full(testrad, npts)/SOL;
  double testinitmom = testmass*initvel/sqrt(1 - initvel*initvel);
  double testchempot = muFn_interp(testrad, npts);

  double test = TT;

  printf("%0.8e\n", test);




  return 0;
}
