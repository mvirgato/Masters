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

#include "constants.h"
#include "cap_funcs.h"
#include "NSinterp.h"

// =========================================================

#define Nrpts 100

double radint[Nrpts];
double cap_full[Nrpts];
double dCdr[Nrpts];

// =========================================================

double FermiVel(double chempot){
  return sqrt(0.5*chempot/NM);
}



// =========================================================

double LambdaPlus(double initvel, double finvel, double DMmass, double chempot){
  return sqrt( FermiVel(chempot) * FermiVel(chempot) + mu(DMmass) * (initvel*initvel - finvel*finvel) );
}

double LambdaMinus(double initvel, double finvel, double DMmass, double chempot){
  double a = sqrt( FermiVel(chempot) * FermiVel(chempot) - mu(DMmass) * (initvel*initvel - finvel*finvel) );

  if (isnan(a) == 1){
 return 0;
 }
 else{
 return a;
 }
}

// =========================================================

double alphaPlus(double initvel, double finvel, double DMmass){
  return (muPlus(DMmass)*finvel + muMinus(DMmass)*initvel);
}

double alphaMinus(double initvel, double finvel, double DMmass){
  return (muPlus(DMmass)*finvel - muMinus(DMmass)*initvel);
}

// =========================================================

double betaPlus(double initvel, double finvel, double DMmass){
  return (muMinus(DMmass)*finvel + muPlus(DMmass)*initvel);
}

double betaMinus(double initvel, double finvel, double DMmass){
  return (muMinus(DMmass)*finvel - muPlus(DMmass)*initvel);
}

// =========================================================

double rhoPlus(double initvel, double finvel, double DMmass){
  return (mu(DMmass)*finvel + muPlus(DMmass)*(finvel + initvel));
}

double rhoMinus(double initvel, double finvel, double DMmass){
  return (mu(DMmass)*finvel - muPlus(DMmass)*(finvel + initvel));
}

// =========================================================

double phiPlus(double initvel, double finvel, double DMmass){
  return (mu(DMmass)*finvel + muPlus(DMmass)*(finvel - initvel));
}

double phiMinus(double initvel, double finvel, double DMmass){
  return (mu(DMmass)*finvel - muPlus(DMmass)*(finvel - initvel));
}

// =========================================================

double step(double x){
  if (x <0){
    return ((double) 0.);
  }
  else{
    return ((double) 1.);
  }
}

// =========================================================


double SixMuI (double initvel, double finvel, double DMmass, double chempot){

  double FV  = FermiVel(chempot);
  double LP  = LambdaPlus(initvel, finvel, DMmass, chempot);
  double LM  = LambdaMinus(initvel, finvel, DMmass, chempot);
  double AP  = alphaPlus(initvel, finvel, DMmass);
  double AM  = alphaMinus(initvel, finvel, DMmass);
  double BP  = betaPlus(initvel, finvel, DMmass);
  double BM  = betaMinus(initvel, finvel, DMmass);
  double RP  = rhoPlus(initvel, finvel, DMmass);
  double RM  = rhoMinus(initvel, finvel, DMmass);
  double PP  = phiPlus(initvel, finvel, DMmass);
  double PM  = phiMinus(initvel, finvel, DMmass);

  double MU  = mu(DMmass);
  double MUP = muPlus(DMmass);
  double MUM = muMinus(DMmass);



  double l1  = (2.*FV + finvel*MU)*(FV - finvel*MU)*(FV - finvel*MU)*( step((RM - FV)/(2.*MUP)) + step( (FV - RP)/(2.*MUP) ) );

  double l2  = 2.*( MU*( 2.*MU + 3.)*finvel*finvel - 6.*MU*MUP*finvel*finvel - FV*MU*finvel + 2.*FV*FV ) * (FV - finvel*MU)*step( (FV - PM)/(2.*MUP) )*step( (RP - FV)/(2.*MUP) );

  double l34  = 3.*MU*MUP*MUP*( (-finvel - initvel)*(finvel*finvel - initvel*initvel)*step( FV*FV - AM*AM ) -(finvel - initvel)*( (finvel*finvel - initvel*initvel)*( step( FV*FV - AP*AP ) - step( FV*FV - BP*BP) )  - (finvel + initvel)*(finvel + initvel)*step( FV*FV - BM*BM ))  );

  double l5  = (FV - initvel*MU)*(FV - initvel*MU)*(2.*FV + initvel*MU)*step( (FV + AM)/(2.*MUP) )*step( (AP - FV)/(2.*MUP) );

  double l6  = (2.*FV - initvel*MU)*(FV + initvel*MU)*(FV + initvel*MU)* step( (AM - FV)/(2.*MUP) )*step( (FV + AM)/(2.*MUP) );

  double l7  = (FV + finvel*MU)*(MU*(2.*MU + 3.)*finvel*finvel - 6.*MU*MUP*finvel*finvel + FV*MU*finvel + 2.*FV*FV)*( step( (FV +RM)/(2.*MUP) ) +step((-FV - RP)/(2.+MUP)) );

  double l8  = (finvel*MU - LP)*(2.*MU*MU*finvel*finvel -6.*MU*MUP*finvel*finvel - MU*LP*finvel + 3.*FV*FV - LP*LP +3.*initvel*initvel*MU)*( step( (RM - LP)/(2.*MUP) ) + step( (LP - RP)/(2.*MUP) ) );

  double l9  = (finvel*MU + LP)*(2.*MU*MU*finvel*finvel -6.*MU*MUP*finvel*finvel + MU*LP*finvel + 3.*FV*FV - LP*LP +3.*initvel*initvel*MU)*( step( (RM + LP)/(2.*MUP) ) + step( (-LP - RP)/(2.*MUP) ) );

  double l10 = (finvel*MU - LP)*( -6.*MU*MUP*finvel*finvel + 2.*FV*FV +MU*( (2.*MU +1)*finvel*finvel - LP*finvel +2.*initvel*initvel ) )*step( (LP - PM)/(2.*MUP) )*step( (RP - LP)/(2.*MUP) );

  double l11 = (finvel*MU + LP)*( -6.*MU*MUP*finvel*finvel + 2.*FV*FV +MU*( (2.*MU +1)*finvel*finvel + LP*finvel +2.*initvel*initvel ) )*step( (-LP - PM)/(2.*MUP) )*step( (RP + LP)/(2.*MUP) );

  double l12 = (initvel*MU - LM)*(2.*FV*FV + MU*( 2.*finvel*finvel + initvel*( 2.*MU*initvel + initvel - LM ) ) - 6.*initvel*initvel*MU*MUP)*step( (AM + LM)/(2.*MUP) )*step( (AP - LM)/(2.*MUP) );

  double l13 = (initvel*MU + LM)*(2.*FV*FV + MU*( 2.*finvel*finvel + initvel*( 2.*MU*initvel + initvel + LM ) ) - 6.*initvel*initvel*MU*MUP)*step( (AM - LM)/(2.*MUP) )*step( (AP + LM)/(2.*MUP) );


  return -l1 - l2 - l34 + l5 - l6 + l7 - l8 - l9 - l10 - l11 + l12 + l13;
}

// =========================================================

struct OmegaParms {double initvel; double chempot; double DMmass;};

//=========================================================

double OmegaIntegrand(double vvel, void *p){

  struct OmegaParms *params = (struct OmegaParms *)p;

  double dm       =  (params->DMmass);
  double chempot  =  (params->chempot);
  double w        =  (params->initvel);

  return (vvel * SixMuI(w, vvel, dm, chempot));

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
  printf("==========================\n");
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

  int i,j;

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  double chempot     = muFn_interp(11.3, npts);
  double initialvel  = esc_vel_full(11.3, npts)/SOL;

  int range = 100;
  double mass_vals[range];
  logspace(-9, 3, range, mass_vals);

  // double test = LambdaMinus(initialvel, 0.9*initialvel, 1e3, chempot);
  // printf("%0.8e\n", test);


  FILE *outfile = fopen("dat_zero_temp.dat", "w");
  //
  for (j =0; j<range;j++){
    // double test = 0;
    // test = 3.*(mass_vals[j]/NM + 1)*(I1Integral(initialvel, initialvel, muFn, mass_vals[j]) + I2Integral(initialvel, initialvel, muFn, mass_vals[j]));
    // fprintf(outfile, "%0.8e\t%0.8E\n",mass_vals[j], test);

    for (i = 0; i<Nrpts; i++){

      radint[i] = rmin + ((double) i)*(rmax-rmin)/(Nrpts-1);

      double initialvel  = esc_vel_full(radint[i], npts)/SOL*(6./8.);
      double nd          = nd_interp(radint[i], npts)* 1.e45 ; // m^-3
      double chempot     = muFn_interp(radint[i], npts);
      double ndfree      = pow(2.*NM*chempot,1.5)/3./M_PI/M_PI/hbarc/hbarc/hbarc * 1.e45 ;


      dCdr[i] = prefactors(mass_vals[j])*constCS()*OmegaIntegral(initialvel, chempot, mass_vals[j]) * nd*nd/ndfree;
      // fprintf(outfile, "%0.10E\t%0.10E\t%0.10E\n", radint[i], dCdr[i], nd*nd/ndfree);

    }
  //
  // double test = capture_rate(rmin, rmax);
  // printf("%0.8e\n", test);

    cap_full[j] = capture_rate(rmin, rmax);
    fprintf(outfile, "%0.10e\t%0.10e\n", mass_vals[j], cap_full[j]);

    for (i = 0; i<Nrpts; i++){
      dCdr[i] = 0;
    }
  }
  //
  fclose(outfile);

  return 0;
}
