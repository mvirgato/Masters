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
#include "NSinterp.h"

// =========================================================

double FermiVel(double chempot){
  return sqrt(0.5*chempot/NM);
}

// =========================================================

double mu(double DMmass){
  return DMmass/NM;
}

double muPlus(double DMmass){
  return (0.5*(DMmass/NM + 1));
}

double muMinus(double DMMass){
  return (0.5*(DMMass/NM - 1));
}

// =========================================================

double LambdaPlus(double initvel, double finvel, double DMmass, double chempot){
  return sqrt( FermiVel(chempot) * FermiVel(chempot) + mu(DMmass) * (initvel*initvel - finvel*finvel) );
}

double LambdaMinus(double initvel, double finvel, double DMmass, double chempot){
  return sqrt( FermiVel(chempot) - mu(DMmass) * (initvel*initvel - finvel*finvel) );
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
  if (x <=0){
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

  double l3  = 3.*MU*MUP*MUP*( (-finvel - initvel)*(finvel*finvel - initvel*initvel)*step( FV*FV - AM*AM ) );

  double l4  = 3.*MU*MUP*MUP*( -(finvel - initvel)*( (finvel*finvel - initvel*initvel)*( step( FV*FV - AP*AP ) - step( FV*FV - BP*BP) ) ) - (finvel + initvel)*(finvel + initvel)*step( FV*FV - BM*BM ) );

  double l5  = (FV - initvel*MU)*(FV - initvel*MU)*(2.*FV + initvel*MU)*step( (FV + AM)/(2.*MUP) )*step( (AP - FV)/(2.*MUP) );

  double l6  = (2.*FV - initvel*MU)*(FV + initvel*MU)*(FV + initvel*MU)* step( (AM - FV)/(2.*MUP) )*step( (FV + AM)/(2.*MUP) );

  double l7  = (FV + finvel*MU)*(MU*(2.*MU + 3.)*finvel*finvel - 6.*MU*MUP*finvel*finvel + FV*MU*finvel + 2.*FV*FV)*( step( (FV +RM)/(2.*MUP) ) +step(-(FV + RP)/(2.+MUP)) );

  double l8  = (finvel*MU - LP)*(2.*MU*MU*finvel*finvel -6.*MU*MUP*finvel*finvel - MUP*LP*finvel + 3.*FV*FV - LP*LP +3.*initvel*initvel*MU)*( step( (RM - LP)/(2.*MUP) ) + step( -(LP - RP)/(2.*MUP) ) );

  double l9  = (finvel*MU + LP)*(2.*MU*MU*finvel*finvel -6.*MU*MUP*finvel*finvel + MUP*LP*finvel + 3.*FV*FV - LP*LP +3.*initvel*initvel*MU)*( step( (RM + LP)/(2.*MUP) ) + step( -(LP + RP)/(2.*MUP) ) );

  double l10 = (finvel*MU - LP)*( -6.*MU*MUP*finvel*finvel + 2.*FV*FV +MU*( (2.*MU +1)*finvel*finvel - LP*finvel +2.*initvel*initvel ) )*step( (LP - PM)/(2.*MUP) )*step( (RP - LP)/(2.*MUP) );

  double l11 = (finvel*MU + LP)*( -6.*MU*MUP*finvel*finvel + 2.*FV*FV +MU*( (2.*MU +1)*finvel*finvel + LP*finvel +2.*initvel*initvel ) )*step( -(LP + PM)/(2.*MUP) )*step( (RP + LP)/(2.*MUP) );

  double l12 = (initvel*MU - LM)*(2.*FV*FV + MU*( 2.*finvel*finvel + initvel*( 2.*MU*initvel + initvel - LM ) ) - 6.*initvel*initvel*MU*MUP)*step( (AM + LM)/(2.*MUP) )*step( (AP - LM)/(2.*MUP) );

  double l13 = (initvel*MU + LM)*(2.*FV*FV + MU*( 2.*finvel*finvel + initvel*( 2.*MU*initvel + initvel + LM ) ) - 6.*initvel*initvel*MU*MUP)*step( (AM - LM)/(2.*MUP) )*step( (AP + LM)/(2.*MUP) );


  return -l1 - l2 - l3 - l4 + l5 - l6 + l7 - l8 - l9 - l10 - l11 + l12 + l13;
}

int main()
{

  int npts;
  npts = readdata("eos_24_lowmass.dat");

  double test_rad   = 11.3;
  double test_mass  = 1.e1;
  double muFn       = muFn_interp(test_rad, npts);
  double initialvel = esc_vel_full(test_rad, npts)/SOL/SOL;


  double test = SixMuI(initialvel, initialvel, test_mass, muFn);
  printf("%0.8E\n", test);



  return 0;
}
