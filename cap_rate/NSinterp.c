#include <stdio.h>
#include <math.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "constants.h"
#include "NSinterp.h"

#define N 1000


double rad[N]; // km
double mass[N]; // Msun
double nb[N];  // Baryon number density (fm^-3)
double Yn[N];
double muFn[N];  // Neutron chemical potential (GeV)
double nd[N]; //neutron number densite


//=========================================================


/*#################################################
                Reading the data
##################################################*/

int readdata(char * filename){

   FILE *datafile;
   int i, npts;
   double Ye,Yp,Ymu;
   double muFp,muFe,muFmu;

   datafile = fopen(filename, "r");

   if (datafile == NULL) {
      // printf("Error: can't open file %s\n",filename);
      return 1;
   }
   else {
      printf("File %s opened successfully.\n",filename);
      // Skip header
      fscanf(datafile, "%*[^\n]\n", NULL);
      i = 0;

      while (!feof(datafile)) {
         fscanf(datafile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &rad[i],&mass[i],&nb[i],&Ye,&Ymu,&Yp,&Yn[i],&muFn[i],&muFp,&muFe,&muFmu);
         muFn[i] /= 1e3; // GeV
         nd[i] = nb[i] * Yn[i];
//         printf("%d\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\n", i,rad[i],mass[i],nb[i],Ye,Ymu,Yp,Yn[i],muFn[i],muFp,muFe,muFmu);
         i++;
      }

   }
   npts=i;
   fclose(datafile);

   rmin = rad[0];
   rmax = 12.387; // NS core
   //rmax = 12.158947177343357; // R inner crust
  //  printf("Npts: %d\n",npts);
   return npts;
}

//=========================================================




// Baryon number density interpolation
double nb_interp(double r, int npts) {
  // printf("called nb_interp\n" );

   double nb_r;

   if (r >= rad[0] && r <= rad[npts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npts);

      gsl_spline_init (spline, rad, nb, npts);

      nb_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return nb_r;
   }
   else
      return 0.;

}

//=========================================================

double Yn_interp(double r, int npts) {
  // printf("called Yn\n");

   double Yn_r;

   if (r >= rad[0] && r <= rad[npts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npts);

      gsl_spline_init (spline, rad, Yn, npts);

      Yn_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return Yn_r;
   }
   else
      return 0.;

}

//=========================================================


double nd_interp(double r, int npts) {
  // printf("called nd\n");

   double nd_r;

   if (r >= rad[0] && r <= rad[npts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npts);

      gsl_spline_init (spline, rad, nd, npts);

      nd_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return nd_r;
   }
   else
      return 0.;

}

//=========================================================

double muFn_interp(double r, int npts) {

  // printf("called muFn\n");

   double muFn_r;

   if (r >= rad[0] && r <= rad[npts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npts);

      gsl_spline_init (spline, rad, muFn, npts);

      muFn_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return muFn_r;
   }
   else
      return 0.;

}

double mass_interp(double r, int npts) {
  // printf("called mass\n");

   double mass_r;

   if (r >= rad[0] && r <= rad[npts-1]) {

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      //Linear splines
/*      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, npts);*/

      // Cubic splines
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npts);

      gsl_spline_init (spline, rad, mass, npts);

      mass_r = gsl_spline_eval (spline, r, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);

      return mass_r;
   }
   else
      return 0.;
    }


//=========================================================

double potn_integrand(double rad, void *p){

  struct potnint_params *params = ( struct potnint_params *)p;

  int np = (params->np);

  return (2 * Grav * mass_interp(rad, np)*2e30)/rad/rad/1e3;

}

//=========================================================

double potnl(double rad1, double rad2, int np){

  size_t calls = 1000;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);


  double res, err;

  struct potnint_params params = {np};

  //    double xl[1] = {1};
  //    double xu[1] = {12}; // r in km

  gsl_function F;
  F.function = &potn_integrand;
  F.params = &params; // {function, dimension, params}


  gsl_integration_qags (&F, rad1 , rad2 , 0.0 , 1e-7, calls, w, &res, &err );


  gsl_integration_workspace_free (w);

  return res;
}

//=========================================================

// double rel_potn(double rad1, double rad2, int np){
//
//   double ev =
//   double I =
//
// }

//=========================================================

double esc_vel(double radius, int npts){ //rad in km
    //return ESCAPE_VEL;
    return sqrt((2.0 * Grav * mass_interp(radius, npts) * 2E30) /radius/1e3);
}

//=========================================================

double B_r(double escvel) {
   return 1.0-escvel*escvel/SOL/SOL;
}


//
// int main(int argc, char** argv)
// {
// //
//      int npts;
//      int i;
//
//      readdata(argv[1]);
//      npts = readdata("eos_24_lowmass.dat");
//
//
// //
//    FILE *outfile = fopen("n_density.dat","w");
//
//    for (i=0;i<=N;i++) {
//       // linear interpolation
//       double radius = rad[0] + (12.1- rad[0])*((double) i)/N;
//       fprintf(outfile,"%.10E\t%.10E\n",radius,nd_interp(radius,npts));
// //      printf("%d\t%.10E\t%.10E\n",i,radius,nb_interp(radius,npts));
//    }
//
//    fclose(outfile);
//
//    FILE *outfile2 = fopen("Ynabund.dat","w");
//
//    for (i=0;i<=N;i++) {
//       // linear interpolation
//       double radius = rad[0] + (12.1- rad[0])*((double) i)/N;
//       fprintf(outfile2,"%.10E\t%.10E\n",radius,Yn_interp(radius,npts));
// //      printf("%d\t%.10E\t%.10E\n",i,radius,nb_interp(radius,npts));
//    }
//
//    fclose(outfile2);
//
//    FILE *outfile3 = fopen("muFnchempot.dat","w");
//
//    for (i=0;i<=N;i++) {
//       // linear interpolation
//       double radius = rad[0] + (12.1- rad[0])*((double) i)/N;
//       fprintf(outfile3,"%.10E\t%.10E\n",radius,muFn_interp(radius,npts));
// //      printf("%d\t%.10E\t%.10E\n",i,radius,nb_interp(radius,npts));
//    }
//
//    fclose(outfile3);
//
//    FILE *outfile4 = fopen("mass_ns.dat","w");
//
//    for (i=0;i<=N;i++) {
//       // linear interpolation
//       double radius = rad[0] + (12.1- rad[0])*((double) i)/N;
//       fprintf(outfile4,"%.10E\t%.10E\n",radius,mass_interp(radius,npts));
// //      printf("%d\t%.10E\t%.10E\n",i,radius,nb_interp(radius,npts));
//    }
//
//    fclose(outfile4);
//
//    return 0;
// }
