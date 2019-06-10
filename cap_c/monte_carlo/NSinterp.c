#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_spline.h>

#define N 1000


double rad[N]; // km
double mass[N]; // Msun
double nb[N];  // Baryon number density (fm^-3)
double Yn[N];
double muFn[N];  // Neutron chemical potential (GeV)





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
//         printf("%d\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\n", i,rad[i],mass[i],nb[i],Ye,Ymu,Yp,Yn[i],muFn[i],muFp,muFe,muFmu);
         i++;
      }

   }
   npts=i;
   fclose(datafile);
  //  printf("Npts: %d\n",npts);
   return npts;
}

// Baryon number density interpolation
double nb_interp(double r, int npts) {

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

double Yn_interp(double r, int npts) {

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

double muFn_interp(double r, int npts) {

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



// int main(int argc, char** argv)
// {
//
//      int npts;
//      int i;
//
//      //readdata(argv[1]);
//      npts = readdata("eos_24_lowmass.dat");


//
//    FILE *outfile = fopen("nbdensity.dat","w");
//
//    for (i=0;i<=N;i++) {
//       // linear interpolation
//       double radius = rad[0] + (12.1- rad[0])*((double) i)/N;
//       fprintf(outfile,"%.10E\t%.10E\n",radius,nb_interp(radius,npts));
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
}
