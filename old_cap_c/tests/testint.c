#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cubature.h"

#define FILE_NAME "integral_test.txt"


#define VERBOSE 0

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif



// SETTING UP INTEGRATION

struct int_params {double alpha; double beta; double gamma; double zcoord;};

//=========================================================

int firstintegrand(unsigned ndim, const double *x, void *p, unsigned fdim, double *fval){

    struct int_params *params = (struct int_params *)p;
    double a = (params->alpha);
    double b = (params->beta);
    double c = (params->gamma);
    double z = (params->zcoord);

    fval[0] = a * x[0] * x[0] + b * x[1] * x[1] + c * z * z; // z = x[0], y = x[1], z = x[2]

    return 0;
}

//=========================================================

double doing_first_integral(double parone, double partwo, double parthree, double z){

    double val, error;

    struct int_params params = {parone, partwo, parthree, z};

    double xmax = 10;
    double ymax = z * z;

    double xl[2] = {0, 0}; // lower bounds for (x, y, z)
    double xu[2] = {ymax, xmax}; // upper bounds for (x, y, z)
    hcubature(1, &firstintegrand, &params, 2, xl, xu, 0, 1e-6, 1e-6, ERROR_INDIVIDUAL, &val, &error);

    return val;
}


// struct second_int_params {double alpha2; double beta2; double gamma2;};

//=========================================================

int secondintegrand(unsigned ndim, const double *x, void *p, unsigned fdim, double *fval){

	struct int_params *params = (struct int_params *)p;

	double a = (params->alpha);
	double b = (params->beta);
	double c = (params->gamma);

	fval[0] = doing_first_integral(a, b, c, x[0]);

	return 0;

}

//=========================================================

double doing_second_integral(double parone, double partwo, double parthree){

    double val, error;

    struct int_params params2 = {parone, partwo, parthree};

    double zmax = 10;

    double xl[1] = {0}; // lower bounds for z
    double xu[1] = {zmax}; // upper bounds for z
    hcubature(1, &secondintegrand, &params2, 1, xl, xu, 0, 1e-6, 1e-6, ERROR_INDIVIDUAL, &val, &error);

    return val;
}




//MAIN
int main( ){


	double output = doing_second_integral(1, 1, 1);
//	printf("%0.10e\n", output);


	// FILE *file_ptr = fopen(FILE_NAME, "w");
	printf(" Integral Evaluation\n %0.8e\n ", output);
	// fclose(file_ptr);

	return 0;
}
