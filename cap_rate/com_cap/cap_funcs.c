#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "cap_funcs.h"




//=========================================================

//GENERAL FUNCTIONS

//=========================================================


//=========================================================

// LOGSPACE

//=========================================================

double* logspace(double a, double b, int n, double u[])
{
    double c;
    int i;

    /* make sure number of points and array are valid */
    if(n < 2 || u == 0)
        return (void*)0;

    /* step size */
    c = (b - a)/(n - 1);

    /* fill vector */
    for(i = 0; i < n - 1; ++i)
        u[i] = pow(10., a + i*c);

    /* fix last entry to 10^b */
    u[n - 1] = pow(10., b);

    /* done */
    return u;
}
