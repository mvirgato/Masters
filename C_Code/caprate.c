#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

//GLOBAL CONSTANTS
#define ESCAPE_VEL sqrt((2 * 6.67408E-11 * 1.4 * 2E30) / (10E3))	//NS escape velocity in m/s. Need to make into function of radius
#define TEMP (1E3 * 1E-9) / (1.16E4)	//NS temp in GeV
#define NM 0.939	//neutron mass in GeV
#define FERMI_ENERGY 0.085 // in GeV
#define FERMI_VEL (2 * FERMI_ENERGY) / NM
#define SOL 299792458 // speed of light in vacuum m/s
#define VELDISP 270e3 //DM velocity dispersion for MB dist
#define NSVEL 200e3 // NS velocity in galactic frame

//GENERAL FUNCTIONS
float mu(float dm){
	/* RATIO OF DARK MATTER TO NEUTRON MASS */

	return dm / NM;
}

float mu_plus(float dm){

	return ( 1 + mu(dm)) / (float)2;
}

//TRIG FUNCTIONS
float initial_cos(float s, float t){
	/* angle between incoming DM and neutron */

	float num = powf(ESCAPE_VEL, 2) - powf(s, 2) - powf(t, 2);
	float den = 2 * s * t;

	return num / den;
}

float final_cos(float s, float t, float v){
	/* angle between outgoing DM and neutron */

	float num = powf(v, 2) - powf(s, 2) - powf(t, 2);
	float den = 2 * s * t;

	return num / den;
}

//STEP FUNCTIONS
float step(float f){
	/* a step function to constrain integration region */
	if ( f > 1 || f < -1){
		return 0;
	}
	else{
		return 1;
	}
}

int heaviside_product(float s, float t, float v){
	/* product of relevant step functions */
	return step(initial_cos(s, t)) * final_cos(s, t, v);
}

//ENERGIES
float initial_vel(float s, float t, float dm){
	/* centre of mass value of initial neutron velocity */
	float a = 2 * mu(dm) * mu_plus(dm) * powf(t, 2) + 2 * mu_plus(dm) * powf(s, 2) - mu(dm) * powf(ESCAPE_VEL, 2);
	return a / powf(SOL, 2);
}

float final_vel(float s, float t, float v, float dm){
	/*centre of mass value of final neutron velocity */
	float a = 2 * mu(dm) * mu_plus(dm) * powf(t, 2) + 2 * mu_plus(dm) * powf(s, 2) - mu(dm) * powf(v, 2);
	return a / powf(SOL, 1);
}

float initial_energy(float s, float t, float dm){
	return 0.5 * NM * initial_vel(s, t, dm);
}

float final_energy(float s, float t, float v, float dm){
	return 0.5 * NM * final_vel(s, t, v, dm);
}
//FERMI-DIRAC DISTRIBUTION FUNCTIONS

float initial_FD(float s, float t, float dm){
	return 1 / (float)( 1 + exp( ( initial_energy(s, t, dm) - FERMI_ENERGY ) / TEMP ) );
}

//MAIN
int main(){
	float a = initial_FD(1, 1, 1);
	printf("%f\n", a );
	return 0;
}
