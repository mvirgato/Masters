#define SOL 299792458.0 // speed of light in vacuum m/s
#define hbarc 0.1973269788/1e3  // GeV fm

#define NM 0.939 *1e9   //neutron mass in GeV
#define Grav 6.67408E-11
#define ESCAPE_VEL sqrt((2 * Grav * 1.4 * 2E30) / (10E3))	//NS escape velocity in m/s. Need to make into function of radius
#define TEMP (1E3 ) / (1.16E4)//NS temp in eV
#define FERMI_ENERGY 0.085 // in GeV (will need to make a variable)
#define VELDISP 270e3 //DM velocity dispersion for MB dist
#define NSVEL 200e3 // NS velocity in galactic frame
