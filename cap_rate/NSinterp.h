double rmin;
double rmax;

int readdata(char * filename);

// Baryon number density interpolation
double nb_interp(double r, int npts);

double Yn_interp(double r, int npts);

double nd_interp(double r, int npts); //neutron number density = nb * Yn

double muFn_interp(double r, int npts);

double mass_interp(double r, int npts);

struct potnint_params {int np;};

double potn_integrand(double rad, void *p);

double potnl(double rad1, double rad2, int np);

double esc_vel(double radius, int npts);

double B_r(double escvel);
