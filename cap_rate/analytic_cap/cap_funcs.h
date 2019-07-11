double prefactors(double dm);

double tbound( double dm, double escvel, double muF, double DMvel);

double sbound( double dm, double escvel, double muF, double DMvel);

struct omega_params {double dm_mass; double muF; double escvel; double DMvel;};

double fvel(double DMvel);

double FD(double s, double t, double vel, double chempot, double dm);

// double OmegaIntegrand(double *x, size_t dim, void *p);

struct DMvelint_params {double dm_mass; double muF; double escvel;};

double DMvel_integrand(double DMvel, void *p);

double constCS();

double mom4CS(double s, double t, double vinit, double vfin);

double* logspace(double a, double b, int n, double u[]);
