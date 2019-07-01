
double tbound( double dm, double escvel, double muF, double DMvel);

double sbound( double dm, double escvel, double muF, double DMvel);

struct omega_params {double dm_mass; double muF; double escvel; double DMvel;};

double w(double escvel, double DMvel);

double fvel(double DMvel);

double OmegaIntegrand(double *x, size_t dim, void *p);

struct DMvelint_params {double wr; double diff_rate;};

double DMvel_integrand(double DMvel, void *p);

double* logspace(double a, double b, int n, double u[]);
