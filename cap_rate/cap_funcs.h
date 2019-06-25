
double tbound( double dm, double escvel, double muF, double DMvel);

double sbound( double dm, double escvel, double muF, double DMvel);

struct int_params {double dm_mass; double muF; double escvel; double DMvel;};

double w(double escvel, double DMvel);

double fvel(double DMvel);

double myintegrand(double *x, size_t dim, void *p);

double* logspace(double a, double b, int n, double u[]);
