double constME(double DMmass);

double prefac(double DMmass);

struct qParams {double initmom; double finmom; double DMmass; double chempot;};

struct kfParams {double initmom; double DMmass; double chempot;};

struct rGammaParams {double DMmass; int npts;};

double qIntegral(double initmom, double finmom, double DMmass, double chempot);

double kfIntegral(double initmom, double DMmass, double chempot);

double Gamma(double initmom, double DMmass, double chempot);

double volumeIntegral();

double volAvgRateIntegral(double DMmass, int npts);

double finalEnergyNumIntegral(double initmom, double DMmass, double chempot);

double nextEnergy(double initmom, double DMmass, double chempot);

double TTintegrand(double rad, double dm, int npts);
