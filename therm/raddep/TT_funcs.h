double constME(double DMmass);

double prefac(double DMmass);

struct qParams {double initmom; double finmom; double DMmass; double chempot;};

struct kfParams {double initmom; double DMmass; double chempot;};

double qIntegral(double initmom, double finmom, double DMmass, double chempot);

double kfIntegral(double initmom, double DMmass, double chempot);
