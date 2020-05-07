
static const gsl_rng_type *rngtype;
static const gsl_rng *rng;

void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2);
void tukey(double *data, double alpha, int N);
void tukey_scale(double *s1, double *s2, double alpha, int N);

double tukeyt(double t, double tr, double Tint);

void instrument_noise(double f, double *SAE, double *SXYZ);
