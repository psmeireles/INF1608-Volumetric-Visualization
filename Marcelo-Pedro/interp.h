#define M_PI_2     1.57079632679489661923
#define M_PI	   2*1.57079632679489661923

double * Chebyshev (int n, double a, double b);
double * NewtonCoef (int n, double * xi, int i, int k, unsigned char *data, double (*f) (int i, double j, int k, unsigned char *data));
double NewtonAval (int n, double *xi, double *bi, double x);