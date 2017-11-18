#include <stdlib.h>
#include <math.h>

double DoubleSimpson (double a, double b, int i, int k, unsigned char *data,  double (*f) (int i, int j, int k, unsigned char *data), double *v);
double AdaptiveSimpson (double a, double b, int i, int k, unsigned char *data, double (*f) (int i, int j, int k, unsigned char *data), double tol);