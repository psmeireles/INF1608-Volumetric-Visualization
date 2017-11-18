#include "simpson.h"

double DoubleSimpson (double a, double b, int i, int k, unsigned char *data,  double (*f) (int i, int j, int k, unsigned char *data), double *v){
	double sac, sab, scb, h = b-a, c = (a+b)/2, err;
	sab = h/6*(f(i, (int) a, k, data) + 4*f(i, (int) (a+b)/2, k, data) + f(i, (int) b, k, data));
	sac = (h/2)/6*(f(i, (int) a, k, data) + 4*f(i, (int) (a+c)/2, k, data) + f(i, (int) c, k, data));
	scb = (h/2)/6*(f(i, (int) c, k, data) + 4*f(i, (int) (c+b)/2, k, data) + f(i, (int) b, k, data));
	err = fabs((sab - (sac + scb))/15);
	*v = sac + scb + err;
	return err;
}

double AdaptiveSimpson (double a, double b, int i, int k, unsigned char *data, double (*f) (int i, int j, int k, unsigned char *data), double tol){
	double c = (a + b)/2, err, v;
	err = DoubleSimpson(a, b, i, k, data, f, &v);
	if(err < tol)
		return v;
	else
		return AdaptiveSimpson(a, c, i, k, data, f, tol/2) + AdaptiveSimpson(c, b, i, k, data, f, tol/2);
}