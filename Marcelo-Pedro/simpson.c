#include "simpson.h"

double DoubleSimpson (double a, double b, double (*f) (double x), double *v){
	double sac, sab, scb, h = b-a, c = (a+b)/2, err;
	sab = h/6*(f(a) + 4*f((a+b)/2) + f(b));
	sac = (h/2)/6*(f(a) + 4*f((a+c)/2) + f(c));
	scb = (h/2)/6*(f(c) + 4*f((c+b)/2) + f(b));
	err = fabs((sab - (sac + scb))/15);
	*v = sac + scb + err;
	return err;
}

double AdaptiveSimpson (double a, double b, double (*f) (double x), double tol){
	double c = (a + b)/2, err, v;
	err = DoubleSimpson(a, b, f, &v);
	if(err < tol)
		return v;
	else
		return AdaptiveSimpson(a, c, f, tol/2) + AdaptiveSimpson(c, b, f, tol/2);
}