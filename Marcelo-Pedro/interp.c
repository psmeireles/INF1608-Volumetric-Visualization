#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interp.h"
#include "matrix.h"

double * Chebyshev (int n, double a, double b){
	double *v = (double *) malloc(n*sizeof(double)), sub = (b-a)/2, add = (a+b)/2;
	int i;
	for(i = 1; i <= 2*n-1; i += 2)
		v[i/2] = sub * cos(i*M_PI_2 / n) + add;
	return v;
}

double * NewtonCoef (int n, double * xi, double (*f) (double)){
	double *b = (double *) malloc(n*sizeof(double)), **M;
	int i, j;

	// Creating cache to store Newton Coefficients
	M = mat_create(n, n);

	for(j = 0; j < n; j++)
		for(i = j; i >= 0; i--){
			if(i == j)
				M[i][j] = f(xi[i]);
			else
				M[i][j] = (M[i+1][j] - M[i][j-1]) / (xi[j] - xi[i]);
			if(i == 0)
				b[j] = M[i][j];
		}

	mat_free(n, M);

	return b;
}

double NewtonAval (int n, double *xi, double *bi, double x){
	double p = 0, aux;
	int i, j;
	for(i = 0; i < n; i++){
		aux = bi[i];
		for(j = 0; j < i; j++){
			aux *= (x - xi[j]);
		}
		p += aux;
	}
	return p;
}