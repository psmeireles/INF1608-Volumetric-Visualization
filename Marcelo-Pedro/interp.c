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

double * NewtonCoef (int n, double * xi, int i, int k, unsigned char *data, double (*f) (int i, double j, int k, unsigned char *data)){
	double *b = (double *) malloc(n*sizeof(double)), **M;
	int cont, j;

	// Creating cache to store Newton Coefficients
	M = mat_create(n, n);

	for(j = 0; j < n; j++)
		for(cont = j; cont >= 0; cont--){
			if(cont == j)
				M[cont][j] = f(xi[cont], j, k, data);
			else
				M[cont][j] = (M[cont+1][j] - M[cont][j-1]) / (xi[j] - xi[cont]);
			if(cont == 0)
				b[j] = M[cont][j];
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