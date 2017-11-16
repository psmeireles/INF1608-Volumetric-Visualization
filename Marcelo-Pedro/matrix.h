#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double ** mat_create (int m, int n);
void mat_free (int m, double **A);
void mat_transpose (int m, int n, double **A, double **T);
void mat_multv (int m, int n, double **A, double *v, double *w);
void mat_multm (int m, int n, int q, double **A, double **B, double **C);
int mat_equals (int m, int n, double **A, double **B, double tol);
void mat_print (int m, int n, double **A, char* format);