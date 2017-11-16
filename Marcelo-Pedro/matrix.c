#include "matrix.h"

double ** mat_create (int m, int n){
	double ** matrix;
	int i;
	matrix = (double **) malloc(m * sizeof( double * ));
	if(matrix == NULL){
		printf("Lack of memory when allocating matrix\n");
		exit(1);
	}
	for(i = 0; i < m; i++){
		matrix[i] = (double *) malloc(n * sizeof( double ));
		if(matrix[i] == NULL){
			printf("Lack of memory when allocating matrix\n");
			exit(1);
		}
	}
	return matrix;
}

void mat_free (int m, double **A){
	int i;
	for(i = 0; i < m; i++)
		free( A[i] );
	free( A );
}

void mat_transpose (int m, int n, double **A, double **T){
	int i, j;
	for(j = 0; j < m; j++)
		for(i = 0; i < n; i++)
			T[i][j] = A[j][i];
}

void mat_multv (int m, int n, double **A, double *v, double *w){
	int i, j;
	for(i = 0; i < m; i++)
		w[i] = 0;

	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			w[i] += A[i][j]*v[j];
}

void mat_multm (int m, int n, int q, double **A, double **B, double **C){
	int i, j, k;
	for(i = 0; i < m; i++)
		for(j = 0; j < q; j++)
			C[i][j] = 0;

	for(i = 0; i < m; i++)
		for(k = 0; k < q; k++)
			for(j = 0; j < n; j++)
				C[i][k] += A[i][j]*B[j][k];
}

int mat_equals (int m, int n, double **A, double **B, double tol){
	int i, j;
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)
			if( tol < fabs(A[i][j] - B[i][j]) )
				return 1;
	return 0;
}

void mat_print (int m, int n, double **A, char* format){
	int i, j;
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			printf(format, A[i][j]);
			printf(" ");
		}
		printf("\n");
	}
}