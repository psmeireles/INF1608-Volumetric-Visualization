/*  Pedro Sousa Meireles		- 1510962
	Marcelo Costalonga Cardoso  - 1421229 */

#include "matrix.h"
#include "simpson.h"
#include "interp.h"

#define CTSIZE 6488064 // 256*256*99
//#define NX 2*128
#define NX 128
#define NY 256
#define NZ 99

double tol;

/*
// Sem interpolação
double transfer_function (int i, int j, int k, unsigned char *data){
	unsigned char c = data[k*NY*2*NX + j*2*NX + i];
	if(c/255. < 0.3)
		return 0;
	else
		return 0.05*(c/255. - 0.3);
}

// Interpolação por reta
double transfer_function (int i, double jf, int k, unsigned char *data){
	
	int j = floor(jf);
	int pos1 = k*NY*2*NX + j*2*NX + i;
	int pos2 = k*NY*2*NX + (j+1)*2*NX + i;
	int pos = k*NY*2*NX + (jf)*2*NX + i;
	unsigned char c1 = data[pos1];
	unsigned char c2 = data[pos2];
	unsigned char c = (c1*(pos2 - pos) + c2*(pos - pos1))/(pos2 - pos1);
	
	if(c/255. < 0.3)
		return 0;
	else
		return 0.05*(c/255. - 0.3);
}
*/

double transfer_function (int i, double j, int k, unsigned char *data){

	unsigned char c;

	if( j - floor(j) > 0.5e-5){
		double x[2], *coef;

		// Interpolation if j isn't integer
		x[0] = data[k*NY*2*NX + (int)floor(j)*2*NX + i];
		x[1] = data[k*NY*2*NX + (int)floor(j+1)*2*NX + i];
		coef = NewtonCoef(2, x, i, k, data, transfer_function);
		return NewtonAval(2, x, coef, j);
	}
	else{
		c = data[k*NY*2*NX + (int)j*2*NX + i];
		if(c/255. < 0.3){
			return 0;
		}
		else
		return 0.05*(c/255. - 0.3);
	}
}

double intensity_function (int i, double j, int k, unsigned char *data){
	double v, integer = 0;
	double L = j, prev = 0;
	double h = 4.5;

	while(prev < L){
		if(h < L){
			//DoubleSimpson((double) prev, (double) h, i, k, data, transfer_function, &v);
			v = AdaptiveSimpson((double) prev, (double) h, i, k, data, transfer_function, tol);
		}
		else{
			//DoubleSimpson((double) prev, (double) L, i, k, data, transfer_function, &v);
			v = AdaptiveSimpson((double) prev, (double) L, i, k, data, transfer_function, tol);
		}
		integer += v;
		prev = h;
		h += 4.5;
	}
	return transfer_function(i, j, k, data)*exp(-integer);
}

double intensity (int i, int k, unsigned char *data){
	double L = 255, prev = 0;
	double v, w, intensity = 0;
	double h = 4.5;

	while(prev < L){
		if(h < L){
			//DoubleSimpson((double) prev, (double) h, 2*i, k, data, intensity_function, &v);
			//DoubleSimpson((double) prev, (double) h, 2*i + 1, k, data, intensity_function, &w);
			v = AdaptiveSimpson((double) prev, (double) h, 2*i, k, data, intensity_function, tol);
			w = AdaptiveSimpson((double) prev, (double) h, 2*i + 1, k, data, intensity_function, tol);
		}
		else{
			//DoubleSimpson((double) prev, (double) L, 2*i, k, data, intensity_function, &v);
			//DoubleSimpson((double) prev, (double) L, 2*i + 1, k, data, intensity_function, &w);
			v = AdaptiveSimpson((double) prev, (double) L, 2*i, k, data, intensity_function, tol);
			w = AdaptiveSimpson((double) prev, (double) L, 2*i + 1, k, data, intensity_function, tol);
		}
		intensity += (v + w)/2;
		prev = h;
		h += 4.5;
	}
	return intensity;
}

void generate_pgm(FILE *file, double *result){
	int i;
	double maior = 0;

	for(i = 0; i < NX*NZ; i++){
		if(maior < result[i])
			maior = result[i];
	}

	for(i = 0; i < NX*NZ; i++){
		result[i] = (unsigned char) (result[i]*255/maior);
	}

	fprintf(file, "P2\n%d %d\n255\n", NX, NZ);

	for(i = 0; i < NX*NZ; i++){
		unsigned char c = (unsigned char) result[i];
		fprintf(file, "%hhu ", c);
		if((i+1) % NZ == 0){
			fprintf(file, "\n");
		}
	}
}

int main(){
	unsigned char *CTscan = (unsigned char*) malloc(CTSIZE*sizeof(unsigned char)), c;
	double *result = (double*) malloc(NX*NZ*sizeof(double));
	FILE *file, *image;
	int i, k;
	double maior = 0;

	tol = 0.5e-2;

	printf("Reading head-8bit.raw\n");

	file = fopen("head-8bit.raw", "rb");
	if(file == NULL){
		printf("Error when opening file head-8bit.raw\n");
		exit(1);
	}

	image = fopen("img.pgm", "wb");
	if(image == NULL){
		printf("Error when opening file image.pgm\n");
		exit(1);
	}

	for(i = 0; i < CTSIZE; i++){
		fread(&CTscan[i], sizeof(unsigned char), 1, file);
	}

	fclose(file);

	printf("Calculating the result...\n");
	for(k = 0; k < NZ; k++){
		for(i = 0; i < NX; i++){
			result[k*NX + i] = intensity(i, k, CTscan);
		}
	}

	printf("Generating pgm...\n");

	generate_pgm(image, result);

	fclose(image);
	free(result);
	free(CTscan);
	return 0;
}