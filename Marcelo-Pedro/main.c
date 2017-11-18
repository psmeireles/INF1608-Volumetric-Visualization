/*  Pedro Sousa Meireles		- 1510962
	Marcelo Costalonga Cardoso  - 1421229 */

#include "matrix.h"
#include "simpson.h"
#include "interp.h"

#define CTSIZE 6488064 // 256*256*99
#define NX 128
#define NY 256
#define NZ 99


double transfer_function (int i, int j, int k, unsigned char *data){
	if(data[k*NY*NX + j*NX + i] < 0.3)
		return 0;
	else
		return 0.05*(data[k*NY*NX + (int)floor((double)j)*NX + i]/255. - 0.3);
}

double intensity_function (int i, int j, int k, unsigned char *data){
	double v;
	DoubleSimpson((double) 0, (double) j, i, k, data, transfer_function, &v);
	return transfer_function(i, j, k, data)*exp(-v);
}

double intensity (int i, int k, unsigned char *data){
	int L = NY - 1;
	double v, w;
	DoubleSimpson((double) 0, (double) L, 2*i, k, data, intensity_function, &v);
	DoubleSimpson((double) 0, (double) L, 2*i + 1, k, data, intensity_function, &w);
	return (v + w)/2;
}

void generate_pgm(FILE *file, double *result){
	int i;
	double maior = 0;

	for(i = 0; i < NX*NZ; i++){
		if(maior < result[i])
			maior = result[i];
	}

	fprintf(file, "P2\n%d %d\n%lf\n", NX, NZ, maior);

	for(i = 0; i < NX*NZ; i++){
		fprintf(file, "%lf ", result[i]);
		if(i % NZ == 0)
			fprintf(file, "\n");
	}
}

int main(){
	unsigned char *CTscan = (unsigned char*) malloc(CTSIZE*sizeof(char));
	double *result = (double*) malloc(NX*NZ*sizeof(double));
	FILE *file, *image;
	int i, k;
	double maior = 0;

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

	printf("head-8bit.raw reading finished\n");

	for(k = 0; k < NZ; k++){
		for(i = 0; i < NX; i++){
			result[k*NX + i] = intensity(i, k, CTscan);
		}
	}

	printf("Generating pgm...\n");

	generate_pgm(image, result);

	fclose(file);
	fclose(image);
	return 0;
}