//Derivatore numerico
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define pi 3.14159


//prototipi
double * linspace(double t0, double tf, int N); //tf-t0 finestra temp, N,n samples
double * cosGen(int dt, int Periods, double amplitude, double * timeLine, double frequency);
double * derivator(int N, double dt, double * value_array_in);
void printFile(double * string1, double * string2, double * string3);

int f0,D;
double dt,T,A;
int dtSpace=20;
int nSamples;
double offset;

int main(){
	printf("Insert cosine frequency[Hz]: ");
	scanf("%d",&f0);
	printf("Inserted cosine frequency: %d\n\nInsert period number: ",f0);
	scanf("%d",&D);
	printf("%d cosine periods,\n\ninsert amplitude: ",D);
	scanf("%lf",&A);
	printf("cosine amplitude: %lf\n\nInsert time offset[ms] ",A);
	scanf("%lf",&offset);
    printf("Insert time offset %lf ms\n",offset);

	T=1000.00/f0; //periodo in [ms]
	nSamples= D*dtSpace;

	printf("Period %f\n", T);
	printf("nSamples %d\n", nSamples);
	//printf("space: %d\n", );

	double * lin_space_vector= linspace(offset, D*T, nSamples);

	//for(int o=0; o<nSamples; o++) printf("%lf ", lin_space_vector[o]); //debug stampo t

	printf("\n\n\n");

	double * sine_Vector = cosGen(dtSpace, D, A,lin_space_vector, f0);

    //for(int o=0; o<nSamples; o++) printf("%lf ", sine_Vector[o]);//debug stampo x

    //printf("before\n");//debug
    double * derivate_Vector = derivator(nSamples, (T/dtSpace*1000), sine_Vector);
    //printf("after\n");//debug

	//for(int i=0; i< nSamples-1; i++)printf("%lf ", derivate_Vector[i]); //debug stampo x'

	printFile(lin_space_vector, sine_Vector, derivate_Vector);


	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////


double* linspace(double ta, double tb, int N){
	if (N <= 1) {     ///non facciamo esplodere il programma
		return 0;
	} else {
		double * vector = (double *) malloc(N * sizeof(double));

		//double c= (tb-ta) / (N-1) ;

		for (int i = 0; i < N; i++) {
			vector[i] = ta + i * ((tb - ta) / N);
		}
		return vector;
	}

}


double* cosGen(int dt, int Periods, double amplitude, double * timeLine, double frequency){
	double* cosine= (double*) malloc(dt * Periods * sizeof(double));
	printf("sineGen elements: %d\n", dt * Periods);

	for (int k = 0; k < Periods; k++) {
		for (int i = 0; i < dt; ++i) {
			cosine[i + (k * dt)] = amplitude * cos(2 * pi * (frequency / 1000) * timeLine[i + (k * dt)]);
		}
	}
	return cosine;
}


double* derivator(int N, double dt, double* value_array_in) {
	double* derivate_array = (double*) malloc(N * sizeof(double));

	for (int i = 0; i < N-1; ++i) {
		derivate_array[i] = (value_array_in[i + 1] - value_array_in[i]) / dt;
	}
	derivate_array[N - 1] = 0.000000;
	return derivate_array;
}

void printFile(double * string1, double * string2, double * string3){
	FILE *file;
	char template[4][50]={"t = [",
						  "];\nx = [",
						  "];\nxdev = [",
						  "];\nplot(t, 100*x);\nhold on, plot(t, 10*xdev);",
	};

	file = fopen("coseno.m", "w");//w lettura

	fprintf(file, "%s", template[0]);

	for (int i = 0; i < nSamples; ++i)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[1]);

	for (int i = 0; i < nSamples; ++i)
		fprintf(file, "%lf ", string2[i]);

	fprintf(file, "%s", template[2]);

	for (int i = 0; i < nSamples; ++i)
		fprintf(file, "%lf ", string3[i]);

	fprintf(file, "%s", template[3]);

	fclose(file);

	printf("done!\n");

	free(string1);
	free(string2);
	free(string3);

}
