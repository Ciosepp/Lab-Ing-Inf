#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.14159

//prototipi
double * linspace(double t0, double tf, int N);
double * cosineGen(int dt, int Periods, double amplitude, double * timeLine, double frequency);
double * trapz(double* x, double dt, int N);
void printFile(double * string1, double * string2, double * string3);

int f0,D,nSamples;
double dt,T,A,offset;
int dtSpace=20;


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
	
	double * lin_space_vector= linspace(offset, D*T, nSamples);

	double * cosine_Vector = cosineGen(dtSpace, D, A,lin_space_vector, f0);

    double * integral_Vector = trapz(cosine_Vector, (T/(1000*dtSpace)),nSamples);

    printFile(lin_space_vector, cosine_Vector, integral_Vector);

	return 0;
}

///////////////////////////////////////////////////////////////////
double* linspace(double ta, double tb, int N){
	if (N <= 1) {     ///non facciamo esplodere il programma
		return 0;
	} else {
		double * vector = (double *) malloc(N * sizeof(double));

		for (int i = 0; i < N; i++) {
			vector[i] = ta + i * ((tb - ta) / (N - 1));
		}
		return vector;
	}
}

double* cosineGen(int dt, int Periods, double amplitude, double * timeLine, double frequency){

	double* cosine= (double*) malloc(dt * Periods * sizeof(double));
	printf("cosineGen elements: %d\n", dt * Periods);

	for (int k = 0; k < Periods; k++) {
		for (int i = 0; i < dt; ++i) {
			cosine[i + (k * dt)] = amplitude * cos(2 * pi * (frequency / 1000) * timeLine[i + (k * dt)] );
		}
	}
	return cosine;
}

double* trapz(double* x, double dt, int N){

	double* I = (double*)malloc(N*sizeof(double));
	I[0]=0.0;
	for (int i = 1; i < N-1; ++i)
	{
		I[i]= ( ( x[i]+x[i+1] )*dt/2 )+I[i-1];
	}
	I[N-1]=0.0;

	return I;
}
void printFile(double * string1, double * string2, double * string3){
	FILE *file;
	char template[4][50]={"t = [",
						  "];\nx = [",
						  "];\nxint = [",
						  "];\nplot(t, x);\nhold on, plot(t, 100*xint);",
	};

	file = fopen("integratore.m", "w");//w lettura

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

/*double * squareGen(double * timeLine, int nSamplesPerPeriod,int nPeriods, double dutyCycle){

	int nElements= nPeriods* nSamplesPerPeriod;

	double* square= (double*) malloc(nElements* sizeof(double));

	int onElements= nSamplesPerPeriod *dutyCycle;

	printf("squareGen elements: %d\n", nElements);

	for (int k = 0; k < nPeriods; k++) {

		for (int i = 0; i < onElements; ++i) {
			square[i + (k * nSamplesPerPeriod)] = 1;
		}
		for (int i = onElements; i < nSamplesPerPeriod; ++i)
		{
				square[i + (k * nSamplesPerPeriod)] = -1;
		}
	}
	return square;
}*/