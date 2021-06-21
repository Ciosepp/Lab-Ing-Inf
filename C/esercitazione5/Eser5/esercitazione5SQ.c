#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.14159

//prototipi
double * linspace(double t0, double tf, int N); //tf-t0 finestra temp, N,n samples
double * squareGen(double * timeLine, int nSamplesPerPeriod,int nPeriods, double dutyCycle);
double * trapz(double* x, double dt, int N);
void printFile(double * string1, double * string2, double * string3);

int D,nSamples,nDivisions;
double dt,T,DC;

int main(){
	printf("Insert period[ms]: ");
	scanf("%lf",&T);
	printf("Inserted period: %lf\n\nInsert period number: ",T);
	scanf("%d",&D);
	printf("%d signal periods,\n\ninsert duty cycle: ",D);
	scanf("%lf",&DC);
	printf("dutyCycle: %lf\n\nInsert sample number ",DC);
	scanf("%d",&nDivisions);
    printf("Insert time offset %lf ms\n",nDivisions);

	nSamples= D*nDivisions;

	printf("Period %f ms\n", T);
	printf("nSamples %d\n", nSamples);
	printf("nSamples %lf\n", T*D);

	double * lin_space_vector= linspace(0.0,T*D, nSamples);

	double * x = squareGen(lin_space_vector, nDivisions ,D,DC);

    double * integral_Vector = trapz(x, (T/(1000*nDivisions)),nSamples);

    printFile(lin_space_vector, x, integral_Vector);

	return 0;
}

///////////////////////////////////////////////////////////////////
double* linspace(double ta, double tb, int N){
	if (N <= 1) {     ///non facciamo esplodere il programma
		return 0;
	} else {
		double * vector = (double *) malloc(N * sizeof(double));

		for (int i = 0; i < N; i++) {
			vector[i] = ta + i * ((tb - ta) / N);
		}
		return vector;
	}
}

double * squareGen(double * timeLine, int nSamplesPerPeriod,int nPeriods, double dutyCycle){

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

	file = fopen("square.m", "w");//w lettura

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

