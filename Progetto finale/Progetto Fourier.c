#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define twoPi 6.28318
#define pi 3.14159

//void readFile();

double* linspace(double t0, double tf, int N);

double* FourierSynthesizer(double* t, int N, double f0, double A,int order);

void printFile(double * string1, double * string2 ,int N);

int samplesPerPeriod = 1000;
double f0,A;
int nSample,O;
double dt,tau,T;

int main()
{
	printf("Fourier Synthesizer\nReading file...\n ");

	//readFile();
	f0 = 1000.00;
	O = 50;
	A = 1.0;

	printf("Signal frequency:\t%lf Hz\n", f0);

	printf("%d orders\n", O);

	if(f0==0.0 || O==0 || A==0.0) {
        printf("unaccepted values");
        return 0;
    }

	T=2000.0/f0; //esprimo in ms così da aumentarne la precisione

	nSample= samplesPerPeriod;

	dt = T/samplesPerPeriod;

	printf("Period %lf\n",T);//DEBUG
	printf("N samples %d\n",nSample);//DEBUG
	printf("dt: %lf",dt);

	double* t = linspace(0, T ,nSample);

	double* x = FourierSynthesizer(t, nSample, (f0/1000.0), A, O);

	printFile(t, x ,nSample);
    printf("return");
	return 0;
}


///////////////////////////////////////////////////////////////////////////
double* linspace(double ta, double tb, int N){
	double * vector = (double *) malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {
		vector[i] = ta + i * ((tb - ta) / N );
	}
	return vector;
	printf("linS Ok");//DEBUG
}



double* FourierSynthesizer(double* t, int N, double f0, double A,int order){

	double* vector= (double*) malloc(N * sizeof(double));
	double omega=pi*f0;

	for (int i = 0; i < N; i++)
	{
	    for(int n=1; n<order+1; n++){

            vector[i]+=((A/(n))*(cos( (n)*omega*t[i] ) +sin( (n)*omega*t[i] )) );
	    }
	    //printf("n:%d  val:%lf\n",i,vector[i]);

	}
	printf("Fourier Ok\n");//DEBUG

	return vector;
}

//void readFile(char*)

void printFile(double* string1, double* string2, int N){

	FILE *file;

	char template[3][50]={"close all\nt = [",
						  "];\nx = [",
						  "];\nplot(t,x);"};

	file = fopen("signal.m", "w");
	if(file== NULL){
		printf("nope\n");
		exit(EXIT_FAILURE);
	}

	fprintf(file, "%s", template[0]);

	for (int i = 0; i < N; i++)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[1]);

	for (int i = 0; i < N; i++)
		fprintf(file, "%lf ", string2[i]);

	fprintf(file, "%s", template[2]);

	fclose(file);

	printf("done!\n");

}

/*double* FourierSynthesizer(double* t, int N, double f0, double A,int order){

	double* vector= (double*) malloc(N * sizeof(double));
	double omega=pi*f0;

	printf("Fourier: \nOmega %lf\n",omega);

	for (int i = 0; i < N; i++)
	{
	    for(int n=1; n<order+1; n++){
            vector[i]=vector[i]+(A/(2*n-1)* sin( (2*n-1)*omega*t[i] ) );
	    }
	    //printf("n:%d  val:%lf\n",i,vector[i]);

	}
	printf("Fourier Ok\n");//DEBUG

	return vector;
}*/
