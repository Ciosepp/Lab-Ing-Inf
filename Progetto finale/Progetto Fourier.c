#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define twoPi 6.28318
#define pi 3.14159

void readFile();

double* linspace(double t0, double tf, int N);

double* FourierSynthesizer(double* t, int N, double f0, double A,int order);

void printFile(double * string1, double * string2, int N);

int samplesPerPeriod = 1000;
double f0,A;
int D,nSample,O;
double dt,tau,T;

int main()
{
	printf("Fourier Synthesizer\nReading file...\n ");

	readFile();

	printf("Signal frequency:\t%lf Hz\n", f0);

	printf("%d orders\n", O);

	printf("Insert number periods:");
	scanf("%d",D);


	if(f0==0.0 || O==0 || A==0.0) {
        printf("unaccepted values");
        return 0;
    }

    D = 1;
	T=1000.0/f0; //esprimo in ms così da aumentare la precisione

	nSample= D*samplesPerPeriod;
	dt=T/samplesPerPeriod;

	printf("Period %lf\n",T);//DEBUG
	printf("Tau %lf\n",tau);//DEBUG
	printf("N samples %d\n",nSample);//DEBUG
	printf("dt: %lf",dt);

	double* t= linspace(0, (T*D) ,nSample);

	double* x= FourierSynthesizer(t, nSample, (f0/1000.0), A, J);

	printFile(t, x ,nSample);
    printf("fine");
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

            vector[i]+=((A/(2*n-1))* sin( (2*n-1)*omega*t[i] ));
	    }
	    //printf("n:%d  val:%lf\n",i,vector[i]);

	}
	printf("Fourier Ok\n");//DEBUG

	return vector;
}

void readFile(char*)

void printFile(double * string1, double * string2, double * string3,int N){

	FILE *file;

	char template[5][50]={"t = [",
						  "];\nx = [",
						  "];\nt2 = [",
						  "];\ny = [",
						  "];\nplot(t,x);\nhold on, plot(t2,y);",
	};

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

	for (int i = 0; i < N-1; i++)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[3]);

	for (int i = 0; i < N-1; i++)
		fprintf(file, "%lf ", string3[i]);

	fprintf(file, "%s", template[4]);

	fclose(file);

	printf("done!\n");

	//free(string1);
	//free(string2);
	//free(string3);

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
