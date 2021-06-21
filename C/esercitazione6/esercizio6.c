#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define twoPi 6.28318

double* linspace(double t0, double tf, int N);
double* signalGen(double* t, int N, double f0, double A);
double* filtroRC(double* x, double dt, double tau, int N);
void printFile(double * string1, double * string2, double * string3,int N);
doubl*e squareGen(double * timeLine, int nSamplesPerPeriod,int nPeriods, double dutyCycle);

const int samplesPerPeriod=100; // numero campioni per periodo
double f0,A,fc;
int D,nSample;
double dt,tau,T;

int main()
{
	printf("Numerical RC filter\nInsert signal frequency [Hz]: ");
	scanf("%lf", &f0);
	printf("frequency: %lf Hz\nInsert period number:", f0);
	scanf("%d",&D);
	printf("%d periods;\nInsert signal amplitude: ", D);
	scanf("%lf",&A);
	printf("Amplitude : %lf\nInsert filter cutoff frequency[Hz]:" , A);
	scanf("%lf",&fc);
	printf("Cutoff frequency: %lf [Hz]\n", fc);

	if(f0==0.0 || fc==0.0) {
            printf("0 Hz unaccepted val");
            return 0;
    }

	T=1000.0/f0; //[ms]esprimo in ms così da aumentare la precisione
	tau= 1000.0/(twoPi*fc);//[ms]
	nSample= D*samplesPerPeriod;//[int]
	dt=T/samplesPerPeriod; //[ms]

	double* t= linspace(0.00, (T*D) ,nSample);

	double* x= signalGen(t, nSample, (f0/1000), A);

	//double* x= squareGen(t, samplesPerPeriod, D, 0.5);

	double* y= filtroRC(x, dt, tau, nSample);

	printFile(t, x, y ,nSample);

	printf("done!\n");

	return 0;
}


///////////////////////////////////////////////////////////////////////////
double* linspace(double ta, double tb, int N){
	if (N < 1) {
		return 0;
	}
	else {
            //printf("linspace N:%d",N);//DEBUG
		double * vector = (double *) malloc(N * sizeof(double));

		for (int i = 0; i < N; i++) {
			vector[i] = ta + i * ((tb - ta) / N );
		}
		return vector;
	}
}

double* signalGen(double* t, int N, double f0, double A){

	double* vector= (double*) malloc(N * sizeof(double));
	double omega= twoPi *f0;

	for (int i = 0; i < N; i++)
	{
		vector[i]=A*( cos(omega * t[i]) + cos(10*omega* t[i]) );
		//printf("signal val %.3f\n",vector[i]);//DEBUG
	}

	return vector;
}

double* filtroRC(double* x, double dt, double tau, int N){

	double* vector= (double*) malloc((N-1) * sizeof(double));

	double K= tau/dt;
	vector[0]= (x[0]/(K+1));
	for (int i = 1; i < N-1 ++i)
	{
		vector[i]= (x[i]/(K+1)) + ((K*vector[i-1])/(K+1));
		//printf("element: %d,filter val %.3f\n",i,vector[i]);//DEBUG
	}
	return vector;
}

void printFile(double * string1, double * string2, double * string3,int N){
	FILE *file;
	char template[5][50]={"t = [",
						  "];\nx = [",
						  "];\nt2 = [",
						  "];\ny = [",
						  "];\nplot(t,x);\nhold on, plot(t2,y);",
	};

	file = fopen("filtratoRC.m", "w");//w lettura

	fprintf(file, "%s", template[0]);

	for (int i = 0; i < N; ++i)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[1]);

	for (int i = 0; i < N; ++i)
		fprintf(file, "%lf ", string2[i]);

	fprintf(file, "%s", template[2]);

	for (int i = 0; i < N-1; ++i)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[3]);

	for (int i = 0; i < N-1; ++i)
		fprintf(file, "%lf ", string3[i]);

	fprintf(file, "%s", template[4]);

	fclose(file);

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
				square[i + (k * nSamplesPerPeriod)] = 0;
		}
	}
	return square;
}
