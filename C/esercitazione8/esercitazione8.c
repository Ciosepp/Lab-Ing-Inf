#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define twoPi 6.283185

char fileName[30]="coeff.dat";

void fileRead(char* fileName, double* fs, double* h, int* M);
double* linspace(int N, int dt);
double* signalGen(double* t, int N, double f0, double A);
double FIR(double xn, double* dlx, double* h, unsigned M);
void printFile(double * string1, double * string2, double * string3,int N);

int M; //ordine del filtro
double fs, dt;//freqeunza campionamento -> dt linea ritardo
double f0, A; //freqeunza e ampiezza segnale da filtrare
int D; // #periodi segnale
FILE file;

int int main()
{
	printf("Numerical FIR\n");
	printf("Reading %s file...",fileName);
	double* fs, h;
	int M; //fir dimension
	fileRead("coeff.dat", fs ,h, M);

	printf("\nInsert signal frequency [Hz]: ");
	scanf("%lf", &f0);
	if(f0==0.0) {
            printf("0 Hz unaccepted val");
            return 0;
    }
	printf("frequency: %lf Hz\nInsert periods' number:", f0);
	scanf("%d",&D);
	printf("%d periods;\nInsert signal amplitude: ", D);
	scanf("%lf",&A);
	printf("Amplitude : %lf\n", A);

	double Ts= 1000/ fs; //fe=requenza campionamento da file
	double T0= 1000/ f0; //frequenza del segnale  
	int nElements= D*fs/f0; //signals' number of elements

	double* t= linspace(0.00, D*T0, Ts, nElements);
	double* x= signalGen(t, nElements, f0, A);


	return 0;
}

double* linspace(int N, int dt){
	int n=D*T/dt
        //printf("linspace N:%d",N);//DEBUG
	double * vector = (double *) malloc((n+1)* sizeof(double));

	for (int i = 0; i < n; i++) {
		vector[i] = i * dt;
	}
	return vector;

}


double* signalGen(double* t, int N, double f0, double A){

	double* vector= (double*) malloc(N * sizeof(double));
	double omega= twoPi *f0;

	for (int i = 0; i < N; i++){
		vector[i]=A*( cos(omega * t[i]) + cos(3*omega* t[i]) );
		//printf("signal val %.3f\n",vector[i]);//DEBUG
	}
	return vector;
}