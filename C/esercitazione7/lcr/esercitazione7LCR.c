#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#define radToDec 57.296
#define twoPi 6.283185

double complex cCast(double x,double y);
void printComplex(double w);

double* linspace(double a, double b, int N);
double* ABS(double* f, double f0, double Q, int N);
double* ARG(double* f, double f0, double Q, int N);

double f0,R,L,C,Q;
int D;
int nElements = 10000;

int main(){
	printf("\nRLC filter Argument and phase plot\nInsert R value [Ohm]:");
	scanf("%lf",&R);
	printf("\nInsert L value [mH]:");
	scanf("%lf",&L);
	printf("\nInsert C value [uF]:");
	scanf("%lf",&C);
	printf("\nInsert plot period number:");
	scanf("%d",&D);

	L = L/1000.0;
	C = C/1000000.0;
    //printf("L:%lf  %lf")
	f0= 1.00/(twoPi* sqrt(L*C));
	Q= twoPi*f0*L/R;
	printf("\nf0=%lf ; Q=%lf\n", f0,Q);

    double* f =   linspace(0.00, (f0*D) , nElements);
    double* Abs = ABS(f,f0,Q, nElements);
    double* Arg = ARG(f,f0,Q, nElements);
    printFile(f, Abs, Arg, nElements);

    printf("Done!!\n");

	return 0;

}
double complex cCast(double x,double y){
    return x+y*I;
}
void printComplex(double w){
    printf("%.3f+%.3fi",creal(w), cimag(w));
}


double* linspace(double a, double b, int N){
	if (N < 1) {
		return 0;
	}
	else{
            printf("linspace N:%d",N);//DEBUG
		double * vector = (double*) malloc(N * sizeof(double));

		for (int i = 0; i < N; i++) {
			vector[i] = a + i * (b - a) / N;
		}
		return vector;
	}
	//printf("linspace OK\n");//debug
}

double* ABS(double* f, double f0, double Q, int N){
	double* vector= (double*)malloc(N*sizeof(double));
	printf("\nModule:");
	vector[0]= 0.0;
	for (int i = 1; i < N; i++)
	{
		vector[i]=1.00/cabs( 1+ I*Q*( (f[i]/f0)-(f0/f[i]) ) );
		//printf("%lf\n",vector[i]);
	}
	//printf("ABS OK\n");//debug
	return vector;
}
double* ARG(double* f, double f0, double Q, int N){
	double* vector= (double*)malloc(N*sizeof(double));
	printf("\nPhase:");
	vector[0]= 90.0;

	for (int i = 1; i < N; i++)
	{
		vector[i]= -carg(1.0+Q*( (f[i]/f0)-(f0/f[i]) )*I )*radToDec;
		//printf("%lf\n",vector[i]);
	}
	//printf("ARG OK\n");//debug
	return vector;
}

void printFile(double * string1, double * string2, double * string3,int N){
	FILE *file;
	char template[4][50]={"f = [",
						  "];\nModulo = [",
						  "];\nFase = [",
						  "];\nloglog(f,Modulo);\nfigure, semilogx(f,Fase);",
	};

	file = fopen("fdtLCR.m", "w");//w lettura

	fprintf(file, "%s", template[0]);

	for (int i = 0; i < N; ++i)
		fprintf(file, "%lf ", string1[i]);

	fprintf(file, "%s", template[1]);

	for (int i = 0; i < N; ++i)
		fprintf(file, "%lf ", string2[i]);

	fprintf(file, "%s", template[2]);

	for (int i = 0; i < N; ++i)
		fprintf(file, "%lf ", string3[i]);

	fprintf(file, "%s", template[3]);

	fclose(file);


}
