#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#define radToDec 57.296

double complex z1,z2;

double complex cCast(double x,double y);
void printComplex(double w);

double complex summ(double complex a, double complex b);
double complex sub(double complex a, double complex b);
double complex prod(double complex a, double complex b);
double complex divi(double complex a, double complex b);

double* linspace(double a, double b, int N);
double* ABS(double* f, double ft, int N);
double* ARG(double* f, double ft, int N);

double ft;
int D;
int nElements = 1000;

int main(){
    printf("Insert z1 real part:\t\t");
	scanf("%lf",&x1);
	printf("Insert z1 immaginary part:\t");
	scanf("%lf",&y1);
	printf("Insert z2 real part:\t\t");
	scanf("%lf",&x2);
	printf("Insert z2 immaginary part:\t");
	scanf("%lf",&y2);

	z1= cCast(x1,y1);
    z2= cCast(x2,y2);

    printf("\nz1= %.3f%+.3fi\n", creal(z1), cimag(z1));
    printf("Module: %.3f,Argument %.3f\n", cabs(z1), carg(z1));
    printf("\nz2: %.3f%+.3fi\n", creal(z2), cimag(z2));
    printf("Module: %.3f,Argument %.3f\n", cabs(z2), carg(z2));
    printf("\nSumm: %f%+fi \nsub: %f%+fi \nprod: %f%+fi \ndiv: %f%+fi\n",
           summ(z1,z2),
           sub(z1,z2),
           prod(z1,z2),
           divi(z1,z2));

	printf("\nRC filter Argument and phase plot\nInsert cut-off frequency[hz]:");
	scanf("%lf",&ft);
	printf("Inserted frequency: %.3f Hz \n\nInsert periods number:",ft);
	scanf("%d",&D);


    double* f =   linspace(0.0,ft*D, nElements);
    double* Abs = ABS(f,ft, nElements);
    double* Arg = ARG(f,ft, nElements);
    printFile(f, Abs, Arg,nElements);

    printf("Done!!\n");

	return 0;

}
double complex cCast(double x,double y){
    return x+y*I;
}
void printComplex(double w){
    printf("%.3f+%.3fi",creal(w), cimag(w));
}
/////////////////OPERATIONS/////////////////////////////////


double complex summ(double complex a, double complex b){
	return a+b;
}
double complex sub(double complex a, double complex b){
	return a-b;
}
double complex prod(double complex a, double complex b){
    return a*b;
}
double complex divi(double complex a, double complex b){
	return a/b;
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

double* ABS(double* f, double ft, int N){
	double* vector= (double*)malloc(N*sizeof(double));
	printf("\nModule:");
	for (int i = 0; i < N; i++)
	{
		vector[i]=1.00/cabs(1+I*f[i]/ft);
		printf("%.4f\n",vector[i]);
	}
	//printf("ABS OK\n");//debug
	return vector;
}
double* ARG(double* f, double ft, int N){
	double* vector= (double*)malloc(N*sizeof(double));
	printf("\nPhase:");
	for (int i = 0; i < N; i++)
	{
		vector[i]= -carg(1+I*f[i]/ft)*radToDec;
		printf("%.4f\n",vector[i]);
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

	file = fopen("fdtRC.m", "w");

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
