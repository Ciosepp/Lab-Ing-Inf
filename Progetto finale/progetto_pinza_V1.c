#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define pi 3.14159

double* linspace(double, double,  int);
double integer(double *, double *, double, int);
double* onda_quadra(int , int);
double* dente_sega(int , int);
double* triangolare(int, int);

int main()
{
    FILE *t;
    double A;
    double T, f, dt;
    int i, k, nSample, N, choice;
    double *asse_x, *funzione, *cos_t, *sin_t, *Cn_real, *Cn_imag;

    printf("Inserisci il tipo di onda che vuoi creare:\n");
    printf("1 per onda quadra\n");
    printf("2 per onda dente di sega\n");
    printf("3 per onda triangolare\n");

    scanf("%d", &choice);

    printf("Inserisci ampiezza dell'onda\n");
    scanf("%lf", &A);

    printf("Inserisci frequenza dell'onda\n");
    scanf("%lf", &f);

    printf("Inserisci numero di samples\n");
    scanf("%d", &nSample);

    T=1/f;
    dt=T/nSample;
    printf("T = %lf\n", T);
    printf("dt = %lf\n", dt);

    asse_x = linspace(0, T, nSample);
    cos_t = (double*)malloc((nSample+1)*sizeof(double));
    sin_t = (double*)malloc((nSample+1)*sizeof(double));

    if(choice==1){
        funzione = onda_quadra (A, nSample);
    }//creazione onda quadra
    if(choice==2){
        funzione = dente_sega (A, nSample);
    }//creazione onda dente di sega
    if(choice==3){
        funzione = triangolare (A, nSample);
    }//creazione onda triangolare


    //inizio creazione text.m
    t=fopen("text.m", "w");
    if (t==NULL){
      printf("text.m non aperto\n");
      exit(EXIT_FAILURE);
    } //verifica apertura
    fprintf(t, "x = [");
    for(i=0;i<nSample+1;i++){
        fprintf(t, "%lf ", asse_x[i]);
    }
    fprintf(t, "];\n\n");
    fprintf(t,"y = [");
    for(i=0;i<nSample+1;i++){
        fprintf(t, "%lf ", funzione[i]);
    }
    fprintf(t, "];\n\n");
    fprintf(t,"plot(x, y);\n\n");
    fclose(t);
    //fine creazione text.m


    //Inizio calcolo coefficienti di fourier
    printf("Quante armoniche consideriamo?\n");
    scanf("%d", &N);
    printf("\n");

    Cn_real = (double*)malloc(N*sizeof(double));
    Cn_imag = (double*)malloc(N*sizeof(double));

    for(i=0; i<N; i++){
        for(k=0;k<nSample+1;k++){
            cos_t[k] = cos(-2*pi*i*asse_x[k]/T);
            sin_t[k] = sin(-2*pi*i*asse_x[k]/T);
        }
        Cn_real[i] = 1/T * integer(funzione, cos_t, dt, nSample);
        Cn_imag[i] = 1/T * integer(funzione, sin_t, dt, nSample);
    }     //calcolo dei Cn

    for(i=0; i<N; i++){
        printf("%d \t %lf + %lfi\n", i, Cn_real[i], Cn_imag[i]);
    }     //stampa in programma dei Cn

    //inizio creazione coeff.dat
    t=fopen("coeff.dat", "w");
    if (t==NULL){
      printf("coeff.dat non aperto\n");
      exit(EXIT_FAILURE);
    } //verifica apertura

    fprintf(t, "%lf\n", f);
    fprintf(t, "%d\n", nSample);
    fprintf(t, "%d\n", N);
    for(i=0;i<N;i++){
        fprintf(t, "%lf \n %lf\n", Cn_real[i], Cn_imag[i]);
    }
    fclose(t);
    //fine creazione coeff.dat

    return 0;
}

double* linspace(double a, double b, int N){
  double *t;
  double dt;
  int i;
  dt=(b-a)/N;
  t=(double*)malloc((N+1)*sizeof(double));

  t[0]=a;
  for (i=1;i<N+1;i++){
      t[i]=dt+t[i-1];
    }
  return t;
}

double integer (double *x, double *y, double dt, int N){
  int i;
  double c=0, a, b;
  for (i=0;i<N;i++){
    a = x[i+1] * y[i+1];
    b = x[i] * y[i];
    c = 0.5*dt*(a+b)+c;
    }
  return c;
}

double* onda_quadra(int A, int nSample){
    int i;
    double *y;
    y = (double*)malloc((nSample+1)*sizeof(double));
    for(i=0;i<nSample+1;i++){
        if(i<nSample/2){
            y[i] = A;
        }
        else{
            y[i] = -A;
        }
    }
    return y;
}

double* dente_sega(int A, int nSample){
    int i;
    double *y, a=A, n=nSample;

    y = (double*)malloc((nSample+1)*sizeof(double));
        for(i=0;i<nSample+1;i++){
            y[i] = ((2*a/n)*i) - A;
        }
    return y;
}

double* triangolare(int A, int nSample){
    int i;
    double *y, a=A, n=nSample;
    y = (double*)malloc((nSample+1)*sizeof(double));

    for(i=0;i<nSample+1;i++){
            if(i<=nSample/2){
                y[i] = A - ((4*a/n)*i);
            }
            else{
                y[i] = ((4*a/n)*i) - 3*A;
            }
    }
    return y;

}
