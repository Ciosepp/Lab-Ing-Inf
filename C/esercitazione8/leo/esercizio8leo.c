#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.14159

    double* linspace(double fs, int N){
      double *t;
      double dt;
      int i;
      dt=1/fs;
      t=(double*)malloc((N+1)*sizeof(double));
      t[0]=0;
      for (i=1;i<N;i++){
          t[i]=dt+t[i-1];
        }
      return t;
    }


    double* coseno(double A, double fo, double *t, int N){
      double *x;
      int i;
      x=(double*)malloc((N+1)*sizeof(double));
      for (i=0;i<N;i++){
          x[i]=(A*cos(2*pi*fo*t[i]))+(A*cos(2*pi*3*fo*t[i]));
        }
      return x;
    }


double FIR (double x, double* dlx, double* h, unsigned M){
    int k;
    double y, c=0;
    for  (k=0; M-1-k>0; k++){
       dlx[M-2-k]= dlx[M-1-k];
    }
    dlx[0]= x;
    for(k=0; k<M; k++){
        y= h[k]*dlx[k]+c;
        c=y;
    }

    return y;
    }

    /*double FIR (double *x, int n, double* dlx, double* h, unsigned M){
    int i, k;
    double y;
    for  (i=0; i<M; i++){
        if(n-i<0) break;
       dlx[i]= x[n-i];
    }
     double c=0;
    for(k=0; k<M; k++){
        y= h[k]*dlx[k]+c;
        c=y;
    }

    return y;
    }*/


    int main(void)
    {
        FILE *f;
        double fs,A,fo;
        int M,i,D;
        double *h;

        f=fopen("coeff.dat", "r");

        if (f==NULL){ //verifica apertura
          printf("coeff.dat non aperto\n");
          exit(EXIT_FAILURE);
        }

        fscanf(f, "%lf", &fs);
        //printf("%lf\n", fs);

        fscanf(f, "%d", &M);
        //printf("%d\n", M);

        h=(double*)malloc(M*sizeof(double));
        for(i=0;i<M;i++){
            fscanf(f, "%lf\n", &h[i]);
            //printf("%lf\n", h[i]);
        }

        fclose(f);

        printf("Inserisci A\t");
        scanf("%lf", &A);
        printf("Inserisci fo\t");
        scanf("%lf", &fo);
        printf("Inserisci D\t");
        scanf("%d", &D);
        printf("\n\n");

        double dt= 1/fs;
        double *dlx;
        double *t;
        double *x;
        int N= (int) fs*D;

         t= linspace(fs, N);
         x= coseno(A, fo, t, N);

        dlx= (double*)malloc(M*sizeof(double));

        for(i=0; i<M; i++){

            dlx[i]= 0;
        }
    int n;
    double* y= (double*)malloc(N*sizeof(double));

    for(n=0; n<N; n++){
        y[n]= FIR(x[n], dlx, h, M);
    }



    FILE *doc;
    doc = fopen("filtratoFIR.m", "w");
    fprintf(doc, "t = [");
    for (i=0;i<N; i++){
        fprintf(doc, "%lf, ",t[i]);
        }
    fprintf(doc, "];\n");

    fprintf(doc, "x = [");
    for (i=0;i<N; i++){
        fprintf(doc, "%lf, ",x[i]);
        }
    fprintf(doc, "];\n");

    fprintf(doc, "t2 = [");
    for (i=0;i<N-1; i++){
        fprintf(doc, "%lf, ",t[i]);
        }
    fprintf(doc, "];\n");

    fprintf(doc, "y = [");
    for (i=0;i<N-1; i++){
        fprintf(doc, "%lf, ",y[i]);
        }
    fprintf(doc, "];\n");

    fprintf(doc, "plot(t, x);\n");
    fprintf(doc, "hold on, plot(t2, y);\n");

    fclose(doc);
    free(t);
    free(x);
    free(y);
    free(dlx);

    return 0;
    }
