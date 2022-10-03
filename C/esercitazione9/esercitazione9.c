#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double* linspace(double, double,  int);
double randn(double, double);
double rande(double);
double* histocounts(double*, int, double*, int);

FILE *f;
FILE *g;
double mu=2, sigma=5;
int N, i, Nbin;
double *a, *b, *mean, *binVector, *cont, sum=0, sum2=0;

int debug=1;

int main()
{

	printf("Inserisci il valore di N\t");
	scanf("%d", &N);

	printf("Inserisci il numero di bin\t");
	scanf("%d", &Nbin);
	printf("\n\n");

	a = (double*)malloc(N*sizeof(double));
	b = (double*)malloc(N*sizeof(double));
	mean = (double*)malloc(Nbin*sizeof(double));

	for(i=0; i<N; i++){

		a[i] = randn(mu, sigma);
		b[i] = rande(mu);
    }

    binVector = linspace((mu - 3.0*sigma), (mu + 3.0*sigma), Nbin);

    cont = histocounts(a, N, binVector, Nbin);

    printf("\n\n------------ Guass ------------\n");

    for(i=0; i<Nbin-1; i++){

    	mean[i] = (binVector[i] + binVector[i+1]) / 2;
        sum += cont[i];
    }

    printf("gauss sum: %f\n",sum);

    for(i=0; i<Nbin-1; i++){

        cont[i] /= sum;
        sum2 += cont[i];
    }
    printf("sum2: %f\n\n", sum2);

    f=fopen("pdfGauss.m", "w");

    if (f==NULL){
    	printf("pdfGauss.m non aperto\n");
    	exit(EXIT_FAILURE);
    }

    fprintf(f,"x = [");

    for(i=0;i<Nbin-1;i++){

    	fprintf(f, "%f, ", mean[i]);
    }

    fprintf(f, "];\n\nf = [");

    for(i=0;i<Nbin-1;i++){
    	fprintf(f, "%f ", cont[i]);
    }

    fprintf(f, "];\n\nfigure;\n\nplot(x, f);print(gcf, 'gauss', '-dpng', '-r300');");
    fclose(f);


    printf("\n\n------------ Exp ------------\n");

    sum=0;
    sum2=0;

    binVector = linspace(0.0, (5.0*mu), Nbin);

    cont = histocounts(b, N, binVector, Nbin);

    for(i=0; i<Nbin-1; i++){

    	mean[i] = (binVector[i] + binVector[i+1]) / 2;
    	sum += cont[i];
    }
    printf("exp sum: %f\n",sum);

    for(i=0; i<Nbin-1; i++){

        .cont[i] /= sum;
    	sum2 += cont[i];
    }
    printf("\n%f \n\n", sum2);

    g=fopen("pdfExp.m", "w");

    if (g==NULL){
    	printf("pdfExp.m non aperto\n");
    	exit(EXIT_FAILURE);
    }
    fprintf(g,"x = [");

    for(i=0;i<Nbin-1;i++){

    	fprintf(g, "%f ", mean[i]);
    }
    fprintf(g, "];\n\nf = [");

    for(i=0;i<Nbin-1;i++){

    	fprintf(g, "%f ", cont[i]);
    }

    fprintf(g,"];\n\nfigure;\n\nplot(x, f);print(gcf, 'exp', '-dpng', '-r300');");

    fclose(g);

    free (a);
    free (b);
    free (mean);
    free (cont);

    return 0;
}

double randn(double mean, double stdDev) {

    static double spare;
    static int hasSpare = 0;

    if (hasSpare) {
        hasSpare = 0;
        return spare * stdDev + mean;
    } else {
        double u, v, s;
        do {
            u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        s = sqrt(-2.0 * log(s) / s);
        spare = v * s;
        hasSpare = 1;
        return mean + stdDev * u * s;
    }
}


double* linspace(double a, double b, int Nbin)
{
    double dn;
	dn=(b-a)/(Nbin-1);
	double *x =(double*)malloc((Nbin+1)*sizeof(double));

	for(int i=0; i<Nbin; i++){

		x[i] = (double)dn*i+a;
	}
	return x;
}

double rande(double mu)
{
	double x,y;
	srand(clock());
	y = (double)rand() / (RAND_MAX+1.0);
	x = -mu * log(1.0-y);
	return x;
}

double* histocounts(double* x, int dimx, double* binVector, int Nintv)
{
	int i,k;
	double *cont= (double*)malloc(Nintv *sizeof(double));

	for(k=0; k < Nintv-1; k++)cont[k]=0;

	for(i=0; i < dimx; i++){

		for(k=0; k < Nintv-1; k++){


            if(x[i] >= binVector[k] && x[i] <binVector[k+1]){

                cont[k]++;
                break;
            }
		}
    }

    return cont;
}
