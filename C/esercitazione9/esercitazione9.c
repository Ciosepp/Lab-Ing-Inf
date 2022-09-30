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
        //printf("%f \n",a[i]);
		b[i] = rande(mu);
    }

    binVector = linspace((mu - 3.0*sigma), (mu + 3.0*sigma), Nbin);

    cont = histocounts(a, N, binVector, Nbin);

    if(debug == 1){

        printf("histocounts result:\n");
        for(int i=0; i< Nbin-1; i++){

            //printf("%f ",cont[i]);
        }
    }
    printf("mean:\n");
    for(i=0; i<Nbin-2; i++){

    	mean[i] = (binVector[i] + binVector[i+1]) / 2;
        sum += mean[i];
    }

    //riscrivo cont normalizandolo
    printf("Stampo vett cont normalizzato\n");

    for(i=0; i<Nbin-1; i++){
    	cont[i] /= sum;
    	//printf("%f  ", cont[i]);
    	sum2 += cont[i];
    	printf("%f  ", cont[i]);
    }

    sum=0;
    sum2=0;

    f=fopen("pdfGauss.m", "w");

    if (f==NULL){ //verifica apertura
    	printf("pdfGauss.m non aperto\n");
    	exit(EXIT_FAILURE);
    }

    //inizio creazione pdfGauss.m
    fprintf(f,"x = [");

    for(i=0;i<Nbin-2;i++){

    	fprintf(f, "%f, ", mean[i]);
    }
    fprintf(f, "];\n\n");    //e' il vettore mean

    fprintf(f,"f = [");
    for(i=0;i<Nbin-2;i++){
    	fprintf(f, "%f ", cont[i]);
    }
    fprintf(f, "];\n\n");    //e' il vettore cont normalizzato
    fprintf(f,"figure;\n\n");
    fprintf(f,"plot(x, f);\n\n");  //fine creazione pdfGauss.m
    fclose(f);

    //Exp
    printf("Exp\n");

    binVector = linspace(0.0, 5.0*mu, Nbin);

    cont = histocounts(b, N, binVector, Nbin);

    printf("Stampo vett mean\n");
    for(i=0; i<Nbin; i++){
    	mean[i] = (binVector[i] + binVector[i+1]) / 2;
    	//printf("%f  ", mean[i]);
    }

    printf("Stampo vett cont\n");

    for(i=0; i<Nbin; i++){
    	//printf("%f  ", (cont[i] / N));
    	sum += cont[i];
    }printf("\n%f \n\n", sum / N);       //stampo cont

    //riscrivo cont normalizandolo
    printf("Stampo vett cont normalizzato\n");
    for(i=0; i<Nbin; i++)
    {
    	cont[i] /= sum;

    	sum2 += cont[i];
    }printf("\n%f \n\n", sum2);       //stampo cont normalizzato

    g=fopen("pdfExp.m", "w");

    if (g==NULL){ //verifica apertura
    	printf("pdfExp.m non aperto\n");
    	exit(EXIT_FAILURE);
    }

    //inizio creazione pdfExp.m
    fprintf(g,"x = [");
    for(i=0;i<Nbin;i++){
    	fprintf(g, "%f ", mean[i]);
    }
    fprintf(g, "];\n\n");    //e' il vettore mean

    fprintf(g,"f = [");
    for(i=0;i<Nbin;i++){
    	fprintf(g, "%f ", cont[i]);
    }
    fprintf(g, "];\n\n");    //e' il vettore cont normalizzato

    fprintf(g,"figure;\n\n");

    fprintf(g,"plot(x, f);\n\n");  //fine creazione pdfExp.m

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
		//printf("l%d:%f \n",i, x[i]);
	}
	return x;
}

/*double randn(double mu, double sigma)
{
	srand(clock());
	double u1, u2, d, f;
	static double x1, x2;
	static int flag=0;
	if(flag==1)
	{
		flag = !flag;
		return (mu+sigma*(double)x2);
	}
	do
	{
		u1 = -1 + ((double)rand()/RAND_MAX)*2;
		u2 = -1 + ((double)rand()/RAND_MAX)*2;
		d = u1*u1+u2*u2;
	}
	while (d>= 1 || d==0);
	f = sqrt ((-2*log(d))/d);
	x1 = u1*f;
	x2 = u2*f;
	flag=!flag;
	return (mu+sigma*(double)x1);
}
*/
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
