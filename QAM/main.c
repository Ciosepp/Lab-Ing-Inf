/*
    Simulazione trasmissione rumorosa M-QAM su canale AWGN

    di Giuseppe Balducci
    versione 2.0 2023, 03, 27
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
time_t t;

FILE *f;
char *fileName = "QAM.m";

const int nPoints= 20;
const double minSNR = 1;
const double maxSNR = 100;

int N=100000;
int L=1;   //L-ASK
int M=1;   //M-QAM
int l=1;   //l bit nibble stringa
int n=1;   //n bit stringa
double sigma;
char *word;
char *An;
double *Yn;
char* Z;
char *received;

double *SNRARRAY;
double *PeSp, *PeTh;
/*
convenzioni : vettore di bit v[i]: i=0 -LSB, i=n -MSB
Segnali I Q: I - stringa LS, Q -stringa MS
*/

double *SNRarrayGen(int nPts, double min, double max){
	double *x = (double*)malloc(sizeof(double)* nPts);
	double c = (max - min)/(nPts-1 );
	x[0] = min;
	for (int i = 1; i < nPts; i++){
		x[i] = x[i-1] + c;
	}
	return x;
}

char *randGen(int nbits){           //genera una stringa di n-bit casuali
	char *o = (char*)malloc(sizeof(char)*nbits);
    //printf("Word:");            //debug
    for(int i=0; i<nbits; i++){     //genero 4 bit
        o[i] = rand()%2;        //bit random

        //printf("%d", o[i] );    //debug
    }
    //printf("\n");               //debug
    return o;
}
int string2int(char *in , int nBits){
	int x=0;
	for (int i = 0; i < nBits; i++){
		x += in[i] * (1<<i);
	}
	return x;
}
//crea mappa livelli con codifica di Grey
char *mapGenerator(int nLevs){
	char *IQMAP = (char*)malloc(sizeof(char)*nLevs);
	IQMAP[0] = -(nLevs-1);
	for(int i=1; i<nLevs; i++){
		IQMAP[i] = IQMAP[i-1]+2;
	}
	return IQMAP;
}


char *IQmapper(int x, char *map,int nLevs, int nbit){ //mappa x in cordinate cartesiane
	char *IQ = (char*)malloc(sizeof(char)*2);
    IQ[0] = map[ x & (nLevs-1)];     //I: è il mapping del valore dei 2 LSB
    IQ[1] = map[ x>>nbit & (nLevs-1)];  //Q: è il mapping del valore dei 2 MSB
    return IQ;
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

double *noiseAdder(char *Ain, double STDEV){
	double  *S = (double  *)malloc(sizeof(double )*2);
	S[0] = Ain[0] + randn(0,STDEV);
	S[1] = Ain[1] + randn(0,STDEV);
	return S;
}

//double thresholds[]= {-2.0, 0.0, 2.0};


double *thrGen(int nThr){              //genera n-1 soglie per la decisione del simbolo
	double *t = (double*)malloc(sizeof(double)*(nThr-1));
	t[0] = -(nThr-2);
	for(int i=1; i<nThr-1; i++){
		t[i] = t[i-1] +2.0;
	}
	return t;
}
char *demod(double *Yn, double *thresholds, char *Map,int h){   //demodulatore segnale QAM
	char *z = (char*)malloc(sizeof(char )*2);

	for(int i = 0; i < 2; i++){
        z[i] = Map[h-1];
		for (int j = 0; j < h-1; j++){

			if(Yn[i] <= thresholds[j]){
				z[i]= Map[j];
				break;
			}
		}
	}
	return z;
}

char *demapper(char *In, char *map,int mlev,int nbit){            //demappatore
	char  *q = (char  *)malloc(sizeof(char)*nbit *2);
	unsigned char val = 0;
    for (int i = 0; i < 2; i++){    //Componenti IQ i=0LSB
        for(int j = 0; j < mlev; j++){ //livelli
            if(In[i] == map[j]){    //comparazione livelli alfabeto
            	val += j<< i*(nbit);
            }
        }
    }
    //printf("\n%d",val);             //debug
    //printf("\n");                   //debug
    for (int i = 0; i < nbit*2; i++)
    {
    	q[i] = 0x01 & val>>i;
        //printf("%d",q[i]);          //debug
    }
    return q;
}

int a;

int Errors =0;
double Pe=0.0;
double PeS =0.0;

char *greyMap;      // mappa livelli ask
double *thresholds; // soglie del decisore

double SNR=0.0;
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
int debug =0;
char go = 0;
int main(){
	srand((unsigned) time(&t));

	printf("%d iterations, ok? y/n  :",N);
	while(go != 'y'&&go != 'n'){
		scanf("%c",&go);
	}
	if(go == 'n'){
		printf("Insert #iteration:");
		scanf("%d",&N);
	}

    if(debug == 1) N= 10;
	//printf("Insert l :");
	//scanf("%d",&l);
	/* while(SNR<=0){
			    	printf("\nInsert SNR:");
			    	scanf("%lf",&SNR);
			    	sigma = (double)sqrt( ((M-1))/(3.0*SNR))*0.707;

			    }*/
/*
    for(int i=0; i < L; i++){    //debug
        printf("%d\n", greyMap[i]);
    }
    for(int i=0; i < L-1; i++){
        printf("%f\n", thresholds[i]);
    }
*/
    //alloco vettore valori SNR
    SNRARRAY = SNRarrayGen(nPoints,minSNR, maxSNR);

    //alloco array delle Pe
    PeSp	= (double*) malloc(sizeof(double) * nPoints *3);//array Pe sperimenali
    PeTh	= (double*) malloc(sizeof(double) * nPoints *3);//array Pe teoriche
    printf("SNR steps:\n");
    for (int i = 0; i < nPoints; i++)
   	{
    	printf("%lf\n",SNRARRAY[i] );
    }

    for(int currentL =0; currentL<3; currentL++){					// varia l


		l = currentL+1;
		n = 2*l;    //numero bit stringa
	    L=1 ,M =1;	//resetto senno esplode
	    for(int i=0; i<l; i++){
	        L *=2;  //numero livelli ASK
	    }
	    //printf("\n%d",L);
	    for(int i=0; i<n; i++){
	        M *=2;  //numero livelli QAM
	    }
		printf("\n%d-QAM\n",M);

	    greyMap = mapGenerator(L);  //generazione mappa

	    thresholds = thrGen(L);     //generazione soglie

    	for (int np = 0; np < nPoints; np++){	//varia SNR

		    sigma = (double)sqrt( ((M-1))/(3.0*SNRARRAY[np]))*0.707;

			printf("\nRun:%d/%d",np+1,nPoints );
		    printf(" -Noise sigma:%.3lf",sigma);

		    Errors =0;

		    for(int z=0; z<N; z++){     //iterazioni Monte Carlo

		        // punto 1: genero un stringa random di 4 bit
		    	word = randGen(n);

		        if(debug == 1){
		            printf("\nword: ");
		            for (int i = 0; i < n; i++){
		                printf("%d",word[n-1-i]);
		            }

		        }

		            // punto 2
		        a= string2int(word,n);//converto in un valore
		        if(debug == 1)printf("\nconv:%d",a);                       //debug

		        An = IQmapper(a,greyMap,L,l);
		        if(debug == 1)printf("\nMAP: %d %d",An[0] , An[1]);   //debug

		            // punto 3
		        Yn = noiseAdder(An, sigma);
		        if(debug == 1)printf("\n+noise: %f %f",Yn[0] ,Yn[1]);  //debug

		            // punto 4: cerco il punto della bitmap più vicino al punto Y
		        double S[] = {Yn[0]-An[0], Yn[1]-An[1]};
		        double ab = sqrt(S[0]*S[0] + S[1]*S[1]);
		        if(debug == 1)printf("\nNoise Mod:%.3f",ab);    //debug


		        Z = demod(Yn, thresholds, greyMap,L);
		        if(debug == 1)printf("\nDecoded: %d %d",Z[0] ,Z[1]);  //debug

		            // punto 5: stringa ricevuta demodulata

		        received = demapper(Z, greyMap,L,l);

		        if(debug == 1){
		            printf("\nreceived: ");
		            for (int i = 0; i < n; i++){
		                printf("%d",received[n-1-i]);
		            }
		            printf("\n");
		        }


		        for (int i = 0; i < n; i++)
		        {
		        	if(word[i] != received[i]){
		        		Errors++;
		                if(debug == 1)printf("E");//debug
		        		//break;
		        	}
		        }

		        free(word);
		        free(An);
		        free(Yn);
		        free(Z);
		        free(received);
		        PeSp[np + nPoints*currentL] = 100.0*(double)Errors / (N*n);
		        PeTh[np + nPoints*currentL] = 100.0*(L-1)*erfc( sqrt(3*SNRARRAY[np]/ (2*(M-1)) ) )/(L*l);
		    }
    	}
	}
    //printf("\nProbabilita' teorica: %.4f %%",PeS);
    //printf("\nProbabilita' errore: %.4f %%",Pe);

    f=fopen(fileName, "w");

    if (f==NULL){
    	printf("%s non aperto\n",fileName);
    	exit(EXIT_FAILURE);
    }
    else{
    	fprintf(f, "sn = [");
    	for (int i = 0; i < nPoints; i++){

    		fprintf(f, "%lf, ", SNRARRAY[i]);
    	}
    	fprintf(f, "];");

    	for (int i = 0; i < 3; i++){

    		fprintf(f, "\nPs%d = [", i);
    		for (int j = 0; j < nPoints; j++){
    			fprintf(f, "%lf, ", PeSp[j+ i*nPoints]);
    		}
    		fprintf(f, "];");

    		fprintf(f, "\nPt%d = [", i);
    		for (int j = 0; j < nPoints; j++){
    			fprintf(f, "%lf, ", PeTh[j+ i*nPoints]);
    		}
			fprintf(f, "];");
    	}

		fprintf(f, "figure;\n plot(sn, Ps0);\nhold on;\nplot(sn,Pt0);\nhold off;\n");
		fprintf(f, "title('4-QAM');legend('sperimentale','teorico');xlabel('SNR');ylabel('Pe%');");

		fprintf(f, "figure;\n plot(sn, Ps1);\nhold on;\nplot(sn,Pt1);\nhold off;\n");
		fprintf(f, "title('16-QAM');legend('sperimentale','teorico');xlabel('SNR');ylabel('Pe%');");

		fprintf(f, "figure;\n plot(sn, Ps2);\nhold on;\nplot(sn,Pt2);\nhold off;\n");
		fprintf(f, "title('64-QAM');legend('sperimentale','teorico');xlabel('SNR');ylabel('Pe%');");

    }

    free(thresholds);
    free(greyMap);
    return 0;
}
