#include <stdio.h>
#include <stdlib.h>

int N=32;
#define c0 48
#define c1 49

char subM[4][4]={
    {'1','1','0','0'},
    {'1','1','1','0'},
    {'0','1','1','1'},
    {'0','0','1','1'}
    };

int main()
{
	//allocazione matrice

	char ** A= (char **)malloc( N* sizeof(char *));
	for(int i=0; i<32; i++){
        A[i]=(char*)malloc(32 * sizeof(char));
	}

	for(int i=0; i<32; i++){     //riempio di "0"la matrice
		for(int j=0; j<32; j++){
			A[i][j]=48;
		}
	}

	for(int i=0; i<32; i++){    //stampo la matrice vuota
		for(int j=0; j<32; j++){
			printf("%c ",A[i][j]);
		}
		printf("\n");
	}

	for(int i=0; i<10; i++) printf("\n");

									//incollo matrice subM nella matrice A(ib,jb)
	for(int ib=0; ib<8; ib++){
		for(int jb=0; jb<8; jb++){

			for(int a=0; a<4; a++){
				for (int b = 0; b<4; ++b){

					A[4*ib+a][4*jb+b]= subM[a][b];
                    //printf("%c ",subM[a][b]);   //debug
				}
				//printf("\n");  //debug
			}
		}
	}
	for(int i=0; i<32; i++){		//stampo la matrice rienpita
		for(int j=0; j<32; j++){

			printf("%c ",A[i][j]);

		}
		printf("\n");
	}

                            //deallocazione memoria
	for(int i=0; i<32; i++){
        free(A[i]);
	}
	free(A);


	return 0;
}
