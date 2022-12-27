#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define abs(a) ((a)>0?(a):-(a))
//#define min(a,b) ((a>b)?(b):(a))
#ifndef SIZEX
#define SIZEX 10000
#endif
#ifndef SIZEY
#define SIZEY 10000
#endif
#ifndef TILEX
#define TILEX 100
#endif
#ifndef TILEY
#define TILEY 100
#endif
#include <omp.h>


void dump_stencil(float *a, FILE* f) {
	int i, j ;
	for (i=0; i < SIZEX; i++) {
	for (j=0; j < SIZEY;j++) {
		fprintf(f, "%lf ", a[i*SIZEY+j]) ;
	}
	fprintf(f, "\n");
	}
}

#include "mipp.h"
void stencil9_vecto(float *a, const float *b) {
	int ii, jj, i, j ;
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int x, y;
	mipp::Reg<float> regA, regij, regip1j, regip2j, regim1j, regim2j, regijp1, regijp2, regijm1, regijm2 ;
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = 2; i < SIZEX-2; i++) {
		for (j = 2; j + nmip-1 < SIZEY-2; j+=nmip) {
			regij.load(b +i*SIZEY+j) ;
			regip1j.load(b+(i+1)*SIZEY+j) ;
			regip2j.load(b+(i+2)*SIZEY+j) ;
			regim1j.load(b+(i-1)*SIZEY+j) ;
			regim2j.load(b+(i-2)*SIZEY+j) ;
			regijp1.load(b+i*SIZEY+j+1) ;
			regijp2.load(b+i*SIZEY+j+2) ;
			regijm1.load(b+i*SIZEY+j-1) ;
			regijm2.load(b+i*SIZEY+j-2) ;
			regA = ((reg8*regij + regim1j + regim2j + regip1j + regip2j + regijm1 + regijm2 + regijp1 + regijp2)/reg9) ;
			regA.store(a + i*SIZEY+j) ;
	        }
		for (; j < SIZEY-2; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
}

void stencil9_tiled(float *a, const float *b) {
	int ii, jj, i, j ;
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int x, y;
	for (x = 0; x < TX; x++)
	for (y = 0; y < TY; y++) {
		for (ii = (x == 0 ? 2 : 0); ii < (x == TX-1 ? TILEX-2 : TILEX); ii++)
		for (jj = (y == 0 ? 2 : 0); jj < (y == TY-1 ? TILEY-2 : TILEY); jj++) {
			i = x*TILEX + ii ; j = y*TILEY +jj ;
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;


		}
	}

}

void stencil9_omp(float *a, const float *b) {
  int i,j;
   
 
#pragma omp parallel for collapse(2) schedule(static)
  for (i=2; i<SIZEX-2; i++){
      for (j=2; j<SIZEY-2; j++) 
      a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;

    }
}

void stencil9(float *a, const float *b) {
  int i,j;
    for (i=2; i<SIZEX-2; i++){
     for (j=2; j<SIZEY-2; j++) 
      a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;

    }
}
float dot1D(float *a,float *b,int n)  {
  int i;
  float s=0;
  for (i=2; i<n-2; i++) 
    s+=a[i]*b[i];
  return s;
}
float dot2D(float *a,float *b)  {
  int i;
  float s=0;
  for (i=2; i<SIZEX-2; i++) 
    s+=dot1D(&a[i*SIZEY],&b[i*SIZEY],SIZEY);
  return s;
}
 int main() {
   int i,j,k;
  float *a,*b;
  float s=0;
  a=(float *)malloc(sizeof(float)*SIZEX*SIZEY);
  b=(float*)malloc(sizeof(float)*SIZEX*SIZEY);
  /* Initialization */
  for (i=0; i<SIZEX; i++)    
  for (j=0; j<SIZEY; j++)    
    a[i*SIZEY+j]=b[i*SIZEY+j]=0;
    for (j=SIZEY/4; j<SIZEY/2; j++) 
  for (i=SIZEX/4; i<SIZEX/2; i++) 
      b[i*SIZEY+j]=a[i*SIZEY+j]=1;
  /* Iterative computation. Iterate while error greater than ERROR */
  
    double st = omp_get_wtime() ;
    for(k=0;k<100;k++) {
    float *c;
    stencil9_omp(a,b);
    s+=dot2D(a,b);
    stencil9_omp(b,a);
    s+=dot2D(a,b);
    fprintf(stderr,".");
  }
    fprintf(stderr,"\n%f\n",s);
    double end = omp_get_wtime() ;
   fprintf(stderr, "It takes %lf secs\n", end - st) ;
   FILE * f = fopen("result_vecto", "w");
   dump_stencil(a, f) ;
   fclose(f) ;
  free(a);
  free(b);
  return 0;
}

