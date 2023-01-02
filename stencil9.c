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
#define TILEX 10 
#endif
#ifndef TILEY
#define TILEY 10
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
void _stencil9_vecto(float *a, const float *b, int sx, int sy, int width, int height) {
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int i, j;
	mipp::Reg<float> regA, regij, regip1j, regip2j, regim1j, regim2j, regijp1, regijp2, regijm1, regijm2 ;
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = sx; i < width ; i++) {
		for (j = sy; j + nmip-1 < height; j+=nmip) {
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
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
}

void _stencil9_vecto_unrolled_j4(float *a, const float *b, int sx, int sy, int width, int height) {
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int i, j;
	mipp::Reg<float> regA, regij, regip1j, regip2j, regim1j, regim2j, regijp1, regijp2, regijm1, regijm2 ;
	mipp::Reg<float> regA2, reg2ij, reg2ip1j, reg2ip2j, reg2im1j, reg2im2j, reg2ijp1, reg2ijp2, reg2ijm1, reg2ijm2 ;
	mipp::Reg<float> regA3, reg3ij, reg3ip1j, reg3ip2j, reg3im1j, reg3im2j, reg3ijp1, reg3ijp2, reg3ijm1, reg3ijm2 ;
	mipp::Reg<float> regA4, reg4ij, reg4ip1j, reg4ip2j, reg4im1j, reg4im2j, reg4ijp1, reg4ijp2, reg4ijm1, reg4ijm2 ;
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = sx; i < width ; i++) {
		for (j = sy; j + 4*nmip-1 < height; j+=4*nmip) {
			regij.load(b +i*SIZEY+j) ;
			reg2ij.load(b+i*SIZEY+j+nmip) ;
			reg3ij.load(b +i*SIZEY+j+2*nmip) ;
			reg4ij.load(b+i*SIZEY+j+3*nmip) ;
			
			regip1j.load(b+(i+1)*SIZEY+j) ;
			reg2ip1j.load(b+(i+1)*SIZEY+j+nmip) ;
			reg3ip1j.load(b+(i+1)*SIZEY+j+2*nmip) ;
			reg4ip1j.load(b+(i+1)*SIZEY+j+3*nmip) ;
			
			regip2j.load(b+(i+2)*SIZEY+j) ;
			reg2ip2j.load(b+(i+2)*SIZEY+j+nmip) ;
			reg3ip2j.load(b+(i+2)*SIZEY+j+2*nmip) ;
			reg4ip2j.load(b+(i+2)*SIZEY+j+3*nmip) ;
			
			regim1j.load(b+(i-1)*SIZEY+j) ;
			reg2im1j.load(b+(i-1)*SIZEY+j+nmip) ;
			reg3im1j.load(b+(i-1)*SIZEY+j+2*nmip) ;
			reg4im1j.load(b+(i-1)*SIZEY+j+3*nmip) ;
			
			regim2j.load(b+(i-2)*SIZEY+j) ;
			reg2im2j.load(b+(i-2)*SIZEY+j+nmip) ;
			reg3im2j.load(b+(i-2)*SIZEY+j+3*nmip) ;
			reg4im2j.load(b+(i-2)*SIZEY+j+4*nmip) ;
			
			regijp1.load(b+i*SIZEY+j+1) ;
			reg2ijp1.load(b+i*SIZEY+j+1+nmip) ;
			reg3ijp1.load(b+i*SIZEY+j+1+2*nmip) ;
			reg4ijp1.load(b+i*SIZEY+j+1+3*nmip) ;
			
			regijp2.load(b+i*SIZEY+j+2) ;
			reg2ijp2.load(b+i*SIZEY+j+2+nmip) ;
			reg3ijp2.load(b+i*SIZEY+j+3+nmip) ;
			reg4ijp2.load(b+i*SIZEY+j+4+nmip) ;
			
			regijm1.load(b+i*SIZEY+j-1) ;
			reg2ijm1.load(b+i*SIZEY+j-1+nmip) ;
			reg3ijm1.load(b+i*SIZEY+j-1+2*nmip) ;
			reg4ijm1.load(b+i*SIZEY+j-1+3*nmip) ;
			
			regijm2.load(b+i*SIZEY+j-2) ;
			reg2ijm2.load(b+i*SIZEY+j-2+nmip) ;
			reg3ijm2.load(b+i*SIZEY+j-2+2*nmip) ;
			reg4ijm2.load(b+i*SIZEY+j-2+3*nmip) ;
			
			regA = ((reg8*regij + regim1j + regim2j + regip1j + regip2j + regijm1 + regijm2 + regijp1 + regijp2)/reg9) ;
			regA2 = ((reg8*reg2ij + reg2im1j + reg2im2j + reg2ip1j + reg2ip2j + reg2ijm1 + reg2ijm2 + reg2ijp1 + reg2ijp2)/reg9) ;
			regA3 = ((reg8*reg3ij + reg3im1j + reg3im2j + reg3ip1j + reg3ip2j + reg3ijm1 + reg3ijm2 + reg3ijp1 + reg3ijp2)/reg9) ;
			regA4 = ((reg8*reg4ij + reg4im1j + reg4im2j + reg4ip1j + reg4ip2j + reg4ijm1 + reg4ijm2 + reg4ijp1 + reg4ijp2)/reg9) ;
			regA.store(a + i*SIZEY+j) ;
			regA2.store(a + i*SIZEY+j+nmip) ;
			regA3.store(a + i*SIZEY+j+2*nmip) ;
			regA4.store(a + i*SIZEY+j+3*nmip) ;
	        }
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
}


void _stencil9_vecto_unrolled_j2(float *a, const float *b, int sx, int sy, int width, int height) {
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int i, j;
	mipp::Reg<float> regA, regij, regip1j, regip2j, regim1j, regim2j, regijp1, regijp2, regijm1, regijm2 ;
	mipp::Reg<float> regA2, reg2ij, reg2ip1j, reg2ip2j, reg2im1j, reg2im2j, reg2ijp1, reg2ijp2, reg2ijm1, reg2ijm2 ;
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = sx; i < width ; i++) {
		for (j = sy; j + 2*nmip-1 < height; j+=2*nmip) {
			regij.load(b +i*SIZEY+j) ;
			reg2ij.load(b+i*SIZEY+j+nmip) ;
			regip1j.load(b+(i+1)*SIZEY+j) ;
			reg2ip1j.load(b+(i+1)*SIZEY+j+nmip) ;
			regip2j.load(b+(i+2)*SIZEY+j) ;
			reg2ip2j.load(b+(i+2)*SIZEY+j+nmip) ;
			regim1j.load(b+(i-1)*SIZEY+j) ;
			reg2im1j.load(b+(i-1)*SIZEY+j+nmip) ;
			regim2j.load(b+(i-2)*SIZEY+j) ;
			reg2im2j.load(b+(i-2)*SIZEY+j+nmip) ;
			regijp1.load(b+i*SIZEY+j+1) ;
			reg2ijp1.load(b+i*SIZEY+j+1+nmip) ;
			regijp2.load(b+i*SIZEY+j+2) ;
			reg2ijp2.load(b+i*SIZEY+j+2+nmip) ;
			regijm1.load(b+i*SIZEY+j-1) ;
			reg2ijm1.load(b+i*SIZEY+j-1+nmip) ;
			regijm2.load(b+i*SIZEY+j-2) ;
			reg2ijm2.load(b+i*SIZEY+j-2+nmip) ;
			regA = ((reg8*regij + regim1j + regim2j + regip1j + regip2j + regijm1 + regijm2 + regijp1 + regijp2)/reg9) ;
			regA2 = ((reg8*reg2ij + reg2im1j + reg2im2j + reg2ip1j + reg2ip2j + reg2ijm1 + reg2ijm2 + reg2ijp1 + reg2ijp2)/reg9) ;
			regA.store(a + i*SIZEY+j) ;
			regA2.store(a + i*SIZEY+j+nmip) ;
	        }
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
}

void _stencil9_vecto_unrolled_i4(float *a, const float *b, int sx, int sy, int width, int height) {
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int i, j;
	mipp::Reg<float> regA, regA2, regA3, regA4, regij, regip1j, regip2j, regip3j, regim1j, regim2j, regijp1, regip1jp1, regijp2, regip1jp2, regijm1, regip1jm1, regijm2, regip1jm2 ;
	mipp::Reg<float> regip4j, regip2jm1, regip2jm2, regip2jp1, regip2jp2, regip5j, regip3jm1, regip3jm2, regip3jp1, regip3jp2 ; 
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = sx; i + 3 < width ; i+=4) {
		for (j = sy; j + nmip-1 < height; j+=nmip) {
			regijm2.load(b+i*SIZEY+j-2) ;
			regijm1.load(b+i*SIZEY+j-1) ;
			regim2j.load(b+(i-2)*SIZEY+j) ;
			regim1j.load(b+(i-1)*SIZEY+j) ;
			
			regij.load(b +i*SIZEY+j) ;
			
			regip1j.load(b+(i+1)*SIZEY+j) ;
			regip2j.load(b+(i+2)*SIZEY+j) ;
			regijp1.load(b+i*SIZEY+j+1) ;
			regijp2.load(b+i*SIZEY+j+2) ;
			regip1jm2.load(b+(i+1)*SIZEY+j-2) ;
			regip1jm1.load(b+(i+1)*SIZEY+j-1) ;
			regip2jm1.load(b+(i+2)*SIZEY+j-1) ;
			regip3jm1.load(b+(i+3)*SIZEY+j-1) ;
			regip2jm2.load(b+(i+2)*SIZEY+j-2) ;
			regip3jm2.load(b+(i+3)*SIZEY+j-2) ;
			regip1jp1.load(b+(i+1)*SIZEY+j+1) ;
			regip2jp1.load(b+(i+2)*SIZEY+j+1) ;
			regip3jp1.load(b+(i+3)*SIZEY+j+1) ;
			regip1jp2.load(b+(i+1)*SIZEY+j+2) ;
			regip2jp2.load(b+(i+2)*SIZEY+j+2) ;
			regip3jp2.load(b+(i+3)*SIZEY+j+2) ;
			regip3j.load(b+(i+3)*SIZEY+j) ;	
			regip4j.load(b+(i+4)*SIZEY+j) ;	
			regip5j.load(b+(i+5)*SIZEY+j) ;	
			regA = ((reg8*regij + regim1j + regim2j + regip1j + regip2j + regijm1 + regijm2 + regijp1 + regijp2)/reg9) ;
			regA2 = ((reg8*regip1j + regij + regim1j + regip2j + regip3j + regip1jm1 + regip1jm2 + regip1jp1 + regip1jp2)/reg9) ;
			regA3 = ((reg8*regip2j + regip1j + regij + regip3j + regip4j + regip2jm1 + regip2jm2 + regip2jp1 + regip2jp2)/reg9) ;
			regA4 = ((reg8*regip3j + regip2j + regip1j + regip4j + regip5j + regip3jm1 + regip3jm2 + regip3jp1 + regip3jp2)/reg9) ;

			regA.store(a + i*SIZEY+j) ;
			regA2.store(a+(i+1)*SIZEY+j) ;
			regA3.store(a+(i+2)*SIZEY+j) ;
			regA4.store(a+(i+3)*SIZEY+j) ;
	        }
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
	for (; i < width ; i++) {
		for (j = sy; j + nmip-1 < height; j+=nmip) {
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
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}

}



void _stencil9_vecto_unrolled_i2(float *a, const float *b, int sx, int sy, int width, int height) {
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int i, j;
	mipp::Reg<float> regA, regA2, regij, regip1j, regip2j, regip3j, regim1j, regim2j, regijp1, regip1jp1, regijp2, regip1jp2, regijm1, regip1jm1, regijm2, regip1jm2 ;
	mipp::Reg<float> reg8 = 8., reg9 = 9. ;
	int nmip = mipp::N<float>() ;
	for (i = sx; i + 1 < width ; i+=2) {
		for (j = sy; j + nmip-1 < height; j+=nmip) {
			regijm2.load(b+i*SIZEY+j-2) ;
			regijm1.load(b+i*SIZEY+j-1) ;
			regim2j.load(b+(i-2)*SIZEY+j) ;
			regim1j.load(b+(i-1)*SIZEY+j) ;
			
			regij.load(b +i*SIZEY+j) ;
			
			regip1j.load(b+(i+1)*SIZEY+j) ;
			regip2j.load(b+(i+2)*SIZEY+j) ;
			regijp1.load(b+i*SIZEY+j+1) ;
			regijp2.load(b+i*SIZEY+j+2) ;
			regip1jm2.load(b+(i+1)*SIZEY+j-2) ;
			regip1jm1.load(b+(i+1)*SIZEY+j-1) ;
			regip1jp1.load(b+(i+1)*SIZEY+j+1) ;
			regip1jp2.load(b+(i+1)*SIZEY+j+2) ;
			regip3j.load(b+(i+3)*SIZEY+j) ;	
			regA = ((reg8*regij + regim1j + regim2j + regip1j + regip2j + regijm1 + regijm2 + regijp1 + regijp2)/reg9) ;
			regA2 = ((reg8*regip1j + regij + regim1j + regip2j + regip3j + regip1jm1 + regip1jm2 + regip1jp1 + regip1jp2)/reg9) ;
			regA.store(a + i*SIZEY+j) ;
			regA2.store(a+(i+1)*SIZEY+j) ;
	        }
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}
	for (; i < width ; i++) {
		for (j = sy; j + nmip-1 < height; j+=nmip) {
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
		for (; j < height; j++) {
      		    a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
		}	
	}

}

void stencil9_vecto(float *a, const float *b) {
	_stencil9_vecto_unrolled_j4(a, b, 2, 2, SIZEX-2, SIZEY-2) ;
}
void stencil9_tiled(float *a, const float *b) {
	int ii, jj, i, j ;
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int x, y;
	for (x = 0; x < TX; x++)
	for (y = 0; y < TY; y++) {
		_stencil9_vecto_unrolled_i4(a, b, x == 0 ? 2 : x*TILEX, y == 0 ? 2 : y*TILEY, x == TX-1 ? SIZEX-2 : (x+1)*TILEX, y == TY-1 ? SIZEY-2 : (y+1)*TILEY) ;
	}	
}
void stencil9_omp(float *a, const float *b) {
	int ii, jj, i, j ;
	int TX = SIZEX / TILEX, TY = SIZEY /TILEY ;
	int x, y;
#pragma omp for collapse(2)
	for (x = 0; x < TX; x++)
	for (y = 0; y < TY; y++) {
		_stencil9_vecto_unrolled_i4(a, b, x == 0 ? 2 : x*TILEX, y == 0 ? 2 : y*TILEY, x == TX-1 ? SIZEX-2 : (x+1)*TILEX, y == TY-1 ? SIZEY-2 : (y+1)*TILEY) ;
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
#pragma omp parallel firstprivate(k)
    for(k=0;k<100;k++) {
    float *c;
    stencil9_omp(a,b);
    //s+=dot2D(a,b);
    stencil9_omp(b,a);
    //s+=dot2D(a,b);
#pragma omp single
    fprintf(stderr,".");
  }
    fprintf(stderr,"\n%f\n",s);
    double end = omp_get_wtime() ;
   fprintf(stderr, "It takes %lf secs\n", end - st) ;
   FILE * f = fopen("result20k_vecto_2", "w");
   //dump_stencil(a, f) ;
   fclose(f) ;
  free(a);
  free(b);
  return 0;
}

