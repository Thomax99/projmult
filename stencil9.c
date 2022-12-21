#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define abs(a) ((a)>0?(a):-(a))
#define min(a,b) ((a>b)?(b):(a))
#define SIZEX 10000
#define SIZEY 10000

#define TILEX 100
#define TILEY 100

void stencil9_tiled(float * restrict a, const float * restrict b, const int tile_x, const int tile_y) {
  int i,j;
  int sizeij = 2*SIZEY+2, sizeip1j, sizeip2j, sizeim1j, sizeim2j, sizeijp1, sizeijm1, sizeijm2, sizeijp2 ;
  sizeip2j = sizeip1j+SIZEY ;
  sizeim1j = sizeij - SIZEY ;
  sizeim2j = sizeim1j - SIZEY ;
  sizeijp1 = sizeij+1 ;
  sizeijp2 = sizeijp1+1;
  sizeijm1 = sizeij-1 ;
  sizeijm2 = sizeijm1-1 ;
  int nb_X = SIZEX / tile_x, nb_Y = SIZEY / tile_y ;
  int tx, ty ;
  for (tx = 0 ; tx < nb_X ; tx++) {
      for (ty = 0 ; ty < nb_Y ; ty++) {
		for (i = (tx == 0)*2; i < tile_x - (tx == (nb_X-1))*2; i++) {
			for (j = (ty == 0)*2; j < tile_y - (ty == (nb_Y-1))*2; j++) {
				sizeij = (tx*tile_x + i)*SIZEX + ty*tile_y + j ;
    				a[sizeij] = (8*b[sizeij] + b[sizeij + SIZEY] + b[sizeij + 2*SIZEY] 
					+ b[sizeij - SIZEY] + b[sizeij - 2*SIZEY] + b[sizeij+1]
					+ b[sizeij - 1] + b[sizeij - 2] + b[sizeij + 2])/9.;
			}
		}
	    }
     }
}

void stencil9(float * restrict a, const float * restrict b) {
  int i,j;
  int sizeij = 2*SIZEY+2, sizeip1j, sizeip2j, sizeim1j, sizeim2j, sizeijp1, sizeijm1, sizeijm2, sizeijp2 ;
  sizeip2j = sizeip1j+SIZEY ;
  sizeim1j = sizeij - SIZEY ;
  sizeim2j = sizeim1j - SIZEY ;
  sizeijp1 = sizeij+1 ;
  sizeijp2 = sizeijp1+1;
  sizeijm1 = sizeij-1 ;
  sizeijm2 = sizeijm1-1 ;	
  for (i=2; i<SIZEX-2; i++){
     sizeij += 4 ;
     sizeip2j += 4;
     sizeip1j += 4;
     sizeim1j += 4;
     sizeim2j += 4;
     sizeijp1 += 4;
     sizeijp2 += 4;
     sizeijm1 += 4;
     sizeijm2 += 4;
     for (j=2; j<SIZEY-2; j++, sizeij++, sizeip1j++, sizeip2j++, sizeim1j++, sizeim2j++, sizeijp1++, sizeijp2++, sizeijm1++, sizeijm2++) {
	          a[sizeij] = (8*b[sizeij] + b[sizeip1j] + b[sizeip2j] 
			+ b[sizeim1j] + b[sizeim2j] + b[sizeijp1]
			+ b[sizeijm1] + b[sizeijm2] + b[sizeijp2])/9.;
      }
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
    for(k=0;k<100;k++) {
    float *c;
    stencil9_tiled(a,b, TILEX, TILEY);
    s+=dot2D(a,b);
    stencil9_tiled(b,a, TILEX, TILEY);
    s+=dot2D(a,b);
    fprintf(stderr,".");
  }
    fprintf(stderr,"\n%f",s);
  free(a);
  free(b);
  return 0;
}

