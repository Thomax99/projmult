#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define abs(a) ((a)>0?(a):-(a))
#define min(a,b) ((a>b)?(b):(a))
#define SIZEX 10000
#define SIZEY 10000

void stencil9(float *a, const float *b) {
  int i,j;
  for (j=2; j<SIZEY-2; j++) 
    for (i=2; i<SIZEX-2; i++) 
      a[i*SIZEY+j] = (8*b[i*SIZEY+j] + b[(i+1)*SIZEY+j] + b[(i+2)*SIZEY+j] 
			+ b[(i-1)*SIZEY+j] + b[(i-2)*SIZEY+j] + b[i*SIZEY+j+1]
			+ b[i*SIZEY+j-1] + b[i*SIZEY+j-2] + b[i*SIZEY+j+2])/9.;
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
    stencil9(a,b);
    s+=dot2D(a,b);
    stencil9(b,a);
    s+=dot2D(a,b);
    fprintf(stderr,".");
  }
    fprintf(stderr,"\n%f",s);
  free(a);
  free(b);
  return 0;
}

