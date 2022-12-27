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
#define TILEX 50
#endif
#ifndef TILEY
#define TILEY 50
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

__global__ void stencil9(float *a, const float *b) {
  int i = threadIdx.x + blockIdx.x * blockDim.x ;
  int j = threadIdx.y + blockIdx.y * blockDim.y ;
  if (i >= 2 && j >= 2 && i < SIZEX-2 && j < SIZEY-2)
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

  float *a_gpu, *b_gpu ;
  cudaMalloc((void**) &a_gpu, SIZEX*SIZEY*sizeof(float)) ;
  cudaMalloc((void**) &b_gpu, SIZEX*SIZEY*sizeof(float)) ;
  cudaMemcpy(a_gpu, a, SIZEX*SIZEY*sizeof(float), cudaMemcpyHostToDevice) ;
  cudaMemcpy(b_gpu, b, SIZEX*SIZEY*sizeof(float), cudaMemcpyHostToDevice) ;
  
  cudaEvent_t start, stop ;
  cudaEventCreate(&start) ;
  cudaEventCreate(&stop) ;
  cudaEventRecord(start);
  dim3 grid(SIZEX/TILEX, SIZEY/TILEY) ;
  dim3 block(TILEX, TILEY) ;
  for(k=0;k<100;k++) {
          
    stencil9<<<grid, block>>>(a_gpu,b_gpu);
    cudaDeviceSynchronize() ;
    stencil9<<<grid, block>>>(b_gpu,a_gpu);
    cudaDeviceSynchronize() ;
    fprintf(stderr,".");
  }
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  fprintf(stderr,"\n%f\n",s);
  float elapsed = 0 ;
  cudaEventElapsedTime(&elapsed, start, stop) ;
  fprintf(stderr, "It takes %lf millisecs\n", elapsed) ;
  cudaMemcpy(a, a_gpu, SIZEX*SIZEY*sizeof(float), cudaMemcpyDeviceToHost) ;
  cudaMemcpy(b, b_gpu, SIZEX*SIZEY*sizeof(float), cudaMemcpyDeviceToHost) ;
   FILE * f = fopen("result_cu", "w");
   dump_stencil(a, f) ;
   fclose(f) ;
  free(a);
  free(b);
  return 0;
}

