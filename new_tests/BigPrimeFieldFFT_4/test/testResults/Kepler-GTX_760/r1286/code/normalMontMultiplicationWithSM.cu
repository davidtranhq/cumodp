/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/


#include <time.h>
#include <stdio.h>
#include <unistd.h>

#include "ct_fft_mont.h"
#include "cudautils.h"
#include "inlines.h"
#include "printing.h"
#include "montmulmod.h"

//BN block number
#define BN 512
//TN thread number in a block
#define TN 256 

//BC one big number decomposed into BC normal int.
#define BC 16

#define P 469762049
#define BASE_1 31

__constant__ fprime_t fp;

static inline void setup_const(const fprime_t * const fpp) {
    cudaMemcpyToSymbol(fp, fpp, sizeof(fprime_t));
}

__global__ void batchNormalMul(sfixn *xs, sfixn *ys, sfixn *us) 
{
  __shared__ sfixn xsm[TN*BC];  
  __shared__ sfixn ysm[TN*BC]; 
  __shared__ sfixn usm[TN*BC]; 
  int tid = blockIdx.x*blockDim.x + threadIdx.x;   
  short i, pos, t;
  sfixn *xd, *yd, *ud, *xm, *ym, *um; 
  pos = threadIdx.x * BC;
  xm = &xsm[pos];
  ym = &ysm[pos];
  um = &usm[pos];
 
  xd = (sfixn*)((char*)xs + tid*sizeof(sfixn)*BC);
  yd = (sfixn*)((char*)ys + tid*sizeof(sfixn)*BC);
  ud = (sfixn*)((char*)us + tid*sizeof(sfixn)*BC);
  for(i=0;i<BC;i++)
  {
    xm[i]=xd[i];
  }  
  __syncthreads();
  
  for(i=0;i<BC;i++)
  {
    ym[i]=yd[i];
  }
  __syncthreads();
  
  for(i=0;i<BC;i++)
  {
    for(t=0;t<100;t++)
    {
      um[i]=fourier_reduction(xm[i], ym[i], &fp); 
    }    
  }
  __syncthreads();
  for(i=0;i<BC;i++)
  {
    ud[i]=um[i];
  }
}


int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1;	 
  int *xs1, *ys1;
  sfixn *xs2, *ys2, *us, *xs_d, *ys_d, *us_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  int i, nc = 0;
  fprime_t fp_h;
  
  
  nc = BC*TN*BN;
  xs1=(int *)malloc((sizeof(int)*nc));
  ys1=(int *)malloc((sizeof(int)*nc));
  xs2=(int *)malloc((sizeof(sfixn)*nc));
  ys2=(int *)malloc((sizeof(sfixn)*nc));
  us=(int *)malloc((sizeof(sfixn)*nc));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "small_number_benchmark.dat");
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  
  memset(xs1,0,sizeof(int)*nc);
  fread(xs1,sizeof(int),nc,fp1);
  memset(ys1,0,sizeof(int)*nc);
  fread(ys1,sizeof(int),nc,fp1);
  
  for(i=0;i<nc;i++)
  {
    xs2[i]=xs1[i];
    ys2[i]=ys1[i];
  }  
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  
  //init gpu memory
  init_fourier_prime(&fp_h, P);
  setup_const(&fp_h);
  
  cudaMalloc((void **)&xs_d, sizeof(sfixn)*nc);
  cudaMemcpy(xs_d, xs2, sizeof(sfixn)*nc, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&ys_d, sizeof(sfixn)*nc);
  cudaMemcpy(ys_d, ys2, sizeof(sfixn)*nc, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&us_d, sizeof(sfixn)*nc);
  cudaMemcpy(us_d, us, sizeof(sfixn)*nc, cudaMemcpyHostToDevice);
  
  batchNormalMul<<<BN,TN>>>(xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(sfixn)*nc, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done 100 times of us = (xs * ys) mod 469762049 without share memory, %d normal number pair.\n", BN*TN*BC);
  printf("the time of gpu is %f ms\n", elapsedTime);
  
  //free
  fclose(fp1);
  cudaFree(xs_d);	
  cudaFree(ys_d);	
  cudaFree(us_d);	
  free(xs1);
  free(ys1);
  free(xs2);
  free(ys2);
  free(us);
}



