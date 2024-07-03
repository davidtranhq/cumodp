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
#include "inlines.h"

//BN block number
#define BN 512
//TN thread number in a block
#define TN 512 
//R  2^63+2^34
#define R 9223372054034644992
//P  near to R
#define P 9223372054034644991
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//ULMAX=2^64-1
#define ULMAX 18446744073709551615

__device__ __host__ __inline__ unsigned long add_mod(unsigned long a, unsigned long b, unsigned long p) {
    int r = a + b;
    r -= p;
    r += (r >> BASE_1) & p;
    return r;
}

__global__ void batchBigAdd(unsigned long r, unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;  
  short i;
  unsigned long num1;
  unsigned long *xd, *yd, *ud; 
  num1=0;
  
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  yd = (unsigned long*)((char*)ys + tid*sizeof(unsigned long)*8);
  ud = (unsigned long*)((char*)us + tid*sizeof(unsigned long)*8);

  
    for(i=0;i<=7;i++)
    {
      //ud[i]=add_mod(xd[i],yd[i],P);
      num1 = xd[i] + yd[i];
      if(num1<xd[i]||num1<yd[i])
      {
        ud[i]=num1+RC;
      }
      else if(num1>=r)
      {  
        ud[i]=num1-r;
      }
      else
      {
        ud[i]=num1;
      }
    }  
    
}

int main(int argc, char *argv[])
{
  int i;
  char fileName[1024];
	FILE *fp1;	 
  unsigned long *xs, *ys, *us, *xs_d, *ys_d, *us_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  ys=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  us=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "addition_benchmark_data.dat");
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  memset(xs,0,sizeof(unsigned long)*TN*BN);
  fread(xs,sizeof(unsigned long),TN*BN,fp1);
  memset(ys,0,sizeof(unsigned long)*TN*BN);
  fread(ys,sizeof(unsigned long),TN*BN,fp1);
  
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*TN*BN);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&ys_d, sizeof(unsigned long)*8*TN*BN);
  cudaMemcpy(ys_d, ys, sizeof(unsigned long)*8*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&us_d, sizeof(unsigned long)*8*TN*BN);
  cudaMemcpy(us_d, us, sizeof(unsigned long)*8*TN*BN, cudaMemcpyHostToDevice);
  
  batchBigAdd<<<BN,TN>>>(R,xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(unsigned long)*8*TN*BN, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done us = xs + ys.\n");
  printf("the xs is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", xs[i]);
  }
  printf("the ys is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", ys[i]);
  }
  printf("the us is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", us[i]);
  }  
  printf("the time of gpu is %f ms\n", elapsedTime);
  
  //free
  fclose(fp1);
  cudaFree(xs_d);	
  cudaFree(ys_d);	
  cudaFree(us_d);	
  free(xs);
  free(ys);
  free(us);
}



