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

//BN block number
#define BN 512
//TN thread number in a block
#define TN 256 

//BC one big number decomposed into BC normal int.
#define BC 16

#define P 469762049
#define BASE_1 31


__device__ __host__ __inline__ int mul_mod(int a, int b, int n) {
    double ninv = 1 / (double)n;
    int q  = (int) ((((double) a) * ((double) b)) * ninv);
    int res = a * b - q * n;
    res += (res >> BASE_1) & n;
    res -= n;
    res += (res >> BASE_1) & n;
    return res;
}

__global__ void batchNormalMul(int *xs, int *ys, int *us) 
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;   
  short i, t;
  int *xd, *yd, *ud; 
  
 
  xd = (int*)((char*)xs + tid*sizeof(int)*BC);
  yd = (int*)((char*)ys + tid*sizeof(int)*BC);
  ud = (int*)((char*)us + tid*sizeof(int)*BC);
  
  for(i=0;i<BC;i++)
  {
    for(t=0;t<100;t++)
    {
      ud[i]=mul_mod(xd[i],yd[i],P);
    }    
  }
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1;	 
  int *xs, *ys, *us, *xs_d, *ys_d, *us_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  
  xs=(int *)malloc((sizeof(int)*BC*TN*BN));
  ys=(int *)malloc((sizeof(int)*BC*TN*BN));
  us=(int *)malloc((sizeof(int)*BC*TN*BN));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "small_number_benchmark.dat");
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  memset(xs,0,sizeof(int)*BC*TN*BN);
  fread(xs,sizeof(int),BC*TN*BN,fp1);
  memset(ys,0,sizeof(int)*BC*TN*BN);
  fread(ys,sizeof(int),BC*TN*BN,fp1);
  
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(int)*BC*TN*BN);
  cudaMemcpy(xs_d, xs, sizeof(int)*BC*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&ys_d, sizeof(int)*BC*TN*BN);
  cudaMemcpy(ys_d, ys, sizeof(int)*BC*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&us_d, sizeof(int)*BC*TN*BN);
  cudaMemcpy(us_d, us, sizeof(int)*BC*TN*BN, cudaMemcpyHostToDevice);
  
  batchNormalMul<<<BN,TN>>>(xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(int)*BC*TN*BN, cudaMemcpyDeviceToHost);
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
  free(xs);
  free(ys);
  free(us);
}



