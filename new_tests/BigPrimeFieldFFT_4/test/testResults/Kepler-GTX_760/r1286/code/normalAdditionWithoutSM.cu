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

#define P 469762049


__global__ void batchNormalAdd(int *xs, int *ys, int *us) 
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;   
  short i, t;
  int num1,num2,num3;
  int *xd, *yd, *ud; 
  
  num1=0;
  
  xd = (int*)((char*)xs + tid*sizeof(int)*8);
  yd = (int*)((char*)ys + tid*sizeof(int)*8);
  ud = (int*)((char*)us + tid*sizeof(int)*8);
  
  for(t=0;t<100;t++)
  {
    for(i=0;i<=7;i++)
    {
      num2=xd[i];
      num3=yd[i];
      num1=num2+num3;
      if(num1<num2||num1<num3) //there is overflow/truncation
      {
        ud[i]=num1;
      }
      else if(num1>=P)
      {  
        ud[i]=num1-P;
      }
      else
      {
        ud[i]=num1;
      }
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
  
  xs=(int *)malloc((sizeof(int)*8*TN*BN));
  ys=(int *)malloc((sizeof(int)*8*TN*BN));
  us=(int *)malloc((sizeof(int)*8*TN*BN));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "small_number_benchmark.dat");
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  memset(xs,0,sizeof(int)*8*TN*BN);
  fread(xs,sizeof(int),8*TN*BN,fp1);
  memset(ys,0,sizeof(int)*8*TN*BN);
  fread(ys,sizeof(int),8*TN*BN,fp1);
  
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(int)*8*TN*BN);
  cudaMemcpy(xs_d, xs, sizeof(int)*8*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&ys_d, sizeof(int)*8*TN*BN);
  cudaMemcpy(ys_d, ys, sizeof(int)*8*TN*BN, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&us_d, sizeof(int)*8*TN*BN);
  cudaMemcpy(us_d, us, sizeof(int)*8*TN*BN, cudaMemcpyHostToDevice);
  
  batchNormalAdd<<<BN,TN>>>(xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(int)*8*TN*BN, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done 100 times of us = (xs + ys) mod 469762049 without share memory, %d normal number pair.\n", BN*TN*8);
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



