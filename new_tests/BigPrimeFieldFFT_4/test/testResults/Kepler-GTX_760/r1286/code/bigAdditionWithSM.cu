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
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//ULMAX=2^64-1
#define ULMAX 18446744073709551615


__global__ void batchBigAdd(unsigned long r, unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
  __shared__ unsigned long xsm[TN*8];  
  __shared__ unsigned long ysm[TN*8]; 
  __shared__ unsigned long usm[TN*8]; 
  int tid = blockIdx.x*blockDim.x + threadIdx.x;   
  short i, pos, pos1, t;
  unsigned short c=0;
  unsigned long num1,num2;
  unsigned long *xd, *yd, *ud, *xm, *ym, *um; 
  
  num1=0;
  num2=r-1;
  
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  yd = (unsigned long*)((char*)ys + tid*sizeof(unsigned long)*8);
  ud = (unsigned long*)((char*)us + tid*sizeof(unsigned long)*8);
  pos = threadIdx.x * 8;
  xm = &xsm[pos];
  ym = &ysm[pos];
  um = &usm[pos];

  for(i=0;i<8;i++)
  {
    xm[i]=xd[i];
  }  
  for(i=0;i<8;i++)
  {
    ym[i]=yd[i];
  }
  __syncthreads();
  
  for(t=0;t<100;t++)
  {
    for(i=0;i<=7;i++)
    {
      num1=xm[i]+ym[i]+c;
      if(num1<xm[i]||num1<ym[i]) //there is overflow/truncation
      {
        c=1;
        um[i]=num1+RC; //2^64-r
      }
      else if(num1>=r)
      {  
        c=1;
        um[i]=num1-r;
      }
      else
      {
        c=0;
        um[i]=num1;
      }
    }
    if(c>0)
    {
      pos1=-1;
      for(i=0;i<8;i++) 
      {
        if(um[i] != 0)
        {
          pos1=i;
          break;
        }
      }
      if(pos1>=0)
      {
        for(i=0;i<pos1;i++)
        {
            um[i] = num2;
        }
        um[pos]--;
      }            
      else
      {
        um[0]=ULMAX;
        for(i=1;i<8;i++)
        {
          um[i]=0;
        }
      }  
    }  
  }
  for(i=0;i<8;i++)
  {
    ud[i]=um[i];
  }
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1;	 
  unsigned long *xs, *ys, *us, *xs_d, *ys_d, *us_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  ys=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  us=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "big_number_benchmark.dat");
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
  
  printf("we have done 100 times of us = xs + ys with shared memory, %d big number pair.\n", BN*TN);
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



