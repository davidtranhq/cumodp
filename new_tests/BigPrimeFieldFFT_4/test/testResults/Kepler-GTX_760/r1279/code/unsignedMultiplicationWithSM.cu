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
#define BN 1
//TN thread number in a block
#define TN 512
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//ULMAX=2^64-1
#define ULMAX 18446744073709551615
//sqrt(R) 3037000502
#define SQRTR 3037000502


// for 0=< x,y < r, compute l, h, c s.t. x*y = l+h*r+c*r^2 with
// 0<= l, h < r and c=0
__device__ void mulTwoNum(unsigned long x, unsigned long y, unsigned long *s)
{
  unsigned int x0,x1,y0,y1;
  unsigned long a0,a1,t,t0;
  s[2]=0;
  x1=x>>32;
  x0=x-(((unsigned long)x1)<<32);
  y1=y>>32;
  y0=y-(((unsigned long)y1)<<32);
  a0=x0*y0;
  a1=x1*y1;
  
  t=x0*y1;
  a1+=(t>>32);
  a0+=((t-(t>>32)<<32)<<32);
  t=x1*y0;
  a1+=(t>>32);
  a0+=((t-(t>>32)<<32)<<32);
  
  a1=(a1<<1);
  a0>=R?s[1]=1:s[1]=0;
  s[0]=(s[1]>0?a0-R:a0);  //a0 finish
  s[1]+=a1;
  //-a1*2^34  0<=a1<2^64
  t=a1>>29;      //d1
  t0=a1-(t<<29); //d0
  if(t0>=t)      //d0>=d1
  {
    a0=t0-t;
    y1=a0>>29;      //e1
    y0=a0-(y1<<29); //e0
    t+=y1;
    if(y0>=y1)
    {
      t0=((unsigned long)(y0-y1))<<34;
    }
    else
    {
      t0=R-(((unsigned long)(y1-y0))<<34);
      t--;
    }
  }
  else
  {
    //d0<d1; x0<t
    a0=t-t0;
    y1=a0>>29;      //e1
    y0=a0-(y1<<29); //e0
    if(y0>=y1)
    {
      t=t-1-y1;
      t0=R-((unsigned long)(y0-y1)<<34);
    }
    else
    {
      t-=y1;
      t0=(unsigned long)(y1-y0)<<34;
    }
  }
  //we calculate a1*2^34 and store into t0+t*r;
  s[1]-=t;
  if(s[0]>=t0)
  {
    s[0]-=t0;
  }
  else
  {
    s[0]=R-t0+s[0];
    s[1]--;
  }
}

// pre-processor for the previous function
__device__ void mulLong(unsigned long x, unsigned long y, unsigned long *s)
{
  if(x==R&&y==R)
  {
    s[0]=0;s[1]=0;s[2]=1;
  }
  else if(x==R&&y!=R)
  {
    s[0]=0;s[1]=y;s[2]=0;
  }
  else if(x!=R&&y!=R)
  {
    s[0]=0;s[1]=x;s[2]=0;
  }
  else if(x<=SQRTR && y<=SQRTR)
  {
    s[0]=x*y;s[1]=0;s[2]=0;
  }
  else
  {
    mulTwoNum(x,y,s);
  }
}



// Test for the above
__global__ void unsignedMul(unsigned long r, unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
  __shared__ unsigned long xsm[TN];  
  __shared__ unsigned long ysm[TN]; 
  __shared__ unsigned long usm[TN*3]; 
  int tid = threadIdx.x;   
  short pos, t;
  
  
 
  xsm[tid]=xs[tid];
  __syncthreads();
  ysm[tid]=ys[tid];
  __syncthreads();
  
  pos=tid*3;
  for(t=0;t<100;t++)
  {
    mulLong(xsm[tid],ysm[tid],&usm[pos]);    
  }
  __syncthreads();
  for(t=0;t<3;t++)
  {
    us[pos]=usm[pos];
    pos++;
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
  
  //read all data
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
  
  unsignedMul<<<BN,TN>>>(R,xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(unsigned long)*8*TN*BN, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done 100 times of us = xs * ys with shared memory, %d big number pair.\n", BN*TN);
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



