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
#define TN 384 
//input number count
#define INC 384
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//ULMAX=2^64-1
#define ULMAX 18446744073709551615

__device__ void bigAdd(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
  
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=xs[i]+ys[i]+c;
    if(num1<xs[i]||num1<ys[i]) //there is overflow/truncation
    {
      c=1;
      us[i]=num1+RC;
    }
    else if(num1>=R)
    {  
      c=1;
      us[i]=num1-R;
    }
    else
    {
      c=0;
      us[i]=num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] != 0)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = num2;
      }
      us[pos] = us[pos] - 1;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }  
}

__device__ void bigSub(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{  
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  num1=0;
  num2=R-1;
  
  for(i=0;i<=7;i++)
  {
    num1=ys[i]+c;
    if(xs[i]<num1) //there is not enough to do subtraction
    {
      c=1;
      us[i]=R-num1+xs[i];
    }
    else
    {
      c=0;
      us[i]=xs[i]-num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(us[i] < num2)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          us[i] = 0;
      }
      us[pos] ++;
    }            
    else
    {
      us[0]=ULMAX;
      for(i=1;i<8;i++)
      {
        us[i]=0;
      }
    }  
  }  
}

__device__ void cyclicShift(unsigned long r, unsigned long *xs, short sn)
{
  short i=0,j=0;
  unsigned long ts[8]={0};
  unsigned long ys[8]={0};
  if(sn<=0)
  {
    return;
  }
  j=8-sn;
  for(i=0;i<sn;i++)
  {
    ts[i]=xs[j++];
  }
  for(i=7-sn;i>=0;i--)
  {
    xs[i+sn]=xs[i];
  }
  for(i=0;i<sn;i++)
  {
    xs[i]=0;
  }
  bigSub(xs,&ts[0],&ys[0]);
  for(i=0;i<8;i++)
  {
    xs[i]=ys[i];
  }
}

__constant__ short ind1[8]={0,4,2,6,1,5,3,7};
__constant__ short ind2[16]={0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
  
__global__ void FFT(unsigned long r, unsigned long *xs) 
{
  //Now, we only use one block to do $24$ FFT16. 
  //int tid = blockIdx.x*blockDim.x + threadIdx.x;  
  short tid = threadIdx.x; 
  short wid = tid >> 5;  //warp no. [0,1,...,23]
  short bwid = tid >> 6; //big warp no.[0,1,...,11]
  short sf = wid - bwid << 1; //sub flag for odd warp[0,1]
  short wpn = (tid-(wid<<5))>>3; //(tid-(wid*32))/8  [0,1,2,3] 
  short wpid = tid-((tid>>3)<<3); //in [0,1,...,7]
  short posid = 0; //in[0,1,...,15]
  short i, pos, pos1;
  unsigned long *xd, *xf, *yf; 
  __shared__ unsigned long xsm[INC*8]; 
  __shared__ unsigned long ysm[INC*8];  
  
  
  //point initialization
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  pos=tid<<3; //tid*8
  for(i=0;i<8;i++)
  {
    xsm[pos++]=xd[i];
  }
  __syncthreads();  //copy the input to shared memory.
  
  //xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  //yd = (unsigned long*)((char*)ys + tid*sizeof(unsigned long)*8);
  xf = (unsigned long*)((char*)xsm + (bwid*64+wpn*16)*sizeof(unsigned long)*8);
  yf = (unsigned long*)((char*)ysm + (bwid*64+wpn*16)*sizeof(unsigned long)*8);
  
  //first round
  if(sf)
  {
    //bigSub
    bigSub(&xf[wpid*8],&xf[(wpid+8)*8],&yf[(wpid+8)*8]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[wpid*8],&xf[(wpid+8)*8],&yf[wpid*8]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  //second round
  posid=wpid>=4?wpid+4:wpid;
  if(posid>=8)
  {
    cyclicShift(r,&xf[(posid+4)*8],4);
  }
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid*8],&xf[(posid+4)*8],&yf[(posid+4)*8]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid*8],&xf[(posid+4)*8],&yf[posid*8]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  //third round
  posid=wpid>=4?(wpid-4)*4+1:wpid*4;
  wpid>=4?cyclicShift(r,&xf[(posid+2)*8],posid>7?posid-7:posid-1):cyclicShift(r,&xf[(posid+2)*8],posid>6?posid-6:posid);
  
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid*8],&xf[(posid+2)*8],&yf[(posid+2)*8]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid*8],&xf[(posid+2)*8],&yf[posid*8]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  //fourth
  posid=wpid>>1;
  cyclicShift(r,&xf[(posid+1)*8],ind1[wpid]);
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid*8],&xf[(posid+1)*8],&yf[(posid+1)*8]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid*8],&xf[(posid+1)*8],&yf[posid*8]);
  }
  __syncthreads();  
  
  if(sf)
  {
    posid=wpid+8;
    pos=posid*8;
    pos1=ind2[posid]*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  else
  {
    posid=wpid;
    pos=posid*8;
    pos1=ind2[posid]*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  __syncthreads();  
  pos=tid<<3; //tid*8
  for(i=0;i<8;i++)
  {
    xd[i]=xsm[pos++];
  }
  __syncthreads();   
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1, *fp2;	 
  unsigned long *xs, *xs_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*INC));
  //ys=(unsigned long *)malloc((sizeof(unsigned long)*8*INC));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_Data_Block.dat"); //INC*8 (unsigned long) data
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  memset(xs,0,sizeof(unsigned long)*8*INC);
  fread(xs,sizeof(unsigned long),8*INC,fp1);
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*INC);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*INC, cudaMemcpyHostToDevice);
  
  //cudaMalloc((void **)&ys_d, sizeof(unsigned long)*8*INC);
  //cudaMemcpy(ys_d, ys, sizeof(unsigned long)*8*INC, cudaMemcpyHostToDevice);
  
  FFT<<<BN,TN>>>(R,xs_d);
  cudaThreadSynchronize();
   
  cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*INC, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done xs = FFT(xs).\n");
  printf("the time of gpu is %f ms\n", elapsedTime);
  printf("we now write the output result.\n");
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_Data_Block_output.dat"); //INC*8 (unsigned long) data
  if((fp2=fopen(fileName,"wb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp2);
  	exit(-1);
  } 
  fwrite(xs,sizeof(unsigned long),8*INC,fp2); 
  
  //free
  fclose(fp1);
  fclose(fp2);
  cudaFree(xs_d);	
  free(xs);
}



