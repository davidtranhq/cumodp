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
#define TN 16 
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

__global__ void FFT(unsigned long r, unsigned long *xs, unsigned long *ys) 
{
  //Now, we only use one thread to do FFT 16. 
  //FFT 16 has four rounds. Each round can be done in parallel. 
  //int tid = blockIdx.x*blockDim.x + threadIdx.x;  
  short i, j;
  unsigned long *xd, *yd;   
  
  //point initialization
  //xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  xd = (unsigned long*)((char*)xs);
  yd = (unsigned long*)((char*)ys);
  
  
  //first round
  //0  8
  bigAdd(&xd[0],&xd[64], &yd[0]);
  bigSub(&xd[0],&xd[64], &yd[64]);
  //4  12
  bigAdd(&xd[32],&xd[96], &yd[32]);
  bigSub(&xd[32],&xd[96], &yd[96]);
  //2  10
  bigAdd(&xd[16],&xd[80], &yd[16]);
  bigSub(&xd[16],&xd[80], &yd[80]);
  //6  14
  bigAdd(&xd[48],&xd[112], &yd[48]);
  bigSub(&xd[48],&xd[112], &yd[112]);
  //1  9
  bigAdd(&xd[8],&xd[72], &yd[8]);
  bigSub(&xd[8],&xd[72], &yd[72]);
  //5  13
  bigAdd(&xd[40],&xd[104], &yd[40]);
  bigSub(&xd[40],&xd[104], &yd[104]);
  //3  11
  bigAdd(&xd[24],&xd[88], &yd[24]);
  bigSub(&xd[24],&xd[88], &yd[88]);
  //7  15
  bigAdd(&xd[56],&xd[120], &yd[56]);
  bigSub(&xd[56],&xd[120], &yd[120]);
  for(i=0;i<128;i++)
  {
    xd[i]=yd[i];
  }
  
  //second round
  //0 4
  bigAdd(&xd[0],&xd[32], &yd[0]);
  bigSub(&xd[0],&xd[32], &yd[32]);
  //8 12
  cyclicShift(r,&xd[96],4);
  bigAdd(&xd[64],&xd[96], &yd[64]);
  bigSub(&xd[64],&xd[96], &yd[96]);
  
  //2  6
  bigAdd(&xd[16],&xd[48], &yd[16]);
  bigSub(&xd[16],&xd[48], &yd[48]);
  //10 14
  cyclicShift(r,&xd[112],4);
  bigAdd(&xd[80],&xd[112], &yd[80]);
  bigSub(&xd[80],&xd[112], &yd[112]);
  
  //1 5
  bigAdd(&xd[8],&xd[40], &yd[8]);
  bigSub(&xd[8],&xd[40], &yd[40]);
  //9 13
  cyclicShift(r,&xd[104],4);
  bigAdd(&xd[72],&xd[104], &yd[72]);
  bigSub(&xd[72],&xd[104], &yd[104]);
  
  //3  7
  bigAdd(&xd[24],&xd[56], &yd[24]);
  bigSub(&xd[24],&xd[56], &yd[56]);  
  //11  15
  cyclicShift(r,&xd[120],4);
  bigAdd(&xd[88],&xd[120], &yd[88]);
  bigSub(&xd[88],&xd[120], &yd[120]);
  for(i=0;i<128;i++)
  {
    xd[i]=yd[i];
  }
  
  //third round
  //0 2
  bigAdd(&xd[0],&xd[16], &yd[0]);
  bigSub(&xd[0],&xd[16], &yd[16]);
  
  //8 10
  cyclicShift(r,&xd[80],2);
  bigAdd(&xd[64],&xd[80], &yd[64]);
  bigSub(&xd[64],&xd[80], &yd[80]);
  
  //4 6
  cyclicShift(r,&xd[48],4);
  bigAdd(&xd[32],&xd[48], &yd[32]);
  bigSub(&xd[32],&xd[48], &yd[48]);  
  
  //12 14
  cyclicShift(r,&xd[112],6);
  bigAdd(&xd[96],&xd[112], &yd[96]);
  bigSub(&xd[96],&xd[112], &yd[112]);
  
  //1 3
  bigAdd(&xd[8],&xd[24], &yd[8]);
  bigSub(&xd[8],&xd[24], &yd[24]);
  
  //9 11
  cyclicShift(r,&xd[88],2);
  bigAdd(&xd[72],&xd[88], &yd[72]);
  bigSub(&xd[72],&xd[88], &yd[88]);
  
  //5 7
  cyclicShift(r,&xd[56],4);
  bigAdd(&xd[40],&xd[56], &yd[40]);
  bigSub(&xd[40],&xd[56], &yd[56]);
  
  //13 15
  cyclicShift(r,&xd[120],6);
  bigAdd(&xd[104],&xd[120], &yd[104]);
  bigSub(&xd[104],&xd[120], &yd[120]);
  for(i=0;i<128;i++)
  {
    xd[i]=yd[i];
  }
  
  //fourth
  //0 1
  bigAdd(&xd[0],&xd[8], &yd[0]);
  bigSub(&xd[0],&xd[8], &yd[8]);
  
  //8 9
  cyclicShift(r,&xd[72],1);
  bigAdd(&xd[64],&xd[72], &yd[64]);
  bigSub(&xd[64],&xd[72], &yd[72]);
  
  //4 5
  cyclicShift(r,&xd[40],2);
  bigAdd(&xd[32],&xd[40], &yd[32]);
  bigSub(&xd[32],&xd[40], &yd[40]);
  
  //12 13
  cyclicShift(r,&xd[104],3);
  bigAdd(&xd[96],&xd[104], &yd[96]);
  bigSub(&xd[96],&xd[104], &yd[104]);
  
  //2 3
  cyclicShift(r,&xd[24],4);
  bigAdd(&xd[16],&xd[24], &yd[16]);
  bigSub(&xd[16],&xd[24], &yd[24]);
  
  //10 11
  cyclicShift(r,&xd[88],5);
  bigAdd(&xd[80],&xd[88], &yd[80]);
  bigSub(&xd[80],&xd[88], &yd[88]);
  
  //6 7
  cyclicShift(r,&xd[56],6);
  bigAdd(&xd[48],&xd[56], &yd[48]);
  bigSub(&xd[48],&xd[56], &yd[56]);
  
  //14 15
  cyclicShift(r,&xd[120],7);
  bigAdd(&xd[112],&xd[120], &yd[112]);
  bigSub(&xd[112],&xd[120], &yd[120]);
  
  
  //0 0
  j=0;
  for(i=0;i<8;i++)
  {
    xd[i]=yd[j++];
  }  
  //1 8
  j=64;
  for(i=8;i<16;i++)
  {
    xd[i]=yd[j++];
  }  
  //2 4
  j=32;
  for(i=16;i<24;i++)
  {
    xd[i]=yd[j++];
  } 
  //3 12
  j=96;
  for(i=24;i<32;i++)
  {
    xd[i]=yd[j++];
  } 
  //4 2
  j=16;
  for(i=32;i<40;i++)
  {
    xd[i]=yd[j++];
  }  
  //5 10
  j=80;
  for(i=40;i<48;i++)
  {
    xd[i]=yd[j++];
  }   
  //6 6 omit
  j=48;
  for(i=48;i<56;i++)
  {
    xd[i]=yd[j++];
  } 
  //7 14
  j=112;
  for(i=56;i<64;i++)
  {
    xd[i]=yd[j++];
  }  
  //8 1
  j=8;
  for(i=64;i<72;i++)
  {
    xd[i]=yd[j++];
  } 
  //9 9 omit
  j=72;
  for(i=72;i<80;i++)
  {
    xd[i]=yd[j++];
  } 
  //10 5
  j=40;
  for(i=80;i<88;i++)
  {
    xd[i]=yd[j++];
  }  
  //11 13
  j=104;
  for(i=88;i<96;i++)
  {
    xd[i]=yd[j++];
  } 
  //12 3
  j=24;
  for(i=96;i<104;i++)
  {
    xd[i]=yd[j++];
  }  
  //13 11
  j=88;
  for(i=104;i<112;i++)
  {
    xd[i]=yd[j++];
  } 
  //14 7
  j=56;
  for(i=112;i<120;i++)
  {
    xd[i]=yd[j++];
  } 
  
  //15 15 
  j=120;
  for(i=120;i<128;i++)
  {
    xd[i]=yd[j++];
  } 
  
}

int main(int argc, char *argv[])
{
  int i;
  char fileName[1024];
	FILE *fp1;	 
  unsigned long *xs, *xs_d, *ys, *ys_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  ys=(unsigned long *)malloc((sizeof(unsigned long)*8*TN*BN));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_Data_16.dat"); //16*8 (unsigned long) data
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  //读取所有的数据。
  memset(xs,0,sizeof(unsigned long)*8*TN*BN);
  fread(xs,sizeof(unsigned long),8*TN*BN,fp1);
  
  printf("the xs is \n");
  for(i=0;i<128;i++)
  {
    printf ("%lu\n", xs[i]);
  }
  
  //begin gpu
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*TN*BN);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*TN*BN, cudaMemcpyHostToDevice);
  
  cudaMalloc((void **)&ys_d, sizeof(unsigned long)*8*TN*BN);
  cudaMemcpy(ys_d, ys, sizeof(unsigned long)*8*TN*BN, cudaMemcpyHostToDevice);
  
  FFT<<<1,1>>>(R,xs_d,ys_d);
  cudaThreadSynchronize();
  cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*TN*BN, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done xs = FFT(xs).\n");
  printf("the calculated xs is \n");
  for(i=0;i<128;i++)
  {
    printf ("%lu\n", xs[i]);
  }
  
  printf("the time of gpu is %f ms\n", elapsedTime);
  
  //free
  fclose(fp1);
  cudaFree(xs_d);	
  free(xs);
}



