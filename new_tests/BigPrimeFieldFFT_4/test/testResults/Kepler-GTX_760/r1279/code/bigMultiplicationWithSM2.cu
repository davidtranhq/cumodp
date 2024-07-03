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
#define BN 256
//TN thread number in a block
#define TN 256 
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//ULMAX=2^64-1
#define ULMAX 18446744073709551615
//sqrt(R) 3037000502
#define SQRTR 3037000502

__device__ void bigAdd2(unsigned long *xs, unsigned long *ys) 
{
  unsigned long us[8];
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
      us[0]=R;
    }  
  }  
  xs[0]=us[0];
  xs[1]=us[1];
  xs[2]=us[2];
  xs[3]=us[3];
  xs[4]=us[4];
  xs[5]=us[5];
  xs[6]=us[6];
  xs[7]=us[7];
}


__device__ void bigSub2(unsigned long *xs, unsigned long *ys) 
{  
  unsigned long us[8];
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
      us[7]=R;
      for(i=0;i<7;i++)
      {
        us[i]=0;
      }
    }  
  }
  xs[0]=us[0];
  xs[1]=us[1];
  xs[2]=us[2];
  xs[3]=us[3];
  xs[4]=us[4];
  xs[5]=us[5];
  xs[6]=us[6];
  xs[7]=us[7];  
}
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

//store in l0, h0, c0
__device__ void smallAdd(unsigned long *l0, unsigned long *h0, short *c0, unsigned long *l1, unsigned long *h1, unsigned long *c1)
{
  short c=0;
  unsigned long s=0;
  s=*l0+*l1;
  s<*l0 || s<*l1?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  *l0=s;
  
  *h1=*h1+c;  //h1<r<2^64-1. This means no overflow
  s=*h0+*h1;
  s<*h0||s<*h1?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  *h0=s;
  
  *c0=*c0+(short)*c1+c;
}

__device__ void bigMul(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
    unsigned long ts1[8];
    unsigned long ts2[8];
    unsigned long rs[3];
    short c0,c1,c2,c3,c4,c5,c6,c7;
    unsigned long l0,l1,l2,l3,l4,l5,l6,l7,h0,h1,h2,h3,h4,h5,h6,h7;
    
    //x0*y0
    mulLong(xs[0],ys[0],rs);
    l0=rs[0];
    h0=rs[1];
    c0=(short)rs[2];
    
    //x0*y1+x1*y0
    mulLong(xs[0],ys[1],rs);    
    l1=rs[0];
    h1=rs[1];
    c1=(short)rs[2];
    mulLong(xs[1],ys[0],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    
    //x0*y2+x1*y1+x2*y0
    mulLong(xs[0],ys[2],rs);    
    l2=rs[0];
    h2=rs[1];
    c2=(short)rs[2];
    mulLong(xs[1],ys[1],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[0],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    
    //x0*y3+x1*y2+x2*y1+x3*y0
    mulLong(xs[0],ys[3],rs);    
    l3=rs[0];
    h3=rs[1];
    c3=(short)rs[2];
    mulLong(xs[1],ys[2],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[1],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[0],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    
    //x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
    mulLong(xs[0],ys[4],rs);    
    l4=rs[0];
    h4=rs[1];
    c4=(short)rs[2];
    mulLong(xs[1],ys[3],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[2],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[1],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[0],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    
    //x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
    mulLong(xs[0],ys[5],rs);    
    l5=rs[0];
    h5=rs[1];
    c5=(short)rs[2];
    mulLong(xs[1],ys[4],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[3],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[2],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[1],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[0],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
 
    //x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
    mulLong(xs[0],ys[6],rs);    
    l6=rs[0];
    h6=rs[1];
    c6=(short)rs[2];
    mulLong(xs[1],ys[5],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[4],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[3],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[2],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[1],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[0],rs);
    smallAdd(&l6,&h6,&c6,&rs[0],&rs[1],&rs[2]);
    
    //x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
    mulLong(xs[0],ys[7],rs);    
    l7=rs[0];
    h7=rs[1];
    c7=(short)rs[2];
    mulLong(xs[1],ys[6],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[2],ys[5],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[4],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[3],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[2],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[1],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[0],rs);
    smallAdd(&l7,&h7,&c7,&rs[0],&rs[1],&rs[2]);
    
    // (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
    ts1[0]=l0;ts1[1]=h0;ts1[2]=c0;ts1[3]=c1;
    ts1[4]=c2;ts1[5]=c3;ts1[6]=c4;ts1[7]=c5;
    ts2[0]=0;ts2[1]=l1;ts2[2]=h1;ts2[3]=h2;
    ts2[4]=h3;ts2[5]=h4;ts2[6]=h5;ts2[7]=h6;
    bigAdd2(ts1,ts2);
    ts2[0]=0;ts2[1]=0;ts2[2]=l2;ts2[3]=l3;
    ts2[4]=l4;ts2[5]=l5;ts2[6]=l6;ts2[7]=l7;
    bigAdd2(ts1,ts2);
    ts2[0]=c6;ts2[1]=c7;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigSub2(ts1,ts2);
    ts2[0]=h7;ts2[1]=0;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigSub2(ts1,ts2);
    
    //(x7*y7)r^6
    mulLong(xs[7],ys[7],rs);
    l6=rs[0];
    h6=rs[1];
    c6=(short)rs[2];
    
    //(x6*y7+x7*y6)r^5
    mulLong(xs[6],ys[7],rs);    
    l5=rs[0];
    h5=rs[1];
    c5=(short)rs[2];
    mulLong(xs[7],ys[6],rs);
    smallAdd(&l5,&h5,&c5,&rs[0],&rs[1],&rs[2]);
    
    //(x5*y7+x6*y6+x7*y5)r^4
    mulLong(xs[5],ys[7],rs);    
    l4=rs[0];
    h4=rs[1];
    c4=(short)rs[2];
    mulLong(xs[6],ys[6],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[5],rs);
    smallAdd(&l4,&h4,&c4,&rs[0],&rs[1],&rs[2]);    
    
    //(x4*y7+x5*y6+x6*y5+x7*y4)r^3
    mulLong(xs[4],ys[7],rs);    
    l3=rs[0];
    h3=rs[1];
    c3=(short)rs[2];
    mulLong(xs[5],ys[6],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    mulLong(xs[6],ys[5],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);    
    mulLong(xs[7],ys[4],rs);
    smallAdd(&l3,&h3,&c3,&rs[0],&rs[1],&rs[2]);
    
    //(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
    mulLong(xs[3],ys[7],rs);    
    l2=rs[0];
    h2=rs[1];
    c2=(short)rs[2];
    mulLong(xs[4],ys[6],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[5],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[4],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[3],rs);
    smallAdd(&l2,&h2,&c2,&rs[0],&rs[1],&rs[2]);
    
    //(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
    mulLong(xs[2],ys[7],rs);    
    l1=rs[0];
    h1=rs[1];
    c1=(short)rs[2];
    mulLong(xs[3],ys[6],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[5],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[4],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[3],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[2],rs);
    smallAdd(&l1,&h1,&c1,&rs[0],&rs[1],&rs[2]);
    
    //(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
    mulLong(xs[1],ys[7],rs);    
    l0=rs[0];
    h0=rs[1];
    c0=(short)rs[2];
    mulLong(xs[2],ys[6],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[3],ys[5],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[4],ys[4],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[5],ys[3],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[6],ys[2],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    mulLong(xs[7],ys[1],rs);
    smallAdd(&l0,&h0,&c0,&rs[0],&rs[1],&rs[2]);
    
    //(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
    ts2[0]=l0;ts2[1]=h0;ts2[2]=c0;ts2[3]=c1;
    ts2[4]=c2;ts2[5]=c3;ts2[6]=c4;ts2[7]=c5;
    bigSub2(ts1,ts2);
    ts2[0]=0;ts2[1]=l1;ts2[2]=h1;ts2[3]=h2;
    ts2[4]=h3;ts2[5]=h4;ts2[6]=h5;ts2[7]=h6;
    bigSub2(ts1,ts2);
    ts2[0]=0;ts2[1]=0;ts2[2]=l2;ts2[3]=l3;
    ts2[4]=l4;ts2[5]=l5;ts2[6]=l6;ts2[7]=0;
    bigSub2(ts1,ts2);
    ts2[0]=c6;ts2[1]=0;ts2[2]=0;ts2[3]=0;
    ts2[4]=0;ts2[5]=0;ts2[6]=0;ts2[7]=0;
    bigAdd2(ts1,ts2);
    us[0]=ts1[0];
    us[1]=ts1[1];
    us[2]=ts1[2];
    us[3]=ts1[3];
    us[4]=ts1[4];
    us[5]=ts1[5];
    us[6]=ts1[6];
    us[7]=ts1[7];
}


__global__ void batchBigMul(unsigned long r, unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
  __shared__ unsigned long xsm[TN*8];  
  __shared__ unsigned long ysm[TN*8]; 
  __shared__ unsigned long usm[TN*8]; 
  int tid = blockIdx.x*blockDim.x + threadIdx.x;   
  short i, pos, t;
  unsigned long *xd, *yd, *ud, *xm, *ym, *um; 
  
  
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  yd = (unsigned long*)((char*)ys + tid*sizeof(unsigned long)*8);
  ud = (unsigned long*)((char*)us + tid*sizeof(unsigned long)*8);
  __syncthreads();
  
  
  pos = threadIdx.x * 8;
  xm = &xsm[pos];
  ym = &ysm[pos];
  um = &usm[pos];
  
  for(i=0;i<8;i++)
  {
    xm[i]=xd[i];
  }  
  __syncthreads();
  for(i=0;i<8;i++)
  {
    ym[i]=yd[i];
  }
  __syncthreads();
  
  for(t=0;t<100;t++)
  {
    bigMul(xm,ym,um);    
  }
  __syncthreads();
  for(i=0;i<8;i++)
  {
    ud[i]=um[i];
  }
  
  /*
  for(t=0;t<100;t++)
  {
    bigMul(xd,yd,ud);    
  }
  */
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
  
  batchBigMul<<<BN,TN>>>(R,xs_d,ys_d,us_d);
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



