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

typedef struct
{
   unsigned long l,h;
   short c;
} prod;

typedef struct
{
   unsigned long u0;
   short u1;
} urform;


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
      us[0]=ULMAX;
      for(i=1;i<8;i++)
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
      us[0]=ULMAX;
      for(i=1;i<8;i++)
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


__device__ void mulTwoNum(unsigned long x, unsigned long y, prod *s)
{
  unsigned int x0,x1,y0,y1;
  unsigned long a0,a1,t,t0;
  s->c=0;
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
  a0>=R?s->h=1:s->h=0;
  s->l=(s->h>0?a0-R:a0);  //a0 finish
  s->h+=a1;
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
  s->h-=t;
  if(s->l>=t0)
  {
    s->l-=t0;
  }
  else
  {
    s->l=R-t0+s->l;
    s->h--;
  }
}

__device__ void mulLong(unsigned long x, unsigned long y, prod *s)
{
  if(x==R&&y==R)
  {
    s->l=0;s->h=0;s->c=1;
  }
  else if(x==R&&y!=R)
  {
    s->l=0;s->h=y;s->c=0;
  }
  else if(x!=R&&y!=R)
  {
    s->l=0;s->h=x;s->c=0;
  }
  else if(x<=SQRTR && y<=SQRTR)
  {
    s->l=x*y;s->h=0;s->c=0;
  }
  else
  {
    mulTwoNum(x,y,s);
  }
}

__device__ void tinyAdd(urform *u, unsigned long a)
{
  short c=0;
  unsigned long s=0;
  
  s=u->u0+a;
  s<u->u0 || s<a?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  u->u0=s;  
  u->u1+=c;
}

__device__ void tinySub(urform *u, unsigned long a)
{
  if(u->u0>=a)
  {
    u->u0-=a;
  }
  else
  {
    u->u0+=(R-a);
    u->u1--;
  }
}

__device__ void smallAdd(prod *p1, prod *p2)
{
  short c=0;
  unsigned long s=0;
  
  s=p1->l+p2->l;
  s<p1->l || s<p2->l?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  p1->l=s;
  
  p2->h=p2->h+c;  //h1<r<2^64-1. This means no overflow
  s=p1->h+p2->h;
  s<p1->h||s<p2->h?c=1:c=0;
  c>0?s=s+RC:s=s;
  if(s>=R)
  {
    s=s-R;
    c=1;
  }
  p1->h=s;
  
  p1->c+=(p2->c+c);
}

__device__ void smallSub(prod *p1, prod *p2)
{
  short c=0;
  
  if(p1->l>=p2->l)
  {
    p1->l-=p2->l;
  }
  else
  {
    c=1;
    p1->l+=(R-p2->l);
  }
  p2->h+=c;
  if(p1->h>=p2->h)
  {
    c=0;
    p1->h-=p2->h;
  }
  else
  {
    c=1;
    p1->h+=(R-p2->h);
  } 
  
  p1->c=p1->c-p2->c-c;
}

__device__ void bigMul(unsigned long *xs, unsigned long *ys, unsigned long *us) 
{
    urform u[8];  //8*8+2*8
    prod  pr1, pr2;  //20+20
    short c;
    
    //x_0y_0-x_1y_7-x_2y_6-x_3y_5-x_4y_4-x_5y_3-x_6y_2-x_7y_1
    mulLong(xs[0],ys[0],&pr1);
    mulLong(xs[1],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[2],ys[6],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[3],ys[5],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[4],ys[4],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[5],ys[3],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[6],ys[2],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[1],&pr2);
    smallSub(&pr1,&pr2); //[l0,h0,c0] stored into pr1
    u[0].u0=pr1.l;
    u[1].u1=pr1.h;
    pr1.c>0?tinyAdd(&u[2],pr1.c):tinySub(&u[2],-pr1.c);
    
    
    //x_0y_1+x_1y_0-x_2y_7-x_3y_6-x_4y_5-x_5y_4-x_6y_3-x_7y_2
    mulLong(xs[0],ys[1],&pr1);    
    mulLong(xs[1],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[3],ys[6],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[4],ys[5],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[5],ys[4],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[6],ys[3],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[2],&pr2);
    smallSub(&pr1,&pr2);//[l1,h1,c1] stored into pr1
    tinyAdd(&u[1],pr1.l);
    tinyAdd(&u[2],pr1.h);
    pr1.c>0?tinyAdd(&u[3],pr1.c):tinySub(&u[3],-pr1.c);
    
    
    //x_0y_2+x_1y_1+x_2y_0-x_3y_7-x_4y_6-x_5y_5-x_6y_4-x_7y_3
    mulLong(xs[0],ys[2],&pr1);    
    mulLong(xs[1],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[4],ys[6],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[5],ys[5],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[6],ys[4],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[3],&pr2);
    smallSub(&pr1,&pr2);//[l2,h2,c2] stored into pr1
    tinyAdd(&u[2],pr1.l);
    tinyAdd(&u[3],pr1.h);
    pr1.c>0?tinyAdd(&u[4],pr1.c):tinySub(&u[4],-pr1.c);
    
    //x_0y_3+x_1y_2+x_2y_1+x_3y_0-x_4y_7-x_5y_6-x_6y_5-x_7y_4
    mulLong(xs[0],ys[3],&pr1);    
    mulLong(xs[1],ys[2],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[4],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[5],ys[6],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[6],ys[5],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[4],&pr2);
    smallSub(&pr1,&pr2);//[l3,h3,c3] stored into pr1
    tinyAdd(&u[3],pr1.l);
    tinyAdd(&u[4],pr1.h);
    pr1.c>0?tinyAdd(&u[5],pr1.c):tinySub(&u[5],-pr1.c);
    
    //x_0y_4+x_1y_3+x_2y_2+x_3y_1+x_4y_0-x_5y_7-x_6y_6-x_7y_5
    mulLong(xs[0],ys[4],&pr1);    
    mulLong(xs[1],ys[3],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[2],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[4],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[5],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[6],ys[6],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[5],&pr2);
    smallSub(&pr1,&pr2);//[l4,h4,c4] stored into pr1
    tinyAdd(&u[4],pr1.l);
    tinyAdd(&u[5],pr1.h);
    pr1.c>0?tinyAdd(&u[6],pr1.c):tinySub(&u[6],-pr1.c);
    
    //x_0y_5+x_1y_4+x_2y_3+x_3y_2+x_4y_1+x_5y_0-x_6y_7-x_7y_6
    mulLong(xs[0],ys[5],&pr1);    
    mulLong(xs[1],ys[4],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[3],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[2],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[4],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[5],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[6],ys[7],&pr2);
    smallSub(&pr1,&pr2);
    mulLong(xs[7],ys[6],&pr2);
    smallSub(&pr1,&pr2);//[l5,h5,c5] stored into pr1
    tinyAdd(&u[5],pr1.l);
    tinyAdd(&u[6],pr1.h);
    pr1.c>0?tinyAdd(&u[7],pr1.c):tinySub(&u[7],-pr1.c);
    
    //x_0y_6+x_1y_5+x_2y_4+x_3y_3+x_4y_2+x_5y_1+x_6y_0-x_7y_7
    mulLong(xs[0],ys[6],&pr1);    
    mulLong(xs[1],ys[5],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[4],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[3],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[4],ys[2],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[5],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[6],ys[0],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[7],ys[7],&pr2);
    smallSub(&pr1,&pr2);//[l6,h6,c6] stored into pr1
    tinyAdd(&u[6],pr1.l);
    tinyAdd(&u[7],pr1.h);
    pr1.c>0?tinySub(&u[0],pr1.c):tinyAdd(&u[0],-pr1.c);
    
    //x_0y_7+x_1y_6+x_2y_5+x_3y_4+x_4y_3+x_5y_2+x_6y_1+x_7y_0
    mulLong(xs[0],ys[7],&pr1);    
    mulLong(xs[1],ys[6],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[2],ys[5],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[3],ys[4],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[4],ys[3],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[5],ys[2],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[6],ys[1],&pr2);
    smallAdd(&pr1,&pr2);
    mulLong(xs[7],ys[0],&pr2);
    smallAdd(&pr1,&pr2);//[l7,h7,c7] stored into pr1
    tinyAdd(&u[7],pr1.l);
    tinySub(&u[0],pr1.h);
    pr1.c>0?tinySub(&u[0],pr1.c):tinyAdd(&u[0],-pr1.c);
    
    c=0;
    if(u[7].u1>0&&u[7].u1>u[0].u0)
    {
      us[0]=R-u[7].u1+u[0].u0;
      u[0].u1--;
    }
    else
    {
      us[0]=u[0].u0-u[7].u1;
      if(us[0]>=R)
      {
        us[0]-=R;
        u[0].u1++;
      }
    }
    if(u[0].u1<0&&(-u[0].u1)>u[0].u0)
    {
      us[1]=R+u[0].u1+u[1].u0;
      u[1].u1--;
    }
    else
    {
      us[1]=u[0].u1+u[1].u0;
      if(us[1]>=R)
      {
        us[1]-=R;
        u[1].u1++;
      }
    }
    if(u[1].u1<0&&(-u[1].u1)>u[2].u0)
    {
      us[2]=R+u[1].u1+u[2].u0;
      u[2].u1--;
    }
    else
    {
      us[2]=u[1].u1+u[2].u0;
      if(us[2]>=R)
      {
        us[2]-=R;
        u[2].u1++;
      }
    }
    if(u[2].u1<0&&(-u[2].u1)>u[3].u0)
    {
      us[3]=R+u[2].u1+u[3].u0;
      u[3].u1--;
    }
    else
    {
      us[3]=u[2].u1+u[3].u0;
      if(us[3]>=R)
      {
        us[3]-=R;
        u[3].u1++;
      }
    }
    if(u[3].u1<0&&(-u[3].u1)>u[4].u0)
    {
      us[4]=R+u[3].u1+u[4].u0;
      u[4].u1--;
    }
    else
    {
      us[4]=u[3].u1+u[4].u0;
      if(us[4]>=R)
      {
        us[4]-=R;
        u[4].u1++;
      }
    }
    if(u[4].u1<0&&(-u[4].u1)>u[5].u0)
    {
      us[5]=R+u[4].u1+u[5].u0;
      u[5].u1--;
    }
    else
    {
      us[5]=u[4].u1+u[5].u0;
      if(us[5]>=R)
      {
        us[5]-=R;
        u[5].u1++;
      }
    }
    if(u[5].u1<0&&(-u[5].u1)>u[6].u0)
    {
      us[6]=R+u[5].u1+u[6].u0;
      u[6].u1--;
    }
    else
    {
      us[6]=u[5].u1+u[6].u0;
      if(us[6]>=R)
      {
        us[6]-=R;
        u[6].u1++;
      }
    }
    if(u[6].u1<0&&(-u[6].u1)>u[7].u0)
    {
      us[7]=R+u[6].u1+u[7].u0;
      u[7].u1--;
    }
    else
    {
      us[7]=u[6].u1+u[7].u0;
      if(us[7]>R)
      {
        us[7]-=R;
        c=1;
      }
    }
    if(c>0)
    {
      for(c=0;c<7;c++)
      {
        if(us[c]!=0)
        {
          us[c]--;
          break;
        }
        else
        {
          us[c]=R-1;
        }
      }
    }
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



