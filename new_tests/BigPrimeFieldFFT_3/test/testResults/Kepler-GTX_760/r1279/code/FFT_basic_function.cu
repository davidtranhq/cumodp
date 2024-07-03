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
#include<math.h>

//N the FFT size, where each element contains eight unsigned integers.
#define N 16777216
//NN the size of unsigned numbers, is equal to N*8 where each element contains eight unsigned integers.
#define NN 134217728
//BN block number
#define BLOCK_DIM 512
//TN thread number in a block
#define THREAD_DIM 512 
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//sqrt(R) 3037000502
#define SQRTR 3037000502
//ULMAX=2^64-1
#define ULMAX 18446744073709551615

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

//These constant variable is the offsets for FFT16

__constant__ short ind2[8]={0,0,0,0,4,4,4,4};
__constant__ short ind3[8]={0,4,2,6,0,4,2,6};
__constant__ short ind4[8]={0,4,2,6,1,5,3,7};
__constant__ short ind5[16]={0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};

/**
 * the function for calculating us=xs+ys
 */
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
      us[pos] --;
    }            
    else
    {
      us[7]=R;
    }  
  }  
}

/**
 * the function for calculating xs=xs+ys
 */
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
      us[pos] --;
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

/**
 * the function for calculating us=xs-ys
 */
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
      //all elements are R-1
      us[7]=R;
      for(i=0;i<7;i++)
      {
        us[i]=0;
      }
    }  
  }  
}

/**
 * the function for calculating xs=xs-ys
 */
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

/**
 * the function for doing right shift $sn$ operation on xs
 */
__device__ void cyclicShift(unsigned long *xs, short sn)
{
  short i=0,j=0;
  unsigned long ts[8]={0};
  if(sn<=0)
  {
    return;
  }
  if(sn>8)
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
  bigSub2(xs,ts);
}


/**
 * the function for multipling x*y into s, where 0<=x,y<r, s is a struct and s=l+h*r+c*r^2
 */
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

/**
 * the function for multipling x*y into s, where 0<=x,y<=r, s is a struct and s=l+h*r+c*r^2. It deals several special cases and calls mulTwoNum function to finish multiplication
 */
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

/**
 * the function for calculating u = u+a,  where u is a struct and u=u0+u1*r. a is an unsigned number. 
 */
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

/**
 * the function for calculating u = u-a,  where u is a struct and u=u0+u1*r;
 */
 
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

/**
 * the function for calculating p1=p1+p2,  where p1,p2 are struct and p1 = l1+h1*r+c1*r^2, p2=l2+h2*r+c2*r^2;
 */

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


/**
 * the function for calculating p1=p1-p2,  where p1,p2 are struct and p1 = l1+h1*r+c1*r^2, p2=l2+h2*r+c2*r^2;
 */
 
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

/**
 * the function for calculating the multiplication of xs and ys, xs=xs*ys,  
 */

__device__ void bigMul(unsigned long *xs, unsigned long *ys) 
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
      xs[0]=R-u[7].u1+u[0].u0;
      u[0].u1--;
    }
    else
    {
      xs[0]=u[0].u0-u[7].u1;
      if(xs[0]>=R)
      {
        xs[0]-=R;
        u[0].u1++;
      }
    }
    if(u[0].u1<0&&(-u[0].u1)>u[0].u0)
    {
      xs[1]=R+u[0].u1+u[1].u0;
      u[1].u1--;
    }
    else
    {
      xs[1]=u[0].u1+u[1].u0;
      if(xs[1]>=R)
      {
        xs[1]-=R;
        u[1].u1++;
      }
    }
    if(u[1].u1<0&&(-u[1].u1)>u[2].u0)
    {
      xs[2]=R+u[1].u1+u[2].u0;
      u[2].u1--;
    }
    else
    {
      xs[2]=u[1].u1+u[2].u0;
      if(xs[2]>=R)
      {
        xs[2]-=R;
        u[2].u1++;
      }
    }
    if(u[2].u1<0&&(-u[2].u1)>u[3].u0)
    {
      xs[3]=R+u[2].u1+u[3].u0;
      u[3].u1--;
    }
    else
    {
      xs[3]=u[2].u1+u[3].u0;
      if(xs[3]>=R)
      {
        xs[3]-=R;
        u[3].u1++;
      }
    }
    if(u[3].u1<0&&(-u[3].u1)>u[4].u0)
    {
      xs[4]=R+u[3].u1+u[4].u0;
      u[4].u1--;
    }
    else
    {
      xs[4]=u[3].u1+u[4].u0;
      if(xs[4]>=R)
      {
        xs[4]-=R;
        u[4].u1++;
      }
    }
    if(u[4].u1<0&&(-u[4].u1)>u[5].u0)
    {
      xs[5]=R+u[4].u1+u[5].u0;
      u[5].u1--;
    }
    else
    {
      xs[5]=u[4].u1+u[5].u0;
      if(xs[5]>=R)
      {
        xs[5]-=R;
        u[5].u1++;
      }
    }
    if(u[5].u1<0&&(-u[5].u1)>u[6].u0)
    {
      xs[6]=R+u[5].u1+u[6].u0;
      u[6].u1--;
    }
    else
    {
      xs[6]=u[5].u1+u[6].u0;
      if(xs[6]>=R)
      {
        xs[6]-=R;
        u[6].u1++;
      }
    }
    if(u[6].u1<0&&(-u[6].u1)>u[7].u0)
    {
      xs[7]=R+u[6].u1+u[7].u0;
      u[7].u1--;
    }
    else
    {
      xs[7]=u[6].u1+u[7].u0;
      if(xs[7]>R)
      {
        xs[7]-=R;
        c=1;
      }
    }
    if(c>0)
    {
      for(c=0;c<7;c++)
      {
        if(xs[c]!=0)
        {
          xs[c]--;
          break;
        }
        else
        {
          xs[c]=R-1;
        }
      }
    }
}

/**
 * This function is to calculate xs=xs*w^pid, and ws is the pointer of w array. The part of w^pid can be executed into cyclic shift.
 */
 
__device__ void mulOmega(unsigned long *xs, int pid, unsigned long *ws)
{
  int rp,wp;
  unsigned long *A;
  if(pid<=0)
  {
    return;
  }
  
  rp=pid>>20; //w^1048576=r  0<=rp<16  caution: depend on N.
  wp=pid-(rp<<20); //wp<1048576
  
  rp>8?cyclicShift(xs,8):cyclicShift(xs,rp);
  rp>8?cyclicShift(xs,rp-8):cyclicShift(xs,0);  
  
  if(wp>0)
  {
    A = (unsigned long*)((char*)ws + (wp<<6)); 
    bigMul(xs,A);
  }
}

/**
 * This function is the kernel to do ysm=FFT16(xsm). It uses one block to do $16$ FFT 16. All data has been copied into shared memory
  
 */

__device__ void FFT16(unsigned long *xsm, unsigned long *ysm) 
{
  short tid = threadIdx.x; 
  short wid = (tid >> 5);  //warp no. [0,1,...,7]
  short bwid = (tid >> 6); //big warp no.[0,1,...,3]
  short sf = wid - (bwid << 1); //sub flag for odd warp[0,1]
  short wpn = ((tid-(wid<<5))>>3); //(tid-(wid*32))/8  [0,1,2,3] 
  short wpid = (tid-((tid>>3)<<3)); //in [0,1,...,7]
  short posid = 0; //in[0,1,...,15]
  short i, pos, pos1;
  unsigned long *xf, *yf;   
  int pos2 = ((bwid<<6)+(wpn<<4))<<3;
  
  xf = (unsigned long*)((char*)xsm + pos2*sizeof(unsigned long)); //
  yf = (unsigned long*)((char*)ysm + pos2*sizeof(unsigned long));
  
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
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //second round
  posid=(wpid>=4?wpid+4:wpid);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+4)<<3],ind2[wpid]);
  }
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+4)<<3],&yf[(posid+4)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+4)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //third round
  posid=(wpid>=4?((wpid-4)<<2)+1:(wpid<<2));
  if(sf>0)
  {
    cyclicShift(&xf[(posid+2)<<3],ind3[wpid]);
  }
  
  __syncthreads(); 
  if(sf>0)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+2)<<3],&yf[(posid+2)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+2)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  if(sf)
  {
    pos=(wpid+8)<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  else
  {
    pos=wpid<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos+i];
    }
  }
  __syncthreads(); 
  
  
  //fourth
  posid=(wpid<<1);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+1)<<3],ind4[wpid]);
  }
  __syncthreads(); 
  if(sf)
  {
    //bigSub
    bigSub(&xf[posid<<3],&xf[(posid+1)<<3],&yf[(posid+1)<<3]);
  }
  else
  {
    //bigAdd
    bigAdd(&xf[posid<<3],&xf[(posid+1)<<3],&yf[posid<<3]);
  }
  __syncthreads();  
  
  
  if(sf)
  {
    posid=wpid+8;
    pos=posid<<3;
    pos1=ind5[posid]<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  else
  {
    posid=wpid;
    pos=posid<<3;
    pos1=ind5[posid]<<3;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  __syncthreads(); 
  
}


/**
 * uniTwiddle is to twiddle a list of matrices xs, then store back to xs.
 * xs: the pointer of $s$ matrices.
 * s: the matrix number
 * dnm: the data number of one matrix. dnm=he*wi
 * he: the height of matrix
 * wi: the width of matrix
 * I_s @ D_{he, wi}
 */
__global__ void uniTwiddle(unsigned long *xs, unsigned long *ws, short sb, short dnmb, short heb, short wib) 
{
    __shared__ unsigned long block[4096];  
    int idx = (blockIdx.x<<9) + threadIdx.x;
    int tid = threadIdx.x;
    short dn = 1;  
    unsigned long *A;
    int i,mno,pos,pos1,nh,nw,j,powind;
    N>=262144?(dn=(N>>18)):dn=1;
    
    for(i=0;i<dn&&idx<N;i++)
    {        
      pos=(tid<<3); //tid*8;
      A = (unsigned long *)((char*)xs + (idx<<6));
      for(j=0;j<8;j++)
      {
        block[pos++]=A[j];
      }
      __syncthreads(); 
      
      mno=idx>>dnmb; //mno=(idx)/dnm; 
      pos1=idx-(mno<<dnmb); 
      nh=pos1>>wib; 
      nw=pos1-(nh<<wib); 
      powind=nh*nw;
      
      if(powind > 0)
      {
        pos=tid<<3; //tid*8;
        mulOmega(&block[pos],powind<<sb,ws);
        for(j=0;j<8;j++)
        {
          A[j]=block[pos++];
        }
      }  
      idx+=262144; //i*262144; 
      __syncthreads(); 
    }    
}

/**
 * uniTranspose is to transpose a list of matrices xs, then store into ys.
 * xs: the pointer of $s$ matrices.
 * s: the matrix number
 * dnmb: log_2^(the data number of one matrix. dnm=he*wi)
 * heb: log_2^(the height of matrix)
 * wib: log_2^(the width of matrix)
 */
__global__ void uniTranspose(unsigned long *xs, unsigned long *ys, int s, short dnmb, short heb, short wib) 
{
    __shared__ unsigned long block[512];  //512*8   4KB
    
    int id = (blockIdx.x<<9) + threadIdx.x;
    int tid = threadIdx.x;
    short dn;  //the number each thread need deal.
    int idx,mno,pos1,h,w;
    short i,j;
    N>=32768?(dn=(N>>15)):dn=1;
    
    for(i=0;i<dn&&id<NN;i++)
    {
      block[tid]=xs[id]; //value
      __syncthreads();
      idx=id>>3;  
      j=id-(idx<<3);
      mno=idx>>dnmb; //mno=idx/dnm; //submatrix no, start from 0
      pos1=idx-(mno<<dnmb); //position in the current submatrix
      h=pos1>>wib; //old height
      w=pos1-(h<<wib); //new width
      idx=(((idx-pos1)+(w<<heb)+h)<<3)+j;
      ys[idx]=block[tid];      
      id+=262144;  
      __syncthreads();
    }    
}

/**
 * the function to copy the array ys to xs. 
 */
__global__ void uniCopy(unsigned long *xs, unsigned long *ys) 
{
    int id = (blockIdx.x<<9) + threadIdx.x;
    short dn,i;  //dn the number each thread need deal.
    N>=32768?(dn=(N>>15)):dn=1;
    
    for(i=0;i<dn&&id<NN;i++)
    {
      xs[id]=ys[id];
      __syncthreads();
      id+=262144;
    }   
}
/**
 * the kernel for calculating FFT16
 */

__global__ void data_fft16_kernel(unsigned long *xs) 
{
  short bid = blockIdx.x;
  short tid = threadIdx.x;
  int needBlock=(N>>8); //N/256;
  short bn = ((needBlock>=BLOCK_DIM)?needBlock/BLOCK_DIM:1);  //the number of block
  short i,j;
  //unsigned long *A, *B;
  unsigned long *A;
  __shared__ unsigned long xsm[2048];  //256*8
  __shared__ unsigned long ysm[2048];  
  int pos,pos1,bno; //tid*8
  
  bno=bid;
  pos=((bno<<8)+tid)<<3;
  for(i=0;i<bn;i++)
  {
    if(bno>=needBlock)
    {
      break;
    }
    pos1 = tid<<3;
    //pos = (bno<<8)+tid;  //each block deal 256 big numbers.
    A = (unsigned long *)((char*)xs + (pos<<3));
    
    for(j=0;j<8;j++)
    {
      xsm[pos1++]=A[j];
    }
    __syncthreads();  
    
    FFT16(xsm,ysm); //ysm is a temp array
    __syncthreads();  
    
    pos1 = (tid<<3);
    for(j=0;j<8;j++)
    {
      A[j]=xsm[pos1++];
    }
    bno+=512;
    pos+=1048576;  //512block*256bignumber*8numbers
  }
}


