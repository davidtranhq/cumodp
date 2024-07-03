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
//BN block number
#define BLOCK_DIM 512
//TN thread number in a block
#define THREAD_DIM 512 
//input number count
#define INC 256
//R  2^63+2^34
#define R 9223372054034644992
//RC R complement  RC=2^64-R
#define RC 9223372019674906624
//sqrt(R) 3037000502
#define SQRTR 3037000502
//ULMAX=2^64-1
#define ULMAX 18446744073709551615

__constant__ short ind2[8]={0,0,0,0,4,4,4,4};
__constant__ short ind3[8]={0,4,2,6,0,4,2,6};
__constant__ short ind4[8]={0,4,2,6,1,5,3,7};
__constant__ short ind5[16]={0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};


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
      us[7]=R;
    }  
  }  
}

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
      us[7]=ULMAX;
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
      us[7]=R;
      for(i=0;i<7;i++)
      {
        us[i]=0;
      }
    }  
  }  
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

//0<=sn<=8
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

//[l,h,c]
__device__ void mulLong(unsigned long x, unsigned long y, unsigned long *s)
{
  short x1,y1;
  unsigned long l,h,c,x0,y0,v2,v5,v9,v10,v11,v14,v15,v16,v17,q,t;
  unsigned long a0,a1,b0,b1,c0,c1,c1prime,d0,d1,d2,e0,e1;
  
  if(x<=SQRTR && y<=SQRTR)
  {
    s[0]=x*y;
    s[1]=0;
    s[2]=0;
    return;
  }
  
  x1=(x>=R?1:0);
  x0=(x1>0?x-R:x);
  y1=(y>=R?1:0);
  y0=(y1>0?y-R:y);
  
  v2=x0*y1; //[0,v2,0];
  v5=x1*y0; //[0,v5,0];
  v9=x1*y1; //[0,0,1];
  
  c=v9;
  l=0;
  h=v5+v2;
  h<v5||h<v2?(c=c+1):(c=c);
  c>v9?(h=h+RC):(h=h);
  
  if(x0<=SQRTR&&y0<=SQRTR)
  {
    s[0]=x0*y0;
    s[1]=h;
    s[2]=c;
    return;
  }
  
  //lhc
  //x0*y0
  a1=x0>>32;
  a0=x0-(a1<<32);
  b1=y0>>32;
  b0=y0-(b1<<32);
  
  c0=0;
  c1=a1*b1;
  
  t=a0*b1;
  q=t>>32;
  t=(t-(q<<32))<<32;
  c1+=q;
  c0+=t;  //safe
  
  t=a1*b0;
  q=t>>32;
  t=(t-(q<<32))<<32;
  c1+=q;
  q=c0+t;               //here, is not related to r.
  q<c0||q<t?(c1++):(c1=c1);  //c0=c0+t and carry, safe
  c0=q;
  
  t=a0*b0;
  q=c0+t;
  q<c0||q<t?(c1++):(c1=c1);  //Now we finish [c0,c1]=x0*y0
  c0=q;
  
  c1prime=c1<<1;
  
  c0>=R?(v11=1):(v11=0);
  v11>0?(v10=c0-R):(v10=c0);
  //v12=0;
  
  q=l+v10;  //[l,h,c] + [v10,v11,0]
  q<l||q<v10?(v11=v11+1):(v11=v11);
  q<l||q<v10?(l=q+RC):(l=q);
  if(l>=R)
  {
    l=l-R;
    v11++;
  }
  q=h+v11;
  q<h||q<v11?(c=c+1):(c=c);
  q<h||q<v11?(h=q+RC):(h=q);
  if(h>=R)
  {
    h=h-R;
    c++;
  }
  //v13=0;
  c1prime>=R?(v15=1):(v15=0);
  v15>0?(v14=c1prime-R):(v14=c1prime); //v13=0;
  
  q=h+v14;  //[l,h,c]+[0,v14,v15]
  q<h||q<v14?(c=c+v15+1):(c=c+v15);
  q<h||q<v14?(h=q+RC):(h=q);
  if(h>=R)
  {
    h=h-R;
    c++;
  }
  //[l,h,c]
  
  d1=c1prime>>29;
  d0=c1prime-(d1<<29);
  if(d0>=d1)
  {
    d2=d0-d1;
    e1=d2>>29;
    e0=d2-(e1<<29);
    e0>=e1?(v16=(e0-e1)<<34):(v16=R-(e1<<34)+(e0<<34));
    e0>=e1?(v17=e1+d1):(v17=e1+d1-1);
    /*
    if(e0>=e1)
    {
      v16=(e0-e1)<<34;
      v17=e1+d1;
    }
    else
    {
      v17=e1+d1-1;
      v16=R-(e1<<34)+(e0<<34);
    }
    */
  }
  else
  {
    //d1>d0
    d2=d1-d0;
    e1=d2>>29;
    e0=d2-(e1<<29);
    e0>=e1?(v16=R-((e0-e1)<<34)):(v16=(e1-e0)<<34);
    e0>=e1?(v17=d1-e1-1):(v17=d1-e1);
    /*
    if(e0>=e1)
    {
      v16=R-((e0-e1)<<34);
      v17=d1-e1-1;
    }
    else
    {
      v16=(e1-e0)<<34;
      v17=d1-e1;
    }
    */
  }
  //[l,h,c]-[v16,v17,0]
  //q
  q=0;
  if(l>=v16)
  {
    l=l-v16;
  }
  else
  {
    l=R-v16+l;
    q=1;
  }
  //t
  if(h<q+v17)
  {
    c=c-1;
    h=R-q-v17+h;
  }
  else
  {
    h=h-q-v17;
  }  
  s[0]=l;
  s[1]=h;
  s[2]=c;
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

__device__ void bigMul(unsigned long *xs, unsigned long *ys) 
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
    xs[0]=ts1[0];
    xs[1]=ts1[1];
    xs[2]=ts1[2];
    xs[3]=ts1[3];
    xs[4]=ts1[4];
    xs[5]=ts1[5];
    xs[6]=ts1[6];
    xs[7]=ts1[7];
}


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
  
  /*
  if(wp<=0)
  {
    return;
  }
  rp=wp>>16;      //0<=rp<16
  wp=wp-(rp<<16); //wp<65536
  if(rp>0)
  {
    bigMul(xs,W8[(rp-1)<<3]);
  }
  rp=wp>>12;
  wp=wp-(rp<<12);
  if(rp>0)
  {
    bigMul(xs,W12[(rp-1)<<3]);
  }
  rp=wp>>8;
  wp=wp-(rp<<8);
  if(rp>0)
  {
    bigMul(xs,W16[(rp-1)<<3]);
  }
  rp=wp>>4;
  wp=wp-(rp<<4);
  if(rp>0)
  {
    bigMul(xs,W20[(rp-1)<<3]);
  }
  if(wp>0)
  {
    bigMul(xs,W24[(rp-1)<<3]);
  }
  */
  if(wp>0)
  {
    //A = (unsigned long*)((char*)ws + (wp<<6)*sizeof(unsigned long)); 
    A = (unsigned long*)((char*)ws + (wp<<6)); 
    //bigMul(xs,&ws[wp<<3]);  
    bigMul(xs,A);
  }
}



__device__ void FFT16(unsigned long *xsm, unsigned long *ysm) 
{
  //Now, we only use one block to do $16$ FFT 16. 
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
 * uniTwiddle is to twiddle a list of matrices xs, then store into ys.
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
      
      mno=idx>>dnmb; //mno=(idx)/dnm; //表示第几个子矩阵 从0开始。
      pos1=idx-(mno<<dnmb); //表示该元素在子矩阵中的位置。
      nh=pos1>>wib; //旧的高度索引
      nw=pos1-(nh<<wib); //旧的宽度索引        
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
 * dnm: the data number of one matrix. dnm=he*wi
 * he: the height of matrix
 * wi: the width of matrix
 */
__global__ void uniTranspose(unsigned long *xs, unsigned long *ys, int s, short dnmb, short heb, short wib) 
{
    __shared__ unsigned long block[4096];  //512*8   4096*8=32KB 
    
    // The size of shared memory in a block is 48kb. 
    // We use 32KB, so the number of big numbers one block deal is 512.
    
    int id = (blockIdx.x<<9) + threadIdx.x;
    int tid = threadIdx.x;
    short dn = 1;  //the number each thread need deal.
    unsigned long *A, *B;
    int idx,mno,pos,pos1,h,w,pos2;
    short i,j;
    //dn=(N+262143)>>18;
    N>=262144?(dn=(N>>18)):dn=1;
    
    idx=id;
    for(i=0;i<dn&&idx<N;i++)
    {
      pos=tid<<3; //tid*8;              
      A = (unsigned long *)((char*)xs + (idx<<6));
      for(j=0;j<8;j++)
      {
        block[pos++]=A[j];
      }
      __syncthreads();
        
      mno=idx>>dnmb; //mno=idx/dnm; //submatrix no, start from 0
      pos1=idx-(mno<<dnmb); //position in the current submatrix
      h=pos1>>wib; //old height
      w=pos1-(h<<wib); //new width
      pos2=(idx-pos1)+(w<<heb)+h;
        
      pos=tid<<3; //tid*8;
      B = (unsigned long *)((char*)ys + (pos2<<6));
      for(j=0;j<8;j++)
      {
        B[j]=block[pos++];
      }
      idx+=262144;  //idx=idx+i*262144;
      __syncthreads();
    }    
}

__global__ void uniCopy(unsigned long *xs, unsigned long *ys) 
{
    int idx = (blockIdx.x<<10) + threadIdx.x;
    short dn = 1;  //the number each thread need deal.
    int pos;
    short i,j;
    unsigned long *A, *B;
    N>=1048576?(dn=(N>>20)):dn=1;
    
    for(i=0;i<dn&&idx<N;i++)
    {
      pos=idx<<6;
      A = (unsigned long *)((char*)xs + pos);
      B = (unsigned long *)((char*)ys + pos);
      for(j=0;j<8;j++)
      {
        A[j]=B[j];
      }
      idx+=1048576;
    }   
}


void data_shuffle(int k, unsigned long *xs, unsigned long *ys)
{
  cudaEvent_t start, stop;
  float elapsedTime, time1;
  short i,j,dnmb,heb,wib;
  int s, he, wi, dnm;  //matrix number, height, width, data number for matrix
  wi=16;wib=4;
  s=1;  
  dnm=N;dnmb=k<<2;  //dnmb=24
  he=N/wi;heb=dnmb-4;
  //s*he*wi is equal to N
  //k=2  s=1  wi=16  he=16
  //k=3  s=1  wi=16  he=256
  //     s=16 wi=16  he=16
  
  j=0;
  elapsedTime=0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  for(i=0;i<=k-2;i++) //0,1,...,k-2
  {
    //xs transpose to ys
    //printf("call uniTranspose...\n");
    //printf("%d submatrix, each %d elements, %d * %d\n",s, dnm, he, wi);
    if(j==0)
    {
      j=1;
      cudaEventRecord(start, 0);
      uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(xs, ys, s, dnmb, heb, wib); 
      cudaThreadSynchronize();
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time1, start, stop);
      elapsedTime+=time1;
      time1=0;
    }
    else
    {
      j=0;
      cudaEventRecord(start, 0);
      uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(ys, xs, s, dnmb, heb, wib); 
      cudaThreadSynchronize();
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time1, start, stop);
      elapsedTime+=time1;
      time1=0;
    }    
    
    he=he>>4;   //he=he/16;
    heb=heb-4;
    s=s<<4;     //s=s*16;
    dnm=dnm>>4; //dnm=dnm/16;
    dnmb=dnmb-4;
  }
  if(j>0)
  {
    //copy ys to xs
    //printf("call uniCopy...\n");
    cudaEventRecord(start, 0);
    uniCopy<<<1024,1024>>>(xs,ys);
    cudaThreadSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time1, start, stop);
    elapsedTime+=time1;
    time1=0; 
  }
  //cudaMemset((void **)&ys, 0, sizeof(unsigned long)*8*N); 
  printf("the transposition in P1 costs %f millseconds\n", elapsedTime);
}

__global__ void data_fft16_kernel(unsigned long *xs) 
{
  short bid = blockIdx.x;
  short tid = threadIdx.x;
  int needBlock=(N>>8); //N/256;
  short bn = ((needBlock>=BLOCK_DIM)?needBlock/BLOCK_DIM:1);  //the number of block
  short i,j;
  //unsigned long *A, *B;
  unsigned long *A;
  __shared__ unsigned long xsm[INC<<3]; 
  __shared__ unsigned long ysm[INC<<3];  
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

void data_fft(unsigned long *xs)
{
  cudaEvent_t start, stop;
  float elapsedTime;
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  //calculate xs in FFT16 batch, and store into ys
  data_fft16_kernel<<<BLOCK_DIM,256>>>(xs);
  cudaThreadSynchronize();
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("the fft time is %f milliseconds\n", elapsedTime);
}


void data_twiddle(int k, unsigned long *xs, unsigned long *ys, unsigned long *ws)
{
  cudaEvent_t start, stop;
  float elapsedTime, twiddleTime, transposeTime;
  short i, sb, dnmb, heb, wib;
  int s;
  //s submatrix number
  
  elapsedTime=0;
  twiddleTime=0;
  transposeTime=0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  for(i=k-2;i>=0;i--)
  {
    //T1  I_{16^i} @ D_{16,16^{k-1-i}}
    //s*dnm=s*he*wi
    s=(int) pow(16,i);
    sb=i<<2;
    //dnm=N/s;
    dnmb=(k-i)<<2; 
    heb=4; //he=16;
    wib=dnmb-heb;
    cudaEventRecord(start, 0);
    uniTwiddle<<<BLOCK_DIM,THREAD_DIM>>>(xs, ws, sb, dnmb, heb, wib);
    cudaThreadSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    twiddleTime+=elapsedTime;
    elapsedTime=0;
    
    //T2  I_{16^i} @ L_{16^{k-1-i}}^{16^{k-i}}
    dnmb=(k-i)<<2; 
    //wi=(int) pow(16,k-1-i);
    wib=(k-1-i)<<2;
    heb=4;
    cudaEventRecord(start, 0);
    uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(xs, ys, s, dnmb, heb, wib); 
    cudaThreadSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    transposeTime+=elapsedTime;
    elapsedTime=0;
    
    
    //T3  I_{16^{k-1}} @ DTT16
    data_fft(ys);  //do fft on ys, store on ys. xs is a temp array
    
    //T4  I_{16^i} @ L_{16}^{16^{k-i}}
    //wi=16;
    wib=4;
    heb=dnmb-wib;
    
    cudaEventRecord(start, 0);
    uniTranspose<<<BLOCK_DIM,THREAD_DIM>>>(ys, xs, s, dnmb, heb, wib); 
    cudaThreadSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    transposeTime+=elapsedTime;
    elapsedTime=0;
  }
  printf("twiddle time in P3 is %f milliseconds\n", twiddleTime);
  printf("transpose time in P3 is %f milliseconds\n", transposeTime);
}

void FFT16777216(unsigned long *xs, unsigned long *ws)
{
  int k=6;
  unsigned long *xs_d, *ys_d, *ws_d;
  
  
  //initiation
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*N);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*N, cudaMemcpyHostToDevice);   cudaMalloc((void **)&ys_d, sizeof(unsigned long)*8*N);
  
  cudaMemset((void **)&ys_d, 0, sizeof(unsigned long)*8*N);
  
  cudaMalloc((void **)&ws_d, sizeof(unsigned long)*8*N/16);
  cudaMemcpy(ws_d, ws, sizeof(unsigned long)*8*N/16, cudaMemcpyHostToDevice);
  
  //P1
  data_shuffle(k, xs_d, ys_d); 
  //printf("P1 end...\n");
  
  //P2
  data_fft(xs_d);
  //printf("P2 end...\n");
  
  
  //P3
  data_twiddle(k, xs_d, ys_d, ws_d);
  //printf("P3 end...\n");
    
  cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*N, cudaMemcpyDeviceToHost);
  
  printf("we have done xs = FFT16777216(xs) in 1 times.\n");
  
  cudaFree(xs_d);	
  cudaFree(ys_d);	
  cudaFree(ws_d);	
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1, *fp2, *fp3;	 
  unsigned long *xs,*ws; 
  
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*N));
  ws=(unsigned long *)malloc((sizeof(unsigned long)*8*(N/16)));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_16777216_input.dat"); 
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  memset(xs,0,sizeof(unsigned long)*8*N);
  fread(xs,sizeof(unsigned long),8*N,fp1);
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "W_16777216.dat"); 
  if((fp3=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp3);
  	exit(-1);
  } 
  memset(ws,0,sizeof(unsigned long)*8*(N/16));
  fread(ws,sizeof(unsigned long),8*N/16,fp3);
  
  
  FFT16777216(xs,ws);
    
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_16777216_output.dat"); //INC*8 (unsigned long) data
  if((fp2=fopen(fileName,"wb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp2);
  	exit(-1);
  } 
  fwrite(xs,sizeof(unsigned long),8*N,fp2); 
  
  //free
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  
  free(xs);
  return 0;
}



