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
#define TN 256 
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
__constant__ unsigned long W[15][8]={4832614680308292603u, 8402438242641148880u, 2486701145004827901u, 1215222897245048005u, 8892305690461566140u, 7297077274619179662u, 5862425953453277502u, 3703360755969114825u, 2302450199823606216u, 27080194135346259u, 5112443829649209337u, 8827039575836218556u, 1287512346408126576u, 4365726449074111407u, 3360903764649200564u, 117722958365041421u, 390324908435919352u, 770891875385868514u, 2999696378145516236u, 3132067064637725081u, 3345465055005335275u, 5514352822140455844u, 6208596698165918434u, 3765508936539829447u, 9210848751572284741u, 1645667176001997922u, 5941128666689277807u, 6768772028034979755u, 4304177272275830547u, 6104712663442007465u, 1954368461670077480u, 6432833914249218714u, 5366875984871234376u, 8713960611431090699u, 7737743610771967836u, 4469224371143296702u, 8872582861821233821u, 3739059486843436842u, 8813835214102132195u, 5147471904605268606u, 2493293057650284961u, 8167589484735262926u, 6404888664753609044u, 2974179143108738507u, 4816710638611137048u, 1700523433557964716u, 5243180851737162921u, 7115155689365786339u, 3506894965669318197u, 7418365113978784395u, 8014352668811394697u, 2170069100293384465u, 4477472493274660183u, 8475805979162198949u, 3965693063579606190u, 7514110372369303684u, 6476678624121558179u, 3216147715999829916u, 4168578045703223305u, 5665582343424371252u, 4802751892420194639u, 7507476765166618649u, 794581407066667206u, 4885142567208165039u, 5149649432978519211u, 5596777036269387984u, 8665828768845738479u, 1575627868412583125u, 6499525351416721532u, 5094347230997088364u, 8715176080936469701u, 1841878500194672339u, 7542379292388838454u, 9056921649810209609u, 8630708080678083361u, 6562740377707965803u, 802590099908464585u, 7056946099352360173u, 7815206358853348244u, 8556126457443025758u, 2735242483360463802u, 1509217343334731913u, 8137216084652364444u, 5194826376948744458u, 8888111600302919097u, 5096574789189059728u, 2602858914188694775u, 1490780946449946389u, 1871173134326118337u, 6362464371896225965u, 8698229383614300186u, 2005452311774545610u, 7556559507024761731u, 5864043929342978154u, 7933436335374648226u, 5827741928560708899u, 1224342742045833269u, 972714779852171763u, 2469969374900767744u, 925611663966851379u, 5253808685097889801u, 7052192790816421625u, 6899658349937940167u, 1281442313926189540u, 5206235348781813121u, 976526785949102141u, 4199848929055484124u, 663415973782467310u, 4667091686862902153u, 5656282773537629911u, 1019143404946265847u, 4695224974402318448u, 6954395578847384962u, 5241574410285476070u, 1626120174922578762u, 2533521736498923254u, 3068873176292580125u, 753173231159502005u, 4879546884467215097u, 6201146281074609078u};


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

__device__ void mulOmegaPower(unsigned long *xs, short po)
{
  short rp,wp;
  if(po<=0)
  {
    return;
  }
  rp=po>>4; //w^16=r
  wp=po-(rp<<4);
  
  rp>8?cyclicShift(xs,8):cyclicShift(xs,rp);
  rp>8?cyclicShift(xs,rp-8):cyclicShift(xs,0);  
  
  if(wp>0)
  {
    bigMul(xs,W[wp-1]);
  }
}
__device__ void transpose(unsigned long *idata, unsigned long *odata)
{
	unsigned long *ip,*op;
	short i,pos,tid,w,h,pos2;
	
	tid = threadIdx.x;	
	h=tid>>4; //h=tid/16
	w=tid-(h<<4); //w=h mod 16
	pos=((w<<4)+h)<<3;
  pos2=tid<<3;
  
  ip = (unsigned long*)((char*)idata + pos*sizeof(unsigned long)); //iregular
	op = (unsigned long*)((char*)odata + pos2*sizeof(unsigned long)); //regular
	
  for(i=0;i<8;i++)
	{
	  op[i]=ip[i];
	}	
	__syncthreads();
}

__device__ void transpose2(unsigned long *idata, unsigned long *odata)
{
	unsigned long *ip,*op;
	short i,pos,tid,w,h,pos2;
	
	tid = threadIdx.x;	
	pos=tid<<3;
	h=tid>>4; //h=tid/16
	w=tid-(h<<4); //w=h mod 16
	pos2=((w<<4)+h)<<3; 
	
	ip = (unsigned long*)((char*)idata + pos*sizeof(unsigned long)); //iregular
	mulOmegaPower(ip,h*w);	
	__syncthreads();
  op = (unsigned long*)((char*)odata + pos2*sizeof(unsigned long)); //regular
  for(i=0;i<8;i++)
	{
	  op[i]=ip[i];
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
  
  xf = (unsigned long*)((char*)xsm + (bwid*64+wpn*16)*sizeof(unsigned long)*8); //
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
  posid=(wpid>=4?wpid+4:wpid);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+4)*8],ind2[wpid]);
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
  posid=(wpid>=4?(wpid-4)*4+1:wpid*4);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+2)*8],ind3[wpid]);
  }
  
  __syncthreads(); 
  if(sf>0)
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
  posid=(wpid<<1);
  if(sf>0)
  {
    cyclicShift(&xf[(posid+1)*8],ind4[wpid]);
  }
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
    pos1=ind5[posid]*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  else
  {
    posid=wpid;
    pos=posid*8;
    pos1=ind5[posid]*8;
    for(i=0;i<8;i++)
    {
      xf[pos+i]=yf[pos1+i];
    }
  }
  __syncthreads(); 
  
}

__global__ void FFT256(unsigned long *xs) 
{
  //1 block 256 threads
  short i,h,w,pos,pos2;
  short tid = threadIdx.x;
  unsigned long *xd;
  __shared__ unsigned long xsm[INC*8]; 
  __shared__ unsigned long ysm[INC*8];  
  
  //0. import data
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  pos=tid<<3; //tid*8
  for(i=0;i<8;i++)
  {
    xsm[pos++]=xd[i];
  }
  __syncthreads();  //copy the input to shared memory.
  
  //1. transposition
  transpose(xsm,ysm);  //transposition xsm into ysm
  __syncthreads();  
  
  
  
  //2. fft16
  FFT16(ysm,xsm);   //calculate FFT16 of ysm. The final result is also stored in ysm. xsm is a temp array.
  __syncthreads();
  
  //3. multiple w^{jk} and transposition
  transpose2(ysm,xsm); //calculate ysm*w^{jk} and store the result into xsm
  __syncthreads();
  
  //4. fft16 again
  FFT16(xsm,ysm); //calculate FFT16 of xsm. The final result is also stored in xsm. ysm is a temp array.
  __syncthreads();
  
  //5. write back, from xsm to xs
  h=tid>>4; //h=tid/16
	w=tid-(h<<4); //w=h mod 16
	pos2=((w<<4)+h)<<3;
	pos=tid<<3;
  for(i=0;i<8;i++)
  {
    xs[pos2++]=xsm[pos++];
  }
  
}

int main(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1, *fp2;	 
  unsigned long *xs, *xs_d; 
  cudaEvent_t start, stop;
  float elapsedTime;
  short i;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*INC));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_256_input.dat"); //INC*8 (unsigned long) data
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  memset(xs,0,sizeof(unsigned long)*8*INC);
  fread(xs,sizeof(unsigned long),8*INC,fp1);
  
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*INC);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*INC, cudaMemcpyHostToDevice);
  
  FFT256<<<BN,TN>>>(xs_d);
  cudaThreadSynchronize();
   
  cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*INC, cudaMemcpyDeviceToHost);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done xs = FFT256(xs) in 1 times.\n");
  printf("the time of gpu is %f ms\n", elapsedTime);
  printf("we now write the output result.\n");
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_256_output.dat"); //INC*8 (unsigned long) data
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
  for(i=0;i<8*INC;i++)
  {
    printf("%lu\n",xs[i]);
  }
  cudaFree(xs_d);	
  free(xs);
  return 0;
}




/*
int main1(int argc, char *argv[])
{
  char fileName[1024];
	FILE *fp1, *fp2;	 
  unsigned long *xs, *xs_d, *ys; 
  cudaEvent_t start, stop;
  float elapsedTime;
  short i;
  
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8*INC));
  ys=(unsigned long *)malloc((sizeof(unsigned long)*8*INC));
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_256_input.dat"); //INC*8 (unsigned long) data
  if((fp1=fopen(fileName,"rb"))==NULL)
  {
  	printf("fail to %s", fileName);
  	fclose(fp1);
  	exit(-1);
  } 
  
  memset(xs,0,sizeof(unsigned long)*8*INC);
  fread(xs,sizeof(unsigned long),8*INC,fp1);
  memcpy(ys,xs,sizeof(unsigned long)*8*INC);
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  cudaMalloc((void **)&xs_d, sizeof(unsigned long)*8*INC);
  for(i=0;i<100;i++)
  {
    cudaMemcpy(xs_d, xs, sizeof(unsigned long)*8*INC, cudaMemcpyHostToDevice);
  
    FFT256<<<BN,TN>>>(xs_d);
    cudaThreadSynchronize();
   
    cudaMemcpy(xs, xs_d, sizeof(unsigned long)*8*INC, cudaMemcpyDeviceToHost);
    memcpy(xs,ys,sizeof(unsigned long)*8*INC);
  }  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  
  printf("we have done xs = FFT256(xs) in 100 times.\n");
  printf("the time of gpu is %f ms\n", elapsedTime);
  printf("we now write the output result.\n");
  
  memset(fileName, 0 , sizeof(char)*1024);
  sprintf(fileName, "FFT_256_output.dat"); //INC*8 (unsigned long) data
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
  for(i=0;i<8*INC;i++)
  {
    printf("%lu\n",xs[i]);
  }
  cudaFree(xs_d);	
  free(xs);
  free(ys);
  return 0;
}

*/

