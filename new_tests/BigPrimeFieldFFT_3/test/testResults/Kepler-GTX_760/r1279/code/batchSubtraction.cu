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

__global__ void batchBigSub(unsigned long r, unsigned long *xs, unsigned long *ys, unsigned long *us) 
{  
  int tid = threadIdx.x;  
  short i, pos;
  unsigned short c=0;
  unsigned long num1,num2;
  unsigned long *xd, *yd, *ud; 
  num1=0;
  num2=r-1;
  
  xd = (unsigned long*)((char*)xs + tid*sizeof(unsigned long)*8);
  yd = (unsigned long*)((char*)ys + tid*sizeof(unsigned long)*8);
  ud = (unsigned long*)((char*)us + tid*sizeof(unsigned long)*8);

  
  for(i=0;i<=7;i++)
  {
    num1=yd[i]+c;
    if(xd[i]<num1) //there is not enough to do subtraction
    {
      c=1;
      ud[i]=r-num1+xd[i];
    }
    else
    {
      c=0;
      ud[i]=xd[i]-num1;
    }
  }
  if(c>0)
  {
    pos=-1;
    for(i=0;i<8;i++) 
    {
      if(ud[i] < num2)
      {
        pos=i;
        break;
      }
    }
    if(pos>=0)
    {
      for(i=0;i<pos;i++)
      {
          ud[i] = 0;
      }
      ud[pos] ++;
    }            
    else
    {
      ud[0]=18446744073709551615u;
      for(i=1;i<8;i++)
      {
        ud[i]=0;
      }
    }  
  }
  
}

int main(int argc, char *argv[])
{
  int i,j;
  int L=16;  //L should be smaller than 1024.
  int pos=0;
  unsigned long r = 9223372054034644992;
  unsigned long *xs, *ys, *us, *xs_d, *ys_d, *us_d; 
  xs=(unsigned long *)malloc((sizeof(unsigned long)*8 * L));
  ys=(unsigned long *)malloc((sizeof(unsigned long)*8 * L));
  us=(unsigned long *)malloc((sizeof(unsigned long)*8 * L));
  //init value
  for(j=0;j<L;j++)
  {
    for(i=0;i<8;i++)
    {
      pos=j*8+i;  // pos++;
      xs[pos] = 1;
      ys[pos] = 9223372054034644990u;
      us[pos] = 0;
    }
  }
  
  //init gpu memory
  cudaMalloc((void **)&xs_d, sizeof(unsigned long) * 8 * L);
  cudaMemcpy(xs_d, xs, sizeof(unsigned long) * 8 * L, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&ys_d, sizeof(unsigned long) * 8 * L);
  cudaMemcpy(ys_d, ys, sizeof(unsigned long) * 8 * L, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&us_d, sizeof(unsigned long) * 8 * L);
  cudaMemcpy(us_d, us, sizeof(unsigned long) * 8 * L, cudaMemcpyHostToDevice);
  
  batchBigSub<<<1,L>>>(r,xs_d,ys_d,us_d);
  cudaThreadSynchronize();
  cudaMemcpy(us, us_d, sizeof(unsigned long) * 8 * L, cudaMemcpyDeviceToHost);
  
  printf("we have done us = xs - ys.\n");
  printf("the xs is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", xs[i]);
  }
  printf("the ys is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", ys[i]);
  }
  printf("the us is \n");
  for(i=0;i<16;i++)
  {
    printf ("%lu\n", us[i]);
  }
  
}

__global__ void addition1(unsigned long *a) 
{
  *a = *a + 1;
}

int main1(int argc, char *argv[])
{	
  unsigned long *a;
  unsigned long *a_d;
  a=(unsigned long *)malloc((sizeof(unsigned long)*1));  
  a[0] = 18446744073709551614u;
  printf("a is %lu\n",a[0]);
  cudaMalloc((void **)&a_d, sizeof(unsigned long) * 1);
  cudaMemcpy(a_d, a, sizeof(unsigned long) * 1, cudaMemcpyHostToDevice);
  addition1<<<1,1>>>(a_d);
  cudaThreadSynchronize();
  cudaMemcpy(a, a_d, sizeof(unsigned long) * 1, cudaMemcpyDeviceToHost);
  printf("a is %lu\n",a[0]);
  return 0;
}

