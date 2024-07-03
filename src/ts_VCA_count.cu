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


#include "taylor_shift_conf.h"
#include "inlines.h"
#include "taylor_shift_cpu.h"
#include "taylor_shift_kernel.h"
#include "taylor_shift.h"
#include "taylor_shift_fft.h"
#include "list_pointwise_mul.h"
#include "list_stockham.h"
#include "ts_chinese_remainder.h"
#include "primesFFT.h"
#include "taylor_comb.h"
#include "ts_VCA_count.h"


// compute the sign of a polynomial in MRR-representation in the GPU
// b is a vector size s which encodes an integer z in MRR representation
// using the s primes of primes_GPU.
// The output is the sign fo z
__device__ sfixn MRR_sign(sfixn *b, sfixn s, sfixn *primes_GPU)
{
  sfixn j = s-1;
  sfixn sign = 0;
  sfixn k = 0;
  sfixn comp;

  // test to check if b=0
  while ((b[k] == 0) && (k<s))
    k++;

  if (k!=s)   // then b!=0 and we can compare coefficients of n with h_j=(primes[j]-1)/2
  {
    while ((sign == 0) && (j>=0))
    {
      comp = b[j] - (primes_GPU[j]-1)/2;
      j--;
      if (comp > 0)
        sign = -1;         // b represents a negative number
      else if (comp < 0)
        sign = 1;          // b represents a positive number
    }

    if (j == -1)           // b represents (M-1)/2
      sign = 1;
  }

  // if (k==s) then b==0 represents the real number 0 so sign = 0

  return sign;
}


// compute the sign of a polynomial in MRR-representation in the CPU
// b is a vector size s which encodes an integer z in MRR representation
// using the first s primes of the array "primes" defined in primesFFT.h
// The output is the sign fo z
sfixn MRR_sign_CPU(sfixn *b, sfixn s)
{
  sfixn j = s-1;
  sfixn sign = 0;
  sfixn k = 0;
  sfixn comp;

  // test to check if b=0
  while ((b[k] == 0) && (k<s))
    k++;

  if (k!=s)   // then b!=0 and we can compare coefficients of n with h_j=(primes[j]-1)/2
  {
    while ((sign == 0) && (j>=0))
    {
      comp = b[j] - (primes[j]-1)/2;
      j--;
      if (comp > 0)
        sign = -1;         // b represents a negative number
      else if (comp < 0)
        sign = 1;          // b represents a positive number
    }

    if (j == -1)           // b represents (M-1)/2
      sign = 1;
  }

  // if (k==s) then b==0 represents the real number 0 so sign = 0

  return sign;
}



// compute locally a number of sign change in X and store data in the arrays count_sg and bound_sg

/* count_sg : each thread computes a local number of sign change
              and stores this number in one position of count_sg
   bound_sg : each thread stores the first and the last non-null
              sign : it's possible we don't count all the change
              as each thread doesn't know what was the first sign
              after the part it deals with and also for the first
              sign before this same part.                         */

// X is a n-times-s array in row major format. As above, s is the number
// of primes in use and n is the size of the polynomial encoded by X
// This input polynomial i sin MRR representation.
// Each thread computes the number of sign changes for a sequence
// of T_CHN (= 32) consecutive coefficients (= rows of X). 
// Moreover, each thread stores the first and the last non-null sign.
// in the array bound_sg.
// By comparing the values stored by bound_sg by all threads, one
// can deduce the "possibly missing sign changes" that were not
// accumated in count_sg

__global__ void sign_change_MRR(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg)
{
  sfixn I = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn Z[T_CHN];
  sfixn i, j, begin, end, base;
  sfixn BOUND_SG_left, BOUND_SG_right;
  sfixn nb_chg_sign = 0;

  // store in Z the signs of the numbers in the rows i = I*T..(I+1)*T-1
  for (j=0; j<T_CHN;j++)
  {
    i = I*T_CHN + j;
    Z[j] = 0;
    if (i<n)
      Z[j] = MRR_sign(X+i*s, s, primes_GPU);
  }

  // parameters to define the length of the loop 
  begin = 0;
  end = T_CHN;

  // find the position of the first sign non-null
  base = 0;
  while ((base==0) && (begin<T_CHN))
  {
    base = Z[begin];
    begin++;
  }

  // adjust the sign bounds and the beginning of the next loop
  BOUND_SG_left  = base; // BOUND_SG_left  known
  BOUND_SG_right = base; // BOUND_SG_right extremely likely modifiable
  // remark : if Z = {0,0,...,0} then bound_sg[2i] = bound_sg[2i+1] = 0

  // conditionnal loop
  for (i=begin; i<end; i++)
  {
    if (Z[i] != 0) // if (Z[i] == 0), nothing has to be done, we can continue the loop
    {
      BOUND_SG_right = Z[i]; // update BOUND_SG_right
      if (Z[i]*base == -1)   // test if there is a sign change
      {
        nb_chg_sign++;       // increment nb_chg_sign
        base = Z[i];         // base become the opposite than previously
      }
    }
  }

  // output of the thread I
  count_sg[I]     =  nb_chg_sign;
  bound_sg[2*I]   =  BOUND_SG_left;
  bound_sg[2*I+1] =  BOUND_SG_right;

}


// same procedure than the previous one but avoid to compute again the signs with an additional parameter Tsign.
// The array Tsign contains the signs of the coefficients of X
// Tsigns resides in the gloabla memory
__global__ void sign_change_MRR2(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg, sfixn *Tsign)
{
  sfixn I = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn Z[T_CHN];
  sfixn i, j, begin, end, base;
  sfixn BOUND_SG_left, BOUND_SG_right;
  sfixn nb_chg_sign = 0;

  // store in Z the signs of the numbers in the rows i = I*T..(I+1)*T-1
  for (j=0; j<T_CHN;j++)
  {
    i = I*T_CHN + j;
    Z[j] = 0;
    if (i<n)
      Z[j] = Tsign[i]; // !!! ONLY THIS LINE IS MODIFIED !!!
  }

  // parameters to define the length of the loop 
  begin = 0;
  end = T_CHN;

  // find the position of the first sign non-null
  base = 0;
  while ((base==0) && (begin<T_CHN))
  {
    base = Z[begin];
    begin++;
  }

  // adjust the sign bounds and the beginning of the next loop
  BOUND_SG_left  = base; // BOUND_SG_left  known
  BOUND_SG_right = base; // BOUND_SG_right extremely likely modifiable
  // remark : if Z = {0,0,...,0} then bound_sg[2i] = bound_sg[2i+1] = 0

  // conditionnal loop
  for (i=begin; i<end; i++)
  {
    if (Z[i] != 0) // if (Z[i] == 0), nothing has to be done, we can continue the loop
    {
      BOUND_SG_right = Z[i]; // update BOUND_SG_right
      if (Z[i]*base == -1)   // test if there is a sign change
      {
        nb_chg_sign++;       // increment nb_chg_sign
        base = Z[i];         // base become the opposite than previously
      }
    }
  }

  // output of the thread I
  count_sg[I]     =  nb_chg_sign;
  bound_sg[2*I]   =  BOUND_SG_left;
  bound_sg[2*I+1] =  BOUND_SG_right;
}


// computes on the CPU the number of sign changes missed by the GPU
sfixn nb_bounds_change_sign(sfixn *bound_sg, sfixn size)
{
  sfixn i, j;
  sfixn nb_chg_sign = 0;
  
  for (i=1; i<size-2; )
  {
    if (bound_sg[i] == 0)
    {
      i += 2;
      continue;
    }

    j = i+1;
    while ((bound_sg[j] == 0) && (j<size))
      j += 2;

    if(j>=size)
      break;
    else // (j<size)
    {
      if (bound_sg[i]*bound_sg[j] == -1)
        nb_chg_sign++;
      i = j+1;
      continue;
    }

  }

  return nb_chg_sign;
}


// computes the maximum in absolute value of an array and store it in MAXI[0]
// T is an array of size n, whose maximum is stored in MAXI[0]
__global__ void max_gpu(sfixn *T, sfixn n, sfixn *MAXI)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn m = (sfixn) ceil( n/(double) NB_THREADS);
  __shared__ sfixn Z[NB_THREADS];
  sfixn j;
  sfixn max = 0;
  Z[i] = 0;

  if (i < m)
  {
    for (j=i*NB_THREADS; (j<(i+1)*NB_THREADS) && (j<n); j++)
    {
      if (abs(T[j]) > max)
        max = T[j];
    }
    __syncthreads();
    Z[i] = max;
  }

  __syncthreads();
  if (i == 0)
  {
    max = 0;
    for (j=0; j<NB_THREADS; j++)
      max += Z[j];
    MAXI[0] = max;
  }
}


// computes the sum of all the non-negative entries of an 
// array and store it in ADD[0]
// IMPORTABT: it is assumed that the sum of these entries 
// is no more that the size of the array.
__global__ void global_add_gpu(sfixn *T, sfixn n, sfixn *ADD)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn m = (sfixn) ceil( n/(double) NB_THREADS);
  __shared__ sfixn Z[NB_THREADS];
  sfixn j;
  sfixn add = 0;
  Z[i] = 0;

  if (i < m)
  {
    for (j=i*NB_THREADS; (j<(i+1)*NB_THREADS) && (j<n); j++)
      add += T[j];
    __syncthreads();
    Z[i] = add;
  }

  __syncthreads();
  if (i == 0)
  {
    add= 0;
    for (j=0; j<NB_THREADS; j++)
      add += Z[j];
    ADD[0] = add;
  }
}


// transform X in MRR representation
/* input  : X in CRT representation
   output : X in MRR representation */
// X as above (thus an n-times-s 2D array in row major layout)
// e is log[2](n)
void MRR_all(sfixn *X, sfixn n, sfixn e, sfixn s, sfixn *primes_GPU)
{
  /* use the Taylor mixed representation by a matrix of modulo to
     obtain coefficients b_i we will use to get the final Taylor shift
     we want                                                           */
  char name_file[100];
  sfixn nb_blocks;
  sfixn *D, *L;
  cudaMalloc( (void **) &D , s * (s-1)/2 * sizeof(sfixn) );
  cudaMalloc( (void **) &L , s * (s-1)/2 * sizeof(sfixn) );
  createDL<<<s*(s-1)/2, NB_THREADS>>>(D, L, primes_GPU, s);

  // save X before doing anything
  /* FOR TESTS */
/*  sprintf(name_file, "Pol%d.X.dat\0", e);
  sfixn *Y;
  Y = (sfixn*) malloc(n*s*sizeof(sfixn));
  cudaMemcpy(Y, X, n*s*sizeof(sfixn), cudaMemcpyDeviceToHost);
  stock_array_in_file(name_file, Y, n*s);
  free(Y); */


  // do operations in the GPU
  sprintf(name_file, "Pol%d.Y.dat", e);
  mixed_radix_phase1<<<n, T_CHN>>>(X, s, n, primes_GPU, D, L);
  nb_blocks = (sfixn) ceil((double) n/(double)T_CHN);
  mixed_radix_phase2<<<nb_blocks, T_CHN>>>(X, s, n, primes_GPU, D, L);


  cudaFree(D);
  cudaFree(L);

/* FOR TESTS */
/*  sfixn *temp2;
  temp2 = (sfixn*) malloc(n*s * sizeof(sfixn));
  cudaMemcpy(temp2, X, s*n*sizeof(sfixn), cudaMemcpyDeviceToHost);
  stock_array_in_file(name_file, temp2, n*s);
  free(temp2); */
}


// function counting the number of sign change in the MRR representation (X) of the polynomial we want

/**/
sfixn count_all(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg, sfixn *ADD, sfixn m)
{
  sfixn Vv[1], V;
  sfixn nb_blocks = (sfixn) ceil(m/(double) T_CHN);

  sign_change_MRR<<<nb_blocks,T_CHN>>>(X, n, s, primes_GPU, count_sg, bound_sg);
  global_add_gpu<<<1,NB_THREADS>>>(count_sg, m, ADD);
  cudaMemcpy(Vv, ADD, 1*sizeof(sfixn), cudaMemcpyDeviceToHost);

  sfixn *bound_sg_cpu;
  bound_sg_cpu = (sfixn *) malloc(2*m * sizeof(sfixn));
  cudaMemcpy(bound_sg_cpu, bound_sg, 2*m*sizeof(sfixn), cudaMemcpyDeviceToHost);

  return Vv[0] + nb_bounds_change_sign(bound_sg_cpu, 2*m);
}

