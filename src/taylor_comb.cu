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
#include "ts_VCA_count.h"
#include "taylor_comb.h"


/* ------------------ Taylor_shift procedure ------------------ */

/** This procedure calls s times the procedure Taylor_shift_Zp.
    It gives the CRT representation of the polynomial
    2^{k*n}*P( (x+a)/2^k ) (called PKC), combination of a
    Taylor_shift by a of a input polynomial P and multiplications for
    each coefficient by a power of 2^k. One column of X
    will be filled modulo p (p prime numbers).

    PARAMETERS :
    ------------
*   @a : integer for the Taylor_shift by a
*   @k : power of 2 considered
*   @n : number of coefficients of the polynomial considered
         (size of the polynomial)
*   @e : e = log2(n) so n = 2^e
*   @f : f = 28, parameter used to compute s
*   @s : number of prime numbers needed
*   @X : array of size n*s containing CRT and MRR representation of
         the output polynomial
*   @Polynomial_device : input polynomial
*   @Factorial_device : array containing the sequence (i!)_{i in [0;n]}
*   @Powers_a : array containing in this code :
                -> firstly : the sequence (a^i)_{i in [0;n]} mod p
                -> secondly : the sequence ((2^k)^i)_{i in [0;n]} mod p
*   @Monomial_shift_device : array containing the sequence of the coeffi-
                             cients of ((x+1)^{2^i})_{i in [0;e-1]} mod p
*   @Mgpu : array used in several intermediate computations mod p
*   @Polynomial_shift_device : array in 2 dimensions containing pairwise
                               intermediate computations to obtain the
                               Taylor_shift by 1 of a polynomial mod p
*   @fft_device : array used for Fast Fourier Transform mod p
**/
void taylor_shift_Z(sfixn a, sfixn k, sfixn n, sfixn e, sfixn f, sfixn s, sfixn *X, sfixn *Polynomial_device, sfixn *Factorial_device, sfixn *Powers_a, sfixn *Monomial_shift_device, sfixn *Mgpu, sfixn **Polynomial_shift_device, sfixn *fft_device)
{
  // declaration of variables
  sfixn p, i;
  double pinv;

  // do the PKC s times in Z/pZ then store the result in X
  for (i=0; i<s; i++)
  {
    p = primes[i];
    pinv = (double) 1 / ((double) p);
    taylor_shift_Zp(n, e, p, pinv, X, i, s, Polynomial_device, a, k, Factorial_device, Powers_a, Monomial_shift_device, Mgpu, Polynomial_shift_device, fft_device);
  }
}


/** Each thread in each thread block computes the maximum in absolute 
absolute value of the positive (resp. negative) numbers in 
a sequence T_CHN (= 32) consecutive coefficients in MRR representation.
   That's why we need 2 numbers as 
    mrrPosNorm and mrrNegNorm. Later we will compute the maximum of the sequence of all
    these numbers. Here each thread computes some maximums according to the sign of the
    coefficients in MRR representation.
* @n : n = d+1, n is the size of the polynomial.p_{k,c}
* @X : array of size n*s storing the MRR representation of
       the polynomial p_{k,c}
* @mrrPosNorm :  largest of positive coefficients of p_{k, c} in MRR
representation. (output)
* @mrrNegNorm : smallest of negative coefficients of p_{k, c} in MRR
representation. (output)
* @s : number of prime numbers needed
* @Tsign : array containing the sign of the coefficients of the polynomial p_{k,c}
* @primes_GPU, input array of prime numbers
// !!! : I'm not sure about the sign of (M-1)/2, positive
**/
__global__ void value_max_gpu(sfixn n, sfixn *X, sfixn* mrrPosNorm, sfixn *mrrNegNorm, sfixn s, sfixn *Tsign, sfixn *primes_GPU)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn j, k, l, t;

  if (threadIdx.x == 0)
    for (j=0; j<NB_THREADS; j++)
    {
      mrrPosNorm[j] = 0;
      mrrNegNorm[j] = primes_GPU[j]-1;
    }
  __syncthreads();


  for (j=0; j<T_CHN; j++)
  {
    l = i*T_CHN+j;

    if (Tsign[l] == 0)
      continue;

    else if (Tsign[l] == 1)
    {
      k = s-1;
      t = 0;

      while ((t == 0) && (k>=0))
      {
        t = X[l*T_CHN+k] - mrrPosNorm[k];
        k--;
      }

      if (k == -1)
        continue;

      if (t>0)
        for (j=0; j<NB_THREADS; j++)
          mrrPosNorm[j] = X[l*T_CHN+j];
    }

    else // (Tsign[l] == -1)
    {
      k = s-1;
      t = 0;

      while ((t == 0) && (k>=0))
      {
        t = X[l*T_CHN+k] - mrrNegNorm[k];
        k--;
      }

      if (k == -1)
        continue;

      if (t<0)
        for (j=0; j<NB_THREADS; j++)
          mrrNegNorm[j] = X[l*T_CHN+j];
    }

  }
}


/** creates the array Tsign given the sign of the number represented in MRR representation in each row of X
* @Tsign : array containing the sign of the coefficients of the input polynomial considered of size n
* @X : array of size n*s containing CRT and MRR representation of
       the output polynomial
* @n : n = d+1, n is the size of the polynomial.
* @s : number of prime numbers needed
* @primes_GPU, input array of prime numbers
**/
__global__ void create_signs_GPU(sfixn *Tsign, sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ sfixn PRIMES[NB_THREADS];

  if (i<n)
  {
    PRIMES[threadIdx.x] = primes_GPU[threadIdx.x];
    __syncthreads();

///////////////////////////
// replace the MRR_sign function, problem with the call
    sfixn j = s-1;
    sfixn sign = 0;
    sfixn k = 0;
    sfixn comp;

    while ((X[i*s+k] == 0) && (k<s))
      k++;

    if (k!=s)
    {
      while ((sign == 0) && (j>=0))
      {
        comp = X[i*s+j] - (PRIMES[j]-1)/2;
        j--;
        if (comp > 0)
          sign = -1;
        else if (comp < 0)
          sign = 1;
      }

      if (j == -1)
      sign = 1;
    }
/////////////

    Tsign[i] = sign;
  }
}


/** Input :
    -------
  (k, c): current node. (c corresponds also to a in the code)
* @d : degree of polynomial.
* @n : n = d+1, n is the size of the polynomial.
* @np : number of extra prime numbers.
* @A : array of images of the input polynomial modulo prime numbers whose size
is np*(d+1), that is, np*n.
* @checkZero: if 1 then check if trailing coefficient of p_{k, c} is 0; 
otherwise do nothing.

    Output :
    --------
* @sc : number of sign changes.
* @isZero : 1 if we checked that the trailing coefficient of p_{k, c} was 0, otherwise 0.
* @mrrPosNorm :  MRR representation of the largest positive coefficient of p_{k, c} 
* @mrrNegNorm : MRR representation of the smallest negative coefficients of p_{k, c} 
**/
// REMAIN : parameters of Changbo's procedure
// extra parameters needed:
//   * primes_GPU
//   * count_sg
//   * bound_sg
//   * ADD
//   * Tsign
//   * s
void desBound_gpu(sfixn k, sfixn c, sfixn n, sfixn np, sfixn *A, sfixn checkZero, sfixn *sc, sfixn *IsZero, sfixn* mrrPosNorm, sfixn *mrrNegNorm, sfixn *primes_GPU, sfixn s)
{
  sfixn sg, nb_blocks;
  sfixn m = (sfixn) ceil(n/(double) T_CHN);
  sfixn Vv[1];

  // mmrPosNorm and mrrNegNorm must be arrays of size s including only zeros at the beginning of the code
  //  sfixn *mmrPosNorm, *mmrNegNorm;
  //  cudaMalloc( (void **) &mmrPosNorm , s * sizeof(sfixn) );
  //  cudaMalloc( (void **) &mmrNegNorm , s * sizeof(sfixn) );

  sfixn *count_sg, *bound_sg, *ADD, *b;
  cudaMalloc( (void **) &count_sg, m * sizeof(sfixn) );
  cudaMalloc( (void **) &bound_sg, 2*m * sizeof(sfixn) );
  cudaMalloc( (void **) &ADD, 1 * sizeof(sfixn) );
  b = (sfixn*) malloc(s*sizeof(sfixn));
  cudaMemcpy(b, A, s*sizeof(sfixn), cudaMemcpyDeviceToHost );
  if (checkZero == 1)
  {
    sg = MRR_sign_CPU(b, s);
    *IsZero = (sfixn) (sg == 0);
  }

  sfixn *Tsign;
  cudaMalloc( (void**) &Tsign, n*sizeof(sfixn) );

  nb_blocks = number_of_blocks(n);
  create_signs_GPU<<<nb_blocks,NB_THREADS>>>(Tsign, A, n, s, primes_GPU);

  nb_blocks = (sfixn) ceil(m/(double) T_CHN);
  sign_change_MRR2<<<nb_blocks,T_CHN>>>(A, n, s, primes_GPU, count_sg, bound_sg, Tsign);
  global_add_gpu<<<1,NB_THREADS>>>(count_sg, m, ADD);
  cudaMemcpy(Vv, ADD, 1*sizeof(sfixn), cudaMemcpyDeviceToHost);
  sfixn *bound_sg_cpu;
  bound_sg_cpu = (sfixn *) malloc(2*m * sizeof(sfixn));
  cudaMemcpy(bound_sg_cpu, bound_sg, 2*m*sizeof(sfixn), cudaMemcpyDeviceToHost);

  *sc = Vv[0] + nb_bounds_change_sign(bound_sg_cpu, 2*m);



  nb_blocks = (sfixn) ceil(m/(double) T_CHN);
//  value_max_gpu<<<nb_blocks,T_CHN>>>(X, n, s, primes_GPU, count_sg, bound_sg);



  cudaFree(count_sg);
  cudaFree(bound_sg);
  cudaFree(ADD);
  free(b);
  cudaFree(Tsign);
  cudaFree(bound_sg_cpu);

}
