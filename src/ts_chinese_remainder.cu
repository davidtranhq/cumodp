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
#include "primesFFT.h"

// gives 'positions' in arrays D and L
// !!! not used any more in the code !!!
/*__device__ void funcPOS(sfixn i, sfixn s, sfixn *out)
{
  out[0] = 0;
  out[1] = 0;

  while (out[0]<=i)
  {
    out[1]++;
    out[0] += s - out[1];
  }

//  if (out[0] != i)
//  {
    out[0] -= s - out[1];
    out[1]--;
//  }
}*/


/* CONVERSION OF MODULAR NUMBERS WITH MIXED RADIX REPRESENTATION BY A MATRIX FORMULA
 
   Given a vector x of size s representing an integer number in CRT rep, 
   the following computes ist MRR rep 
         b = (...((x A_1)A_2)...)A_{s-1} using A_k such that :

           [               |                                            ]
           [               |                                            ]
           [     I_{k-1}   |                   0                        ]
           [               |                                            ]
           [    ___________|_________________________________________   ]
           [               |                                            ]
           [               |  1   n_{k,k+1}  n_{k,k+2}  ...  n_{k,s}    ]
   A_k =   [               |  0   m_{k,k+1}      0      ...     0       ]
           [               |  0       0      m_{k,k+2}          :       ]
           [               |  :                  0   .          :       ]
           [       0       |  :                    .   .        :       ]
           [               |  :        ...           .   .      :       ]
           [               |  :                        .   .    :       ]
           [               |  0          ...             0   m_{k,s}    ]
           [               |                                            ]


    We don't really need to do matrix multiplication, but just to consider diagonal elements
    and the k-th lines of these matrices. We need to use vectors D and L.
*/


/** procedure to compute the elements of vectors D and L
* @D, output array containing diagonal elements (m_{i,j}) of all the matrices
      A_k such as the example given as follow
* @L, output array containing k_th rows of all the matrices A_k (n_{i,j}) such
      as the example given as follow
* @primes_GPU, input array of prime numbers
* @s, number of such prime numbers
**/
__global__ void createDL(sfixn *D, sfixn *L, sfixn *primes_GPU, sfixn s)
{
  sfixn I = blockIdx.x * blockDim.x + threadIdx.x;

/* ------------------------------------------------------------------------------------

                              STRUCTURE OF THE ARRAYS L & D

Example for L :
---------------
                                                       Q
             PART 0                      PART 1        |                 PART s-1
    _____________________________ _____________________V_______________ ___________
   |                             |                     |               |           |
   |   n_1,2  n_1,3  ...  n_1,s  |  n_2,3  ...  n_2,s  |  ...........  |  n_s-1,s  |
   |_____________________________|_____________________|_______________|___________|
(we want)      L[1]                        L[2]                             L[s]
(we store)   ~L'[0]                      ~L'[1]                          ~L'[s-1]

   ------------------------------------------------------------------------------------

  m_i * m_i,j  =  1  mod m_j    with 0 < m_i,j < m_j
        n_i,j  = m_j - m_i,j

  L_k = {n_k,j | k+1 <= j <= s}
  D_k = {m_k,j | k+1 <= j <= s}

  Q gives the "absolute" positions at the beginning of each PART
  alpha = i
  alpha is the index of the current part plus one.
*/


  if (I < s*(s-1)/2)
  {
    sfixn Q = 0;
    sfixn alpha = 1;
    sfixn Q_next = s-1;

    // compute Q & alpha for a thread I so as to identify which m_{alpha,j} or n_{alpha,j} the thread I has to compute
    while (I >= Q_next)
    {
      Q = Q_next;
      alpha++;
      Q_next += (s-alpha);
    }
    alpha--;

    sfixn pos = I - Q;
    sfixn j = alpha + 1 + pos;
    sfixn a = inv_mod(primes_GPU[alpha], primes_GPU[j]);

    D[I] = a;
    L[I] = primes_GPU[j] - a;
  }
}


/** procedure copying the elements of the Taylor shift by 1 in Z/pZ in X
* @X, output matrix containing the coefficients of the Polynomial shifted modulo prime numbers
* @s, number of prime numbers used in the program
* @Polynomial_shift, input polynomial shifted by 1 modulo a prime number (primes[s])
* @n, size of Polynomial_shift
* @step, integer allowing to know what is the column of X we fill for this step of rank 'step'
**/
__global__ void copyX(sfixn *X, sfixn s, sfixn *Polynomial_shift, sfixn n, sfixn step)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n)
    X[i * s + step] = Polynomial_shift[i];
}


/** phase 1 procedure computing (...(X A_{k0})...)A_{s-T}
* Phase 1 is used when rows of X have at least 33 elements.
* It works until there is at most 32 active elements.
* @X, input matrix containing the coefficients of the Polynomial shifted modulo prime numbers ;
      output matrix such that X := (...(X A_{k0})...)A_{s-T}
* @s, number of prime numbers used in the program
* @n, size of the input polynomial
* @primes_GPU, input array of prime numbers
* @D, input array containing diagonal elements (m_{i,j}) of all the matrices A_k
* @L, input array containing diagonal elements (n_{i,j}) of all the matrices A_k
**/
__global__ void mixed_radix_phase1(sfixn* X, sfixn s, sfixn n, sfixn* primes_GPU, sfixn* D, sfixn* L)
{
  __shared__ sfixn Z[512]; // ideally Z[s] but impossible as s is not a const int !!!
  sfixn B = ceil(s / (double) T_CHN); // B is the shared part
  sfixn b = blockIdx.x;
  sfixn i = threadIdx.x;
  sfixn j, Q, k, temp1, temp2, p0;

  if(b<n)
  {
    // store X in the shared memory
    for (j=B*i; (j < B*(i+1)) && (j<s); j++)
      Z[j] = X[b*s+j];
    __syncthreads();


    // phase 1 : computes (...(X A_1)...)X_{s-T} with the strategy '1 block = 1 row'
    // (one block computes a row of the matrix product)
    /*
       Let's reduce the product (...(X A_{k0})...)A_{s-T} considering a loop containing
       product like 'Y = X A_k' :

       For all j in [1,k-1] :  Y_j = X_j
       For all j in [k,s]   :  Y_j = ( X_j * m_k,j + X_k * n_k,j )  mod m_j

       This is what we do in the following code, in the second 'for' loop :

       p0 := m_j
       temp1 := X_k * n_k,j  mod m_j
       temp2 := X_j * m_k,j  mod m_j
       
       then Z[j] = temp1 + temp 2  mod m_j
       so   Z[j] = ( X_j * m_k,j + X_k * n_k,j )  mod m_j
    */
    for (k=0; k < s-T_CHN; k++)
    {
      B  = ceil((s-k) / (double) T_CHN);
      Q  = k*(s-1) - (k-1)*k/2;

      for (j=B*i+k+1; (j < B*(i+1)+k+1) && (j<s); j++) // one  thread computes B positions
      {
        p0 = primes_GPU[j];
        temp1 = mul_mod(Z[k], L[Q+j-k-1], p0);
        temp2 = mul_mod(Z[j], D[Q+j-k-1], p0);
        Z[j]  = add_mod(temp1, temp2, p0);
      }
      __syncthreads();
    }

    // store Z in the vector variable X (global memory)
    B = ceil(s / (double) T_CHN);
    for (j=B*i; (j < B*(i+1)) && (j<s); j++)
      X[b*s+j] = Z[j];
  }
}


/** procedure computing (...(X A_{s-T+1})...)A_{s-1}
* Phase 2 assumes that there is at most 32 active elements per row
* @X, input matrix containing the output matrix X of the procedure mixed_radix_phase1, this phase continue the phase1 of the computation which is not over ;
      output matrix corresponds to X := X (...(X A_{s-T+1})...)A_{s-1}
* @s, number of prime numbers used in the program
* @n, size of the input polynomial
* @primes_GPU, input array of prime numbers
* @D, input array containing diagonal elements (m_{i,j}) of all the matrices A_k
* @L, input array containing diagonal elements (n_{i,j}) of all the matrices A_k
**/
__global__ void mixed_radix_phase2(sfixn* X, sfixn s, sfixn n, sfixn* primes_GPU, sfixn* D, sfixn* L)
{
  sfixn Z[T_CHN];
  sfixn i = blockIdx.x*blockDim.x + threadIdx.x;
  sfixn j, Q, k, temp1, temp2, p0;
  sfixn M, N, K;
  	
  // According to s, computations are different, indeed, if s is not big enough
  // then there's no need to do the phase1 and then the phase 2 is the only phase to do
  // so as to have the mixed radix representation of X.
  if (s-T_CHN<0)
  {
    M = s;
    N = 0;
    K = 0;
  }
  else
  {
    M = T_CHN;
    N = s-T_CHN;
    K = s- T_CHN;
  }

  if(i<n)
  {
    // store the last coefficients of the row that the thread will modify in the local memory
    for (j=0; j<M; j++)
      Z[j] = X[i*s + N + j];

    // phase 2 : computes (...(X A_{s-T+1})...)X_{s-1} with the strategy 1 'thread = 1 row'
    // same idea than for the procedure mixed_radix_phase1
    for (k=0; k<M; k++, ++K)
    {
      Q  = K*(s-1) - (K-1)*K/2;

      for (j=k+1; j<M; j++)
      {
        p0 = primes_GPU[N+j];
        temp1 = mul_mod(Z[k], L[Q+j-k-1], p0);
        temp2 = mul_mod(Z[j], D[Q+j-k-1], p0);
        Z[j]  = add_mod(temp1, temp2, p0);
      }
    }

    // store Z in the vector variable X (global memory)
    for (j=1; j<M; j++)
      X[i*s + N + j] = Z[j];
  }
}


/** procedure computing (...(X A_1)...)A_{s-1}
* @X, input matrix containing the coefficients of the Polynomial shifted modulo prime numbers ;
      output matrix corresponds to (considering X before using mixed_radix_phase1)
      X := X (...(X A_{1})...)A_{s-1}
* @s, number of prime numbers used in the program
* @n, size of the input polynomial
* @primes_GPU, input array of prime numbers
**/
void recombine(sfixn* X, sfixn s, sfixn n, sfixn *primes_GPU)
{
  sfixn nb_blocks;
  sfixn *D, *L;
  cudaMalloc( (void **) &D, s*(s-1)/2 * sizeof(sfixn) );
  cudaMalloc( (void **) &L, s*(s-1)/2 * sizeof(sfixn) );

  nb_blocks = number_of_blocks(s*(s-1)/2);
  createDL<<<nb_blocks, NB_THREADS>>>(D, L, primes_GPU, s);

  // computes (...(X A_{k0})...)A_{s-T}
  mixed_radix_phase1<<<n, T_CHN>>>(X, s, n, primes_GPU, D, L);

  // computes (...(X A_{s-T+1})...)A_{s-1}
  nb_blocks = ceil(n / (double) T_CHN);
  mixed_radix_phase2<<<nb_blocks, T_CHN>>>(X, s, n, primes_GPU, D, L);
}


/** procedure computing (...(X A_1)...)A_{s-1} in the CPU
* @X, input matrix containing the coefficients of the Polynomial shifted modulo prime numbers ;
      output matrix corresponds to (considering X before using mixed_radix_phase1)
      X := X (...(X A_{1})...)A_{s-1}
* @s, number of prime numbers used in the program
* @n, size of the input polynomial
* @D, input array containing diagonal elements (m_{i,j}) of all the matrices A_k
* @L, input array containing diagonal elements (n_{i,j}) of all the matrices A_k
**/
void mixed_radix_phase_cpu(sfixn* X, sfixn s, sfixn n, sfixn* D, sfixn* L)
{
  sfixn b, j, Q, k, temp1, temp2, p0;


 //////////// TEST //////////////
/*if (n == 8){
 printf("\nL : ");
 for (i=0;i<s*(s-1)/2; i++)
   printf("%d ", L[i]);
 printf("\nD : ");
 for (i=0;i<s*(s-1)/2; i++)
   printf("%d ", D[i]);
  sfixn X2[n*s];
  for (i=0; i<n*s; i++)
  {
    if (i%4 == 1)
      X2[i] = 2;
    else
      X2[i] = 0;
  }
 printf("\nX2 : ");
 for (i=0;i<s*n; i++)
   printf("%d ", X2[i]);
 printf("\n");
  for (k=0; k < s-1; k++)
  {
    Q  = k*(s-1) - (k-1)*k/2;
//    printf("Q = %d\n", Q);

    for (b=0; b<n; b++)
      for (j=k+1; j<s; j++)
      {
        printf("b=%d, j=%d, Q=%d, Q+j=%d\n", b, j, Q, Q+j);
        p0 = primes[j];
        temp1 = mul_mod(X2[b*s+k], L[Q+j-k-1], p0);
        temp2 = mul_mod(X2[b*s+j], D[Q+j-k-1], p0);
        X2[b*s+j] = add_mod(temp1, temp2, p0);
      }
  }

 printf("\nX2 : ");
 for (i=0;i<s*n; i++)
   printf("%d ", X2[i]);
 printf("\n");
}*/
 ////////////////////////////////



  // phase 1 : computes (...(X A_1)...)X_{s-T} with the strategy 1 row = 1 block
  for (k=0; k < s-1; k++)
   {
    Q  = k*(s-1) - (k-1)*k/2;

    for (b=0; b<n; b++)
      for (j=k+1; j<s; j++)
      {
        p0 = primes[j];
        temp1 = mul_mod(X[b*s+k], L[Q+j-k-1], p0);
        temp2 = mul_mod(X[b*s+j], D[Q+j-k-1], p0);
        X[b*s+j] = add_mod(temp1, temp2, p0);
      }
  }
}
