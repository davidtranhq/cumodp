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



#include "taylor_shift.h"
#include "taylor_shift_cpu.h"
#include "taylor_shift_conf.h"
#include "taylor_shift_kernel.h"
#include "ts_chinese_remainder.h"
#include "taylor_comb.h"
#include "ts_VCA_count.h"
#include "ts_VCA_all.h"
#include "primesFFT.h"

/** display of an array
* @T, input array
* @size, size of T
* This procedure displays signs : '+', '-' and '0'.
**/
void display_array2(sfixn *T, sfixn size)
{
  sfixn k, i;
  sfixn M = (sfixn) ceil(size/ (double) T_CHN);

  for (i=0; i<M; i++)
  {
    for (k=i*T_CHN; (k<(i+1)*T_CHN) && (k<size); k++)
    {
      if (T[k] == 1)
        printf("+ ");
      else if (T[k] == -1)
        printf("- ");
      else
        printf("%d ", T[k]);
    }
    printf("\n");
  }
}


/**
* @file, name of the file containing the coefficients of the polynomial considered, say p
* n-1, degree of p
* @a, integer to do the Taylor shift by a
* @k, integer in p_{k, a} := 2^{kn}p((x+a)/2^k).
* This procedure counts the number of sign changes of p_{k, a}
after doing several computations
  with functions and procedures detailled in all the files.
* ! : n must be such that n = 2^e with e in [3,13].
      Indeed, for e>13, we need more prime numbers than 512.
**/
sfixn main_prog(sfixn a, sfixn k, char* file)
{
  // temporal data
//  float total_time;
//  cudaEvent_t start, stop;     /* Initial and final time */

  // TIME
/*  cudaEventCreate(&start); 
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);*/


  // declaration of variables
  sfixn n, e, s, f, i;//,S;
  sfixn nb_blocks;
  n = size_file(file);
  e = (sfixn) log2((double) n);
  //
  f = 28;
  // s should the number of prime numbers sufficient to represent the coefficients of the input polynomial in CRT representation
  s = (sfixn) ceil( (n + 2) / (double) f);
  // S is the number of primes numbers stored in the file /include/primesFFT.h
  //S = NS; this line is made comment by Sardar Haque on 6th December as S is never used afterwards. And it creates compiler warning

  // arrays and initialization
  // X should be removed when it will be a parameter
  // X is an array representing a matrix of size n*s where each line gives the CRT representation of each coefficient of the polynomial in 'input' 
  // Polynomial and Polynomial_device won't ne needed when X will be a parameter
  sfixn *X, *Polynomial, *Polynomial_device;
  cudaMalloc( (void **) &X , n * s * sizeof(sfixn) );
  Polynomial = (sfixn*) malloc(n * sizeof(sfixn));
  cudaMalloc( (void **) &Polynomial_device , n * sizeof(sfixn) );
  stock_file_in_array(file, n, Polynomial);
  cudaMemcpy(Polynomial_device, Polynomial, n*sizeof(sfixn), cudaMemcpyHostToDevice);
  free(Polynomial);

  /* ----------------------------------------------------------- */

  // declaration of the arrays used in the GPU in the procedure taylor_shift_Z
  sfixn *Factorial_device, *Powers_a;
  sfixn *Monomial_shift_device;
  sfixn *Mgpu, *fft_device;
  sfixn *Polynomial_shift_device[2];


  // allocation of arrays we use for Taylor_shift_Z
  cudaMalloc( (void **) &Factorial_device, (n+1) * sizeof(sfixn) );
  cudaMalloc( (void **) &Powers_a, (n+1) * sizeof(sfixn) );
  cudaMalloc( (void **) &Monomial_shift_device , n * sizeof(sfixn) );
  cudaMalloc( (void **) &Polynomial_shift_device[0], n * sizeof(sfixn) );
  cudaMalloc( (void **) &Polynomial_shift_device[1], n * sizeof(sfixn) );
  cudaMalloc((void **)&Mgpu, n * sizeof(sfixn));
  cudaMalloc( (void **) &fft_device, 2 * n * sizeof(sfixn) );


  //  initialization of X
  // !!! : should be removed when X will be a parameter of the function main_prog
  nb_blocks = number_of_blocks(s*n);
  Zeros_GPU<<<nb_blocks, NB_THREADS>>>(X, n*s);


  // returns CRT representation of P_{k,a}
  taylor_shift_Z(a, k, n, e, f, s, X, Polynomial_device, Factorial_device, Powers_a, Monomial_shift_device, Mgpu, Polynomial_shift_device, fft_device);


  // deallocation
  cudaFree(Polynomial_device);
  cudaFree(Factorial_device);
  cudaFree(Powers_a);
  cudaFree(Monomial_shift_device);
  cudaFree(Mgpu);
  cudaFree(fft_device);
  for (i=0; i<2; i++)
    cudaFree(Polynomial_shift_device[i]);

  /* ----------------------------------------------------------- */

  // converting CRT in MRR representation
  sfixn *primes_GPU;
  cudaMalloc( (void **) &primes_GPU , s * sizeof(sfixn) );
  cudaMemcpy(primes_GPU, primes, s*sizeof(sfixn), cudaMemcpyHostToDevice);

  MRR_all(X, n, e, s, primes_GPU);

  /* ----------------------------------------------------------- */


  // counting
  sfixn *count_sg, *bound_sg, *ADD;
  sfixn V;
  sfixn m = (sfixn) ceil(n/(double) T_CHN);
  cudaMalloc( (void **) &count_sg, m * sizeof(sfixn) );
  cudaMalloc( (void **) &bound_sg, 2*m * sizeof(sfixn) );
  cudaMalloc( (void **) &ADD, 1 * sizeof(sfixn) );

  // count the number of sign changes of X in MRR representation
  V = count_all(X, n, s, primes_GPU, count_sg, bound_sg, ADD, m);

  printf("\nV = %d\n", V);

/*  sfixn *B, *C;
  C = (sfixn*) malloc(m * sizeof(sfixn));
  B = (sfixn*) malloc(2*m * sizeof(sfixn));
  cudaMemcpy(C, count_sg, m * sizeof(sfixn), cudaMemcpyDeviceToHost);
  cudaMemcpy(B, bound_sg, 2*m * sizeof(sfixn), cudaMemcpyDeviceToHost);

  printf("\ncount_sg :\n");
  display_array(C, m);
  printf("\n-----------------------------------------\n");
  printf("\nbound_sg :\n");
  display_array2(B, 2*m);

  sfixn NB = nb_bounds_change_sign(B, 2*m);
  printf("\nnb_chg_sign in bounds = %d\n\n", NB);*/

/////////////////////////
  cudaFree(count_sg);
  cudaFree(bound_sg);
  cudaFree(primes_GPU);

  // X contient les b_i pour chaque coeff, maintenant les associer pour avoir a coeff avec la representation donnee
  cudaFree(X);


  // TIME
/*  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&total_time, start, stop);
  cudaEventDestroy(stop);
  total_time /= 1000.0;
  printf("execution time = %.6lf\n\n", total_time);*/

  return 0;
}
