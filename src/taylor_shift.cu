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

/* Important to notice :

  n  : number of coefficients of the polynomial considered
       (size of the polynomial)
 n-1 : degree of the polynomial considered
  p  : prime number, it must be greater than n
*/


/* ------------- Taylor_shift procedure modulo Zp ------------- */

/** This procedure is called s times by the CPU with s different prime
    numbers. It gives the CRT representation of the polynomial
    2^{k*n}*P( (x+a)/2^k ) modulo p (called PKC), combination of a
    Taylor_shift by a of a input polynomial P and multiplications for
    each coefficient by a power of 2^k. In each call, one column of X
    will be filled modulo p.

    PARAMETERS :
    ------------
*   @n : number of coefficients of the polynomial considered
         (size of the polynomial)
*   @e : e = log2(n) so n = 2^e
*   @p : prime number
*   @pinv : invert of p
*   @X : array of size n*s containing CRT and MRR representation of
         the output polynomial
*   @step : indicates what column of X we will fill
*   @s : number of prime numbers needed
*   @Polynomial_device : input polynomial
*   @a : integer for the Taylor_shift by a
*   @k : power of 2 considered
*   @Factorial_device : array containing the sequence (i!)_{i in [0;n]}
*   @Powers_a : array containing in this code :
                -> firstly : the sequence (a^i)_{i in [0;n]} mod p
                -> secondly : the sequence ((2^k)^i)_{i in [0;n]} mod p
*   @Monomial_shift_device : array containing the sequence of the coeffi-
                             cients of ((x+1)^{2^i})_{i in [0;e-1]} mod p
*   @Mgpu : array used in several intermediate computations mod p
*           this array has size n
*   @Polynomial_shift_device : array in 2 dimensions containing pairwise
                               intermediate computations to obtain the
                               Taylor_shift by 1 of a polynomial mod p
             this array has 2 rows of n elements
*   @fft_device : array used for Fast Fourier Transform mod p      
                  The size of this array is 2 n
    */

void taylor_shift_Zp(sfixn n, sfixn e, sfixn p, double pinv, sfixn *X, sfixn step, sfixn s, sfixn *Polynomial_device, sfixn a, sfixn k, sfixn *Factorial_device, sfixn *Powers_a, sfixn *Monomial_shift_device, sfixn *Mgpu, sfixn **Polynomial_shift_device, sfixn *fft_device)
{
  // declaration of variables
  sfixn i, nb_blocks;
  sfixn local_n = 2;

  float cpu_time, gpu_time, outerTime;
  cudaEvent_t start, stop;     /* Initial and final time */
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);


  // TIME
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&cpu_time, start, stop);
  cudaEventDestroy(stop);
  cpu_time /= 1000.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);


  // Create the array Factorial containing (i!)_{i in [0;n]}
  nb_blocks = number_of_blocks(n+1);
  cudaThreadSynchronize();
  identity_GPU<<<nb_blocks, NB_THREADS>>>(Factorial_device, n+1);
  cudaThreadSynchronize();

  nb_blocks = number_of_blocks(n/2);
  create_factorial_step0_GPU<<<nb_blocks, NB_THREADS>>>(Factorial_device + 1, n, e, p, pinv);
  cudaThreadSynchronize();
  sfixn L = 1;
  for (i=1; i<e; i++)
  {
    L *= 2;
    create_factorial_stepi_GPU<<<nb_blocks, NB_THREADS>>>(Factorial_device + 1, n, e, p, pinv, L);
    cudaThreadSynchronize();
  }

    // display Factorial_device
//    temp = (sfixn*) calloc(n, sizeof(sfixn));
//    cudaMemcpy( temp, Factorial_device, n*sizeof(sfixn), cudaMemcpyDeviceToHost );
//    printf("\nFactorial_device AFTER IDENTITY_GPU :\n");
//    display_array(temp, n);
//    free(temp);


  // Create the array Powers_a containing (a^i)_{i in [0;n]}
  nb_blocks = number_of_blocks(n+1);
  cudaThreadSynchronize();
  identity_a_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a, n+1, a);
  cudaThreadSynchronize();

  nb_blocks = number_of_blocks(n/2);
  create_powers_a_step0_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a + 1, n, e, p, pinv);
  cudaThreadSynchronize();
  L = 1;
  for (i=1; i<e; i++)
  {
    L *= 2;
    create_powers_a_stepi_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a + 1, n, e, p, pinv, L);
    cudaThreadSynchronize();
  }


    // display Powers_a
    /*  if (step == 0){
      temp = (sfixn*) calloc(n+1, sizeof(sfixn));
      cudaMemcpy( temp, Powers_a, (n+1)*sizeof(sfixn), cudaMemcpyDeviceToHost );
      printf("\nPowers_a :\n");
      display_array(temp, n+1);
      free(temp);
    } */

  // Create the array of the ((x+1)^{2^i})_{i in [0;e-1]}
  cudaThreadSynchronize();
  nb_blocks = number_of_blocks(n);
  develop_xshift_GPU<<<nb_blocks, NB_THREADS>>>(Monomial_shift_device, n, Factorial_device, p, pinv);
  cudaThreadSynchronize();


    // display Factorial_device
//    temp = (sfixn*) calloc(n, sizeof(sfixn));
//    cudaMemcpy( temp, Factorial_device, n*sizeof(sfixn), cudaMemcpyDeviceToHost );
//    printf("\nFactorial_device AFTER DEVELOP_XSHIFT_GPU :\n");
//    display_array(temp, n);
//    free(temp);

    // display Monomial_shift_device
//    temp = (sfixn*) calloc(n, sizeof(sfixn));
//    cudaMemcpy( temp, Monomial_shift_device, n*sizeof(sfixn), cudaMemcpyDeviceToHost );
//    printf("\nMonomial_shift_device :\n");
//    display_array(temp, n);
//    free(temp);


  /* ************************************************************

                     0th step : divide by Powers_a

     ************************************************************ */

  // initialize polynomial_shift dividing by a^i each coefficient
  nb_blocks = number_of_blocks(n);
  divide_powers_a<<<nb_blocks, NB_THREADS>>>(Polynomial_device, Powers_a, n, p, pinv);
  cudaThreadSynchronize();


  /* ************************************************************

                     1st step : initialization 

     ************************************************************ */

  // do the first step of the tree
  nb_blocks = number_of_blocks(n/2);
  init_polynomial_shift_GPU<<<nb_blocks, NB_THREADS>>>(Polynomial_device, Polynomial_shift_device[0], n, p);
  cudaThreadSynchronize();


  /* ************************************************************

                          next steps (i<10)

     ************************************************************ */

  sfixn polyOnLayerCurrent = n/2;
  sfixn mulInThreadBlock;
  sfixn I = 9;
  if (e < 9)
    I = e;


  // LOOP
  for (i=1; i<I; i++)
  {
    // transfer the polynomials which will be computed
    nb_blocks = number_of_blocks(n);
    transfert_array_GPU<<<nb_blocks, NB_THREADS>>>(Mgpu, Polynomial_shift_device[(i+1)%2], Monomial_shift_device, n, local_n, p, pinv);
    cudaThreadSynchronize();


    // Compute the product of the polynomials in Mgpu ('P2 * Bin' with Bin the array of binomials) and store them in Polynomial_shift_device[i%2] shifted at the right for the multiplication by x so do [( (x+1)^i - 1 ) / x ] * P2(x+1), then multiply it by x so we have [(x+1)^i - 1] * P2(x+1)
    mulInThreadBlock = (sfixn) floor((double) NB_THREADS / (double) (2*local_n));
    nb_blocks = (sfixn) ceil(((double) polyOnLayerCurrent / (double) mulInThreadBlock) * 0.5);
    listPlainMulGpu_and_right_shift_GPU<<<nb_blocks, NB_THREADS>>>(Mgpu, Polynomial_shift_device[i%2], local_n, polyOnLayerCurrent, 2*local_n, mulInThreadBlock, p, pinv);
    cudaThreadSynchronize();


    // add [(x+1)^i - 1] * P2(x+1) with P2(x+1) then we get (x+1)^i * P2(x+1) then do P1(x+1) + (x+1)^i*P2(x+1)
    nb_blocks = number_of_blocks(n/2);
    semi_add_GPU<<<nb_blocks, NB_THREADS>>>(Polynomial_shift_device[i%2], Mgpu, Polynomial_shift_device[(i+1)%2], n, local_n, p);
    cudaThreadSynchronize();


    // for the next step
    polyOnLayerCurrent /= 2;
    local_n *= 2;
  }



  /* ************************************************************

                       next steps : FFT (i >= 10)

     ************************************************************ */

  sfixn J = e;
  if (e < 9)
    J = 9;
  sfixn w;


  // LOOP
  for (i=9; i<J; i++)
  {

    // transfer the polynomials which will be FFTed and Mgpu
    nb_blocks = number_of_blocks(n);
    transfert_array_GPU<<<nb_blocks, NB_THREADS>>>(Mgpu, Polynomial_shift_device[(i+1)%2], Monomial_shift_device, n, local_n, p, pinv);
    cudaThreadSynchronize();
    nb_blocks = number_of_blocks(2*n);
    transfert_array_fft_GPU<<<nb_blocks, NB_THREADS>>>(fft_device, Mgpu, n, local_n);
    cudaThreadSynchronize();


    // Convert the polynomials in the FFT world
    w = primitive_root(i+1, p);
    list_stockham_dev(fft_device, polyOnLayerCurrent, i+1, w, p);
    cudaThreadSynchronize();


    // same operation than for ListPlainMul but in the FFT world
    nb_blocks = number_of_blocks(2*n);
    list_pointwise_mul<<<nb_blocks, NB_THREADS>>>(fft_device, 2*local_n, p, pinv, 2*n);
    cudaThreadSynchronize();


    // return to the real world
    w = inv_mod(w, p);
    list_stockham_dev(fft_device, polyOnLayerCurrent, i+1, w, p);
    cudaThreadSynchronize();


    // adjust the real coefficients : we need to multiplicate by the following w to have to correct size
    w = inv_mod(2*local_n, p);
    nb_blocks = number_of_blocks(n);
    mult_adjust_GPU<<<nb_blocks, NB_THREADS>>>(Polynomial_shift_device[i%2], fft_device, n, local_n, w, p, pinv);
    cudaThreadSynchronize();


    // semi_add
    nb_blocks = number_of_blocks(n/2);
    semi_add_GPU<<<nb_blocks, NB_THREADS>>>(Polynomial_shift_device[i%2], Mgpu, Polynomial_shift_device[(i+1)%2], n, local_n, p);
    cudaThreadSynchronize();


    // for the next steps
    polyOnLayerCurrent /= 2;
    local_n *= 2;
  }


  /* ************************************************************

                     last step : divide by Powers_a

     ************************************************************ */

  // multiply by Powers_a
  nb_blocks = number_of_blocks(n);
  multiply_powers_a<<<nb_blocks, NB_THREADS>>>(Polynomial_shift_device[(e-1)%2], Powers_a, n, p, pinv);
  cudaThreadSynchronize();


  /* ************************************************************

                         last step bis : PKC

     ************************************************************ */

  // Create the array Powers_a containing (a^i)_{i in [0;n]}
  sfixn pow2 = (sfixn) pow(2.0, k);
  nb_blocks = number_of_blocks(n+1);
  identity_a_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a, n+1, pow2);
  cudaThreadSynchronize();
  nb_blocks = number_of_blocks(n/2);
  create_powers_a_step0_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a + 1, n, e, p, pinv);
  cudaThreadSynchronize();
  L = 1;
  for (i=1; i<e; i++)
  {
    L *= 2;
    create_powers_a_stepi_GPU<<<nb_blocks, NB_THREADS>>>(Powers_a + 1, n, e, p, pinv, L);
    cudaThreadSynchronize();
  }


  // do the PKC to get P_{k,c}(x)
  nb_blocks = number_of_blocks(n);
  PKC_powers_2k<<<nb_blocks, NB_THREADS>>>(Polynomial_shift_device[(e-1)%2], Powers_a, n, p, pinv);
  cudaThreadSynchronize();


  /* ************************************************************

                         end : results

     ************************************************************ */


  // Copy the last array containing the Taylor shift by 1 of the input polynomial
/*  temp = (sfixn*) malloc(n * sizeof(sfixn));
  cudaMemcpy( temp, Polynomial_shift_device[(e-1)%2], n*sizeof(sfixn), cudaMemcpyDeviceToHost );
  cudaThreadSynchronize();*/

  // stockes the array of Newton's coefficients in a file
//  char name_file[100];
//  sprintf(name_file, "Pol%d.shiftGPU_%d.dat\0", e, p);
//  stock_array_in_file(name_file, temp, n);
//  free(temp);


  // TIME
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&gpu_time, start, stop);
  cudaEventDestroy(stop);
  gpu_time /= 1000.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);


  // add
  nb_blocks = number_of_blocks(n);
  copyX<<<nb_blocks, NB_THREADS>>>(X, s, Polynomial_shift_device[(e-1)%2], n, step);


  // TIME
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&outerTime, start, stop);
  cudaEventDestroy(stop);
  outerTime /= 1000.0;
  cpu_time += outerTime;


  // execution time
//  printf("  * cpu_time = %.6f s\n", cpu_time);
//  printf("  * gpu_time = %.6f s\n", gpu_time);
}
