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


# ifndef _TAYLOR_SHIFT_KERNEL_H_
# define _TAYLOR_SHIFT_KERNEL_H_

// fast multiplication of two polynomials, created by Sardar Haque, I modify just a line to use it in my code
__global__ void listPlainMulGpu_and_right_shift_GPU(sfixn *Mgpu1, sfixn *Mgpu2 , sfixn length_poly, sfixn poly_on_layer, sfixn threadsForAmul, sfixn mulInThreadBlock, sfixn p, double pinv);


// create array identity (initialization of the array Fact)
__global__ void identity_GPU(sfixn *T, sfixn size);


// create array Powers_a (initialization of the array Powers_a)
__global__ void identity_a_GPU(sfixn *T, sfixn n, sfixn a);


// create all the elements of Factorial (mod p)
// !!! not used any more in the code !!!
__global__ void create_factorial_GPU(sfixn *Fact, sfixn n, sfixn e, sfixn p, double pinv);


// create all the elements of Factorial (mod p) first step
__global__ void create_factorial_step0_GPU(sfixn *Fact, sfixn n, sfixn e, sfixn p, double pinv);
// warning : n+1 is the size of Fact but we will just full the n last element, not the first one


// create all the elements of Factorial (mod p) other steps
// (after e-1 steps, Fact contains the factorial sequence)
__global__ void create_factorial_stepi_GPU(sfixn *Fact, sfixn n, sfixn e, sfixn p, double pinv, sfixn L);


// create all the elements of Powers_a (mod p)
__global__ void create_powers_a_step0_GPU(sfixn *Powers_a, sfixn n, sfixn e, sfixn p, double pinv);


// create all the elements of Powers_a (mod p) other steps
// (after e-1 steps, Powers_a contains the sequence (a^i)_{i in [0,n]})
__global__ void create_powers_a_stepi_GPU(sfixn *Powers_a, sfixn n, sfixn e, sfixn p, double pinv, sfixn L);


// create an array of the inverse numbers in Z/pZ
// !!! not used any more in the code !!!
__global__ void inverse_p_GPU(sfixn *T, sfixn p, double pinv);


// create the inverse of a number in Z/pZ
// !!! not used any more in the code !!!
__device__ sfixn inverse_GPU(sfixn k, sfixn p, double pinv);


// creates an array of the Newton's Binomials until n modulo p (! size of the array = n+1)
__device__ sfixn create_binomial_GPU(sfixn *Factorial, sfixn *Inverse_p, sfixn n, sfixn p, double pinv, sfixn id);


// create the Newton's Binomial coefficient "n choose id" modulo p
// return "n choose id" = n! / [id!(n-id)!] mod p
__device__ sfixn create_binomial2_GPU(sfixn *Factorial, sfixn n, sfixn p, double pinv, sfixn id);



// create the array of the coefficients of (x+1)^k for k in [1,2^(e-1)]
__global__ void develop_xshift_GPU(sfixn *T, sfixn n, sfixn *Factorial, sfixn p, double pinv);



// create the product of two arrays representing polynomials
// !!! not used any more in the code !!!
__device__ void conv_prod_GPU(sfixn *res, sfixn *T1, sfixn *T2, sfixn m, sfixn p, sfixn local_n);


// addition of two arrays : res = T1 + T2
__global__ void add_arrays_GPU(sfixn *res, sfixn *T1, sfixn *T2, sfixn size, sfixn p);


// creates an array of zeros
__global__ void Zeros_GPU(sfixn *T, sfixn n);


// creates an array of zeros
__global__ void Ones_GPU(sfixn *T, sfixn n);


// do the first step of the tree to compute the Taylor_shift by 1
__global__ void init_polynomial_shift_GPU(sfixn *Polynomial, sfixn *Polynomial_shift, sfixn n, sfixn p);
/* EXAMPLE for n=8 :
   after this procedure, Polynomial_shift = [f0+f1, f1, f2+f3, f3, f4+f5, f5, f6+f7, f7] */


// transfer at each step the polynomials which need to be multiplicated 
__global__ void transfert_array_GPU(sfixn *Mgpu, sfixn *Polynomial_shift, sfixn *Monomial_shift, sfixn n, sfixn local_n, sfixn p, double pinv);

  /*                   EXAMPLE

  --------------------------------------------------------

     ARRAY Polynomial_shift_device[i-1] considered
       _________ _________ _________ _________
      |         |         |         |         |
      |         |    X    |         |    Y    |
      |_________|_________|_________|_________|
        part=0    part=1    part=2    part=3
     
     local_n = size of a part

  --------------------------------------------------------

     ARRAY Mgpu[i] considered
       ___________________ ___________________
      |                   |                   |
      |   X      (x+1)^m  |   Y      (x+1)^m  |
      |___________________|___________________|
              PART=0               PART=1
     
     B = 2 * local_n = size of a PART
     m = local_n

We want to fill the array Mgpu[i] like this : the polynomials
which need to be multiplicated by (x+1)^m are of odd part and
we store them at the beginning of each PART of Mgpu[i]. The end
of each part doesn't really contain (x+1)^m as we need arrays 
to be multiplicated, so we avoid the multiplication by 1.
Thus the end of each PART contains exactly :

 [(x+1)^m - 1] / x = m + ... + x^(m-1)    {m elements}          */


// do a right shift of an array
__global__ void right_shift_GPU(sfixn *T, sfixn n);


// add parts of three arrays between them
__global__ void semi_add_GPU(sfixn *NewPol, sfixn *PrevPol1, sfixn *PrevPol2, sfixn n, sfixn local_n, sfixn p);


// divide each coefficient of T by its corresponding power of 2 to do the PKC
__global__ void divide_powers_a(sfixn *T, sfixn *Powers_a, sfixn n, sfixn p, double pinv);


// multiply each coefficient of T by its corresponding power of a to do the Taylor_shift by a
__global__ void multiply_powers_a(sfixn *T, sfixn *Powers_a, sfixn n, sfixn p, double pinv);


// multiply each coefficient of T by its corresponding power of 2^k to do the PKC
__global__ void PKC_powers_2k(sfixn *T, sfixn *Powers_2k, sfixn n, sfixn p, double pinv);


// invert the coefficients of an array on the GPU
__global__ void invert_matrix_GPU(sfixn *X, sfixn n, sfixn s);
  /* Example :
     If X = { 0,  1,  2,
              3,  4,  5,
              6,  7,  8,
              9, 10, 11 },

     after using invert_array :
        X = { 9, 10, 11,
              6,  7,  8,
              3,  4,  5,
              0,  1,  0 }       */


// computes (x+1)^n * P_{k,c}( 1/(x+1) )
/* input  : X given the PKC in MRR representation of P_{k,c}(x+1)
   output : X with the invert order of rows, so given the MRR
            representation of (x+1)^n * P_{k,c}( 1/(x+1) )        */
void invert_matrix_X(sfixn *X, sfixn n, sfixn s);


# endif // _TAYLOR_SHIFT_KERNEL_H_
