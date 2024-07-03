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



#ifndef _TS_VCA_COUNT_H_
#define _TS_VCA_COUNT_H_
#include "taylor_shift_conf.h"


// compute the sign of a polynomial in MRR-representation in the GPU
__device__ sfixn MRR_sign(sfixn *b, sfixn s, sfixn *primes_GPU);


// compute the sign of a polynomial in MRR-representation in the CPU
sfixn MRR_sign_CPU(sfixn *b, sfixn s);


// compute locally a number of sign change in X and store data in the arrays count_sg and bound_sg
/* count_sg : each thread computes a local number of sign change
              and stores this number in one position of count_sg
   bound_sg : each thread stores the first and the last non-null
              sign : it's possible we don't count all the change
              as each thread doesn't know what was the first sign
              after the part it deals with and also for the first
              sign before this same part.
   Tsign : array containing the sign of each coefficient i.e. the
           sign of each row of X.                                 */
__global__ void sign_change_MRR(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg);
__global__ void sign_change_MRR2(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg, sfixn *Tsign);


// computes on the CPU the number of sign changes missed by the GPU
sfixn nb_bounds_change_sign(sfixn *bound_sg, sfixn size);


// computes the maximum in absolute value of an array and store it in MAXI[0]
__global__ void max_gpu(sfixn *T, sfixn n, sfixn *MAXI);


// computes the sum of all the positions of an array and store it in ADD[0]
__global__ void global_add_gpu(sfixn *T, sfixn n, sfixn *ADD);


// transform X in MRR representation
/* input  : X in CRT representation
   output : X in MRR representation */
void MRR_all(sfixn *X, sfixn n, sfixn e, sfixn s, sfixn *primes_GPU);


// function counting the number of sign change in the MRR representation (X) of the polynomial we want
sfixn count_all(sfixn *X, sfixn n, sfixn s, sfixn *primes_GPU, sfixn *count_sg, sfixn *bound_sg, sfixn *ADD, sfixn m);


#endif // _TS_VCA_COUNT_H_
