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


#ifndef _LIST_OPT_POLY_MUL_H_
#define _LIST_OPT_POLY_MUL_H_

#include <iostream>
#include <ctime>
#include <cmath>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

//#include "../include/naive_poly_mul.h"
#include "inlines.h"

/*
 IMPORTANT
 n >= m
 That is, degree(A) >= degree(B)
 where A and B are in the input of the kernel mul and mulCUDA
 */

/*********************************
 *    Tx*BLOCK <= 2^11             *
 *    where Tx is the number of    *
 *    threads per block and BLOCK  *
 *    is the number of coefficients*
 *    per thread.                  *
 *                                 *
 *    (n*m)/(Tx*BLOCK^2) <  2^16   *
 *    Note that this expression    *
 *    gives the number of thread   *
 *    blocks necessary for that    *
 *    muktiplicaiton.              *
 *********************************/
//const sfixn BASE_1 = 31;
const sfixn BLOCK = 4;
const sfixn Tx = 256;
const sfixn T = Tx * BLOCK + BLOCK - 1;

void
mulCUDA (sfixn *A, sfixn *B, sfixn *C, sfixn n, sfixn m, sfixn Prime,
				 sfixn verify);

__global__ void
zeroC (sfixn *C, sfixn n);

__global__ void
copyA (sfixn *A, sfixn *C, sfixn n);

__global__ void
addFinal (sfixn *Cgpu, sfixn *CgpuFinal, sfixn c, sfixn n, sfixn w, sfixn P);

__global__ void
add (sfixn *Cgpu, sfixn *CgpuFinal, sfixn c, sfixn n, sfixn w, sfixn i, sfixn j,
		 sfixn k, sfixn P);

__global__ void
mul (sfixn *A, sfixn *B, sfixn *C, sfixn n, sfixn m, sfixn c, sfixn P,
		 sfixn unitBlocks);

#endif

