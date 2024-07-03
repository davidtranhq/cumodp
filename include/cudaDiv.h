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



#ifndef _OPT_POLY_DIV_H_
#define _OPT_POLY_DIV_H_

#include <iostream>
#include <ctime>
#include <cmath>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"

#include "cudautils.h"
#include "types.h"
#include "list_inv_stockham.h"
#include "subproduct_tree.h"
#include "list_stockham.h"

/*
T is the number of threads in a thread block.
s is the length of the common part (see  HPCS 2012 paper)
*/

const sfixn T = 512;
const sfixn s = T;

/*
The number t is the ratio of the length of the private 
part (A^{+} or B^{+} according to HPCS 2012 paper)
to the comon (= shared) part, within  each thread block.
*/
const sfixn t = 2;


		
__global__ void statusUpdate(sfixn *A,  sfixn *status);
__global__ void zero(sfixn *A, sfixn n);
__global__ void	copyRem(sfixn *A, sfixn *R, sfixn *status);
__global__ void	reduceM(sfixn *B,  sfixn *status);
__global__ void divGPU(sfixn *A, sfixn *B, sfixn *status, sfixn *Q, sfixn *R);
float divPrimeField(sfixn *A, sfixn *B, sfixn *R, sfixn *Q, sfixn n, sfixn m, sfixn p);

#endif
