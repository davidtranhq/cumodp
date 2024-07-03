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



#ifndef _OPT_POLY_GCD_H_
#define _OPT_POLY_GCD_H_

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
*/

const sfixn T = 480;

__global__ void	shrinkA(sfixn *A, sfixn *status);
__global__ void	shrinkB(sfixn *B, sfixn *status);
__global__ void MonicConversionA( sfixn *A, sfixn *status, sfixn p);
__global__ void MonicConversionB( sfixn *B, sfixn *status, sfixn p);
__global__ void	reduceMgcd(sfixn *B,  sfixn *status);
__global__ void	reduceNgcd(sfixn *A,  sfixn *status);

__global__ void	status3(sfixn *status);
__global__ void	status4(sfixn *status);
__global__ void	status5(sfixn *status);

__global__ void gcdGPU(sfixn *A, sfixn *B, sfixn *status);



float gcdPrimeField(sfixn *Dgcd, sfixn *A, sfixn *B, sfixn *G, sfixn n, sfixn m, sfixn p, sfixn optionNormalize);

#endif
