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



#ifndef _NAIVE_POLY_MUL_H_
#define _NAIVE_POLY_MUL_H_

#include <stdio.h>
#include <stdlib.h>

#include "inlines.h"
#include "defines.h"
#include "types.h"
#include "printing.h"
#include "cudautils.h"
#include "rdr_poly.h"

__global__ void mul_eff_ker(sfixn *A, sfixn *B, sfixn n, sfixn *C, sfixn p, double pinv);

void mul_dev(sfixn *A, sfixn *B, sfixn n, sfixn *C, sfixn p);

__host__ __device__ void mul_ker(sfixn *A, sfixn *B, sfixn n, sfixn *C, sfixn p, double pinv);

#endif
