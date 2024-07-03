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



#ifndef _CONDENSATION_FF_H_
#define _CONDENSATION_FF_H_

#include "types.h"
#include "inlines.h"

#include <cstdio>
#include <cmath>

#include <stdio.h>

//sfixn detFF;
float detFFtime;

//RULES OF THUMB FOR THIS CODE:
//NEVER NEVER NEVER 
//MAKE BLOCK >=	n

/************************************************/
void
egcdCPUFF (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo, int P);

/************************************************/
sfixn
inv_modCPUFF (sfixn n, int P);

/************************************************/
sfixn
mul_modCPUFF (sfixn a, sfixn b, int P);

/************************************************/
sfixn
sub_modCPUFF (sfixn a, sfixn b, int P);

/************************************************/
__global__ void
condensationFF (int *Aa, int *Bb, int *L, size_t n, int P, int BK);
/************************************************/
__global__ void
findingLFF (int *A, int *B, int *C, size_t n, int P);
/************************************************/
int
detFF (int *A, int n, int P, int BK, int TH, int BS, int verbose);
/************************************************/
#endif
