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


#ifndef _CONDENSATION_REAL_H_
#define _CONDENSATION_REAL_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <time.h>
#include "types.h"

const sfixn THno = 128;

typedef double realV;

long double detDirect, detCon;
double timeDet;
//RULES OF THUMB FOR THIS CODE:
//NEVER NEVER NEVER 
//MAKE BK >=	n

long double
mulSeq (int n, long double *A);

void
condensation (realV *Aa, realV *Bb, int n, int l, int BKReq, int BK, int BKID,
							int THid);
void
det (realV *A, int n, int BS, int BK);

void
direct (realV *A, int n);

__global__ void
find (realV *A, int *ls, realV *lvalues, int n, int temp);

__global__ void
divideToOne (realV *A, int *B, realV *C, int n, int temp);

__global__ void
condensGPU (realV *Aa, realV *Bb, int n, int *l, int BK);

void
deter (realV *A, int n, int BS, int BK);
#endif
