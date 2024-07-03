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



#ifndef _SMALL_FUNC_SOLVE2_H_
#define _SMALL_FUNC_SOLVE2_H_

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
#include "printing.h"
#include "cudautils.h"
#include "types.h"
#include "list_inv_stockham.h"
#include "subproduct_tree.h"
#include "list_stockham.h"

/*
T is the number of threads in a thread block.
*/

const int T = 512;


__global__ void diffUniDevPar(sfixn *A, sfixn n, sfixn p);
		
void lcBiDev(sfixn *dest, sfixn *source, sfixn *dgs);
void tailBiDev(sfixn *dest, sfixn *source, sfixn *dgs);
void diffUniDev(sfixn *dest, sfixn *source, sfixn n, sfixn p);

#endif
