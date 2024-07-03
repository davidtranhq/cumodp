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



#ifndef _LIST_INV_STOCKHAM_H_
#define _LIST_INV_STOCKHAM_H_

#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "list_stockham.h"
#include "printing.h"

void list_inv_stockham_host(sfixn *X, sfixn m, sfixn k, sfixn w, sfixn p);

__global__ void list_inv_mul_ker(sfixn *X, sfixn invn, sfixn m, sfixn p, double pinv, sfixn length_layer);

__global__ void inv_mul_ker(sfixn *X, sfixn invn, sfixn m, sfixn, sfixn p, double pinv);

__device__ __host__ __inline__ sfixn get_num_of_blocks_for_inv(sfixn m, sfixn k)
{
	sfixn n = (1L << k);
	if (m < 512)
	{
		if (n > m)
			return (n/m + 1);
		else
			return n;
	}
	else
	{
		return ((n * m)/512);
	}
}
__device__ __host__ __inline__ sfixn get_num_of_threads_for_inv(sfixn m, sfixn k)
{
	if (m < 512)
	{
			return m;
	}
	else return 512;
}

#endif
