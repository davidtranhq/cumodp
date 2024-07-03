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



#ifndef _SUBPRODUCT_TREE_H
#define _SUBPRODUCT_TREE_H
/**
 * this part is for building the subproduct tree
 *
 * @author: Jiajian Yang
 * */
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "printing.h"
#include "cudautils.h"
#include "fft_aux.h"
#include "naive_poly_mul.h"
#include "stockham.h"
#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "list_inv_stockham.h"

void subproduct_tree_host(sfixn *X, sfixn *M, sfixn k, sfixn p);

__host__ __device__ __inline__ sfixn get_subtree_size(sfixn k)
{
	sfixn sum = 0;
	for(sfixn i = 1; i <= k-1; ++i)
		sum += ( (1L <<(i-1)) +1 ) * ( 1L << (k-i)) ;
	return sum;
}


__host__ __device__ __inline__ sfixn get_layer_size(sfixn k, sfixn i)
{
	sfixn deg = (1L << (i-1));
	return (deg + 1)*(1L<<(k-i));
}

__host__ __device__ __inline__ sfixn get_polylength_on_layer(sfixn i)
{
	return ((1L << (i-1)) +1);
}

__host__ __device__ __inline__ sfixn get_subtree_depth(sfixn k)
{
	return k-1;
}

__host__ __device__ __inline__ sfixn get_layer_offset(sfixn k, sfixn i)
{
	sfixn sum = 0;
	for (sfixn l = 1; l < i; l++)
		sum += get_layer_size(k,l);
	return sum;
}

void subproduct_tree_host_bchmk(sfixn *X, sfixn *M, sfixn k, sfixn p);
void subproduct_tree_bench(sfixn p, sfixn k);

void subproduct_tree_dev(sfixn *M_dev, sfixn k, sfixn p, double pinv);


#endif

