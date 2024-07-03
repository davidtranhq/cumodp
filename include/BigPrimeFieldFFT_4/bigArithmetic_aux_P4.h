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


#ifndef BIG_ARITHMETIC_AUX_H_P4
#define BIG_ARITHMETIC_AUX_H_P4

#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"
/****************************************************/
//__global__ void
//kernel_p4_fft2_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
//;
/****************************************************/

__global__ void
kernel_p4_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters);

/************************************************/

__global__ void
kernel_p4_transposition_v0 (usfixn64* xs, usfixn64 * ys,
														usfixn64 nPermutationBlocks, usfixn64 blockStride,
														usfixn64* parameters);
/********************************************************/
__global__ void
kernel_p4_transposition (usfixn64* xs, usfixn64 * ys,
												 usfixn64 nPermutationBlocks, usfixn64 blockStride,
												 usfixn64* parameters);

/********************************************************/

__global__ void
kernel_p4_cpy_d2d_v0 (usfixn64* __restrict__ dest,
											const usfixn64 * __restrict__ src, const usfixn64 n,
											const usfixn64* __restrict__ parameters);
/********************************************************/

__global__ void
kernel_p4_cpy_d2d (usfixn64* __restrict__ dest,
									 const usfixn64 * __restrict__ src, const usfixn64 n,
									 const usfixn64* __restrict__ parameters);

/********************************************************/
void
host_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd);

#endif
