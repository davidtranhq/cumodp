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


#ifndef BIG_ARITHMETIC_AUX_DEVICE_H_P4
#define BIG_ARITHMETIC_AUX_DEVICE_H_P4

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"
/****************************************************/
//__global__ void
//kernel_p4_fft2_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
//{
//	short operation = parameters[0];
//	usfixn64 permutationStride = parameters[5];
////	short shuffle = parameters[6];
//	short padding = 0;	//parameters[?]
//	usfixn64 idx;
//	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);
//
//	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
//	//following indexing is slightly faster than above indexing
//
//	idx = tid;
//	if (padding == 0)
//		fft_base2_permutated (&xs[idx], &ys[idx], permutationStride);
//}

#endif
