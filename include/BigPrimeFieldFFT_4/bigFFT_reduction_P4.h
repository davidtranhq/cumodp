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


// for 16-word prime
#ifndef BIG_FFT_REDUCTION_H_
#define BIG_FFT_REDUCTION_H_

#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"
/***************************************************/
__global__ void
kernel_p4_reduction_16_v0 (usfixn64 * xs, usfixn64 * ys, usfixn64 prime,
												usfixn64 * precomputed_r_array, short index,
												usfixn64 * parameters)
;

/***************************************************/
__global__ void
kernel_p4_reduction_16 (const usfixn64 * __restrict__ xs, usfixn64 * __restrict__ ys,
										 const usfixn64 * __restrict__ primeList, const usfixn64 * __restrict__ precomputed_r_array,
										 const usfixn64 * __restrict__ parameters)
;
/***************************************************/
__global__ void
kernel_p4_reduction_16_v2 (usfixn64 * xs, usfixn64 * ys1, usfixn64* ys2,
												usfixn64 *primeList, usfixn64 * precomputed_r_array,
												usfixn64 * parameters)
;

#endif /* end of BIG_FFT_REDUCTION_H_ */
