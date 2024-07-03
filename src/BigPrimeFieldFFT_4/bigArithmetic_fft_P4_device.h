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
#ifndef BIG_ARITHMETIC_FFT_DEVICE_H_
#define BIG_ARITHMETIC_FFT_DEVICE_H_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"
#include "bigArithmetic_addition_P4_device.h"
#include "bigArithmetic_multiplication_P4_device.h"
#include "bigArithmetic_subtraction_P4_device.h"
#include "bigArithmetic_cyclicShift_P4_device.h"

/************************************************/

__device__ void inline
device_p4_swap_permutated (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
	const usfixn64 permutationStride,
	const usfixn64 coefficientSize)
{
	short i = 0;
	usfixn64 tmp = 0;
	usfixn64 offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		tmp = xm[offset];
		xm[offset] = ym[offset];
		ym[offset] = tmp;
		offset += permutationStride;
	}
}
/************************************************/

__device__ void __inline__
device_p4_swap_permutated_2 (usfixn64 * __restrict__ xm,
	usfixn64 * __restrict__ ym,
	const usfixn64 & permutationStride)
{
	short i = 0;
	usfixn64 tmp = 0;
	usfixn64 offset = 0;

	usfixn64 y[COEFFICIENT_SIZE];
//	           , y[COEFFICIENT_SIZE];
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		tmp = xm[offset];
//		xm[offset] = ym[offset];
//		ym[offset] = tmp;
//		offset += permutationStride;
//	}
//
//	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		x[i] = xm[offset];
//		offset += permutationStride;
//	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		y[i] = ym[offset];
		offset += permutationStride;
	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		ym[offset] = xm[offset];  //x[i];
		offset += permutationStride;
	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		xm[offset] = y[i];
		offset += permutationStride;
	}
}
/************************************************/
/************************************************/

//new code for performing fft32, 16 threads for fft-32
__device__ __inline__ void
device_base_fft32 (usfixn64 *xm, const usfixn64 permutationStride,
	const usfixn64 tid)
{
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32; //stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r4[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r4[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 5
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r5[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r5[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	idx = (tid / 16) * 32;
	i = tid % (16);
	idx = idx + i;
	//round 6 permutation
	shNo = shNo_FFT_32_r6[i];
	if (shNo > 0)
		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
	shNo = shNo_FFT_32_r6[i + 16];
	idx += 16;
	if (shNo > 0)
		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
//	if(shNo<0)
//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx-shNo], permutationStride);
}

#endif
