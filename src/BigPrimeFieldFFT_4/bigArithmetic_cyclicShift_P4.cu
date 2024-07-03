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
#ifndef BIG_ARITHMETIC_CYCLIC_SHIFT_CU_
#define BIG_ARITHMETIC_CYCLIC_SHIFT_CU_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"
//#include "BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.h"
#include "bigArithmetic_cyclicShift_P4_device.h"

//__device__ void cyclicShift(usfixn64 *xs, short sn);
//__device__ void cyclicShift(usfixn64 *xs, short sn);

#include <stdio.h>

/****************************************************/
__global__ void
kernel_p4_cyclicShift_plain (usfixn64 * xs, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
{
	short operation = parameters[0];
	short iterations = parameters[1];
	short paddingMethod = parameters[2];
	short dynamicMemSize = parameters[3];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	short shiftNo = 3; //parameters ?
	short i;
//	, pos, pos1, t;
	unsigned short c = 0;
//	usfixn64 num1, num2;
//	usfixn64 *xd, *yd, *ud, *xm, *ym, *um;
//	usfixn64 *xd, *yd, *xm, *ym;
	usfixn64 *xd, *xm;
	shiftNo = tid & 0x7;
	shiftNo = 3;
//	num1 = 0;
//	num2 = R - 1;

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * 8);

	device_p4_cyclicShift_plain (xd, shiftNo);
}

/**********************************************/

__global__ void
kernel_p4_cyclicShift_permutated (usfixn64 * xs, usfixn64 *parameters)
{
	short operation = parameters[0];
	usfixn64 permutationStride = parameters[5];
	short shuffle = parameters[6];
	short padding = parameters[2];
	short nShift = 3;	//parameters [?]
	usfixn64 idx;
	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);

	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
	//following indexing is slightly faster than above indexing

	nShift = tid & 0x7;
	nShift = padding;
	idx = tid;

//		device_cyclicShift_permutated(&xs[idx], nShift, permutationStride);
//	if(padding==3)
//		device_cyclicShift_permutated_3(&xs[idx], nShift, permutationStride);

	//cyclic5 is a bad function for random input, IPC of 0.2 and ins_replay of 2
//	device_cyclicShift_permutated_5(&xs[idx], nShift, permutationStride);
//	device_cyclicShift_permutated_6(&xs[idx], nShift, permutationStride);
//	device_cyclicShift_permutated_7(&xs[idx], nShift, permutationStride);
//	device_cyclicShift_permutated_8(&xs[idx], nShift, permutationStride);
//	device_cyclicShift_permutated_9(&xs[idx], nShift, permutationStride);
//	device_cyclicShift_permutated_10(&xs[idx], nShift, permutationStride);

	//cyclic11 is best function with random input, IPC of 0.5 and ins_replay of 1.0

//	usfixn64 n = parameters[7];
//	for (tid; tid < n; tid += blockDim.x * gridDim.x)
	{
		idx = tid;
		if (nShift > COEFFICIENT_SIZE)
		{
			device_p4_cyclicShift_permutated_11 (&xs[idx], COEFFICIENT_SIZE,
																				permutationStride);
			nShift -= COEFFICIENT_SIZE;
		}

		device_p4_cyclicShift_permutated_11 (&xs[idx], nShift, permutationStride);
	}
}

/**********************************************/

#endif
