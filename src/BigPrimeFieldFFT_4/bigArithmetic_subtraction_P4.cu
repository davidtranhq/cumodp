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
#ifndef BIG_ARITHMETIC_SUBTRACTION_CU_
#define BIG_ARITHMETIC_SUBTRACTION_CU_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"

#include "bigArithmetic_subtraction_P4_device.h"

/**********************************************/
__global__ void
kernel_p4_subtraction_plain (usfixn64 * xs, usfixn64 * ys, usfixn64 *parameters)
{
	//
	short operation = parameters[0];
	short iterations = parameters[1];
	short paddingMethod = parameters[2];
	short dynamicMemSize = parameters[3];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	short i;
//	, pos, pos1, t;
	unsigned short c = 0;
//	usfixn64 num1, num2;
//	usfixn64 *xd, *yd, *ud, *xm, *ym, *um;
	usfixn64 *xd, *yd, *xm, *ym;

//	num1 = 0;
//	num2 = R - 1;

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * COEFFICIENT_SIZE);
	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	ud = (usfixn64*) ((char*) us + tid * sizeof(usfixn64) * 8);

//	device_bigSub_plain(xd, yd);
//	device_bigPrimeSub_plain_inPlace(xd,yd);
//	device_p4_bigPrimeSub_permutated_2(xd,yd, 1);
	device_p4_bigPrimeSub_plain_ptx_v0 (xd, yd);
	//second solution, going over every 2 digits, checking the carry at the end
//	i = 0;
//	device_smallSub2_plain(&xd[i], &yd[i], c);
//	i += 2;
//	device_smallSub2_plain(&xd[i], &yd[i], c);
//	i += 2;
//	device_smallSub2_plain(&xd[i], &yd[i], c);
//	i += 2;
//	device_smallSub2_plain(&xd[i], &yd[i], c);
//	bigPrimeSub_check(xd,c);
}

/**********************************************/
__global__ void
kernel_p4_subtraction_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
{
	short operation = parameters[0];
	usfixn64 permutationStride = parameters[5];
	short shuffle = parameters[6];
	short padding = 0;	//parameters[?]
	usfixn64 idx;
	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);

	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
	//following indexing is slightly faster than above indexing

	idx = tid;
	if (padding == 0)
		device_p4_bigPrimeSub_permutated_2 (&xs[idx], &ys[idx], permutationStride);
}

/**********************************************/

#endif
