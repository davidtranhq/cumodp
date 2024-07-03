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
#ifndef BIG_ARITHMETIC_SUBTRACTION_DEVICE_H_
#define BIG_ARITHMETIC_SUBTRACTION_DEVICE_H_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"

/**********************************************/
// In place subtraction in Z/(R^8 + 1)Z (xs = xs - ys)
__device__ __forceinline__ void
device_p4_bigSub (usfixn64 *xm, usfixn64 *ym, usfixn64 *um)
{
//	usfixn64 um[8];
	short i, pos;
	unsigned short c = 0;
	usfixn64 num1, num2;
	num1 = 0;
	num2 = R - 1;

	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			um[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			um[i] = xm[i] - num1;
		}
	}
	if (c > 0)
	{
		pos = -1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (um[i] <= num2)
			{
				pos = i;
				break;
			}
		}
		if (pos >= 0)
		{
			for (i = 0; i < pos; i++)
			{
				um[i] = 0;
			}
			um[pos]++;
		}
		else
		{
			um[0] = ULMAX;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
				um[i] = 0;
			}
		}
	}
//	xm[0] = um[0];
//	xm[1] = um[1];
//	xm[2] = um[2];
//	xm[3] = um[3];
//	xm[4] = um[4];
//	xm[5] = um[5];
//	xm[6] = um[6];
//	xm[7] = um[7];
}

/**********************************************/
/**********************************************/

__device__ __forceinline__ void
device_p4_bigSub_plain (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym)
{
//	usfixn64 um[8];
	short i, pos;
	unsigned short c = 0;
	usfixn64 num1, num2;
	num1 = 0;
	num2 = R - 1;

	for (i = 0; i < 8; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			xm[i] = xm[i] - num1;
		}
	}
	if (c > 0)
	{
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[i] <= num2)
			{
				pos = i;
				break;
			}
//			printf("pos = %d \n",pos);
		}
		if (pos >= 0)
		{
			for (i = 0; i < pos; i++)
			{
				xm[i] = 0;
			}
			xm[pos]++;
		}
		else
		{
			xm[0] = ULMAX;
			for (i = 1; i < 8; i++)
			{
				xm[i] = 0;
			}
		}
	}
//	xm[0] = um[0];
//	xm[1] = um[1];
//	xm[2] = um[2];
//	xm[3] = um[3];
//	xm[4] = um[4];
//	xm[5] = um[5];
//	xm[6] = um[6];
//	xm[7] = um[7];
}

/*************************************************************/
__device__ inline void
device_p4_bigSub_plain_1 (usfixn64 * __restrict__ xm, const usfixn64 * __restrict
ym)
{
//	usfixn64 um[8];
	short i, posSub = -1;
	unsigned short c = 0;
	usfixn64 num1;
//	, num2;
	num1 = 0;
//	num2 = R - 1;

//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
//			xm[i] = R - num1 + xm[i];
			num1 = R - num1 + xm[i];
		}
		else
		{
			c = 0;
//			xm[i] = xm[i] - num1;
			num1 = xm[i] - num1;
		}
		xm[i] = num1;
		posSub = (posSub == -1 && num1 <= R_MINUS_ONE) * i;
	}

//	if (c > 0)
//	{
//		pos = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[i] < R_MINUS_ONE)
//			{
//				pos = i;
//				break;
//			}
//		}
	if (c > 0 && posSub >= 0)
	{
		for (i = 0; i < posSub; i++)
		{
			xm[i] = 0;
		}
		xm[posSub]++;
	}
	if (c > 0 && posSub < 0)
	{
		xm[0] = ULMAX;
//#pragma unroll (COEFFICIENT_SIZE-1)
		for (i = 1; i < 8; i++)
		{
			xm[i] = 0;
		}
	}
//	xm[0] = um[0];
//	xm[1] = um[1];
//	xm[2] = um[2];
//	xm[3] = um[3];
//	xm[4] = um[4];
//	xm[5] = um[5];
//	xm[6] = um[6];
//	xm[7] = um[7];
}

/*************************************************************/
__device__ __forceinline__ void
device_p4_bigSub_uVector_plain_2 (uvector8 & __restrict__ xm,
															 const uvector8 & __restrict__
															 ym)
{
//	usfixn64 um[8];
	short i, posSub = -1;
	unsigned short c = 0;
	usfixn64 num1;
//	, num2;
	num1 = 0;
//	num2 = R - 1;

//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
//		num1 = ym[i] + c;
		num1 = getUvector8Element (ym, i) + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
		if (getUvector8Element (xm, i) < num1) //there is not enough to do subtraction
		{
			c = 1;
//			xm[i] = R - num1 + xm[i];
//			num1 = R - num1 + xm[i];
			num1 = R - num1 + getUvector8Element (xm, i);
		}
		else
		{
			c = 0;
//			xm[i] = xm[i] - num1;
//			num1 = xm[i] - num1;
			num1 = getUvector8Element (xm, i) - num1;
		}
//		xm[i] = num1;
		setUvector8Element (xm, i, num1);
		posSub = (posSub == -1 && num1 <= R_MINUS_ONE) * i;
	}

	if (c > 0 && posSub >= 0)
	{
		for (i = 0; i < posSub; i++)
		{
//			xm[i] = 0;
			setUvector8Element (xm, i, 0);
		}
//		memset(xm, 0x00, posSub*sizeof(usfixn64));
//		xm[posSub]++;
		setUvector8Element (xm, posSub, getUvector8Element (xm, posSub) + 1);
	}
	if (c > 0 && posSub < 0)
	{
//		xm[0] = ULMAX;
		setUvector8Element (xm, 0, ULMAX);
//#pragma unroll (COEFFICIENT_SIZE-1)
		for (i = 1; i < 8; i++)
		{
//			xm[i] = 0;
			setUvector8Element (xm, i, 0);
		}
	}

}

///*************************************************************/
//__device__ __forceinline__ void  smallSub3(usfixn64 *xm, usfixn64 *ym, usfixn64 *um, short step)
//{
////	usfixn64 um[8];
//	short i, pos;
//	unsigned short c = 0;
//	usfixn64 num1, num2;
//	num1 = 0;
//	num2 = R - 1;
//
//	for (i = 0; i < step; i++)
//	{
//		num1 = ym[i] + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
//		{
//			c = 1;
//			um[i] = R - num1 + xm[i];
//		}
//		else
//		{
//			c = 0;
//			um[i] = xm[i] - num1;
//		}
//	}
//	if (c > 0)
//	{
//		pos = -1;
//		for (i = 0; i < step; i++)
//		{
//			if (um[i] < num2)
//			{
//				pos = i;
//				break;
//			}
//		}
//		if (pos >= 0)
//		{
//			for (i = 0; i < pos; i++)
//			{
//				um[i] = 0;
//			}
//			um[pos]++;
//		}
//		else
//		{
//			um[0] = ULMAX;
//			for (i = 1; i < step; i++)
//			{
//				um[i] = 0;
//			}
//		}
//	}
////	xm[0] = um[0];
////	xm[1] = um[1];
////	xm[2] = um[2];
////	xm[3] = um[3];
////	xm[4] = um[4];
////	xm[5] = um[5];
////	xm[6] = um[6];
////	xm[7] = um[7];
//}

/**********************************************/

__device__ __forceinline__ void
device_p4_bigPrimeSub_permutated (usfixn64 *xm, usfixn64 * ym, usfixn64 *um,
												unsigned short idx, short permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1, num2;

	num1 = 0;
	num2 = R - 1;
	short i;
	usfixn64 pos = threadIdx.x;

	short offset = 0;
	short pos1;
	//	usfixn64 tmp = 0;

	offset = 0;
	for (i = 0; i < 8; i++)
	{
		num1 = ym[idx + offset] + c;
		if (xm[idx + offset] < num1) //there is not enough to do subtraction
		{
			c = 1;
			um[idx + offset] = R - num1 + xm[idx + offset];
		}
		else
		{
			c = 0;
			um[idx + offset] = xm[idx + offset] - num1;
		}
		offset += permutationStride;
	}

	if (c > 0)
	{
		offset = 0;
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (um[idx + offset] <= num2)
			{
				pos = i;
				break;
			}
			offset += permutationStride;
		}

		if (pos >= 0)
		{
			offset = 0;
			for (i = 0; i < pos; i++)
			{
				um[idx + offset] = 0;
				offset += permutationStride;
			}
			offset = pos * permutationStride;
//			um[pos]++;
			um[idx + offset]++;
		}
		else
		{

//			um[0] = ULMAX;
			um[idx] = ULMAX;
			offset = 0;
			for (i = 1; i < 8; i++)
			{
//				um[i] = 0;
				um[idx + offset];
				offset += permutationStride;
			}
		}
	}

//	offset=0;
//	for(i=0; i<8; i++)
//	{
//		um[idx+offset]=;
//		offset+=permutationStride;
//	}

}
/**********************************************/
__device__ __forceinline__ void
device_p4_bigPrimeSub_permutated_2 (usfixn64 *xm, usfixn64 * ym,
																 usfixn64 permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1, num2;

	num1 = 0;
	num2 = R - 1;
	short i;
	usfixn64 pos = threadIdx.x;

//	short offset = 0;
	usfixn64 offset = 0;
	short pos1;
	//	usfixn64 tmp = 0;

	offset = 0;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		num1 = ym[idx + offset] + c;
		num1 = ym[offset] + c;
//		if (xm[idx + offset] < num1) //there is not enough to do subtraction
//		{
//			c = 1;
//			um[idx + offset] = R - num1 + xm[idx + offset];
//		}
//		else
//		{
//			c = 0;
//			um[idx + offset] = xm[idx + offset] - num1;
//		}
		if (xm[offset] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[offset] = R - num1 + xm[offset];
		}
		else
		{
			c = 0;
			xm[offset] = xm[offset] - num1;
		}
		offset += permutationStride;
	}

	if (c > 0)
	{
		offset = 0;
		pos = -1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (xm[offset] <= num2)
			{
				pos = i;
				break;
			}
			offset += permutationStride;
		}

		if (pos >= 0)
		{
			offset = 0;
			for (i = 0; i < pos; i++)
			{
				xm[offset] = 0;
				offset += permutationStride;
			}
			offset = pos * permutationStride;
//			um[pos]++;
			xm[offset]++;
		}
		else
		{

//			um[0] = ULMAX;
			xm[0] = ULMAX;
			offset = permutationStride;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
//				um[i] = 0;
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
	}
}

/**********************************************/
__device__ __inline__ void
device_p4_bigPrimeSub_permutated_2_bound (usfixn64 * __restrict__ xm,
																			 const usfixn64 * __restrict__ ym,
																			 usfixn64 permutationStride,
																			 usfixn64 bound)
{
	unsigned short c = 0;
	usfixn64 num1, num2;

	num1 = 0;
	num2 = R - 1;
	short i;
	usfixn64 pos = threadIdx.x;

//	short offset = 0;
	usfixn64 offset = 0;
	short pos1;
	//	usfixn64 tmp = 0;

	offset = 0;
	for (i = 0; i < bound; i++)
	{
//		num1 = ym[idx + offset] + c;
		num1 = ym[offset] + c;
//		if (xm[idx + offset] < num1) //there is not enough to do subtraction
//		{
//			c = 1;
//			um[idx + offset] = R - num1 + xm[idx + offset];
//		}
//		else
//		{
//			c = 0;
//			um[idx + offset] = xm[idx + offset] - num1;
//		}
		if (xm[offset] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[offset] = R - num1 + xm[offset];
		}
		else
		{
			c = 0;
			xm[offset] = xm[offset] - num1;
		}
		offset += permutationStride;
	}

	for (i = bound; i < COEFFICIENT_SIZE; i++)
	{
		//		num1 = ym[idx + offset] + c;
		num1 = c;
		//		if (xm[idx + offset] < num1) //there is not enough to do subtraction
		//		{
		//			c = 1;
		//			um[idx + offset] = R - num1 + xm[idx + offset];
		//		}
		//		else
		//		{
		//			c = 0;
		//			um[idx + offset] = xm[idx + offset] - num1;
		//		}
		if (xm[offset] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[offset] = R - num1 + xm[offset];
		}
		else
		{
			c = 0;
			xm[offset] = xm[offset] - num1;
		}
		offset += permutationStride;
	}

	if (c > 0)
	{
		offset = 0;
		pos = -1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (xm[offset] <= num2)
			{
				pos = i;
				break;
			}
			offset += permutationStride;
		}

		if (pos >= 0)
		{
			offset = 0;
			for (i = 0; i < pos; i++)
			{
				xm[offset] = 0;
				offset += permutationStride;
			}
			offset = pos * permutationStride;
//			um[pos]++;
			xm[offset]++;
		}
		else
		{

//			um[0] = ULMAX;
			xm[0] = ULMAX;
			offset = permutationStride;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
//				um[i] = 0;
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
	}
}

/***********************************************************/
//
////permutation stride for ym is = 1
////permtuation stride for xm is = permutationStride
//__device__ __forceinline__ void  device_bigPrimeSub_permutated_3(usfixn64 *xm, usfixn64 * ym,
//		usfixn64 permutationStride)
//{
//	unsigned short c = 0;
//	usfixn64 num1, num2;
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
//	usfixn64 pos = threadIdx.x;
//
//	short offset = 0;
//	short pos1;
//	//	usfixn64 tmp = 0;
//
//	offset = 0;
//	for (i = 0; i < 8; i++)
//	{
////		num1 = ym[idx + offset] + c;
//		num1 = ym[i] + c;
////		if (xm[idx + offset] < num1) //there is not enough to do subtraction
////		{
////			c = 1;
////			um[idx + offset] = R - num1 + xm[idx + offset];
////		}
////		else
////		{
////			c = 0;
////			um[idx + offset] = xm[idx + offset] - num1;
////		}
//		if (xm[offset] < num1) //there is not enough to do subtraction
//		{
//			c = 1;
//			xm[offset] = R - num1 + xm[offset];
//		}
//		else
//		{
//			c = 0;
//			xm[offset] = xm[offset] - num1;
//		}
//		offset += permutationStride;
//	}
//
//	if (c > 0)
//	{
//		offset = 0;
//		pos = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[offset] < num2)
//			{
//				pos = i;
//				break;
//			}
//			offset += permutationStride;
//		}
//
//		if (pos >= 0)
//		{
//			offset = 0;
//			for (i = 0; i < pos; i++)
//			{
//				xm[offset] = 0;
//				offset += permutationStride;
//			}
//			offset = pos * permutationStride;
////			um[pos]++;
//			xm[offset]++;
//		}
//		else
//		{
//
////			um[0] = ULMAX;
//			xm[0] = ULMAX;
//			offset = permutationStride;
//			for (i = 1; i < 8; i++)
//			{
////				um[i] = 0;
//				xm[offset] = 0;
//				offset += permutationStride;
//			}
//		}
//	}
//}
/**********************************************/
//permutation stride for xm=1
//permutation stride for ym=permutationStride
__device__ inline void
device_p4_bigPrimeSub_permutated_4 (usfixn64 * __restrict__ xm,
																 usfixn64 * __restrict__ ym,
																 usfixn64 permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1, num2;

	num1 = 0;
	num2 = R - 1;
	short i;
	usfixn64 pos = 0;

	usfixn64 offset = 0;
	usfixn64 pos1;
	//	usfixn64 tmp = 0;

	offset = 0;
	for (i = 0; i < 8; i++)
	{
//		num1 = ym[idx + offset] + c;
		num1 = ym[offset] + c;
//		if (xm[idx + offset] < num1) //there is not enough to do subtraction
//		{
//			c = 1;
//			um[idx + offset] = R - num1 + xm[idx + offset];
//		}
//		else
//		{
//			c = 0;
//			um[idx + offset] = xm[idx + offset] - num1;
//		}
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			xm[i] = xm[i] - num1;
		}
		offset += permutationStride;
	}

	if (c > 0)
	{
//		offset = 0;
		pos = -1;
		for (i = 0; i < 8; i++)
		{
//			if (xm[i] < num2)
			if (xm[i] <= R_MINUS_ONE)
			{
				pos = i;
				break;
			}
//			offset += permutationStride;
		}

		if (pos >= 0)
		{
//			offset = 0;
			for (i = 0; i < pos; i++)
			{
				xm[i] = 0;
//				offset += permutationStride;
			}
//			offset = pos * permutationStride;
//			um[pos]++;
			xm[pos]++;
		}
		else
		{

//			um[0] = ULMAX;
			xm[0] = ULMAX;
//			offset = permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
//				um[i] = 0;
				xm[i] = 0;
//				offset += permutationStride;
			}
		}
	}
}

/**********************************************/
__device__ __inline__ void
device_p4_bigPrimeSub_plain_inPlace (usfixn64 * __restrict__ xm,
																	usfixn64 * __restrict__ ym)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1;

//	num1 = 0;
//	num2 = R - 1;
	short i;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			xm[i] = xm[i] - num1;
		}
	}

	if (c > 0)
	{
		pos1 = -1;
		if (xm[0] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 0;
			xm[0]++;
			return;
		}
		if (xm[1] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 1;
			xm[0] = 0;
			xm[1]++;
			return;
		}
		if (xm[2] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 2;
			xm[0] = 0;
			xm[1] = 0;
			xm[2]++;
			return;
		}
		if (xm[3] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 3;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3]++;
			return;
		}
		if (xm[4] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 4;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4]++;
			return;
		}
		if (xm[5] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 5;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5]++;
			return;
		}
		if (xm[6] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 6;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6]++;
			return;
		}

		if (xm[7] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 7;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7]++;
			return;
		}

		if (xm[8] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 8;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8]++;
			return;
		}

		if (xm[9] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 9;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9]++;
			return;
		}

		if (xm[10] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 10;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10]++;
			return;
		}
		if (xm[11] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 11;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11]++;
			return;
		}

		if (xm[12] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 12;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11] = 0;
			xm[12]++;
			return;
		}

		if (xm[13] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 13;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11] = 0;
			xm[12] = 0;
			xm[13]++;
			return;
		}

		if (xm[14] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 14;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11] = 0;
			xm[12] = 0;
			xm[13] = 0;
			xm[14]++;
			return;
		}

		if (xm[15] <= R_MINUS_ONE && pos1 == -1)
		{
			pos1 = 15;
			xm[0] = 0;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11] = 0;
			xm[12] = 0;
			xm[13] = 0;
			xm[14] = 0;
			xm[15]++;
			return;
		}

//else (c>0) but (pos ==-1)
		{
			xm[0] = ULMAX;
			xm[1] = 0;
			xm[2] = 0;
			xm[3] = 0;
			xm[4] = 0;
			xm[5] = 0;
			xm[6] = 0;
			xm[7] = 0;
			xm[8] = 0;
			xm[9] = 0;
			xm[10] = 0;
			xm[11] = 0;
			xm[12] = 0;
			xm[13] = 0;
			xm[14] = 0;
			xm[15] = 0;
		}
	}
}
/**********************************************/

/**********************************************/
__device__ __inline__ void
device_p4_bigPrimeSub_plain_ptx_v0 (usfixn64 * __restrict__ xm,
																 usfixn64 * __restrict__ ym)
{
	unsigned short c = 0;

//	short pos1;
	usfixn32 pos1;
	usfixn64 num1;

	usfixn32 bitsFlag = 0, bit = 0;

//	num1 = 0;
//	num2 = R - 1;
	short i;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			xm[i] = xm[i] - num1;
		}
		bit = (xm[i] <= R_MINUS_ONE);
		bitsFlag |= (bit << i);
	}

	asm("brev.b32 %0,%0;":"+r"(bitsFlag):);
	asm("bfind.u32 %0, %1;":"=r"(pos1):"r"(bitsFlag));
//	if (pos1 == 0xFFFFFFFF)
//		pos1 = -1;
//	else
//		pos1 = 31 - pos1;

	pos1 = 31 - pos1;
//
//	if (c > 0)
//	{
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//				xm[i] = 0;
//			xm[i]++;
//		}
//		else
//		{
//			xm[0] = ULMAX;
//			for (i = 1; i < 8; i++)
//				xm[i] = 0;
//		}
//	}

	if (c > 0 && pos1 < 32)
	{
		for (i = 0; i < pos1; i++)
			xm[i] = 0;
		xm[i]++;
	}
	if (c > 0 && pos1 == 32)
	{
		xm[0] = ULMAX;
		for (i = 1; i < COEFFICIENT_SIZE; i++)
			xm[i] = 0;
	}
}
/**********************************************/

__device__ __forceinline__ void
device_p4_bigSubZero_3 (usfixn64 &l, usfixn64 & h, usfixn64& carry)
{
//	short i, pos;
	unsigned short c = 0;
	usfixn64 num1;
	num1 = 0;
//	num2 = R - 1;

	usfixn64 s = 0;
	//step0
	num1 = l - 1;
//		if (xm[i] < num1) //there is not enough to do subtraction
	if (0 < num1)
	{
		c = 1;
		s = R - num1;
	}
	else
	{
		c = 0;
		s = 0 - num1;
	}
	l = s;

	num1 = h + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
	if (0 < num1)
	{
		c = 1;
		s = R - num1;
	}
	else
	{
		c = 0;
		s = 0 - num1;
	}
	h = s;

	num1 = carry + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
	if (0 < num1)
	{
		c = 1;
		s = R - num1;
	}
	else
	{
		c = 0;
		s = 0 - num1;
	}
	carry = s;
}
/**********************************************/

#endif
