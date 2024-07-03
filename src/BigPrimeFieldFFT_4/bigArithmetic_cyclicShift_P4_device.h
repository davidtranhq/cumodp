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
#ifndef BIG_ARITHMETIC_CYCLIC_SHIFT_DEVICE_H_
#define BIG_ARITHMETIC_CYCLIC_SHIFT_DEVICE_H_

#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"
//#include "BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.h"
#include "bigArithmetic_subtraction_P4_device.h"

//__device__ __forceinline__ void cyclicShift(usfixn64 *xs, short sn);
//__device__ __forceinline__ void cyclicShift(usfixn64 *xs, short sn);

#include <stdio.h>

/**********************************************/
//__device__ __inline__ void
//device_oneShiftRight16 (usfixn64 * xs)
//{
//	usfixn64 tmp;
//	short i;
//
////	for (j = 0; j < shiftNo; j++)
//	{
//		tmp = xs[15];
////#pragma unroll 7
////		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
////			xs[i] = xs[i - 1];
//		xs[15] = xs[14];
//		xs[14] = xs[13];
//		xs[13] = xs[12];
//		xs[12] = xs[11];
//		xs[11] = xs[10];
//		xs[10] = xs[9];
//		xs[9] = xs[8];
//		xs[8] = xs[7];
//		xs[7] = xs[6];
//		xs[6] = xs[5];
//		xs[5] = xs[4];
//		xs[4] = xs[3];
//		xs[3] = xs[2];
//		xs[2] = xs[1];
//		xs[1] = xs[0];
////		xs[0] = device_negate_plain(tmp);
//		xs[0] = tmp;
//	}
//}
/**********************************************/
__device__ __forceinline__ void
device_p4_oneShiftRight_uvector16 (uvector16 &xs)
{
	usfixn64 tmp;
	short i;

//	for (j = 0; j < shiftNo; j++)
	{
		tmp = xs.i15;
//#pragma unroll 7
//		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
//			xs[i] = xs[i - 1];
		xs.i15 = xs.i14;
		xs.i14 = xs.i13;
		xs.i13 = xs.i12;
		xs.i12 = xs.i11;
		xs.i11 = xs.i10;
		xs.i10 = xs.i9;
		xs.i9 = xs.i8;
		xs.i8 = xs.i7;
		xs.i7 = xs.i6;
		xs.i6 = xs.i5;
		xs.i5 = xs.i4;
		xs.i4 = xs.i3;
		xs.i3 = xs.i2;
		xs.i2 = xs.i1;
		xs.i1 = xs.i0;
//		xs[0] = device_negate_plain(tmp);
		xs.i0 = tmp;
	}
}

/**********************************************/
__device__ __forceinline__ void
device_p4_oneShiftLeft_uvector16 (uvector16& xs)
{
	usfixn64 tmp;
	short i;

//	for (j = 0; j < shiftNo; j++)
	{
		tmp = xs.i0;
//#pragma unroll 7
//		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
//			xs[i] = xs[i - 1];
		xs.i0 = xs.i1;
		xs.i1 = xs.i2;
		xs.i2 = xs.i3;
		xs.i3 = xs.i4;
		xs.i4 = xs.i5;
		xs.i5 = xs.i6;
		xs.i6 = xs.i7;
		xs.i7 = xs.i8;
		xs.i8 = xs.i9;
		xs.i9 = xs.i10;
		xs.i10 = xs.i11;
		xs.i11 = xs.i12;
		xs.i12 = xs.i13;
		xs.i13 = xs.i14;
		xs.i14 = xs.i15;
//		xs[0] = device_negate_plain(tmp);
		xs.i15 = tmp;
	}
}

/**********************************************/
__device__ __forceinline__ void
device_p4_oneShiftLeft (usfixn64 * xs)
{
	usfixn64 tmp;
//	short i;

//	for (j = 0; j < shiftNo; j++)
	{
		tmp = xs[0];
//#pragma unroll 7
//		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
//			xs[i] = xs[i - 1];
		xs[0] = xs[1];
		xs[1] = xs[2];
		xs[2] = xs[3];
		xs[3] = xs[4];
		xs[4] = xs[5];
		xs[5] = xs[6];
		xs[6] = xs[7];

		xs[7] = xs[8];
		xs[8] = xs[9];
		xs[9] = xs[10];
		xs[10] = xs[11];
		xs[11] = xs[12];
		xs[12] = xs[13];
		xs[13] = xs[14];
		xs[14] = xs[15];
//		xs[0] = device_negate_plain(tmp);
		xs[15] = tmp;
	}
}

/**********************************************/
__device__ __forceinline__ void
device_p4_cyclicShift (usfixn64 *xs, short sn)
{
	short i = 0, j = 0;
	usfixn64 ts[COEFFICIENT_SIZE] =
		{ 0 };
	if (sn <= 0)
	{
		return;
	}
	if (sn > COEFFICIENT_SIZE)
	{
		return;
	}
	j = COEFFICIENT_SIZE - sn;
	for (i = 0; i < sn; i++)
	{
		ts[i] = xs[j++];
	}
	for (i = COEFFICIENT_SIZE - 1 - sn; i >= 0; i--)
	{
		xs[i + sn] = xs[i];
	}
	for (i = 0; i < sn; i++)
	{
		xs[i] = 0;
	}
	device_p4_bigSub (xs, ts, xs);
}

/**********************************************/
__device__ __forceinline__ void
device_p4_cyclicShift_plain (usfixn64 *xs, short sn)
{
	short i = 0, j = 0;
	usfixn64 ts[COEFFICIENT_SIZE] =
		{ 0 };
	if (sn <= 0)
	{
		return;
	}
	if (sn > COEFFICIENT_SIZE)
	{
		return;
	}
	j = COEFFICIENT_SIZE - sn;
	for (i = 0; i < sn; i++)
	{
		ts[i] = xs[j++];
	}
	for (i = COEFFICIENT_SIZE - 1 - sn; i >= 0; i--)
	{
		xs[i + sn] = xs[i];
	}
	for (i = 0; i < sn; i++)
	{
		xs[i] = 0;
	}
	device_p4_bigSub (xs, ts, xs);
}

/**********************************************/

__device__ __inline__ void
device_p4_negate_1 (usfixn64 * __restrict__ xs, const usfixn64 * __restrict__ ys,
								 short * __restrict__ c)
{
	usfixn64 num1;
//	num1 = getUvector8Element(ym,i) + c;
	num1 = *ys + *c;
	//		if (xm[i] < num1) //there is not enough to do subtraction
//			if (getUvector8Element(xm, i) < num1) //there is not enough to do subtraction
	if (*xs < num1) //there is not enough to do subtraction
	{
		*c = 1;
		//			xm[i] = R - num1 + xm[i];
		//			num1 = R - num1 + xm[i];
//				num1 = R - num1 + getUvector8Element(xm,i);
		num1 = R - num1 + *xs;
	}
	else
	{
		*c = 0;
		//			xm[i] = xm[i] - num1;
		//			num1 = xm[i] - num1;
//				num1 = getUvector8Element(xm,i) - num1;
		num1 = *xs - num1;
	}
	//		xm[i] = num1;
//			setUvector8Element(xm, i, num1);
	*xs = num1;
}
/***********************************/
__device__ __inline__ void
device_p4_negate_2 (usfixn64 &xs, const usfixn64 &ys, short & c)
{
	usfixn64 num1;
//	num1 = getUvector8Element(ym,i) + c;
	num1 = ys + c;
	//		if (xm[i] < num1) //there is not enough to do subtraction
//			if (getUvector8Element(xm, i) < num1) //there is not enough to do subtraction
	if (xs < num1) //there is not enough to do subtraction
	{
		c = 1;
		//			xm[i] = R - num1 + xm[i];
		//			num1 = R - num1 + xm[i];
//				num1 = R - num1 + getUvector8Element(xm,i);
//		num1 = R - num1 + xs;
		xs = R - num1 + xs;
		return;
	}
//	else
	{
		c = 0;
		//			xm[i] = xm[i] - num1;
		//			num1 = xm[i] - num1;
//				num1 = getUvector8Element(xm,i) - num1;
//		num1 = xs - num1;
		xs -= num1;
		return;
	}
	//		xm[i] = num1;
//			setUvector8Element(xm, i, num1);
//	*xs = num1;
}
/***********************************/
//using uvector8 on register which makes it register intensive (a lot, up to 63)
//this device_cyclicShift_twice_permutated_1 computes two cyclic shifts if sn>8 and one if sn<=8
//this is cyclicShift over 8 digits, not 16
//e.g. sn=12; first does cyclicShift(sn) then computes cyclicShift(sn-8)
__device__ __inline__ void
device_p4_cyclicShift_twice_permutated_1 (usfixn64 * __restrict__ xs,
																			 const short & snInput,
																			 const usfixn64 & permutationStride)
{
//	if (sn <= 0 || sn > 8)
	short sn, sn2 = 0;
	sn = snInput;
	if (sn <= 0)
	{
		return;
	}
//	short i = 0;
	short j = 0;
	short c = 0;
	short pos = -1;

	if (sn > 8)
	{
		sn2 = sn - 8;
		sn = 8;
	}
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
//	uConstArray8_align8 ts, ys;
//	uvector8 ts, ys;
	uvector8 ys;

	usfixn64 offset = 0;
	usfixn64 tmpX = 0;
	usfixn64 tmpY = 0;

	ys.i0 = xs[offset];
	offset += permutationStride;
	ys.i1 = xs[offset];
	offset += permutationStride;
	ys.i2 = xs[offset];
	offset += permutationStride;
	ys.i3 = xs[offset];
	offset += permutationStride;
	ys.i4 = xs[offset];
	offset += permutationStride;
	ys.i5 = xs[offset];
	offset += permutationStride;
	ys.i6 = xs[offset];
	offset += permutationStride;
	ys.i7 = xs[offset];

//	for(j=0;j<8;j++)
//	{
//		setUvector8Element(ys,j, xs[offset]);
//		offset+=permutationStride;
//	}

	//######################################################
	//doing the subtraction
	for (j = 8 - sn; j < 8; j++)
	{
		tmpX = 0;
		tmpY = getUvector8Element (ys, j);
		device_p4_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

	tmpY = 0;
	for (j = 0; j < 8 - sn; j++)
	{
		tmpX = getUvector8Element (ys, j);
		device_p4_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

//######################################################

	if (sn <= 8)
	{
		//######################################################
		//Finding pos
		for (j = 8 - sn; j < 8; j++)
		{
			if (getUvector8Element (ys, j) < R_MINUS_ONE)
			{
				pos = j;
				break;
			}
		}

		if (pos == -1)
		{
			for (j = 0; j < 8 - sn; j++)
			{
				if (getUvector8Element (ys, j) < R_MINUS_ONE)
				{
					pos = j;
					break;
				}
			}
		}

		//######################################################
		//Making changes based on value of pos
		if (c > 0 && pos > -1)
		{
			if (pos >= (8 - sn))	//in first group
			{
				for (j = 8 - sn; j < pos; j++)
				{
					setUvector8Element (ys, j, 0);
				}
				incUvector8Element (ys, pos, 1);
			}
			else
			{
				for (j = 8 - sn; j < 8; j++)
				{
					setUvector8Element (ys, j, 0);
				}
				for (j = 0; j < pos; j++)
				{
					setUvector8Element (ys, j, 0);
				}
				incUvector8Element (ys, pos, 1);
			}
		}
		if (c > 0 && pos == -1)
		{
			setUvector8Element (ys, 8 - sn, ULMAX);
			for (j = 8 - sn + 1; j < 8; j++)
				setUvector8Element (ys, j, 0);
			for (j = 0; j < 8 - sn; j++)
				setUvector8Element (ys, j, 0);
		}
	}

	if (sn2 <= 0)
	{
		// copying results back to global mem
		offset = 0;

		for (j = 8 - sn; j < 8; j++)
		{
			xs[offset] = getUvector8Element (ys, j);
			offset += permutationStride;
		}

		for (j = 0; j < 8 - sn; j++)
		{
			xs[offset] = getUvector8Element (ys, j);
			offset += permutationStride;
		}

		return;
	}
	pos = -1;
	c = 0;
	//######################################################
//the same for for sn2
	//######################################################
	//doing the subtraction
	for (j = 8 - sn2; j < 8; j++)
	{
		tmpX = 0;
		tmpY = getUvector8Element (ys, j);
		device_p4_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

	tmpY = 0;
	for (j = 0; j < 8 - sn2; j++)
	{
		tmpX = getUvector8Element (ys, j);
		device_p4_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

//######################################################

	//######################################################
	//Finding pos
	for (j = 8 - sn2; j < 8; j++)
	{
		if (getUvector8Element (ys, j) < R_MINUS_ONE)
		{
			pos = j;
			break;
		}
	}

	if (pos == -1)
	{
		for (j = 0; j < 8 - sn2; j++)
		{
			if (getUvector8Element (ys, j) < R_MINUS_ONE)
			{
				pos = j;
				break;
			}
		}
	}

	//######################################################
	//Making changes based on value of pos
	if (c > 0 && pos > -1)
	{
		if (pos >= (8 - sn2))	//in first group
		{
			for (j = 8 - sn2; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
			}
			incUvector8Element (ys, pos, 1);
		}
		else
		{
			for (j = 8 - sn2; j < 8; j++)
			{
				setUvector8Element (ys, j, 0);
			}
			for (j = 0; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
			}
			incUvector8Element (ys, pos, 1);
		}
	}
	if (c > 0 && pos == -1)
	{
		setUvector8Element (ys, 8 - sn2, ULMAX);
		for (j = 8 - sn2 + 1; j < 8; j++)
			setUvector8Element (ys, j, 0);
		for (j = 0; j < 8 - sn2; j++)
			setUvector8Element (ys, j, 0);
	}
	//######################################################
	// copying results back to global mem
	offset = 0;

	for (j = 8 - sn2; j < 8; j++)
	{
		xs[offset] = getUvector8Element (ys, j);
		offset += permutationStride;
	}

	for (j = 0; j < 8 - sn2; j++)
	{
		xs[offset] = getUvector8Element (ys, j);
		offset += permutationStride;
	}
}
/**********************************************/
/**********************************************/
//using uConstArray on lmem
__device__ __inline__ void
device_p4_cyclicShift_permutated_11 (usfixn64 * xs, const short & sn,
																	const usfixn64 & permutationStride)
{
	if (sn <= 0 || sn > COEFFICIENT_SIZE)
	{
		return;
	}
//	short i = 0;
	short j = 0;
	short c = 0;
	usfixn32 pos = -1;
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
//	uConstArray8_align8 ys;
	uvector16 ys;
	usfixn64 offset = 0;
	usfixn64 tmpX = 0;
	usfixn64 tmpY = 0;

	usfixn32 bitFlag = 0;
	ys.i0 = xs[offset];
	offset += permutationStride;
	ys.i1 = xs[offset];
	offset += permutationStride;
	ys.i2 = xs[offset];
	offset += permutationStride;
	ys.i3 = xs[offset];
	offset += permutationStride;
	ys.i4 = xs[offset];
	offset += permutationStride;
	ys.i5 = xs[offset];
	offset += permutationStride;
	ys.i6 = xs[offset];
	offset += permutationStride;
	ys.i7 = xs[offset];
	offset += permutationStride;
	ys.i8 = xs[offset];
	offset += permutationStride;
	ys.i9 = xs[offset];
	offset += permutationStride;
	ys.i10 = xs[offset];
	offset += permutationStride;
	ys.i11 = xs[offset];
	offset += permutationStride;
	ys.i12 = xs[offset];
	offset += permutationStride;
	ys.i13 = xs[offset];
	offset += permutationStride;
	ys.i14 = xs[offset];
	offset += permutationStride;
	ys.i15 = xs[offset];

//#pragma unroll (COEFFICIENT_SIZE)
//	for (j = 0; j < 8; j++)
//	{
//		setUvector8Element(ys,j, xs[offset]);
////		ys.i[j] = xs[offset];
//		offset += permutationStride;
//	}

	//######################################################
	//doing the subtraction
//	for (j = 8 - sn; j < 8; j++)
//	{
//		tmpX = 0;
////		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
////		device_p4_negate_1(&tmpX, &tmpY, &c);
//		device_p4_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//	}

//	for (j = 8 - sn; j < 8; j++)
//		{
//			tmpX = 0;
//	//		tmpY = getUvector8Element(ys, j);
//			tmpY = ys.i[j];
//	//		device_p4_negate_1(&tmpX, &tmpY, &c);
//			device_p4_negate_2(tmpX, tmpY, c);
//	//		setUvector8Element(ys, j, tmpX);
//			ys.i[j] = tmpX;
//			bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		}

	for (j = 0; j < sn; j++)
		device_p4_oneShiftRight_uvector16 (ys);

	j = 0;
//	j = 8 - sn;
//	if (j < 8)
	if (sn > 0)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i0;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	if (j < 8)
	if (sn > 1)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i1;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}
//	if (j < 8)
	if (sn > 2)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i2;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 3)
//	if (j < 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i3;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}
//	if (j < 8)
	if (sn > 4)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i4;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 5)
//	if (j < 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i5;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	if (j < 8)
	if (sn > 6)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i6;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 7)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i7;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i8;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i8 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 9)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i9;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i9 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 10)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i10;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i10 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 11)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i11;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i11 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 12)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i12;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i12 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 13)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i13;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i13 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 14)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i14;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i14 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 15)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i15;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i15 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	tmpY = 0;
//	for (j = 0; j < 8 - sn; j++)
//	{
////		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
////		device_p4_negate_1(&tmpX, &tmpY, &c);
//		device_p4_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
//	}

//	for(j=0;j<sn;j++)
//	device_oneShiftLeft(ys.i);

	for (j = 0; j < sn; j++)
		device_p4_oneShiftLeft_uvector16 (ys);

	tmpY = 0;
	j = 0;
//		for (j = 0; j < 8 - sn; j++)
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i0;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 1)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i1;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 2)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i2;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 3)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i3;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 4)
//	if (j < 8 - sn)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i4;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		j++;
	}
//	if (j < 8 - sn)

	if (sn < COEFFICIENT_SIZE - 5)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i5;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
//	if (sn < 2)
	if (sn < COEFFICIENT_SIZE - 6)
	{
//				tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i6;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
//	if (sn < 1)
	if (sn < COEFFICIENT_SIZE - 7)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i7;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 8)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i8;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i8 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 9)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i9;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i9 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 10)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i10;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i10 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 11)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i11;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i11 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 12)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i12;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i12 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 13)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i13;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i13 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 14)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i14;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i14 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 15)
	{
		//		tmpX = getUvector8Element(ys, j);
		//		tmpX = ys.i[j];
		tmpX = ys.i15;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		//		ys.i[j] = tmpX;
		ys.i15 = tmpX;
		j++;
	}

//######################################################

	//######################################################
	//Finding pos
//	for (j = 8 - sn; j < 8; j++)
//	{
////		if (getUvector8Element(ys, j) < R_MINUS_ONE)
//		if (ys.i[j] < R_MINUS_ONE)
//		{
//			pos = j;
//			break;
//		}
//	}
//
//	if (pos == -1)
//	{
//		for (j = 0; j < 8 - sn; j++)
//		{
////			if (getUvector8Element(ys, j) < R_MINUS_ONE)
//			if (ys.i[j] < R_MINUS_ONE)
//			{
//				pos = j;
//				break;
//			}
//		}
//	}

	asm("brev.b32 %0,%0;":"+r"(bitFlag):);
	asm("bfind.u32 %0, %1;":"=r"(pos):"r"(bitFlag));
	if (pos == 0xFFFFFFFF)
		pos = -1;
	else
		pos = 31 - pos;

//		pos1 = 31-pos1;

	//######################################################
	//Making changes based on value of pos
	if (c > 0 && pos > -1)
	{
		if (pos >= (COEFFICIENT_SIZE - sn))	//in first group
		{

			for (j = COEFFICIENT_SIZE - sn; j < pos; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;

			}
			incUvector16Element (ys, pos, 1);
//			ys.i[pos]++;

		}
		else
		{
			for (j = COEFFICIENT_SIZE - sn; j < COEFFICIENT_SIZE; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			for (j = 0; j < pos; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			incUvector16Element (ys, pos, 1);
//			ys.i[pos]++;
		}
	}

//	if (c > 0 && pos == -1)
//	{
//		setUvector8Element(ys, 8 - sn, ULMAX);
////		ys.i[8 - sn] = ULMAX;
//		for (j = 8 - sn + 1; j < 8; j++)
////			ys.i[j] = 0;
//			setUvector8Element(ys, j, 0);
//		for (j = 0; j < 8 - sn; j++)
////			ys.i[j] = 0;
//			setUvector8Element(ys, j, 0);
//	}

	for (j = 0; j < sn; j++)
		device_p4_oneShiftRight_uvector16 (ys);

	if (c > 0 && pos == -1)
	{
//			setUvector8Element(ys, 8 - sn, ULMAX);
//	//		ys.i[8 - sn] = ULMAX;
//			for (j = 8 - sn + 1; j < 8; j++)
//	//			ys.i[j] = 0;
//				setUvector8Element(ys, j, 0);
//			for (j = 0; j < 8 - sn; j++)
//	//			ys.i[j] = 0;
//				setUvector8Element(ys, j, 0);
		setUvector16Element (ys, 0, ULMAX);
		ys.i1 = 0x0;
		ys.i2 = 0x0;
		ys.i3 = 0x0;
		ys.i4 = 0x0;
		ys.i5 = 0x0;
		ys.i6 = 0x0;
		ys.i7 = 0x0;
		ys.i8 = 0x0;
		ys.i9 = 0x0;
		ys.i10 = 0x0;
		ys.i11 = 0x0;
		ys.i12 = 0x0;
		ys.i13 = 0x0;
		ys.i14 = 0x0;
		ys.i15 = 0x0;
	}

	//######################################################
	// copying results back to global mem
	offset = 0;

//	for(j=0;j<sn;j++)
//		device_oneShiftRight(ys.i);
//	for (j = 0; j < sn; j++)
//		device_oneShiftRight_uvector8(ys);

//	for (j = 8 - sn; j < 8; j++)
//	{
////		xs[offset] = getUvector8Element(ys, j);
//		xs[offset] = ys.i[j];
//		offset += permutationStride;
//	}

//	for (j = 0; j < 8 - sn; j++)
//	for (j = 0; j < 8; j++)
//	{
//		xs[offset] = getUvector8Element(ys, j);
////		xs[offset] = ys.i[j];
//		offset += permutationStride;
//	}
	if (sn == COEFFICIENT_SIZE)
		ys.i0++;
	{
		xs[offset] = ys.i0;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i1;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i2;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i3;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i4;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i5;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i6;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i7;
		offset += permutationStride;

		xs[offset] = ys.i8;
		offset += permutationStride;

		xs[offset] = ys.i9;
		offset += permutationStride;

		xs[offset] = ys.i10;
		offset += permutationStride;

		xs[offset] = ys.i11;
		offset += permutationStride;

		xs[offset] = ys.i12;
		offset += permutationStride;

		xs[offset] = ys.i13;
		offset += permutationStride;

		xs[offset] = ys.i14;
		offset += permutationStride;

		xs[offset] = ys.i15;
		offset += permutationStride;
	}
}

/**********************************************/
__device__ __inline__ void
device_p4_cyclicShift_permutated_12 (usfixn64 * xs)
{

	const short sn = 1;
	const usfixn64 permutationStride = 1;
	if (sn <= 0 || sn > COEFFICIENT_SIZE)
	{
		return;
	}
//	short i = 0;
	short j = 0;
	short c = 0;
	usfixn32 pos = -1;
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
//	uConstArray8_align8 ys;
	uvector16 ys;
	usfixn64 offset = 0;
	usfixn64 tmpX = 0;
	usfixn64 tmpY = 0;

	usfixn32 bitFlag = 0;
	ys.i0 = xs[offset];
	offset += permutationStride;
	ys.i1 = xs[offset];
	offset += permutationStride;
	ys.i2 = xs[offset];
	offset += permutationStride;
	ys.i3 = xs[offset];
	offset += permutationStride;
	ys.i4 = xs[offset];
	offset += permutationStride;
	ys.i5 = xs[offset];
	offset += permutationStride;
	ys.i6 = xs[offset];
	offset += permutationStride;
	ys.i7 = xs[offset];
	offset += permutationStride;
	ys.i8 = xs[offset];
	offset += permutationStride;
	ys.i9 = xs[offset];
	offset += permutationStride;
	ys.i10 = xs[offset];
	offset += permutationStride;
	ys.i11 = xs[offset];
	offset += permutationStride;
	ys.i12 = xs[offset];
	offset += permutationStride;
	ys.i13 = xs[offset];
	offset += permutationStride;
	ys.i14 = xs[offset];
	offset += permutationStride;
	ys.i15 = xs[offset];

//#pragma unroll (COEFFICIENT_SIZE)
//	for (j = 0; j < 8; j++)
//	{
//		setUvector8Element(ys,j, xs[offset]);
////		ys.i[j] = xs[offset];
//		offset += permutationStride;
//	}

	//######################################################
	//doing the subtraction
//	for (j = 8 - sn; j < 8; j++)
//	{
//		tmpX = 0;
////		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
////		device_p4_negate_1(&tmpX, &tmpY, &c);
//		device_p4_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//	}

//	for (j = 8 - sn; j < 8; j++)
//		{
//			tmpX = 0;
//	//		tmpY = getUvector8Element(ys, j);
//			tmpY = ys.i[j];
//	//		device_p4_negate_1(&tmpX, &tmpY, &c);
//			device_p4_negate_2(tmpX, tmpY, c);
//	//		setUvector8Element(ys, j, tmpX);
//			ys.i[j] = tmpX;
//			bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		}

	for (j = 0; j < sn; j++)
		device_p4_oneShiftRight_uvector16 (ys);

	j = 0;
//	j = 8 - sn;
//	if (j < 8)
	if (sn > 0)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i0;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	if (j < 8)
	if (sn > 1)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i1;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}
//	if (j < 8)
	if (sn > 2)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i2;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 3)
//	if (j < 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i3;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}
//	if (j < 8)
	if (sn > 4)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i4;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 5)
//	if (j < 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i5;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	if (j < 8)
	if (sn > 6)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i6;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 7)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i7;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 8)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i8;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i8 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 9)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i9;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i9 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 10)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i10;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i10 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 11)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i11;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i11 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 12)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i12;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i12 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 13)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i13;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i13 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

	if (sn > 14)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
		//		tmpY = ys.i[j];
		tmpY = ys.i14;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		//		ys.i[j] = tmpX;
		ys.i14 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	tmpY = 0;
//	for (j = 0; j < 8 - sn; j++)
//	{
////		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
////		device_p4_negate_1(&tmpX, &tmpY, &c);
//		device_p4_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
//	}

//	for(j=0;j<sn;j++)
//	device_oneShiftLeft(ys.i);

	for (j = 0; j < sn; j++)
		device_p4_oneShiftLeft_uvector16 (ys);

	tmpY = 0;
	j = 0;
//		for (j = 0; j < 8 - sn; j++)
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i0;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 1)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i1;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 2)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i2;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < COEFFICIENT_SIZE - 3)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i3;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 4)
//	if (j < 8 - sn)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i4;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		j++;
	}
//	if (j < 8 - sn)

	if (sn < COEFFICIENT_SIZE - 5)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i5;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
//	if (sn < 2)
	if (sn < COEFFICIENT_SIZE - 6)
	{
//				tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i6;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
//	if (sn < 1)
	if (sn < COEFFICIENT_SIZE - 7)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i7;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 8)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i8;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i8 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 9)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i9;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i9 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 10)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i10;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i10 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 11)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i11;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i11 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 12)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i12;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i12 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 13)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i13;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i13 = tmpX;
		j++;
	}

	if (sn < COEFFICIENT_SIZE - 14)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i14;
		//		device_p4_negate_1(&tmpX, &tmpY, &c);
		device_p4_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i14 = tmpX;
		j++;
	}

//######################################################

	//######################################################
	//Finding pos
//	for (j = 8 - sn; j < 8; j++)
//	{
////		if (getUvector8Element(ys, j) < R_MINUS_ONE)
//		if (ys.i[j] < R_MINUS_ONE)
//		{
//			pos = j;
//			break;
//		}
//	}
//
//	if (pos == -1)
//	{
//		for (j = 0; j < 8 - sn; j++)
//		{
////			if (getUvector8Element(ys, j) < R_MINUS_ONE)
//			if (ys.i[j] < R_MINUS_ONE)
//			{
//				pos = j;
//				break;
//			}
//		}
//	}

	asm("brev.b32 %0,%0;":"+r"(bitFlag):);
	asm("bfind.u32 %0, %1;":"=r"(pos):"r"(bitFlag));
	if (pos == 0xFFFFFFFF)
		pos = -1;
	else
		pos = 31 - pos;

//		pos1 = 31-pos1;

	//######################################################
	//Making changes based on value of pos
	if (c > 0 && pos > -1)
	{
		if (pos >= (COEFFICIENT_SIZE - sn))	//in first group
		{

			for (j = COEFFICIENT_SIZE - sn; j < pos; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;

			}
			incUvector16Element (ys, pos, 1);
//			ys.i[pos]++;

		}
		else
		{
			for (j = COEFFICIENT_SIZE - sn; j < COEFFICIENT_SIZE; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			for (j = 0; j < pos; j++)
			{
				setUvector16Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			incUvector16Element (ys, pos, 1);
//			ys.i[pos]++;
		}
	}

//	if (c > 0 && pos == -1)
//	{
//		setUvector8Element(ys, 8 - sn, ULMAX);
////		ys.i[8 - sn] = ULMAX;
//		for (j = 8 - sn + 1; j < 8; j++)
////			ys.i[j] = 0;
//			setUvector8Element(ys, j, 0);
//		for (j = 0; j < 8 - sn; j++)
////			ys.i[j] = 0;
//			setUvector8Element(ys, j, 0);
//	}

	for (j = 0; j < sn; j++)
		device_p4_oneShiftRight_uvector16 (ys);

	if (c > 0 && pos == -1)
	{
//			setUvector8Element(ys, 8 - sn, ULMAX);
//	//		ys.i[8 - sn] = ULMAX;
//			for (j = 8 - sn + 1; j < 8; j++)
//	//			ys.i[j] = 0;
//				setUvector8Element(ys, j, 0);
//			for (j = 0; j < 8 - sn; j++)
//	//			ys.i[j] = 0;
//				setUvector8Element(ys, j, 0);
		setUvector16Element (ys, 0, ULMAX);
		ys.i1 = 0x0;
		ys.i2 = 0x0;
		ys.i3 = 0x0;
		ys.i4 = 0x0;
		ys.i5 = 0x0;
		ys.i6 = 0x0;
		ys.i7 = 0x0;
		ys.i8 = 0x0;
		ys.i9 = 0x0;
		ys.i10 = 0x0;
		ys.i11 = 0x0;
		ys.i12 = 0x0;
		ys.i13 = 0x0;
		ys.i14 = 0x0;
		ys.i15 = 0x0;
	}

	//######################################################
	// copying results back to global mem
	offset = 0;

//	for(j=0;j<sn;j++)
//		device_oneShiftRight(ys.i);
//	for (j = 0; j < sn; j++)
//		device_oneShiftRight_uvector8(ys);

//	for (j = 8 - sn; j < 8; j++)
//	{
////		xs[offset] = getUvector8Element(ys, j);
//		xs[offset] = ys.i[j];
//		offset += permutationStride;
//	}

//	for (j = 0; j < 8 - sn; j++)
//	for (j = 0; j < 8; j++)
//	{
//		xs[offset] = getUvector8Element(ys, j);
////		xs[offset] = ys.i[j];
//		offset += permutationStride;
//	}
	{
		xs[offset] = ys.i0;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i1;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i2;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i3;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i4;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i5;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i6;
		offset += permutationStride;
//		__syncthreads();
		xs[offset] = ys.i7;
		offset += permutationStride;

		xs[offset] = ys.i8;
		offset += permutationStride;

		xs[offset] = ys.i9;
		offset += permutationStride;

		xs[offset] = ys.i10;
		offset += permutationStride;

		xs[offset] = ys.i11;
		offset += permutationStride;

		xs[offset] = ys.i12;
		offset += permutationStride;

		xs[offset] = ys.i13;
		offset += permutationStride;

		xs[offset] = ys.i14;
		offset += permutationStride;
		xs[offset] = ys.i15;
		offset += permutationStride;
	}
}

/**********************************************/
__device__ __inline__ void
device_p4_cyclicShift_lhc_complement_v0 (usfixn64 * __restrict__ xs, short i)
{
	short j = 0;
	xs[0] = 1;
	if ((i + 3) < COEFFICIENT_SIZE)
		xs[i + 3] = R - 2;
	else
		xs[(i + 3) & 0xF]++; //%COEFFICIENT_SIZE

	for (j = 1; j < i; j++)
		xs[j] = 0;

//	for(j=i;j<COEFFICIENT_SIZE;j++)
//		xs[j]=R-1;
	if (i + 3 < COEFFICIENT_SIZE)
		for (j = i; j < i + 3; j++)
			xs[j] = R - 1;
	else
		for (j = i; j < COEFFICIENT_SIZE; j++)
			xs[j] = R - 1;
}
/**********************************************/
__device__ __inline__ void
device_p4_cyclicShift_lhc_complement (usfixn64 * xs, short i)
{
//	if (i==0)
//	{
//		xs[0]=0;
//		xs[1]=1;
//		xs[2]
//	}
	short j = 0;
	xs[0] = 1;
	if ((i + 3) < COEFFICIENT_SIZE)
		xs[i + 3] = R - 2;
	else
		xs[(i + 3) & 0xF]++; //%COEFFICIENT_SIZE

	for (j = 1; j < i; j++)
		xs[j] = 0;

//	for(j=i;j<COEFFICIENT_SIZE;j++)
//		xs[j]=R-1;

	xs[i] = R - 1;
	if (i + 1 < COEFFICIENT_SIZE)
		xs[i + 1] = R - 1;
	if (i + 2 < COEFFICIENT_SIZE)
		xs[i + 2] = R - 1;
//	if (i + 3 < COEFFICIENT_SIZE)
//		xs[i + 3] = R - 1;

}
/**********************************************/
__device__ __inline__ usfixn64
device_p4_one_step_addition (usfixn64 x, usfixn64 y, usfixn64 &c)
{
	usfixn64 num1;
	num1 = x + y + c;
	if (num1 < x || num1 < y) //there is overflow/truncation
	{
		num1 = num1 + RC;
		c = 1;
	}
	else if (num1 >= R)
	{
		c = 1;
		num1 = num1 - R;
	}
	else
	{
		c = 0;
	}
	return num1;
}
/**********************************************/
__device__ __inline__ void
device_p4_cyclicShift_lhc_complement_and_add (usfixn64 * xm, short sn)
{
//	if (i==0)
//	{
//		xs[0]=0;
//		xs[1]=1;
//		xs[2]
//	}

	usfixn64 num1 = 0, c = 0, pos1 = -1;
	num1 = R_MINUS_ONE - 1;
	usfixn64 min;
	short i = 0;
	short j;

	if (sn == 0)
	{
		for (i = 3; i < COEFFICIENT_SIZE; i++)
		{
			xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
		}
	}

//	if (sn>0 && sn < COEFFICIENT_SIZE - 3)
	for (j = 1; j < COEFFICIENT_SIZE - 3; j++)
	{
		if (j != sn)
			continue;
		xm[0] = device_p4_one_step_addition (xm[0], 1, c);

		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (i < j)
				xm[i] = device_p4_one_step_addition (xm[i], 0, c);
			else if (i < j + 3)
				xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
			else if (i == j + 3)
				xm[i] = device_p4_one_step_addition (xm[i], num1, c);
			else
//		for (i = sn + 4; i < COEFFICIENT_SIZE; i++)
				xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
		}
	}
//	if (sn>=COEFFICIENT_SIZE-3)
	for (j = COEFFICIENT_SIZE - 3; j < COEFFICIENT_SIZE; j++)
	{
		if (j != sn)
			continue;
//		min = (j+ 3) & 0xF;
//				min = (j-COEFFICIENT_SIZE+3);
		xm[0] = device_p4_one_step_addition (xm[0], 1, c);
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (i < j - COEFFICIENT_SIZE + 3)
				xm[i] = device_p4_one_step_addition (xm[i], 0, c);
			else if (i == j - COEFFICIENT_SIZE + 3)
//				xm[min] = device_p4_one_step_addition (xm[min], 1, c);
				xm[i] = device_p4_one_step_addition (xm[i], 1, c);
			else if (i < j)
//		for (i = min + 1; i < sn; i++)
				xm[i] = device_p4_one_step_addition (xm[i], 0, c);

//		for (i = sn; i < COEFFICIENT_SIZE; i++)
			else
				xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
		}
	}
//	if(threadIdx.x==0 && blockIdx.x==0)
//	printf("c=%lu \n",c);
	if (c > 0)
	{
		pos1 = -1;
		if (xm[0] != 0 && pos1 == -1)
		{
			pos1 = 0;
			xm[0]--;
			return;
		}
		if (xm[1] != 0 && pos1 == -1)
		{
			pos1 = 1;
			xm[0] = R_MINUS_ONE;
			xm[1]--;
			return;
		}
		if (xm[2] != 0 && pos1 == -1)
		{
			pos1 = 2;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2]--;
			return;
		}
		if (xm[3] != 0 && pos1 == -1)
		{
			pos1 = 3;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3]--;
			return;
		}
		if (xm[4] != 0 && pos1 == -1)
		{
			pos1 = 4;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4]--;
			return;
		}
		if (xm[5] != 0 && pos1 == -1)
		{
			pos1 = 5;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5]--;
			return;
		}
		if (xm[6] != 0 && pos1 == -1)
		{
			pos1 = 6;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6]--;
			return;
		}

		if (xm[7] != 0 && pos1 == -1)
		{
			pos1 = 7;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7]--;
			return;
		}

		if (xm[8] != 0 && pos1 == -1)
		{
			pos1 = 8;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8]--;
			return;
		}

		if (xm[9] != 0 && pos1 == -1)
		{
			pos1 = 9;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9]--;
			return;
		}

		if (xm[10] != 0 && pos1 == -1)
		{
			pos1 = 10;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10]--;
			return;
		}

		if (xm[11] != 0 && pos1 == -1)
		{
			pos1 = 11;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11]--;
			return;
		}

		if (xm[12] != 0 && pos1 == -1)
		{
			pos1 = 12;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12]--;
			return;
		}

		if (xm[13] != 0 && pos1 == -1)
		{
			pos1 = 13;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13]--;
			return;
		}

		if (xm[14] != 0 && pos1 == -1)
		{
			pos1 = 14;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13] = R_MINUS_ONE;
			xm[14]--;
			return;
		}

		if (xm[15] != 0 && pos1 == -1)
		{
			pos1 = 15;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13] = R_MINUS_ONE;
			xm[14] = R_MINUS_ONE;
			xm[15]--;
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
//	short j = 0;
//	xs[0] = 1;
//	if ((i + 3) < COEFFICIENT_SIZE)
//		xs[i + 3] = R - 2;
//	else
//		xs[(i + 3) & 0xF]++; //%COEFFICIENT_SIZE
//
//	for (j = 1; j < i; j++)
//		xs[j] = 0;
//
////	for(j=i;j<COEFFICIENT_SIZE;j++)
////		xs[j]=R-1;
//
//	xs[i] = R - 1;
//	if (i + 1 < COEFFICIENT_SIZE)
//		xs[i + 1] = R - 1;
//	if (i + 2 < COEFFICIENT_SIZE)
//		xs[i + 2] = R - 1;
//	if (i + 3 < COEFFICIENT_SIZE)
//		xs[i + 3] = R - 1;

}
/**********************************************/
__device__ __inline__ void
device_p4_cyclicShift_lhc_complement_and_add_v0 (usfixn64 * xm, short sn)
{
//	if (i==0)
//	{
//		xs[0]=0;
//		xs[1]=1;
//		xs[2]
//	}

	usfixn64 num1 = 0, c = 0, pos1 = -1;
	num1 = R_MINUS_ONE - 1;
	usfixn64 min;
	short i = 0;

	if (sn == 0)
	{
		for (i = 3; i < COEFFICIENT_SIZE; i++)
		{
			xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
		}
	}
	else if (sn < COEFFICIENT_SIZE - 3)
	{
		xm[0] = device_p4_one_step_addition (xm[0], 1, c);
		for (i = 1; i < sn; i++)
			xm[i] = device_p4_one_step_addition (xm[i], 0, c);
		for (i = sn; i < sn + 3; i++)
			xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);

		xm[sn + 3] = device_p4_one_step_addition (xm[sn + 3], num1, c);

		for (i = sn + 4; i < COEFFICIENT_SIZE; i++)
			xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
	}
	else
	{
		min = (sn + 3) & 0xF;
		xm[0] = device_p4_one_step_addition (xm[0], 1, c);
		for (i = 1; i < min; i++)
			xm[i] = device_p4_one_step_addition (xm[i], 0, c);

		xm[min] = device_p4_one_step_addition (xm[min], 1, c);

		for (i = min + 1; i < sn; i++)
			xm[i] = device_p4_one_step_addition (xm[i], 0, c);

		for (i = sn; i < COEFFICIENT_SIZE; i++)
			xm[i] = device_p4_one_step_addition (xm[i], R_MINUS_ONE, c);
	}
//	if(threadIdx.x==0 && blockIdx.x==0)
//	printf("c=%lu \n",c);
	if (c > 0)
	{
		pos1 = -1;
		if (xm[0] != 0 && pos1 == -1)
		{
			pos1 = 0;
			xm[0]--;
			return;
		}
		if (xm[1] != 0 && pos1 == -1)
		{
			pos1 = 1;
			xm[0] = R_MINUS_ONE;
			xm[1]--;
			return;
		}
		if (xm[2] != 0 && pos1 == -1)
		{
			pos1 = 2;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2]--;
			return;
		}
		if (xm[3] != 0 && pos1 == -1)
		{
			pos1 = 3;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3]--;
			return;
		}
		if (xm[4] != 0 && pos1 == -1)
		{
			pos1 = 4;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4]--;
			return;
		}
		if (xm[5] != 0 && pos1 == -1)
		{
			pos1 = 5;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5]--;
			return;
		}
		if (xm[6] != 0 && pos1 == -1)
		{
			pos1 = 6;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6]--;
			return;
		}

		if (xm[7] != 0 && pos1 == -1)
		{
			pos1 = 7;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7]--;
			return;
		}

		if (xm[8] != 0 && pos1 == -1)
		{
			pos1 = 8;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8]--;
			return;
		}

		if (xm[9] != 0 && pos1 == -1)
		{
			pos1 = 9;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9]--;
			return;
		}

		if (xm[10] != 0 && pos1 == -1)
		{
			pos1 = 10;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10]--;
			return;
		}

		if (xm[11] != 0 && pos1 == -1)
		{
			pos1 = 11;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11]--;
			return;
		}

		if (xm[12] != 0 && pos1 == -1)
		{
			pos1 = 12;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12]--;
			return;
		}

		if (xm[13] != 0 && pos1 == -1)
		{
			pos1 = 13;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13]--;
			return;
		}

		if (xm[14] != 0 && pos1 == -1)
		{
			pos1 = 14;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13] = R_MINUS_ONE;
			xm[14]--;
			return;
		}

		if (xm[15] != 0 && pos1 == -1)
		{
			pos1 = 15;
			xm[0] = R_MINUS_ONE;
			xm[1] = R_MINUS_ONE;
			xm[2] = R_MINUS_ONE;
			xm[3] = R_MINUS_ONE;
			xm[4] = R_MINUS_ONE;
			xm[5] = R_MINUS_ONE;
			xm[6] = R_MINUS_ONE;
			xm[7] = R_MINUS_ONE;
			xm[8] = R_MINUS_ONE;
			xm[9] = R_MINUS_ONE;
			xm[10] = R_MINUS_ONE;
			xm[11] = R_MINUS_ONE;
			xm[12] = R_MINUS_ONE;
			xm[13] = R_MINUS_ONE;
			xm[14] = R_MINUS_ONE;
			xm[15]--;
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
//	short j = 0;
//	xs[0] = 1;
//	if ((i + 3) < COEFFICIENT_SIZE)
//		xs[i + 3] = R - 2;
//	else
//		xs[(i + 3) & 0xF]++; //%COEFFICIENT_SIZE
//
//	for (j = 1; j < i; j++)
//		xs[j] = 0;
//
////	for(j=i;j<COEFFICIENT_SIZE;j++)
////		xs[j]=R-1;
//
//	xs[i] = R - 1;
//	if (i + 1 < COEFFICIENT_SIZE)
//		xs[i + 1] = R - 1;
//	if (i + 2 < COEFFICIENT_SIZE)
//		xs[i + 2] = R - 1;
//	if (i + 3 < COEFFICIENT_SIZE)
//		xs[i + 3] = R - 1;

}
/**********************************************/

#endif
