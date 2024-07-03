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


#ifndef BIG_ARITHMETIC_CYCLIC_SHIFT_DEVICE_H_
#define BIG_ARITHMETIC_CYCLIC_SHIFT_DEVICE_H_

#include <stdio.h>
//#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "../../include/BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "bigArithmetic_subtraction_P3_device.h"

/**********************************************/
__device__ __inline__ void
device_p3_oneShiftRight (usfixn64 * xs)
{
	usfixn64 tmp;
	short i;
	{
		tmp = xs[COEFFICIENT_SIZE-1];
		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
			xs[i] = xs[i - 1];
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_oneShiftRight_uvector8 (uvector8 &xs)
{
	usfixn64 tmp;
	short i;

//	for (j = 0; j < shiftNo; j++)
	{
		tmp = xs.i7;
//#pragma unroll 7
//		for (i = COEFFICIENT_SIZE - 1; i > 0; i--)
//			xs[i] = xs[i - 1];
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
__device__ __inline__ void
device_p3_oneShiftLeft_uvector8 (uvector8& xs)
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
//		xs[0] = device_negate_plain(tmp);
		xs.i7 = tmp;
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_oneShiftLeft (usfixn64 * xs)
{
	usfixn64 tmp;
	short i;

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
//		xs[0] = device_negate_plain(tmp);
		xs[7] = tmp;
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_cyclicShift (usfixn64 *xs, short sn)
{
	short i = 0, j = 0;
	usfixn64 ts[8] =
		{ 0 };
	if (sn <= 0)
	{
		return;
	}
	if (sn > 8)
	{
		return;
	}
	j = 8 - sn;
	for (i = 0; i < sn; i++)
	{
		ts[i] = xs[j++];
	}
	for (i = 7 - sn; i >= 0; i--)
	{
		xs[i + sn] = xs[i];
	}
	for (i = 0; i < sn; i++)
	{
		xs[i] = 0;
	}
	device_p3_bigSub (xs, ts, xs);
}

/**********************************************/
__device__ __inline__ void
device_p3_cyclicShift_plain (usfixn64 *xs, short sn)
{
	short i = 0, j = 0;
	usfixn64 ts[8] =
		{ 0 };
	if (sn <= 0)
	{
		return;
	}
	if (sn > 8)
	{
		return;
	}
	j = 8 - sn;
	for (i = 0; i < sn; i++)
	{
		ts[i] = xs[j++];
	}
	for (i = 7 - sn; i >= 0; i--)
	{
		xs[i + sn] = xs[i];
	}
	for (i = 0; i < sn; i++)
	{
		xs[i] = 0;
	}
	device_p3_bigSub (xs, ts, xs);
}

/**********************************************/
//working and faster than permutated_3 and permutated_4
__device__ __inline__ void
device_p3_cyclicShift_permutated_5 (usfixn64 * __restrict__ xs, const short sn,
																 const usfixn64 permutationStride)
{

	if (sn <= 0 || sn > 8)
	{
		return;
	}
	short i = 0, j = 0;
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
	uConstArray8_align8 ts, ys;

//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
		ts.i[i] = 0;
		ys.i[i] = 0;
	}

//	memset(ts.i,0x00, COEFFICIENT_SIZE*sizeof(usfixn64));
//	memset(ys.i,0x00, COEFFICIENT_SIZE*sizeof(usfixn64));

//	memset(ts.i,0x00,64);
//		memset(ys.i,0x00,64);

	j = 8 - sn;
	usfixn64 offset = j * permutationStride;
//	for (i = 0; i < sn; i++)
//	{
//		ts[i] = xs[j++];
//	}
//	for (i = 0; i < sn; i++)
//	{
////		ts[i] = xs[offset];
//		ts.i[i] = xs[offset];
//		offset += permutationStride;
//	}

	for (i = 0; i < sn; i++)
	{
		//		ts[i] = xs[offset];
		ts.i[i] = xs[offset];
		offset += permutationStride;
	}

	offset = (7 - sn) * permutationStride;
	j = 7;
	for (i = 7 - sn; i >= 0; i--)
	{
//		xs[i + sn] = xs[i];
//		ys[i + sn] = xs[offset];
//		ys.i[i + sn] = xs[offset];
		ys.i[j--] = xs[offset];
		offset -= permutationStride;
	}
//	for (i = 0; i < sn; i++)
//	{
//		xs[i] = 0;
//	}
//	device_p3_bigSub(ys, ts, xs);
//	device_bigSub_plain(&ys.i[0], &ts.i[0]);
//	device_bigSub_plain_1(&ys.i[0], &ts.i[0]);
	device_p3_bigPrimeSub_plain_ptx_v0 (ys.i, ts.i);
//	device_bigPrimeSub_permutated_3(xs,ts, permutationStride);
	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
//		xs[offset] = ys[i];
//		xs[offset] = ys.i[i];
		xs[offset] = ys.i[i];
		offset += permutationStride;
//		xs[offset8[i]] = ys.i[i];
	}
}

/**********************************************/
//computing cyclic shift with vectors of size 8 on register
__device__ __inline__  void
device_p3_cyclicShift_permutated_7 (usfixn64 * __restrict__ xs, const short & sn,
																 const usfixn64 & permutationStride)
{
	if (sn <= 0 || sn > 8)
	{
		return;
	}
	short i = 0, j = 0;
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
//	uConstArray8_align8 ts, ys;
	uvector8 ts, ys;
	usfixn64 offset = 0;
//#pragma unroll COEFFICIENT_SIZE
//	for(i=0;i<COEFFICIENT_SIZE;i++)
//	{
//		ys.i[i]=xs[offset];
//		ts.i[i]=0;
//		offset+=permutationStride;
//	}

//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
////		ys.i[i]=xs[offset];
//		setUvector8Element(ys, i, xs[offset]);
////		ts.i[i]=0;
//		setUvector8Element(ts, i, 0);
//		offset += permutationStride;
//	}

//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		ys.i[i]=xs[offset];
//
//		setUvector8Element(ys, i, xs[offset]);
//		ts.i[i]=0;
//		setUvector8Element(ts, i, 0);

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

		ts.i0 = 0;
		ts.i1 = 0;
		ts.i2 = 0;
		ts.i3 = 0;
		ts.i4 = 0;
		ts.i5 = 0;
		ts.i6 = 0;
		ts.i7 = 0;
	}

	j = 8 - sn;
//	usfixn64 offset = j * (permutationStride);
//	usfixn64 offset = (8-sn) * (permutationStride);
//	for (i = 0; i < sn; i++)
//	{
//		ts[i] = xs[j++];
//	}
	for (i = 0; i < sn; i++)
	{
//		ts[i] = xs[offset];
//		ts.i[i] = xs[offset];

//		ts.i[i] = ys.i[j];
		setUvector8Element (ts, i, getUvector8Element (ys, j));
		j++;
//		offset += permutationStride;
	}

//	offset = (7 - sn) * (permutationStride);
	for (i = 7 - sn; i >= 0; i--)
	{
//		xs[i + sn] = xs[i];
//		ys[i + sn] = xs[offset];
//		ys.i[i + sn] = xs[offset];

//		ys.i[i + sn] = ys.i[i];
		j = i + sn;
		setUvector8Element (ys, j, getUvector8Element (ys, i));

//		offset -= permutationStride;
	}
	for (i = 0; i < sn; i++)
	{
////		xs[i] = 0;
//		ys.i[i] = 0;
		setUvector8Element (ys, i, 0);
	}
//	memset(ys.i, 0x00, sn*sizeof(usfixn64));
//	device_p3_bigSub(ys, ts, xs);
//	device_bigSub_plain(&ys.i[0], &ts.i[0]);
//	device_bigSub_plain_1(&ys.i[0], &ts.i[0]);
//	device_bigSub_plain_1(ys.i, ts.i);
//	device_bigSub_plain_1(ys.i, ts.i);
	device_p3_bigSub_uVector_plain_2 (ys, ts);
//	device_bigPrimeSub_permutated_3(xs,ts, permutationStride);
	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
////		xs[offset] = ys[i];
////		xs[offset] = ys.i[i];
//		xs[offset] = getUvector8Element(ys, i);
//		offset += (permutationStride);
//	}

//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
	{
//		xs[offset] = ys[i];
//		xs[offset] = ys.i[i];
		xs[offset] = ys.i0;
		offset += (permutationStride);
		xs[offset] = ys.i1;
		offset += (permutationStride);
		xs[offset] = ys.i2;
		offset += (permutationStride);
		xs[offset] = ys.i3;
		offset += (permutationStride);
		xs[offset] = ys.i4;
		offset += (permutationStride);
		xs[offset] = ys.i5;
		offset += (permutationStride);
		xs[offset] = ys.i6;
		offset += (permutationStride);
		xs[offset] = ys.i7;
//														offset += (permutationStride);
	}
}

/**********************************************/
__device__ __inline__  void
device_p3_negate_1 (usfixn64 * __restrict__ xs, const usfixn64 * __restrict__ ys,
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

/**********************************************/
__device__ __inline__  void
device_p3_negate_2 (usfixn64 &xs, const usfixn64 &ys, short & c)
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

/**********************************************/
//using uvector8 on register which makes it register intensive (a lot, up to 63)
__device__ __inline__  void
device_p3_cyclicShift_permutated_9 (usfixn64 * __restrict__ xs, const short & sn,
																 const usfixn64 & permutationStride)
{
	if (sn <= 0 || sn > 8)
	{
		return;
	}
//	short i = 0;
	short j = 0;
	short c = 0;
	short pos = -1;
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
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

	tmpY = 0;
	for (j = 0; j < 8 - sn; j++)
	{
		tmpX = getUvector8Element (ys, j);
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

//######################################################

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
	//######################################################
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
}

/**********************************************/
//using uvector8 on register which makes it register intensive (a lot, up to 63)
//this device_cyclicShift_twice_permutated_1 computes two cyclic shifts if sn>8 and one if sn<=8
//this is cyclicShift over 8 digits, not 16
//e.g. sn=12; first does cyclicShift(sn) then computes cyclicShift(sn-8)
__device__ __inline__  void
device_p3_cyclicShift_twice_permutated_1 (usfixn64 * __restrict__ xs,
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

	if (sn > COEFFICIENT_SIZE)
	{
		sn2 = sn - COEFFICIENT_SIZE;
		sn = COEFFICIENT_SIZE;
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
	offset += permutationStride;

//	ys.i8 = xs[offset];
//	offset += permutationStride;
//	ys.i9 = xs[offset];
//	offset += permutationStride;
//	ys.i10 = xs[offset];
//	offset += permutationStride;
//	ys.i11 = xs[offset];
//	offset += permutationStride;
//	ys.i12 = xs[offset];
//	offset += permutationStride;
//	ys.i13 = xs[offset];
//	offset += permutationStride;
//	ys.i14 = xs[offset];
//	offset += permutationStride;
//	ys.i15 = xs[offset];
//	offset += permutationStride;

//	for(j=0;j<8;j++)
//	{
//		setUvector16Element(ys,j, xs[offset]);
//		offset+=permutationStride;
//	}

	//######################################################
	//doing the subtraction
	for (j = COEFFICIENT_SIZE - sn; j < COEFFICIENT_SIZE; j++)
	{
		tmpX = 0;
		tmpY = getUvector8Element (ys, j);
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

	tmpY = 0;
	for (j = 0; j < COEFFICIENT_SIZE- sn; j++)
	{
		tmpX = getUvector8Element (ys, j);
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

//######################################################

	if (sn <= COEFFICIENT_SIZE)
	{
		//######################################################
		//Finding pos
		for (j = COEFFICIENT_SIZE - sn; j < COEFFICIENT_SIZE; j++)
		{
			if (getUvector8Element (ys, j) < R_MINUS_ONE)
			{
				pos = j;
				break;
			}
		}

		if (pos == -1)
		{
			for (j = 0; j < COEFFICIENT_SIZE - sn; j++)
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
			if (pos >= (COEFFICIENT_SIZE - sn))	//in first group
			{
				for (j = COEFFICIENT_SIZE - sn; j < pos; j++)
				{
					setUvector8Element (ys, j, 0);
				}
				incUvector8Element (ys, pos, 1);
			}
			else
			{
				for (j = 8 - sn; j < COEFFICIENT_SIZE; j++)
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
			setUvector8Element (ys, COEFFICIENT_SIZE - sn, ULMAX);
			for (j = COEFFICIENT_SIZE- sn + 1; j < COEFFICIENT_SIZE; j++)
				setUvector8Element (ys, j, 0);
			for (j = 0; j < COEFFICIENT_SIZE - sn; j++)
				setUvector8Element (ys, j, 0);
		}
	}

	if (sn2 <= 0)
	{
		// copying results back to global mem
		offset = 0;

		for (j = COEFFICIENT_SIZE - sn; j < COEFFICIENT_SIZE; j++)
		{
			xs[offset] = getUvector8Element (ys, j);
			offset += permutationStride;
		}

		for (j = 0; j < COEFFICIENT_SIZE - sn; j++)
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
	for (j = COEFFICIENT_SIZE - sn2; j < COEFFICIENT_SIZE; j++)
	{
		tmpX = 0;
		tmpY = getUvector8Element (ys, j);
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

	tmpY = 0;
	for (j = 0; j < COEFFICIENT_SIZE - sn2; j++)
	{
		tmpX = getUvector8Element (ys, j);
		device_p3_negate_1 (&tmpX, &tmpY, &c);
		setUvector8Element (ys, j, tmpX);
	}

//######################################################

	//######################################################
	//Finding pos
	for (j = COEFFICIENT_SIZE - sn2; j < COEFFICIENT_SIZE; j++)
	{
		if (getUvector8Element (ys, j) < R_MINUS_ONE)
		{
			pos = j;
			break;
		}
	}

	if (pos == -1)
	{
		for (j = 0; j < COEFFICIENT_SIZE - sn2; j++)
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
		if (pos >= (COEFFICIENT_SIZE- sn2))	//in first group
		{
			for (j = COEFFICIENT_SIZE - sn2; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
			}
			incUvector8Element (ys, pos, 1);
		}
		else
		{
			for (j = COEFFICIENT_SIZE - sn2; j < COEFFICIENT_SIZE; j++)
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
		setUvector8Element (ys, COEFFICIENT_SIZE - sn2, ULMAX);
		for (j = COEFFICIENT_SIZE - sn2 + 1; j < COEFFICIENT_SIZE; j++)
			setUvector8Element (ys, j, 0);
		for (j = 0; j < COEFFICIENT_SIZE - sn2; j++)
			setUvector8Element (ys, j, 0);
	}
	//######################################################
	// copying results back to global mem
	offset = 0;

	for (j = COEFFICIENT_SIZE - sn2; j < COEFFICIENT_SIZE; j++)
	{
		xs[offset] = getUvector8Element (ys, j);
		offset += permutationStride;
	}

	for (j = 0; j < COEFFICIENT_SIZE - sn2; j++)
	{
		xs[offset] = getUvector8Element (ys, j);
		offset += permutationStride;
	}
}

/**********************************************/
//using uConstArray on lmem
__device__ __inline__  void
device_p3_cyclicShift_permutated_11 (usfixn64 * xs, const short & sn,
																	const usfixn64 & permutationStride)
{
	if (sn <= 0 || sn > 8)
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
	uvector8 ys;
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
////		device_negate_1(&tmpX, &tmpY, &c);
//		device_p3_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//	}

//	for (j = 8 - sn; j < 8; j++)
//		{
//			tmpX = 0;
//	//		tmpY = getUvector8Element(ys, j);
//			tmpY = ys.i[j];
//	//		device_negate_1(&tmpX, &tmpY, &c);
//			device_p3_negate_2(tmpX, tmpY, c);
//	//		setUvector8Element(ys, j, tmpX);
//			ys.i[j] = tmpX;
//			bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		}

	for (j = 0; j < sn; j++)
		device_p3_oneShiftRight_uvector8 (ys);

	j = 0;
//	j = 8 - sn;
//	if (j < 8)
	if (sn > 0)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i0;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	tmpY = 0;
//	for (j = 0; j < 8 - sn; j++)
//	{
////		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
////		device_negate_1(&tmpX, &tmpY, &c);
//		device_p3_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
//	}

//	for(j=0;j<sn;j++)
//	device_oneShiftLeft(ys.i);

	for (j = 0; j < sn; j++)
		device_p3_oneShiftLeft_uvector8 (ys);

	tmpY = 0;
	j = 0;
//		for (j = 0; j < 8 - sn; j++)
//	if (j < 8 - sn)
	if (sn < 8)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i0;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 7)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i1;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 6)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i2;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 5)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i3;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		j++;
	}
	if (sn < 4)
//	if (j < 8 - sn)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i4;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 3)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i5;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 2)
	{
//				tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i6;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 1)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i7;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
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
		if (pos >= (8 - sn))	//in first group
		{

			for (j = 8 - sn; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;

			}
			incUvector8Element (ys, pos, 1);
//			ys.i[pos]++;

		}
		else
		{
			for (j = 8 - sn; j < 8; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			for (j = 0; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			incUvector8Element (ys, pos, 1);
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
		device_p3_oneShiftRight_uvector8 (ys);

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
		setUvector8Element (ys, 0, ULMAX);
		ys.i1 = 0x0;
		ys.i2 = 0x0;
		ys.i3 = 0x0;
		ys.i4 = 0x0;
		ys.i5 = 0x0;
		ys.i6 = 0x0;
		ys.i7 = 0x0;

	}

	//######################################################
	// copying results back to global mem
	offset = 0;

//	for(j=0;j<sn;j++)
//		device_oneShiftRight(ys.i);
//	for (j = 0; j < sn; j++)
//		device_p3_oneShiftRight_uvector8(ys);

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
	}
}

/**********************************************/
__device__ __inline__  void
device_p3_cyclicShift_permutated_12 (usfixn64 * xs)
{
//	if (sn <= 0 || sn > 8)
//	{
//		return;
//	}
//	short i = 0;
	usfixn64 permutationStride = 1;
	const short sn = 1;
	short j = 0;
	short c = 0;
	usfixn32 pos = -1;
//	usfixn64 ts[8] = { 0 };
//	usfixn64 ys[8] = { 0 };
//	uConstArray8_align8 ys;
	uvector8 ys;
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
////		device_negate_1(&tmpX, &tmpY, &c);
//		device_p3_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//	}

//	for (j = 8 - sn; j < 8; j++)
//		{
//			tmpX = 0;
//	//		tmpY = getUvector8Element(ys, j);
//			tmpY = ys.i[j];
//	//		device_negate_1(&tmpX, &tmpY, &c);
//			device_p3_negate_2(tmpX, tmpY, c);
//	//		setUvector8Element(ys, j, tmpX);
//			ys.i[j] = tmpX;
//			bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		}

	for (j = 0; j < sn; j++)
		device_p3_oneShiftRight_uvector8 (ys);

	j = 0;
//	j = 8 - sn;
//	if (j < 8)
	if (sn > 0)
	{
		tmpX = 0;
		//		tmpY = getUvector8Element(ys, j);
//		tmpY = ys.i[j];
		tmpY = ys.i0;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
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
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
		j++;
	}

//	tmpY = 0;
//	for (j = 0; j < 8 - sn; j++)
//	{
////		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
////		device_negate_1(&tmpX, &tmpY, &c);
//		device_p3_negate_2(tmpX, tmpY, c);
////		setUvector8Element(ys, j, tmpX);
//		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
//	}

//	for(j=0;j<sn;j++)
//	device_oneShiftLeft(ys.i);

	for (j = 0; j < sn; j++)
		device_p3_oneShiftLeft_uvector8 (ys);

	tmpY = 0;
	j = 0;
//		for (j = 0; j < 8 - sn; j++)
//	if (j < 8 - sn)
	if (sn < 8)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i0;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i0 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 7)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i1;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i1 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 6)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i2;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i2 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 5)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i3;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i3 = tmpX;
		j++;
	}
	if (sn < 4)
//	if (j < 8 - sn)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i4;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i4 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 3)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i5;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i5 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 2)
	{
//				tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i6;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i6 = tmpX;
		j++;
	}
//	if (j < 8 - sn)
	if (sn < 1)
	{
		//		tmpX = getUvector8Element(ys, j);
//		tmpX = ys.i[j];
		tmpX = ys.i7;
		//		device_negate_1(&tmpX, &tmpY, &c);
		device_p3_negate_2 (tmpX, tmpY, c);
		//		setUvector8Element(ys, j, tmpX);
		bitFlag |= ((tmpX < R_MINUS_ONE) << j);
//		ys.i[j] = tmpX;
		ys.i7 = tmpX;
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
		if (pos >= (8 - sn))	//in first group
		{

			for (j = 8 - sn; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;

			}
			incUvector8Element (ys, pos, 1);
//			ys.i[pos]++;

		}
		else
		{
			for (j = 8 - sn; j < 8; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			for (j = 0; j < pos; j++)
			{
				setUvector8Element (ys, j, 0);
//				ys.i[j] = 0;
			}
			incUvector8Element (ys, pos, 1);
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
		device_p3_oneShiftRight_uvector8 (ys);

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
		setUvector8Element (ys, 0, ULMAX);
		ys.i1 = 0x0;
		ys.i2 = 0x0;
		ys.i3 = 0x0;
		ys.i4 = 0x0;
		ys.i5 = 0x0;
		ys.i6 = 0x0;
		ys.i7 = 0x0;

	}

	//######################################################
	// copying results back to global mem
	offset = 0;

//	for(j=0;j<sn;j++)
//		device_oneShiftRight(ys.i);
//	for (j = 0; j < sn; j++)
//		device_p3_oneShiftRight_uvector8(ys);

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
	}
}

/**********************************************/

#endif
