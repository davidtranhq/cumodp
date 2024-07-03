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
#ifndef BIG_ARITHMETIC_ADDITION_DEVICE_H_
#define BIG_ARITHMETIC_ADDITION_DEVICE_H_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"
#include "bigArithmetic_subtraction_P4_device.h"
#include "bigArithmetic_cyclicShift_P4_device.h"

/**********************************************/
__device__ __forceinline__ void
device_p4_bigPrimeAdd_correct (usfixn64 *xm, usfixn64 *ym, usfixn64 *um)
{
	unsigned short c = 0;
	short pos1;
	usfixn64 num1, num2;

//	usfixn64 um[8];

	num1 = 0;
	num2 = R - 1;
	short i;
//	usfixn64 tmp = 0;

	for (i = 0; i < 8; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
			um[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			um[i] = num1 - R;
		}
		else
		{
			um[i] = num1;
			c = 0;
		}
	}
	if (c > 0)
	{
		pos1 = -1;
		for (i = 0; i < 8; i++)
		{
			if (um[i] != 0)
			{
				pos1 = i;
				break;
			}
		}
		if (pos1 >= 0)
		{
			for (i = 0; i < pos1; i++)
			{
				um[i] = num2;
			}
			um[pos1]--;
		}
		else
		{
			um[0] = ULMAX;
			for (i = 1; i < 8; i++)
			{
				um[i] = 0;
			}
		}
	}
}

/******************************************************************************/
//
//__device__ __forceinline__ void  bigPrimeAdd_correct_improved1(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, usfixn64 * __restrict__ um)
//{
//	unsigned short c = 0;
//	short pos1;
//	usfixn64 num1;
////	, num2;
//
////	usfixn64 um[8];
//
//	num1 = 0;
////	num2 = R - 1;
//	short i;
////	usfixn64 tmp = 0;
//
//	for (i = 0; i < 8; i++)
//	{
//		num1 = xm[i] + ym[i] + c;
//		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
//		{
//			um[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			um[i] = num1 - R;
//		}
//		else
//		{
//			um[i] = num1;
//			c = 0;
//		}
//	}
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (um[i] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//		}
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//			{
//				um[i] = R - 1;
//			}
//			um[pos1]--;
//		}
//		else
//		{
//			um[0] = ULMAX;
//			for (i = 1; i < 8; i++)
//			{
//				um[i] = 0;
//			}
//		}
//	}
//}
/******************************************************************************/
//
//__device__ __forceinline__ void  bigPrimeAddStepped_correct(usfixn64 *xm, usfixn64 *ym,
//		usfixn64 *um, short step)
//{
//	unsigned short c = 0;
//	short pos1;
//	usfixn64 num1, num2;
//
////	usfixn64 um[8];
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
////	usfixn64 tmp = 0;
//
//	for (i = 0; i < step; i++)
//	{
//		num1 = xm[i] + ym[i] + c;
//		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
//		{
//			um[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			um[i] = num1 - R;
//		}
//		else
//		{
//			um[i] = num1;
//			c = 0;
//		}
//	}
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < step; i++)
//		{
//			if (um[i] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//		}
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//			{
//				um[i] = num2;
//			}
//			um[pos1]--;
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
//}
/******************************************************************************/
//
//__device__ __forceinline__ void bigPrimeAdd2const_correct(
//		const usfixn64 * __restrict__ x, const usfixn64 * __restrict__ y,
//		usfixn64 * __restrict__ u)
//{
//
//	*u = *u + *x + *y;
////	unsigned short c = 0;
////	short pos1;
////	usfixn64 num1, num2;
////	usfixn64 uTmp = *u;
//////	usfixn64 um[8];
////
////	num1 = 0;
////	num2 = R - 1;
////	short i;
//////	usfixn64 tmp = 0;
////
//////	for (i = 0; i < step; i++)
//////		num1 = xm[i] + ym[i] + c;
////	//		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
////	//		{
////	//			um[i] = num1 + RC;
////	//			c = 1;
////	//		}
////	//		else if (num1 >= R)
////	//		{
////	//			c = 1;
////	//			um[i] = num1 - R;
////	//		}
////	//		else
////	//		{
////	//			um[i] = num1;
////	//			c = 0;
////	//		}
////	num1= *x + *y + c;
////
////	if (num1 < *x || num1 < *y)//there is overflow/truncation
////	{
////		*u = num1 + RC;
////		c = 1;
////	}
////	else if (num1 >= R)
////	{
////		c = 1;
////		*u = num1 - R;
////	}
////	else
////	{
////		*u = num1;
////		c = 0;
////	}
////
////	num1= *u + uTmp + c;
////
////	if (num1 < *u || num1 < uTmp)		//there is overflow/truncation
////	{
////		*u = num1 + RC;
////		c = 1;
////	}
////	else if (num1 >= R)
////	{
////		c = 1;
////		*u = num1 - R;
////	}
////	else
////	{
////		*u = num1;
////		c = 0;
////	}
//
////	if (c > 0)
////	{
////		pos1 = -1;
////		for (i = 0; i < step; i++)
////		{
////			if (um[i] != 0)
////			{
////				pos1 = i;
////				break;
////			}
////		}
////		if (pos1 >= 0)
////		{
////			for (i = 0; i < pos1; i++)
////			{
////				um[i] = num2;
////			}
////			um[pos1]--;
////		}
////		else
////		{
////			um[0] = ULMAX;
////			for (i = 1; i < step; i++)
////			{
////				um[i] = 0;
////			}
////		}
//}
/******************************************************************************/

__device__ __inline__ void
device_bigPrimeAdd_plain (usfixn64 *__restrict__ xm, usfixn64 *__restrict__ ym)
{
	unsigned short c = 0;
	short pos1;
	usfixn64 num1, num2;

//	usfixn64 um[8];

	num1 = 0;
	num2 = R - 1;
	short i;
//	usfixn64 tmp = 0;

	for (i = 0; i < 8; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
	}
	if (c > 0)
	{
		pos1 = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[i] != 0)
			{
				pos1 = i;
				break;
			}
		}
		if (pos1 >= 0)
		{
			for (i = 0; i < pos1; i++)
			{
				xm[i] = num2;
			}
			xm[pos1]--;
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
}
/******************************************************************************/
//__device__ __inline__ void device_bigPrimeAdd_plain_uvector8(
//		usfixn64 * __restrict__ xm, uvector8 &ymVector)
//{
//	unsigned short c = 0;
//
//	short pos1;
//	usfixn64 num1, num2;
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
//	usfixn64 ym;
////	for (i = 0; i < 8; i++)
//	{
//
//		i = 0;
//		ym = ymVector.i0;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 1;
//		ym = ymVector.i1;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 2;
//		ym = ymVector.i2;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 3;
//		ym = ymVector.i3;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 4;
//		ym = ymVector.i4;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 5;
//		ym = ymVector.i5;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 6;
//		ym = ymVector.i6;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//
//		i = 7;
//		ym = ymVector.i7;
//		num1 = xm[i] + ym + c;
//		if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//	}
//
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
////			if (xm[i] != 0)
//			if (xm[i] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//		}
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//			{
//				xm[i] = num2;
//			}
//			xm[pos1]--;
//		}
//		else
//		{
//			xm[0] = ULMAX;
//			for (i = 1; i < 8; i++)
//			{
//				xm[i] = 0;
//			}
//		}
//	}
//}
/******************************************************************************/
/******************************************************************************/
//__device__ __inline__ void device_bigPrimeAdd_plain_uvector8_step(
//		usfixn64 & __restrict__ xm, usfixn64& ym, unsigned short & c, short i)
//{
////	unsigned short c = 0;
//
////	short pos1;
//	usfixn64 num1;
////	, num2;
//
//	num1 = 0;
////	num2 = R - 1;
////	short i;
////	usfixn64 ym;
////	for (i = 0; i < 8; i++)
//
////	num1 = xm[i] + ym + c;
//	num1 = xm + ym + c;
////	if (num1 < xm[i] || num1 < ym) //there is overflow/truncation
//	if (num1 < xm || num1 < ym) //there is overflow/truncation
//	{
////		xm[i] = num1 + RC;
//		xm = num1 + RC;
//		c = 1;
//	}
//	else if (num1 >= R)
//	{
//		c = 1;
////		xm[i] = num1 - R;
//		xm = num1 - R;
//	}
//	else
//	{
////		xm[i] = num1;
//		xm = num1;
//		c = 0;
//	}
//
//}
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_plain_inPlace (usfixn64 * __restrict__ xm,
																	usfixn64 * __restrict__ ym)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
	}

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
}
/******************************************************************************/
//ym is permutated
__device__ __inline__ void
device_bigPrimeAdd_plain_inPlace_steps (usfixn64 * __restrict__ xm,
																				const usfixn64 * __restrict__ ym,
																				usfixn64 permutationStride,
																				usfixn64 steps)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;

	usfixn64 offset = steps * permutationStride;

//	for (i = steps; i < COEFFICIENT_SIZE; i++)
//	{
//		num1 = xm[i] + ym[offset] + c;
//		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//		offset += permutationStride;
//	}

	i = 1;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

	if (steps <= i)
	{
		num1 = xm[i] + ym[offset] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}
	i++;

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
}
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_plain_inPlace_lhc_complement (usfixn64 * __restrict__ xm,
																								 usfixn64 * ym, short j)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
	}

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
}
/******************************************************************************/
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_plain_x_permutated_y_inPlace (usfixn64 * __restrict__ xm,
																								 usfixn64 * __restrict__ ym,
																								 usfixn64 permutationStride)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

	usfixn64 offset = 0;
//	num1 = 0;
//	num2 = R - 1;
	short i;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[offset]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		offset += permutationStride;
	}

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
}
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_plain_ptx_v0 (usfixn64 * __restrict__ xm,
																 usfixn64 * __restrict__ ym)
{
//	unsigned short c = 0;

	usfixn16 c = 0;
//	short pos1;
	usfixn32 pos1;
	usfixn64 num1, num2;
	usfixn32 bitsFlag = 0x0;
//	num1 = 0;
//	num2 = R - 1;
	short i;
	short bit;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
			xm[i] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[i] = num1 - R;
		}
		else
		{
			xm[i] = num1;
			c = 0;
		}
		bit = (xm[i] > 0);
		bitsFlag |= (bit << i);
	}

	asm("brev.b32 %0,%0;":"+r"(bitsFlag):);
	asm("bfind.u32 %0, %1;":"=r"(pos1):"r"(bitsFlag));
	if (pos1 == 0xFFFFFFFF)
		pos1 = -1;
	else
		pos1 = 31 - pos1;

	if (c > 0)
	{
		if (pos1 == -1)
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
		else
		{
			xm[pos1]--;
			for (i = 0; i < pos1; i++)
				xm[i] = R_MINUS_ONE;
		}
	}
}
/******************************************************************************/

//__device__ __forceinline__ void  bigPrimeAdd_permutated(usfixn64 *xm, usfixn64 *ym, usfixn64 *um)
//{
//
////	short pos1;
////	usfixn64 num1, num2;
////	unsigned short c = 0;
//	short pos1;
//	unsigned short c = 0;
//	usfixn64 num1, num2;
//
//	//	usfixn64 um[8];
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
//	usfixn64 pos = threadIdx.x;
//
//	//	usfixn64 tmp = 0;
//
//	for (i = 0; i < 8; i++)
//	{
////		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
//
////		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
////		{
////			um[i] = num1 + RC;
////			c = 1;
////		}
////		else if (num1 >= R)
////		{
////			c = 1;
////			um[i] = num1 - R;
////		}
////		else
////		{
////			um[i] = num1;
////			c = 0;
////		}
//		if (num1 < xm[pos] || num1 < ym[pos]) //there is overflow/truncation
//		{
//			um[pos] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			um[pos] = num1 - R;
//		}
//		else
//		{
//			um[pos] = num1;
//			c = 0;
//		}
//	}
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (um[i] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//		}
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//			{
//				um[i] = num2;
//			}
//			um[pos1]--;
//		}
//		else
//		{
//			um[0] = ULMAX;
//			for (i = 1; i < 8; i++)
//			{
//				um[i] = 0;
//			}
//		}
//	}
//
//}
/******************************************************************************/

__device__ __forceinline__ void
device_p4_bigPrimeAdd_permutated (const usfixn64 *xm, const usfixn64 * ym, usfixn64 *um,
												const usfixn64 idx, const short permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1, num2;

	num1 = 0;
	num2 = R - 1;
	short i;
//	usfixn64 pos = threadIdx.x;

	short offset = 0;
	short pos1;
//	usfixn64 tmp = 0;

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
		//		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
		num1 = xm[idx + offset] + ym[idx + offset] + c;
		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])	//there is overflow/truncation
		{
			um[idx + offset] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			um[idx + offset] = num1 - R;
		}
		else
		{
			um[idx + offset] = num1;
			c = 0;
		}
		offset += permutationStride;
	}

	offset = 0;
	if (c > 0)
	{
		pos1 = -1;
		for (i = 0; i < 8; i++)
		{
			if (um[idx + offset] != 0)
			{
				pos1 = i;
				break;
			}
			offset += permutationStride;
		}
		if (pos1 >= 0)	// shouldn't it be >0?
		{
			offset = 0;
			for (i = 0; i < pos1; i++)
			{
				um[idx + offset] = num2;
				offset += permutationStride;
			}
			offset = pos1 * permutationStride;
			um[idx + offset]--;
//			um[pos1*permutationStride+idx]--;
		}
		else
		{
//			um[0] = ULMAX;
			um[idx] = ULMAX;
			offset = permutationStride;
			for (i = 1; i < 8; i++)
			{
//				um[i] = 0;
				um[idx + offset] = 0;
			}
		}
	}

//	offset=0;
////	permutationStride=32;
//	for(i=0; i<8; i++)
//	{
////		um[idx+offset]=idx;
////		um[idx+offset]=blockIdx.x;
////		um[idx+offset]=111;
//		offset+=permutationStride;
//	}
}
/******************************************************************************/
//__device__ __forceinline__ void  bigPrimeAdd_permutated_2(const usfixn64 __restrict__ * xm,
//		const usfixn64 * __restrict__ ym, usfixn64 *__restrict__ um,
//		const usfixn64 idx, const short permutationStride)
//{
//	unsigned short c = 0;
//	usfixn64 num1;
////	, num2;
//
//	num1 = 0;
////	num2 = R - 1;
//	short i;
////	usfixn64 pos = threadIdx.x;
//
//	short offset = 0;
//	short pos1;
////	usfixn64 tmp = 0;
//
//	offset = 0;
//
//#pragma unroll 8
//	for (i = 0; i < 8; i++)
//	{
//		//		num1 = xm[i] + ym[i] + c;
////		num1 = xm[pos] + ym[pos] + c;
//		num1 = xm[idx + offset] + ym[idx + offset] + c;
//		c = 0;
////		um[idx + offset] = num1;
////		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])		//there is overflow/truncation
//		{
////			um[idx + offset] = num1 + RC;
//			um[idx + offset] = num1
//					+ RC * (num1 < xm[idx + offset] || num1 < ym[idx + offset]);
//			c = (num1 < xm[idx + offset] || num1 < ym[idx + offset]);
//		}
////		else if (num1 >= R)
//		{
////			c = 1;
//			c = (num1 >= R);
//
////			um[idx + offset] = num1 - R;
//			um[idx + offset] = num1 - R * (num1 >= R);
//		}
////		else
//		{
////			um[idx + offset] = num1;
////			c = 0;
//		}
//		offset += permutationStride;
//	}
//
//	offset = 0;
//	if (c > 0)
//	{
//		pos1 = -1;
//#pragma unroll 8
//		for (i = 0; i < 8; i++)
//		{
//			if (um[idx + offset] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//			offset += permutationStride;
//		}
//		if (pos1 >= 0)
//		{
//			offset = 0;
//			for (i = 0; i < pos1; i++)
//			{
//				um[idx + offset] = R - 1;
////				um[idx + offset] = num2;
//				offset += permutationStride;
//			}
//			offset = pos1 * permutationStride;
//			um[idx + offset]--;
////			um[pos1*permutationStride+idx]--;
//		}
//		else
//		{
////			um[0] = ULMAX;
//			um[idx] = ULMAX;
//			offset = permutationStride;
//#pragma unroll 7
//			for (i = 1; i < 8; i++)
//			{
////				um[i] = 0;
//				um[idx + offset] = 0;
////				memset(&um[idx + offset],0x00, 7*sizeof(usfixn64));
//			}
//		}
//	}
//
////	offset=0;
//////	permutationStride=32;
////	for(i=0; i<8; i++)
////	{
//////		um[idx+offset]=idx;
//////		um[idx+offset]=blockIdx.x;
//////		um[idx+offset]=111;
////		offset+=permutationStride;
////	}
//}
/******************************************************************************/
//__device__ __forceinline__ void  bigPrimeAdd_permutated_3(const usfixn64 __restrict__ * xm,
//		const usfixn64 * __restrict__ ym, usfixn64 *__restrict__ um,
//		const usfixn64 idx, const short permutationStride)
//{
//	uConstArray8_align8 x, y, u;
//	short i = 0;
//	int offset = 0;
//#pragma unroll 8
//	for (i = 0; i < 8; i++)
//	{
//		x.i[i] = xm[idx + offset];
//		y.i[i] = ym[idx + offset];
//		offset += permutationStride;
//	}
//
//	bigPrimeAdd2const_correct(&x.i[0], &y.i[0], &u.i[0]);
//
//	offset = 0;
//#pragma unroll 8
//	for (i = 0; i < 8; i++)
//	{
//		um[idx + offset] = u.i[i];
//		offset += permutationStride;
//	}
//}
/******************************************************************************/
__device__ __forceinline__ void
device_bigPrimeAdd_permutated_ptx_v0 (usfixn64 * xm, usfixn64 *ym,
																			const usfixn64 permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1 = 0;
//	, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
//	usfixn64 pos = threadIdx.x;

	usfixn64 offset = 0;
//	short posAdd = -1;
	usfixn32 posAdd = -1;
	usfixn32 bit;
	usfixn32 bitsFlag;
//	usfixn64 tmp = 0;

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
		num1 = xm[offset] + ym[offset] + c;
		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
		{
//			xm[offset] = num1 + RC;
			num1 += RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
//			xm[offset] = num1 - R;
			num1 -= R;
		}
		else
		{
//			xm[offset] = num1;
			c = 0;
		}
//		posAdd = (posAdd == -1 && num1 > 0) * i;

		xm[offset] = num1;
		bit = (num1 > 0);
		bitsFlag |= (bit << i);
		offset += permutationStride;
	}

	asm("brev.b32 %0,%0;":"+r"(bitsFlag):);
	asm("bfind.u32 %0, %1;":"=r"(posAdd):"r"(bitsFlag));
	if (posAdd == 0xFFFFFFFF)
		posAdd = -1;
	else
		posAdd = 31 - posAdd;

	offset = 0;

	offset = 0;
	if (c > 0)
	{

		if (posAdd == -1)
		{
			xm[offset] = ULMAX;
			offset += permutationStride;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
		else
		{
			for (i = 0; i < posAdd; i++)
			{
				xm[offset] = R_MINUS_ONE;
				offset += permutationStride;
			}
			xm[offset]--;
		}
	}
}
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_permutated_ptx_plain_y (usfixn64 * __restrict__ xm,
																					 const usfixn64 * __restrict__ ym,
																					 const usfixn64 permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1 = 0;
//	, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
//	usfixn64 pos = threadIdx.x;

	usfixn64 offset = 0;
//	short posAdd = -1;
	usfixn32 posAdd = -1;
	usfixn32 bit;
	usfixn32 bitsFlag;
//	usfixn64 tmp = 0;

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
		num1 = xm[offset] + ym[i] + c;
		if (num1 < xm[offset] || num1 < ym[i])	//there is overflow/truncation
		{
//			xm[offset] = num1 + RC;
			num1 += RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
//			xm[offset] = num1 - R;
			num1 -= R;
		}
		else
		{
//			xm[offset] = num1;
			c = 0;
		}
		posAdd = (posAdd == -1 && num1 > 0) * i;

		xm[offset] = num1;
//		bit = (num1 > 0);
//		bitsFlag |= (bit << i);
		offset += permutationStride;
	}

//	asm("brev.b32 %0,%0;":"+r"(bitsFlag):);
//	asm("bfind.u32 %0, %1;":"=r"(posAdd):"r"(bitsFlag));
//	if (posAdd == 0xFFFFFFFF)
//		posAdd = -1;
//	else
//		posAdd = 31 - posAdd;
//
//	offset = 0;

	offset = 0;
	if (c > 0)
	{

		if (posAdd == -1)
		{
			xm[offset] = ULMAX;
			offset += permutationStride;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
		else
		{
			for (i = 0; i < posAdd; i++)
			{
				xm[offset] = R_MINUS_ONE;
				offset += permutationStride;
			}
			xm[offset]--;
		}
	}
}
/******************************************************************************/
__device__ __forceinline__ void
device_bigPrimeAdd_permutated_ptx_steps (usfixn64 * __restrict__ xm,
																				 const usfixn64 * __restrict__ ym,
																				 const usfixn64 permutationStride,
																				 usfixn64 steps)
{
	unsigned short c = 0;
	usfixn64 num1 = 0;
//	, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
//	usfixn64 pos = threadIdx.x;

	usfixn64 offset = 0;
//	short posAdd = -1;
	usfixn32 posAdd = -1;
	usfixn32 bit;
	usfixn32 bitsFlag = 0;
//	usfixn64 tmp = 0;

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < steps; i++)
	{
		num1 = xm[offset];
		bit = (num1 > 0);
		bitsFlag |= (bit << i);
		offset += permutationStride;
	}

	for (i = steps; i < COEFFICIENT_SIZE; i++)
	{
//		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
		num1 = xm[offset] + ym[offset] + c;
		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
		{
//			xm[offset] = num1 + RC;
			num1 += RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
//			xm[offset] = num1 - R;
			num1 -= R;
		}
		else
		{
//			xm[offset] = num1;
			c = 0;
		}
//		posAdd = (posAdd == -1 && num1 > 0) * i;

		xm[offset] = num1;
//		bit = (num1 > 0);
//		bitsFlag |= (bit << i);
		offset += permutationStride;
		posAdd = (posAdd == -1 && num1 > 0) * i;
	}

//	for (i = steps; i < COEFFICIENT_SIZE; i++)
//	{
////		num1 = xm[i] + ym[i] + c;
////		num1 = xm[pos] + ym[pos] + c;
//		num1 = xm[offset] + 0 + c;
//		if (num1 < xm[offset] )	//there is overflow/truncation
//		{
////			xm[offset] = num1 + RC;
//			num1 += RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
////			xm[offset] = num1 - R;
//			num1 -= R;
//		}
//		else
//		{
////			xm[offset] = num1;
//			c = 0;
//		}
////		posAdd = (posAdd == -1 && num1 > 0) * i;
//
//		xm[offset] = num1;
//		bit = (num1 > 0);
//		bitsFlag |= (bit << i);
//		offset += permutationStride;
//	}

//	asm("brev.b32 %0,%0;":"+r"(bitsFlag):);
//	asm("bfind.u32 %0, %1;":"=r"(posAdd):"r"(bitsFlag));
//	if (posAdd == 0xFFFFFFFF)
//		posAdd = -1;
//	else
//		posAdd = 31 - posAdd;

	offset = 0;

	offset = 0;
	if (c > 0)
	{

		if (posAdd == -1)
		{
			xm[offset] = ULMAX;
			offset += permutationStride;
			for (i = 1; i < COEFFICIENT_SIZE; i++)
			{
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
		else
		{
			for (i = 0; i < posAdd; i++)
			{
				xm[offset] = R_MINUS_ONE;
				offset += permutationStride;
			}
			xm[offset]--;
		}
	}
}

/******************************************************************************/
/******************************************************************************/
__device__ __inline__ void
device_bigPrimeAdd_permutated_x_plain_y_inPlace (usfixn64 * __restrict__ xm,
																								 const usfixn64 * __restrict__ ym,
																								 usfixn64 permutationStride)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

	usfixn64 offset = 0;
//	num1 = 0;
//	num2 = R - 1;
	short i;
	usfixn64 offset2 = permutationStride << 1;
	usfixn64 offset3 = permutationStride + offset2;
	usfixn64 offset4 = permutationStride << 2;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		num1 = xm[offset] + ym[i] + c;
		if (num1 < xm[offset] || num1 < ym[i]) //there is overflow/truncation
		{
			xm[offset] = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
			xm[offset] = num1 - R;
		}
		else
		{
			xm[offset] = num1;
			c = 0;
		}
		offset += permutationStride;
	}

	offset = 0;
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
//			xm[1]--;
			xm[permutationStride]--;
		return;
	}
	if (xm[2] != 0 && pos1 == -1)
	{
		pos1 = 2;
		xm[0] = R_MINUS_ONE;

		//xm[1]
		xm[permutationStride] = R_MINUS_ONE; //1

		//xm[2]
		xm[offset2]--;
		return;
	}
	if (xm[3] != 0 && pos1 == -1)
	{
		pos1 = 3;
		xm[0] = R_MINUS_ONE;

//			xm[1] = R_MINUS_ONE;
		xm[permutationStride] = R_MINUS_ONE;

//			xm[2] = R_MINUS_ONE;
		xm[offset2] = R_MINUS_ONE;

//			xm[3]--;
		xm[offset3]--;
		return;
	}
	if (xm[4] != 0 && pos1 == -1)
	{
		pos1 = 4;
		//xm[0]
		xm[0] = R_MINUS_ONE;
		//xm[1]
		xm[permutationStride] = R_MINUS_ONE;
		//xm[2]
		xm[offset2] = R_MINUS_ONE;
		//xm[3]
		xm[offset3] = R_MINUS_ONE;
		//xm[4]
		xm[offset4]--;
		return;
	}
	if (xm[5] != 0 && pos1 == -1)
	{
		pos1 = 5;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride]--;
		return;
	}
	if (xm[6] != 0 && pos1 == -1)
	{
		pos1 = 6;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2]--;
		return;
	}

	if (xm[7] != 0 && pos1 == -1)
	{
		pos1 = 7;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3]--;
		return;
	}

	if (xm[8] != 0 && pos1 == -1)
	{
		pos1 = 8;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3]--;
		return;
	}

	if (xm[9] != 0 && pos1 == -1)
	{
		pos1 = 9;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)]--;
		return;
	}

	if (xm[10] != 0 && pos1 == -1)
	{
		pos1 = 10;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)]--;
		return;
	}

	if (xm[11] != 0 && pos1 == -1)
	{
		pos1 = 11;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)] = R_MINUS_ONE;
		//11
		xm[offset3 + (permutationStride << 3)]--;
		return;
	}

	if (xm[12] != 0 && pos1 == -1)
	{
		pos1 = 12;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)] = R_MINUS_ONE;
		//11
		xm[offset3 + (permutationStride << 3)] = R_MINUS_ONE;
		//12
		xm[offset3 << 2]--;

		return;
	}

	if (xm[13] != 0 && pos1 == -1)
	{
		pos1 = 13;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)] = R_MINUS_ONE;
		//11
		xm[offset3 + (permutationStride << 3)] = R_MINUS_ONE;
		//12
		xm[offset3 << 2] = R_MINUS_ONE;
		//13
		xm[permutationStride + (offset3 << 2)]--;
		return;
	}

	if (xm[14] != 0 && pos1 == -1)
	{
		pos1 = 14;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)] = R_MINUS_ONE;
		//11
		xm[offset3 + (permutationStride << 3)] = R_MINUS_ONE;
		//12
		xm[offset3 << 2] = R_MINUS_ONE;
		//13
		xm[permutationStride + (offset3 << 2)] = R_MINUS_ONE;
		//14
		xm[offset2 + (offset3 << 2)]--;
		return;
	}

	if (xm[15] != 0 && pos1 == -1)
	{
		pos1 = 15;
		//0
		xm[0] = R_MINUS_ONE;
		//1
		xm[permutationStride] = R_MINUS_ONE;
		//2
		xm[offset2] = R_MINUS_ONE;
		//3
		xm[offset3] = R_MINUS_ONE;
		//4
		xm[offset4] = R_MINUS_ONE;
		//5
		xm[offset4 + permutationStride] = R_MINUS_ONE;
		//6
		xm[offset4 + offset2] = R_MINUS_ONE;
		//7
		xm[offset4 + offset3] = R_MINUS_ONE;
		//8
		xm[permutationStride << 3] = R_MINUS_ONE;
		//9
		xm[permutationStride + (permutationStride << 3)] = R_MINUS_ONE;
		//10
		xm[offset2 + (permutationStride << 3)] = R_MINUS_ONE;
		//11
		xm[offset3 + (permutationStride << 3)] = R_MINUS_ONE;
		//12
		xm[offset3 << 2] = R_MINUS_ONE;
		//13
		xm[permutationStride + (offset3 << 2)] = R_MINUS_ONE;
		//14
		xm[offset2 + (offset3 << 2)] = R_MINUS_ONE;
		//15
		xm[offset3 + (offset3 << 2)]--;
		return;
	}
//else (c>0) but (pos ==-1)
	{
		//0
		xm[0] = ULMAX;
		//1
		xm[permutationStride] = 0;
		//2
		xm[offset2] = 0;
		//3
		xm[offset3] = 0;
		//4
		xm[offset4] = 0;
		//5
		xm[offset4 + permutationStride] = 0;
		//6
		xm[offset4 + offset2] = 0;
		//7
		xm[offset4 + offset3] = 0;
		//8
		xm[permutationStride << 3] = 0;
		//9
		xm[permutationStride + (permutationStride << 3)] = 0;
		//10
		xm[offset2 + (permutationStride << 3)] = 0;
		//11
		xm[offset3 + (permutationStride << 3)] = 0;
		//12
		xm[offset3 << 2] = 0;
		//13
		xm[permutationStride + (offset3 << 2)] = 0;
		//14
		xm[offset2 + (offset3 << 2)] = 0;
		//15
		xm[offset3 + (offset3 << 2)] = 0;
	}
}
}
/******************************************************************************/
__device__ inline void
device_bigPrimeAdd_permutated_5 (usfixn64 * __restrict__ xm,
																 usfixn64 * __restrict__ ym,
																 usfixn64 * __restrict__ old_xm,
																 const usfixn64 permutationStride)
{
unsigned short c = 0;
usfixn64 num1 = 0;
//	, num2;

//	num1 = 0;
//	num2 = R - 1;
short i;
//	usfixn64 pos = threadIdx.x;

usfixn64 offset = 0;
short posAdd = -1;
//	usfixn64 tmp = 0;

offset = 0;
//#pragma unroll COEFFICIENT_SIZE
for (i = 0; i < COEFFICIENT_SIZE; i++)
{
//		num1 = xm[i] + ym[i] + c;
//		num1 = xm[pos] + ym[pos] + c;
	old_xm[i] = xm[offset];
	num1 = xm[offset] + ym[offset] + c;
	if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
	{
//			xm[offset] = num1 + RC;
		num1 += RC;
		c = 1;
	}
	else if (num1 >= R)
	{
		c = 1;
//			xm[offset] = num1 - R;
		num1 -= R;
	}
	else
	{
//			xm[offset] = num1;
		c = 0;
	}
	posAdd = (posAdd == -1 && num1 > 0) * i;

	xm[offset] = num1;
	offset += permutationStride;
}

offset = 0;
if (c > 0 && posAdd >= 0)

{
//		posAdd = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[offset] != 0)
//			{
//				posAdd = i;
//				break;
//			}
//		posAdd = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[offset] != 0)
//			{
//				posAdd = i;
//				break;
//			}
//			offset += permutationStride;

	// shouldn't it be >0?
	{
		offset = 0;
		for (i = 0; i < posAdd; i++)
		{
//				xm[offset] = num2;
			xm[offset] = R_MINUS_ONE;
			offset += permutationStride;
		}
//			offset = posAdd * permutationStride;
		xm[offset]--;
	}
}
if (c > 0 && posAdd < 0)
{
//			um[0] = ULMAX;
	xm[0] = ULMAX;
	offset = permutationStride;
//#pragma unroll COEFFICIENT_SIZE-1
	for (i = 1; i < 8; i++)
	{
//				um[i] = 0;
		xm[offset] = 0;
	}
}
}
/******************************************************************************/
//xm = xm + ym
//ym = xm - ym
//__device__ inline void fft_base2_permutated(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 xIdx,const usfixn64 yIdx,
//		const usfixn64 permutationStride)
__device__ __inline__ void
device_p4_fft_base2_permutated (usfixn64 * xm, usfixn64 * ym,
											const usfixn64 & permutationStride)
{
short c = 0;
usfixn64 num1, num2;
unsigned short c2 = 0;
usfixn64 num3 = 0;

num1 = 0;
num2 = R - 1;
short i;
//	usfixn64 pos = threadIdx.x;

usfixn64 offset = 0;
short pos1, pos = 0;
//	usfixn64 tmp = 0;
usfixn64 tAdd = 0;
offset = 0;

//#pragma unroll COEFFICIENT_SIZE
for (i = 0; i < COEFFICIENT_SIZE; i++)
{
	num1 = xm[offset] + ym[offset] + c;
	num3 = ym[offset] + c2;

	//addition part
	if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
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

	//subtraction part
	if (xm[offset] < num3) //there is not enough to do subtraction
	{
		c2 = 1;
		num3 = R - num3 + xm[offset];
	}
	else
	{
		c2 = 0;
		num3 = xm[offset] - num3;
	}
	xm[offset] = num1;
	ym[offset] = num3;
	offset += permutationStride;
	__syncthreads ();
}
//	return;

offset = 0;
if (c > 0)
{
	pos1 = -1;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		if (xm[offset] != 0)
		{
			pos1 = i;
			break;
		}
		offset += permutationStride;
	}
	if (pos1 >= 0)	// shouldn't it be >0?
	{
		offset = 0;
		for (i = 0; i < pos1; i++)
		{
			xm[offset] = num2;
			offset += permutationStride;
			__syncthreads ();
		}
		offset = pos1 * permutationStride;
		xm[offset]--;
		//xm[pos1*permutationStride+idx]--;
	}
	else
	{
		//			xm[0] = ULMAX;
		xm[offset] = ULMAX;
		offset = permutationStride;
		for (i = 1; i < COEFFICIENT_SIZE; i++)
		{
			//				xm[i] = 0;
			xm[offset] = 0;
			__syncthreads ();
		}
	}
}

if (c2 > 0)
{
	offset = 0;
	pos = -1;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		if (ym[offset] <= num2)
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
			ym[offset] = 0;
			offset += permutationStride;
			__syncthreads ();
		}
		offset = pos * permutationStride;
		//			um[pos]++;
		ym[offset]++;
		__syncthreads ();
	}
	else
	{

		offset = 0;
		//			um[0] = ULMAX;
		ym[offset] = ULMAX;
		offset = permutationStride;
		for (i = 1; i < COEFFICIENT_SIZE; i++)
		{
			//				um[i] = 0;
			ym[offset] = 0;
			offset += permutationStride;
			__syncthreads ();
		}
	}
}
}
/******************************************************************************/
/******************************************************************************/
//xm = xm + ym
//ym = xm - ym
//__device__ inline void fft_base2_permutated(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 xIdx,const usfixn64 yIdx,
//		const usfixn64 permutationStride)
__device__ __forceinline__ void
device_p4_fft_base2_permutated_v1 (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
												 const usfixn64 & permutationStride)
{
usfixn32 c = 0;
usfixn64 num1, num2;
//	unsigned short c2 = 0;
usfixn32 c2 = 0;
usfixn64 num3 = 0;

num1 = 0;
num2 = R - 1;
short i;
//	usfixn64 pos = threadIdx.x;

usfixn64 offset = 0;
usfixn32 pos1, pos = 0;
//	usfixn64 tmp = 0;
usfixn64 tAdd = 0;
offset = 0;

usfixn32 bitmapAdd = 0;
usfixn32 bitmapSub = 0;

usfixn64 x[COEFFICIENT_SIZE];
offset = 0;
for (i = 0; i < COEFFICIENT_SIZE; i++)
{
	x[i] = xm[offset];
	offset += permutationStride;
}

offset = 0;
//	for (i = 0; i < 8; i++)
//	{
//
//		num1 = c;
////		num1 = xm[offset] + ym[offset] + c;
//		asm("{\n\t"
//				"add.u64 %0,%0,%2;\n\t"
//				"add.cc.u64 %0,%0,%3;\n\t"
//				"addc.u32 %1,0,0;\n\t"
//				"}"
//				:"+l"(num1),"=r"(c):"l"(xm[offset]),"l"(ym[offset]));
//
//		//addition part
////		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
////		{
////			num1 = num1 + RC;
////			c = 1;
////		}
////		else if (num1 >= R)
////		{
////			c = 1;
////			num1 = num1 - R;
////		}
////		else
////		{
////			c = 0;
////		}
//
//		if (c)	//there is overflow/truncation
//		{
//			num1 = num1 + RC;
////					c = 1;
//		}
//		else if (num1 >= R)
//		{
////					c = 1;
//			num1 = num1 - R;
//		}
//		//		else
//		//		{
//		//			c = 0;
//		//		}
//
//		num3 = ym[offset] + c2;
//
////		num3 = R-ym[offset]-c2;
////		num3=c2;
////		asm("{\n\t"
////				"add.cc.u64 %0,%0,%2;\n\t"
////				"addc.u32 %1,0,0;\n\t"
////				"}":"+l"(num3),"=r"(c2):"l"(xm[offset]),"l"(ym[offset]));
//		//subtraction part
////		if (xm[offset] < num3) //there is not enough to do subtraction
////		{
////			c2 = 1;
////			num3 = R - num3 + xm[offset];
////		}
////		else
////		{
////			c2 = 0;
////			num3 = xm[offset] - num3;
////		}
//
//		if (x[i] < num3) //there is not enough to do subtraction
//				{
//					c2 = 1;
//					num3 = R - num3 + x[i];
//				}
//				else
//				{
//					c2 = 0;
//					num3 = x[i] - num3;
//				}
//
////		printf("c2=%d\n",c2);
//
////		if (xm[offset] < num3) //there is not enough to do subtraction
////				{
////					c2 = 1;
////					num3 = R - num3 + xm[offset];
////				}
////				else
////				{
////					c2 = 0;
////					num3 = xm[offset] - num3;
////				}
//		bitmapAdd |= ((num1 > 0) << i);
//		bitmapSub |= ((num3 < R - 1) << i);
//		xm[offset] = num1;
//		ym[offset] = num3;
//		offset += permutationStride;
//	}

offset = 0;
for (i = 0; i < 8; i++)
{

	num1 = c;

	asm("{\n\t"
			"add.u64 %0,%0,%2;\n\t"
			"add.cc.u64 %0,%0,%3;\n\t"
			"addc.u32 %1,0,0;\n\t"
			"}"
			:"+l"(num1),"=r"(c):"l"(x[i]),"l"(ym[offset]));

	//addition part

//				num1 = xm[offset] + ym[offset] + c;
//				if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//				{
//					num1 = num1 + RC;
//					c = 1;
//				}
//				else if (num1 >= R)
//				{
//					c = 1;
//					num1 = num1 - R;
//				}
//				else
//				{
//					c = 0;
//				}

	if (c)	//there is overflow/truncation
	{
		num1 = num1 + RC;
		//					c = 1;
	}
	else if (num1 >= R)
	{
		//					c = 1;
		num1 = num1 - R;
	}
	bitmapAdd |= ((num1 > 0) << i);
//			bitmapSub |= ((num3 < R - 1) << i);
	xm[offset] = num1;
//			ym[offset] = num3;
	offset += permutationStride;
}

offset = 0;
for (i = 0; i < 8; i++)
{

//			num1 = c;
	//		num1 = xm[offset] + ym[offset] + c;
//			asm("{\n\t"
//					"add.u64 %0,%0,%2;\n\t"
//					"add.cc.u64 %0,%0,%3;\n\t"
//					"addc.u32 %1,0,0;\n\t"
//					"}"
//					:"+l"(num1),"=r"(c):"l"(xm[offset]),"l"(ym[offset]));

	//addition part
	//		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
	//		{
	//			num1 = num1 + RC;
	//			c = 1;
	//		}
	//		else if (num1 >= R)
	//		{
	//			c = 1;
	//			num1 = num1 - R;
	//		}
	//		else
	//		{
	//			c = 0;
	//		}

//			if (c)	//there is overflow/truncation
//			{
//				num1 = num1 + RC;
//	//					c = 1;
//			}
//			else if (num1 >= R)
//			{
//	//					c = 1;
//				num1 = num1 - R;
//			}
	//		else
	//		{
	//			c = 0;
	//		}

	num3 = ym[offset] + c2;

	//		num3 = R-ym[offset]-c2;
	//		num3=c2;
	//		asm("{\n\t"
	//				"add.cc.u64 %0,%0,%2;\n\t"
	//				"addc.u32 %1,0,0;\n\t"
	//				"}":"+l"(num3),"=r"(c2):"l"(xm[offset]),"l"(ym[offset]));
	//subtraction part
	//		if (xm[offset] < num3) //there is not enough to do subtraction
	//		{
	//			c2 = 1;
	//			num3 = R - num3 + xm[offset];
	//		}
	//		else
	//		{
	//			c2 = 0;
	//			num3 = xm[offset] - num3;
	//		}

	if (x[i] < num3) //there is not enough to do subtraction
	{
		c2 = 1;
		num3 = R - num3 + x[i];
	}
	else
	{
		c2 = 0;
		num3 = x[i] - num3;
	}

	//		printf("c2=%d\n",c2);

	//		if (xm[offset] < num3) //there is not enough to do subtraction
	//				{
	//					c2 = 1;
	//					num3 = R - num3 + xm[offset];
	//				}
	//				else
	//				{
	//					c2 = 0;
	//					num3 = xm[offset] - num3;
	//				}
//			bitmapAdd |= ((num1 > 0) << i);
	bitmapSub |= ((num3 < (R_MINUS_ONE)) << i);
//			xm[offset] = num1;
	ym[offset] = num3;
	offset += permutationStride;
}

asm("brev.b32 %0,%0;":"+r"(bitmapAdd):);
asm("brev.b32 %0,%0;":"+r"(bitmapSub):);
asm("bfind.u32 %0, %1;":"=r"(pos1):"r"(bitmapAdd));
asm("bfind.u32 %0, %1;":"=r"(pos):"r"(bitmapSub));
if (pos1 == 0xFFFFFFFF)
	pos1 = -1;
else
	pos1 = 31 - pos1;
//	pos = 31 - pos;

if (pos == 0xFFFFFFFF)
	pos = -1;
else
	pos = 31 - pos;

offset = 0;
if (c > 0)
{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[offset] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//			offset += permutationStride;
//		}
//		if (pos1 < 32)	// shouldn't it be >0?
	if (pos1 >= 0)
	{
		offset = 0;
		for (i = 0; i < pos1; i++)
		{
			xm[offset] = num2;
			offset += permutationStride;
		}
//			offset = 0;
//			if (pos1 > 0)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 1)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 2)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 3)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 4)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 5)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			if (pos1 > 6)
//			{
//				xm[offset] = R - 1;
//				offset += permutationStride;
//			}
//			offset = pos1 * permutationStride;
		xm[offset]--;
		//xm[pos1*permutationStride+idx]--;
	}
	else
	{
		//			xm[0] = ULMAX;
		xm[offset] = ULMAX;
		offset = permutationStride;
		for (i = 1; i < 8; i++)
		{
			//				xm[i] = 0;
			xm[offset] = 0;
			offset += permutationStride;
		}
	}
}

if (c2 > 0)
{
//		offset = 0;
//		pos = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (ym[offset] < num2)
//			{
//				pos = i;
//				break;
//			}
//			offset += permutationStride;
//		}

	if (pos >= 0)
//		if (pos < 32)
	{
		offset = 0;
		for (i = 0; i < pos; i++)
		{
			ym[offset] = 0;
			offset += permutationStride;
		}

//			if (pos > 0)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if (pos > 1)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if (pos > 2)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if (pos > 3)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if (pos > 4)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if (pos > 5)
//			{
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//			if(pos>6)
//						{
//							ym[offset]=0;
//							offset+=permutationStride;
//						}
//			offset = pos * permutationStride;
		//			um[pos]++;
		ym[offset]++;
	}
	else
	{

		offset = 0;
		//			um[0] = ULMAX;
		ym[offset] = ULMAX;
		offset = permutationStride;

		for (i = 1; i < 8; i++)
		{
			//				um[i] = 0;
			ym[offset] = 0;
			offset += permutationStride;
		}
	}
}
}
/******************************************************************************/
__device__ __inline__ void
device_p4_fft_base2_plain_inPlace (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
												 const usfixn64 & permutationStride)
{
unsigned short c = 0;
usfixn64 num1, num2;
unsigned short c2 = 0;
usfixn64 num3 = 0;

num1 = 0;
//	num2 = R - 1;
short i;
//	usfixn64 pos = threadIdx.x;

usfixn64 offset = 0;
short pos = 0;
//	usfixn64 tmp = 0;
//	usfixn64 tAdd = 0;
offset = 0;
uConstArray8_align8 xs, ys;

//#pragma unroll COEFFICIENT_SIZE
for (i = 0; i < 8; i++)
{
	xs.i[i] = xm[offset];
	ys.i[i] = ym[offset];
//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;
	num1 = xs.i[i] + ys.i[i] + c;
	num3 = ys.i[i] + c2;

	//addition part
//		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
	if (num1 < xs.i[i] || num1 < ys.i[i])	//there is overflow/truncation
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

	//subtraction part
//		if (xm[offset] < num3) //there is not enough to do subtraction
	if (xs.i[i] < num3) //there is not enough to do subtraction
	{
		c2 = 1;
//			num3 = R - num3 + xm[offset];
		num3 = R - num3 + xs.i[i];
	}
	else
	{
		c2 = 0;
//			num3 = xm[offset] - num3;
		num3 = xs.i[i] - num3;
	}
//		xm[offset] = num1;
//		ym[offset] = num3;
	xs.i[i] = num1;
	ys.i[i] = num3;
	offset += permutationStride;
}

offset = 0;
if (c > 0)
{
	pos = -1;
//		i=0;
	if (xs.i[0] != 0 && pos == -1)
	{
		pos = 0;
		xs.i[0]--;
	}

	//		i=1;
	if (xs.i[1] != 0 && pos == -1)
	{
		pos = 1;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1]--;
	}
	//		i=2;
	if (xs.i[2] != 0 && pos == -1)
	{
		pos = 2;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2]--;
	}
	//		i=3;
	if (xs.i[3] != 0 && pos == -1)
	{
		pos = 3;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2] = R_MINUS_ONE;
		xs.i[3]--;
	}
	//		i=4;
	if (xs.i[4] != 0 && pos == -1)
	{
		pos = 3;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2] = R_MINUS_ONE;
		xs.i[3] = R_MINUS_ONE;
		xs.i[4]--;
	}
	//		i=5;
	if (xs.i[5] != 0 && pos == -1)
	{
		pos = 5;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2] = R_MINUS_ONE;
		xs.i[3] = R_MINUS_ONE;
		xs.i[4] = R_MINUS_ONE;
		xs.i[5]--;
	}
	//		i=6;
	if (xs.i[6] != 0 && pos == -1)
	{
		pos = 6;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2] = R_MINUS_ONE;
		xs.i[3] = R_MINUS_ONE;
		xs.i[4] = R_MINUS_ONE;
		xs.i[5] = R_MINUS_ONE;
		xs.i[6]--;
	}
	//		i=7;
	if (xs.i[7] != 0 && pos == -1)
	{
		pos = 7;
		xs.i[0] = R_MINUS_ONE;
		xs.i[1] = R_MINUS_ONE;
		xs.i[2] = R_MINUS_ONE;
		xs.i[3] = R_MINUS_ONE;
		xs.i[4] = R_MINUS_ONE;
		xs.i[5] = R_MINUS_ONE;
		xs.i[6] = R_MINUS_ONE;
		xs.i[7]--;
	}
	if (pos == -1)
	{
		xs.i[0] = ULMAX;
		xs.i[1] = 0;
		xs.i[2] = 0;
		xs.i[3] = 0;
		xs.i[4] = 0;
		xs.i[5] = 0;
		xs.i[6] = 0;
		xs.i[7] = 0;
	}
}

if (c2 > 0)
{
	pos = -1;
	//		i=0;
	if (ys.i[0] != 0 && pos == -1)
	{
		pos = 0;
		ys.i[0]++;
	}

	//		i=1;
	if (ys.i[1] != 0 && pos == -1)
	{
		pos = 1;
		ys.i[0] = 0;
		ys.i[1]++;
	}
	//		i=2;
	if (ys.i[2] != 0 && pos == -1)
	{
		pos = 2;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2]++;
	}
	//		i=3;
	if (ys.i[3] != 0 && pos == -1)
	{
		pos = 3;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3]++;
	}
	//		i=4;
	if (ys.i[4] != 0 && pos == -1)
	{
		pos = 4;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3] = 0;
		ys.i[4]++;
	}
	//		i=5;
	if (ys.i[5] != 0 && pos == -1)
	{
		pos = 5;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3] = 0;
		ys.i[4] = 0;
		ys.i[5]++;
	}
	//		i=6;
	if (ys.i[6] != 0 && pos == -1)
	{
		pos = 6;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3] = 0;
		ys.i[4] = 0;
		ys.i[5] = 0;
		ys.i[6]++;
	}
	//		i=7;
	if (ys.i[7] != 0 && pos == -1)
	{
		pos = 7;
		ys.i[0] = 0;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3] = 0;
		ys.i[4] = 0;
		ys.i[5] = 0;
		ys.i[6] = 0;
		ys.i[7]++;
	}
	if (pos == -1)
	{
		ys.i[0] = ULMAX;
		ys.i[1] = 0;
		ys.i[2] = 0;
		ys.i[3] = 0;
		ys.i[4] = 0;
		ys.i[5] = 0;
		ys.i[6] = 0;
		ys.i[7] = 0;
	}
}

offset = 0;
for (i = 0; i < 8; i++)
{
	xm[offset] = xs.i[i];
	ym[offset] = ys.i[i];
	offset += permutationStride;
}
}
/******************************************************************************/

/******************************************************************************/
//xm = xm + ym
//ym = xm - ym
//__device__ inline void fft_base2_permutated(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 xIdx,const usfixn64 yIdx,
//		const usfixn64 permutationStride)
//__device__ __inline__ void fft_base2_vector8_permutated(
//		usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
//		const usfixn64 & permutationStride)
//{
//	short c = 0;
//	usfixn64 num1, num2;
//	unsigned short c2 = 0;
//	usfixn64 num3 = 0;
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
////	usfixn64 pos = threadIdx.x;
//
//	usfixn64 offset = 0;
//	short pos1, pos = 0;
////	usfixn64 tmp = 0;
//	usfixn64 tAdd = 0;
//	offset = 0;
//
//	uvector8 ys, xs;
//
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////
////		setUvector8Element(xs, i, xm[offset]);
////		setUvector8Element(ys, i, ym[offset]);
////	}
////
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////
//////		setUvector8Element(xs, i, xm[offset]);
////		setUvector8Element(ys, i, ym[offset]);
////	}
//
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
//
//		setUvector8Element(xs, i, xm[offset]);
//		setUvector8Element(ys, i, ym[offset]);
//
////		num1 = xm[offset] + ym[offset] + c;
////		num3 = ym[offset] + c2;
//
//		num1 = getUvector8Element(xs, i) + getUvector8Element(ys, i) + c;
//		num3 = getUvector8Element(ys, i) + c2;
//		//addition part
////		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//		if (num1 < getUvector8Element(xs, i) || num1 < getUvector8Element(ys, i))	//there is overflow/truncation
//		{
//			num1 = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			num1 = num1 - R;
//		}
//		else
//		{
//			c = 0;
//		}
//
//		//subtraction part
//		if (getUvector8Element(xs, i) < num3) //there is not enough to do subtraction
//		{
//			c2 = 1;
////			num3 = R - num3 + xm[offset];
//			num3 = R - num3 + getUvector8Element(xs, i);
//		}
//		else
//		{
//			c2 = 0;
////			num3 = xm[offset] - num3;
//			num3 = getUvector8Element(xs, i) - num3;
//		}
//
////		xm[offset] = num1;
////		ym[offset] = num3;
//		setUvector8Element(xs, i, num1);
//		setUvector8Element(ys, i, num3);
//		offset += permutationStride;
//	}
//
//	offset = 0;
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
////			if (xm[offset] != 0)
//			if (getUvector8Element(xs, i) != 0)
//			{
//				pos1 = i;
//				break;
//			}
////			offset += permutationStride;
//		}
//		if (pos1 >= 0)	// shouldn't it be >0?
//		{
////			offset = 0;
//			for (i = 0; i < pos1; i++)
//			{
////				xm[offset] = num2;
//				setUvector8Element(xs, i, num2);
////				offset += permutationStride;
//			}
////			offset = pos1 * permutationStride;
////			xm[offset]--;
////			setUvector8Element(xs, i, getUvector8Element(xs, pos1)-1);
//			decUvector8Element(xs, i, 1);
//			//xm[pos1*permutationStride+idx]--;
//		}
//		else
//		{
//			//			xm[0] = ULMAX;
////			xm[offset] = ULMAX;
//			setUvector8Element(xs, 0, ULMAX);
////			offset = permutationStride;
//			for (i = 1; i < 8; i++)
//			{
//				//				xm[i] = 0;
////				xm[offset] = 0;
//				setUvector8Element(xs, i, ZERO);
//			}
//		}
//	}
//
//	if (c2 > 0)
//	{
//		offset = 0;
//		pos = -1;
//		for (i = 0; i < 8; i++)
//		{
////			if (ym[offset] < num2)
//			if (getUvector8Element(xs, i) < num2)
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
////				ym[offset] = 0;
//				setUvector8Element(ys, i, ZERO);
//				offset += permutationStride;
//			}
//			offset = pos * permutationStride;
//			//			um[pos]++;
////			ym[offset]++;
////			setUvector8Element(ys,pos, getUvector8Element(ys,pos)+1);
//			incUvector8Element(ys, pos, 1);
//		}
//		else
//		{
////			offset = 0;
//			//			um[0] = ULMAX;
////			ym[offset] = ULMAX;
//			setUvector8Element(ys, 0, ULMAX);
////			offset = permutationStride;
//			for (i = 1; i < 8; i++)
//			{
//				//				um[i] = 0;
////				ym[offset] = 0;
//				setUvector8Element(ys, i, ZERO);
////				offset += permutationStride;
//			}
//		}
//	}
//
//	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		xm[offset] = getUvector8Element(xs, i);
//		ym[offset] = getUvector8Element(ys, i);
//		offset += permutationStride;
//	}
//}
/******************************************************************************/
//__device__ inline void fft_base2_permutated_2(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 idx,
//		const usfixn64 permutationStride)
//{
//	unsigned short c = 0;
//	usfixn64 num1, num2;
//	unsigned short c2 = 0;
//	usfixn64 num3 = 0;
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
////	usfixn64 pos = threadIdx.x;
//
//	usfixn64 offset = 0;
//	short pos = 0;
////	usfixn64 tmp = 0;
//	offset = idx;
//
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
////		num1 = xm[idx + offset] + ym[idx + offset] + c;
////		num3 = ym[idx + offset] + c2;
//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;
//
//		//addition part
////		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])	//there is overflow/truncation
//		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//		{
//			num1 = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			num1 = num1 - R;
//		}
//		else
//		{
//			c = 0;
//		}
//
//		//subtraction part
////		if (xm[idx + offset] < num3) //there is not enough to do subtraction
//		if (xm[offset] < num3) //there is not enough to do subtraction
//		{
//			c2 = 1;
////			num3 = R - num3 + xm[idx + offset];
//			num3 = R - num3 + xm[offset];
//		}
//		else
//		{
//			c2 = 0;
////			num3 = xm[idx + offset] - num3;
//			num3 = xm[offset] - num3;
//		}
////		xm[idx + offset] = num1;
////		ym[idx + offset] = num1;
//		xm[offset] = num1;
//		ym[offset] = num3;
//		offset += permutationStride;
//	}
//
////	offset = idx;
////	if (c > 0)
////	{
////		pos = -1;
////		for (i = 0; i < 8; i++)
////		{
//////			if (xm[idx+offset] != 0)
////			if (xm[offset] != 0)
////			{
////				pos = i;
////				break;
////			}
////			offset += permutationStride;
////		}
////		if (pos >= 0)	// shouldn't it be >0?
////		{
//////			offset = 0;
////			offset = idx;
////			for (i = 0; i < pos; i++)
////			{
//////				xm[idx + offset] = num2;
////				xm[offset] = num2;
////				offset += permutationStride;
////			}
//////			offset = pos * permutationStride;
//////			xm[idx + offset]--;
////			xm[offset]--;
////			//xm[pos*permutationStride+idx]--;
////		}
////		else
////		{
////			//			xm[0] = ULMAX;
////			xm[idx] = ULMAX;
//////			offset = permutationStride;
////			offset = idx+permutationStride;
////			for (i = 1; i < 8; i++)
////			{
////				//				xm[i] = 0;
//////				xm[idx + offset] = 0;
////				xm[offset] = 0;
////				offset+=permutationStride;
////			}
////		}
////	}
//
////	if (c2 > 0)
////	{
//////		offset = 0;
////		offset = idx;
////		pos = -1;
////		for (i = 0; i < 8; i++)
////		{
////			if (ym[idx + offset] < num2)
////			{
////				pos = i;
////				break;
////			}
////			offset += permutationStride;
////		}
////
////		if (pos >= 0)
////		{
////			offset = 0;
////			for (i = 0; i < pos; i++)
////			{
////				ym[idx + offset] = 0;
////				offset += permutationStride;
////			}
////			offset = pos * permutationStride;
////			//			um[pos]++;
////			ym[idx + offset]++;
////		}
////		else
////		{
////
////			//			um[0] = ULMAX;
////			ym[idx] = ULMAX;
////			offset = 0;
////			for (i = 1; i < 8; i++)
////			{
////				//				um[i] = 0;
////				ym[idx + offset];
////				offset += permutationStride;
////			}
////		}
////	}
//}
/******************************************************************************/

//using vectorized data structure + offset is computed more efficiently
//__device__ inline void fft_base2_permutated_3(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 idx,
//		const usfixn64 permutationStride)
//{
//	short c = 0;
//	usfixn64 num1, num2;
//	unsigned short c2 = 0;
//	usfixn64 num3 = 0;
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
////	usfixn64 pos = threadIdx.x;
//
//	usfixn64 offset = 0;
//	short posAdd = -1, posSub = -1;
////	usfixn64 tmp = 0;
//	offset = idx;
//
//	uConstArray8_align8 xConst, yConst;
//	offset = idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		xConst.i[i] = xm[offset];
////		offset += permutationStride;
////	}
////
////	offset = idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		yConst.i[i] = ym[offset];
////		offset += permutationStride;
////	}
//
////	offset=idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		xConst.i[i] = xm[offset];
//////	offset += permutationStride;
////		offset += STRIDE;
////	}
//
////	offset=idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		yConst.i[i] = ym[offset];
//////	offset += permutationStride;
////		offset += STRIDE;
////	}
//
//	offset = idx;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
////		num1 = xm[idx + offset] + ym[idx + offset] + c;
////		num3 = ym[idx + offset] + c2;
////		num1 = xm[offset] + ym[offset] + c;
////		num3 = ym[offset] + c2;
//		xConst.i[i] = xm[offset];
//		yConst.i[i] = ym[offset];
//		num1 = xConst.i[i] + yConst.i[i] + c;
//		num3 = yConst.i[i] + c2;
//		offset += STRIDE;
//		//addition part
////		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])	//there is overflow/truncation
////		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//		if (num1 < xConst.i[i] || num1 < yConst.i[i])	//there is overflow/truncation
//		{
//			num1 = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			num1 = num1 - R;
//		}
//		else
//		{
//			c = 0;
//		}
//
//		//subtraction part
////		if (xm[idx + offset] < num3) //there is not enough to do subtraction
////		if (xm[offset] < num3) //there is not enough to do subtraction
////		{
////			c2 = 1;
//////			num3 = R - num3 + xm[idx + offset];
////			num3 = R - num3 + xm[offset];
////		}
////		else
////		{
////			c2 = 0;
//////			num3 = xm[idx + offset] - num3;
////			num3 = xm[offset] - num3;
////		}
//
//		if (xConst.i[i] < num3) //there is not enough to do subtraction
//		{
//			c2 = 1;
//			//			num3 = R - num3 + xm[idx + offset];
//			num3 = R - num3 + xConst.i[i];
//		}
//		else
//		{
//			c2 = 0;
//			//			num3 = xm[idx + offset] - num3;
//			num3 = xConst.i[i] - num3;
//		}
////		xm[idx + offset] = num1;
////		ym[idx + offset] = num1;
////		xm[offset] = num1;
////		ym[offset] = num1;
////		offset += permutationStride;
//		xConst.i[i] = num1;
//		yConst.i[i] = num3;
//		posAdd = (posAdd == -1 && num1 > 0) * i;
//		posSub = (posSub == -1 && num2 > 0) * i;
//	}
//
////	offset = idx;
//	if (c > 0)
//	{
////		posAdd = -1;
////		for (i = 0; i < 8; i++)
//		{
////			if (xm[idx+offset] != 0)
////			if (xm[offset] != 0)
////			if (xConst.i[i] != 0)
////			{
////				posAdd = i;
////				break;
////			}
////			offset += permutationStride;
//		}
//		if (posAdd >= 0)	// shouldn't it be >0?
//		{
////			offset = 0;
////			offset = idx;
//			for (i = 0; i < posAdd; i++)
//			{
////				xm[idx + offset] = num2;
////				xm[offset] = num2;
//				xConst.i[i] = num2;
////				offset += permutationStride;
//			}
////			offset = pos * permutationStride;
////			xm[idx + offset]--;
////			xm[offset]--;
//			xConst.i[posAdd]--;
//			//xm[pos*permutationStride+idx]--;
//		}
//		else
//		{
//			//			xm[0] = ULMAX;
////			xm[idx] = ULMAX;
//			xConst.i[0] = ULMAX;
////			offset = permutationStride;
////			offset = idx+permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
//			for (i = 1; i < 8; i++)
//			{
//				//				xm[i] = 0;
////				xm[idx + offset] = 0;
////				xm[offset] = 0;
//				xConst.i[i] = 0;
////				offset+=permutationStride;
//			}
//		}
//	}
////
//	if (c2 > 0)
//	{
////		offset = 0;
////		offset = idx;
////		posSub = -1;
////		for (i = 0; i < 8; i++)
////		{
//////			if (ym[idx + offset] < num2)
////			if (yConst.i[i] < num2)
////			{
////				posSub = i;
////				break;
////			}
//////			offset += permutationStride;
////		}
//
//		if (posSub >= 0)
//		{
////			offset = 0;
//			for (i = 0; i < posSub; i++)
//			{
////				ym[idx + offset] = 0;
////				offset += permutationStride;
//				yConst.i[i] = 0;
//			}
////			offset = pos * permutationStride;
//			//			um[pos]++;
////			ym[idx + offset]++;
//			yConst.i[posSub]++;
//		}
//		else
//		{
//			//			um[0] = ULMAX;
////			ym[idx] = ULMAX;
//			yConst.i[0] = ULMAX;
////			offset = 0;
//#pragma unroll (COEFFICIENT_SIZE-1)
//			for (i = 1; i < 8; i++)
//			{
//				//				um[i] = 0;
////				ym[idx + offset]=0;
//				yConst.i[i] = 0;
////				offset += permutationStride;
//			}
//		}
//	}
//
////	offset = idx;
////	#pragma unroll COEFFICIENT_SIZE
////		for (i = 0; i < 8; i++)
////		{
////			xm[offset] = xConst.i[i];
//////			ym[offset] = yConst.i[i];
//////			offset+=permutationStride;
////			offset+=STRIDE;
////		}
//
////		offset = idx;
////	#pragma unroll COEFFICIENT_SIZE
////		for (i = 0; i < 8; i++)
////		{
////			ym[offset] = yConst.i[i];
//////			offset+=permutationStride;
////			offset+=STRIDE;
////		}
//
//}
/******************************************************************************/

//using vectorized data structure + offset is computed more efficiently
__device__ inline void
deviec_p4_fft_base2_permutated_4 (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
												const usfixn64 & permutationStride)
{
short c = 0;
usfixn64 num1;
//	num2;
short c2 = 0;
usfixn64 num3 = 0;

num1 = 0;
//	num2 = R - 1;
short i;
//	usfixn64 pos = threadIdx.x;

usfixn64 offset = 0;
short posAdd = -1, posSub = -1;
//	usfixn64 tmp = 0;
offset = 0;

//	usfixn64 xConst, yConst;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
//		xConst.i[i] = xm[offset];
//		offset += permutationStride;
//	}
//
//	offset = idx;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
//		yConst.i[i] = ym[offset];
//		offset += permutationStride;
//	}

//	offset=idx;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
//		xConst.i[i] = xm[offset];
////	offset += permutationStride;
//		offset += STRIDE;
//	}

//	offset=idx;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
//		yConst.i[i] = ym[offset];
////	offset += permutationStride;
//		offset += STRIDE;
//	}

offset = 0;
//#pragma unroll COEFFICIENT_SIZE
for (i = 0; i < 8; i++)
{
//		num1 = xm[idx + offset] + ym[idx + offset] + c;
//		num3 = ym[idx + offset] + c2;
//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;
//		xConst = xm[offset];
//		yConst = ym[offset];

//		num1 = xConst + yConst + c;
//		num3 = yConst + c2;

	num1 = xm[offset] + ym[offset] + c;
	num3 = ym[offset] + c2;
	//addition part
//		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])	//there is overflow/truncation
//		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
	if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
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

	//subtraction part
//		if (xm[idx + offset] < num3) //there is not enough to do subtraction
//		if (xm[offset] < num3) //there is not enough to do subtraction
//		{
//			c2 = 1;
////			num3 = R - num3 + xm[idx + offset];
//			num3 = R - num3 + xm[offset];
//		}
//		else
//		{
//			c2 = 0;
////			num3 = xm[idx + offset] - num3;
//			num3 = xm[offset] - num3;
//		}

	if (xm[offset] < num3) //there is not enough to do subtraction
	{
		c2 = 1;
		//			num3 = R - num3 + xm[idx + offset];
		num3 = R - num3 + xm[offset];
	}
	else
	{
		c2 = 0;
		//			num3 = xm[idx + offset] - num3;
		num3 = xm[offset] - num3;
	}
//		xm[idx + offset] = num1;
//		ym[idx + offset] = num1;
//		xm[offset] = num1;
//		ym[offset] = num1;
//		offset += permutationStride;
//		xConst=num1;
//		yConst=num3;
	posAdd = (posAdd == -1 && num1 > 0) * i;
	posSub = (posSub == -1 && num3 > 0) * i;
	xm[offset] = num1;
	ym[offset] = num3;
	offset += permutationStride;
}

//	offset = idx;
//	if (c > 0)
//	{
////		posAdd = -1;
////		for (i = 0; i < 8; i++)
//		{
////			if (xm[idx+offset] != 0)
////			if (xm[offset] != 0)
////			if (xConst.i[i] != 0)
////			{
////				posAdd = i;
////				break;
////			}
////			offset += permutationStride;
//		}
//		offset = idx;
//		if (posAdd >= 0)	// shouldn't it be >0?
//		{
////			offset = 0;
//
//			for (i = 0; i < posAdd; i++)
//			{
////				xm[idx + offset] = num2;
//				xm[offset] = num2;
////				xConst = num2;
//				offset += permutationStride;
//			}
////			offset = pos * permutationStride;
////			xm[idx + offset]--;
//			xm[offset]--;
////			xConst.i[posAdd]--;
//			//xm[pos*permutationStride+idx]--;
//		}
//		else
//		{
//			//			xm[0] = ULMAX;
////			xm[idx] = ULMAX;
////			xConst.i[0] = ULMAX;
//			xm[offset] = ULMAX;
//			offset += permutationStride;
////			offset = idx+permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
//			for (i = 1; i < 8; i++)
//			{
//				//				xm[i] = 0;
////				xm[idx + offset] = 0;
//				xm[offset] = 0;
////				xConst.i[i] = 0;
//				offset+=permutationStride;
//			}
//		}
//	}

offset = 0;
if (c > 0 && posAdd >= 0)
{
	for (i = 0; i < posAdd; i++)
	{
//			xm[offset] = R_MINUS_ONE;
		xm[offset] = R - 1;
		offset += permutationStride;
	}
	xm[offset]--;

}
else if (c > 0 && posAdd == -1)
{
	xm[offset] = ULMAX;
	offset += permutationStride;

//#pragma unroll (COEFFICIENT_SIZE-1)
	for (i = 1; i < 8; i++)
	{
		xm[offset] = 0;
		offset += permutationStride;
	}
}
//
//	if (c2 > 0)
//	{
////		offset = 0;
////		offset = idx;
////		posSub = -1;
////		for (i = 0; i < 8; i++)
////		{
//////			if (ym[idx + offset] < num2)
////			if (yConst.i[i] < num2)
////			{
////				posSub = i;
////				break;
////			}
//////			offset += permutationStride;
////		}
//
//		offset = idx;
//		if (posSub >= 0)
//		{
////			offset = 0;
//			for (i = 0; i < posSub; i++)
//			{
////				ym[idx + offset] = 0;
//				ym[offset] = 0;
//				offset += permutationStride;
////				yConst.i[i] = 0;
//			}
////			offset = pos * permutationStride;
//			//			um[pos]++;
////			ym[idx + offset]++;
////			yConst.i[posSub]++;
//			ym[offset]++;
//		}
//		else
//		{
//			//			um[0] = ULMAX;
////			ym[idx] = ULMAX;
////			yConst.i[0] = ULMAX;
//			ym[offset] = ULMAX;
//			offset += permutationStride;
////			offset = 0;
//#pragma unroll (COEFFICIENT_SIZE-1)
//			for (i = 1; i < 8; i++)
//			{
//				//				um[i] = 0;
////				ym[idx + offset]=0;
////				yConst.i[i] = 0;
//				ym[offset] = 0;
//				offset += permutationStride;
//			}
//		}
//	}

//		offset = 0;
//		offset = idx;
//		posSub = -1;
//		for (i = 0; i < 8; i++)
//		{
////			if (ym[idx + offset] < num2)
//			if (yConst.i[i] < num2)
//			{
//				posSub = i;
//				break;
//			}
////			offset += permutationStride;
//		}

offset = 0;
if (c2 > 0 && posSub >= 0)
{
//			offset = 0;
	for (i = 0; i < posSub; i++)
	{
//				ym[idx + offset] = 0;
		ym[offset] = 0;
		offset += permutationStride;
//				yConst.i[i] = 0;
	}
//			offset = pos * permutationStride;
	//			um[pos]++;
//			ym[idx + offset]++;
//			yConst.i[posSub]++;
	ym[offset]++;
}
else if (c2 > 0 && posSub == -1)
{
	//			um[0] = ULMAX;
//			ym[idx] = ULMAX;
//			yConst.i[0] = ULMAX;
	ym[offset] = ULMAX;
	offset += permutationStride;
//			offset = 0;
//#pragma unroll (COEFFICIENT_SIZE-1)
	for (i = 1; i < 8; i++)
	{
		//				um[i] = 0;
//				ym[idx + offset]=0;
//				yConst.i[i] = 0;
		ym[offset] = 0;
		offset += permutationStride;
	}
}

//	offset = idx;
//	#pragma unroll COEFFICIENT_SIZE
//		for (i = 0; i < 8; i++)
//		{
//			xm[offset] = xConst.i[i];
////			ym[offset] = yConst.i[i];
////			offset+=permutationStride;
//			offset+=STRIDE;
//		}

//		offset = idx;
//	#pragma unroll COEFFICIENT_SIZE
//		for (i = 0; i < 8; i++)
//		{
//			ym[offset] = yConst.i[i];
////			offset+=permutationStride;
//			offset+=STRIDE;
//		}

}
/******************************************************************************/
//
////using vectorized data structure + offset is computed more efficiently
//__device__ __forceinline__ void  device_fft2_permutated_5(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 permutationStride)
//{
//	usfixn64 idx = 0; //should fix idx
//
//	short c = 0;
//	usfixn64 num1;
////	num2;
//	unsigned short c2 = 0;
//	usfixn64 num3 = 0;
//
//	num1 = 0;
////	num2 = R - 1;
//	short i;
////	usfixn64 pos = threadIdx.x;
//
//	usfixn64 offset = 0;
//	short posAdd = -1, posSub = -1;
////	usfixn64 tmp = 0;
//	offset = idx;
//
//	usfixn64 xConst, yConst;
//	offset = idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		xConst.i[i] = xm[offset];
////		offset += permutationStride;
////	}
////
////	offset = idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		yConst.i[i] = ym[offset];
////		offset += permutationStride;
////	}
//
////	offset=idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		xConst.i[i] = xm[offset];
//////	offset += permutationStride;
////		offset += STRIDE;
////	}
//
////	offset=idx;
////#pragma unroll COEFFICIENT_SIZE
////	for (i = 0; i < 8; i++)
////	{
////		yConst.i[i] = ym[offset];
//////	offset += permutationStride;
////		offset += STRIDE;
////	}
//
//	offset = idx;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < 8; i++)
//	{
////		num1 = xm[idx + offset] + ym[idx + offset] + c;
////		num3 = ym[idx + offset] + c2;
////		num1 = xm[offset] + ym[offset] + c;
////		num3 = ym[offset] + c2;
////		xConst = xm[offset];
////		yConst = ym[offset];
//
////		num1 = xConst + yConst + c;
////		num3 = yConst + c2;
//
//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;
//		//addition part
////		if (num1 < xm[idx + offset] || num1 < ym[idx + offset])	//there is overflow/truncation
////		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//		if (num1 < xm[offset] || num1 < ym[offset])	//there is overflow/truncation
//		{
//			num1 = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			num1 = num1 - R;
//		}
//		else
//		{
//			c = 0;
//		}
//
//		//subtraction part
////		if (xm[idx + offset] < num3) //there is not enough to do subtraction
////		if (xm[offset] < num3) //there is not enough to do subtraction
////		{
////			c2 = 1;
//////			num3 = R - num3 + xm[idx + offset];
////			num3 = R - num3 + xm[offset];
////		}
////		else
////		{
////			c2 = 0;
//////			num3 = xm[idx + offset] - num3;
////			num3 = xm[offset] - num3;
////		}
//
//		if (xm[offset] < num3) //there is not enough to do subtraction
//		{
//			c2 = 1;
//			//			num3 = R - num3 + xm[idx + offset];
//			num3 = R - num3 + xm[offset];
//		}
//		else
//		{
//			c2 = 0;
//			//			num3 = xm[idx + offset] - num3;
//			num3 = xm[offset] - num3;
//		}
////		xm[idx + offset] = num1;
////		ym[idx + offset] = num1;
////		xm[offset] = num1;
////		ym[offset] = num1;
////		offset += permutationStride;
////		xConst=num1;
////		yConst=num3;
//		posAdd = (posAdd == -1 && num1 > 0) * i;
//		posSub = (posSub == -1 && num3 > 0) * i;
//		xm[offset] = num1;
//		ym[offset] = num3;
//		offset += STRIDE;
//	}
//
////	offset = idx;
////	if (c > 0)
////	{
//////		posAdd = -1;
//////		for (i = 0; i < 8; i++)
////		{
//////			if (xm[idx+offset] != 0)
//////			if (xm[offset] != 0)
//////			if (xConst.i[i] != 0)
//////			{
//////				posAdd = i;
//////				break;
//////			}
//////			offset += permutationStride;
////		}
////		offset = idx;
////		if (posAdd >= 0)	// shouldn't it be >0?
////		{
//////			offset = 0;
////
////			for (i = 0; i < posAdd; i++)
////			{
//////				xm[idx + offset] = num2;
////				xm[offset] = num2;
//////				xConst = num2;
////				offset += permutationStride;
////			}
//////			offset = pos * permutationStride;
//////			xm[idx + offset]--;
////			xm[offset]--;
//////			xConst.i[posAdd]--;
////			//xm[pos*permutationStride+idx]--;
////		}
////		else
////		{
////			//			xm[0] = ULMAX;
//////			xm[idx] = ULMAX;
//////			xConst.i[0] = ULMAX;
////			xm[offset] = ULMAX;
////			offset += permutationStride;
//////			offset = idx+permutationStride;
////#pragma unroll (COEFFICIENT_SIZE-1)
////			for (i = 1; i < 8; i++)
////			{
////				//				xm[i] = 0;
//////				xm[idx + offset] = 0;
////				xm[offset] = 0;
//////				xConst.i[i] = 0;
////				offset+=permutationStride;
////			}
////		}
////	}
//
//	offset = idx;
//	if (c > 0 && posAdd >= 0)
//	{
//		for (i = 0; i < posAdd; i++)
//		{
//			xm[offset] = R_MINUS_ONE;
//			offset += permutationStride;
//		}
//		xm[offset]--;
//
//	}
//	if (c > 0 && posAdd == -1)
//	{
//		xm[offset] = ULMAX;
//		offset += permutationStride;
//
//#pragma unroll (COEFFICIENT_SIZE-1)
//		for (i = 1; i < 8; i++)
//		{
//			xm[offset] = 0;
//			offset += permutationStride;
//		}
//	}
////
////	if (c2 > 0)
////	{
//////		offset = 0;
//////		offset = idx;
//////		posSub = -1;
//////		for (i = 0; i < 8; i++)
//////		{
////////			if (ym[idx + offset] < num2)
//////			if (yConst.i[i] < num2)
//////			{
//////				posSub = i;
//////				break;
//////			}
////////			offset += permutationStride;
//////		}
////
////		offset = idx;
////		if (posSub >= 0)
////		{
//////			offset = 0;
////			for (i = 0; i < posSub; i++)
////			{
//////				ym[idx + offset] = 0;
////				ym[offset] = 0;
////				offset += permutationStride;
//////				yConst.i[i] = 0;
////			}
//////			offset = pos * permutationStride;
////			//			um[pos]++;
//////			ym[idx + offset]++;
//////			yConst.i[posSub]++;
////			ym[offset]++;
////		}
////		else
////		{
////			//			um[0] = ULMAX;
//////			ym[idx] = ULMAX;
//////			yConst.i[0] = ULMAX;
////			ym[offset] = ULMAX;
////			offset += permutationStride;
//////			offset = 0;
////#pragma unroll (COEFFICIENT_SIZE-1)
////			for (i = 1; i < 8; i++)
////			{
////				//				um[i] = 0;
//////				ym[idx + offset]=0;
//////				yConst.i[i] = 0;
////				ym[offset] = 0;
////				offset += permutationStride;
////			}
////		}
////	}
//
////		offset = 0;
////		offset = idx;
////		posSub = -1;
////		for (i = 0; i < 8; i++)
////		{
//////			if (ym[idx + offset] < num2)
////			if (yConst.i[i] < num2)
////			{
////				posSub = i;
////				break;
////			}
//////			offset += permutationStride;
////		}
//
//	offset = idx;
//	if (c2 > 0 && posSub >= 0)
//	{
////			offset = 0;
//		for (i = 0; i < posSub; i++)
//		{
////				ym[idx + offset] = 0;
//			ym[offset] = 0;
//			offset += permutationStride;
////				yConst.i[i] = 0;
//		}
////			offset = pos * permutationStride;
//		//			um[pos]++;
////			ym[idx + offset]++;
////			yConst.i[posSub]++;
//		ym[offset]++;
//	}
//	if (c2 > 0 && posSub == -1)
//	{
//		//			um[0] = ULMAX;
////			ym[idx] = ULMAX;
////			yConst.i[0] = ULMAX;
//		ym[offset] = ULMAX;
//		offset += permutationStride;
////			offset = 0;
//#pragma unroll (COEFFICIENT_SIZE-1)
//		for (i = 1; i < 8; i++)
//		{
//			//				um[i] = 0;
////				ym[idx + offset]=0;
////				yConst.i[i] = 0;
//			ym[offset] = 0;
//			offset += permutationStride;
//		}
//	}
//
////	offset = idx;
////	#pragma unroll COEFFICIENT_SIZE
////		for (i = 0; i < 8; i++)
////		{
////			xm[offset] = xConst.i[i];
//////			ym[offset] = yConst.i[i];
//////			offset+=permutationStride;
////			offset+=STRIDE;
////		}
//
////		offset = idx;
////	#pragma unroll COEFFICIENT_SIZE
////		for (i = 0; i < 8; i++)
////		{
////			ym[offset] = yConst.i[i];
//////			offset+=permutationStride;
////			offset+=STRIDE;
////		}
//
//}
/******************************************************************************/
//
//__device__ inline void fft_base2_permutated_6(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 permutationStride)
//{
//
//	short i = 0;
//	uConstArray8_align8 old_xm;
////	uConstArray8_align8 old_ym;
////	uConstArray8_align8 old_um;
//
//	usfixn64 offset = 0;
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
//////		old_xm.i[i]=xm[i*permutationStride];
//////		old_um.i[i]=xm[i*permutationStride];
//////		old_ym.i[i]=ym[i*permutationStride];
////		old_xm.i[i] = xm[offset];
////		old_um.i[i] = xm[offset];
////		old_ym.i[i] = ym[offset];
////		offset += permutationStride;
////	}
//
////	device_bigPrimeAdd_plain(old_xm.i, old_ym.i);
////	device_bigSub_plain(old_um.i, old_ym.i);
//
//	offset = 0;
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
//////		xm[i*permutationStride] = old_xm.i[i];
//////		ym[i*permutationStride] = old_um.i[i];
////		xm[offset] = old_xm.i[i];
////		ym[offset] = old_um.i[i];
////		offset+=permutationStride;
////	}
//	device_bigPrimeAdd_permutated_5(xm, ym, old_xm.i, permutationStride);
////	__syncthreads();
////	device_bigPrimeSub_permutated_2(xm, ym, permutationStride);
//	device_bigPrimeSub_permutated_4(old_xm.i, ym, permutationStride);
//////
////__syncthreads();
////
////	usfixn64 offset=0;
//
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
////	ym[i*permutationStride]=old_xm.i[i];
//		ym[offset] = old_xm.i[i];
//		offset += permutationStride;
//	}
//}

#endif
