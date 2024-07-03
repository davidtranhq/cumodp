#ifndef BIG_ARITHMETIC_SUBTRACTION_H_
#define BIG_ARITHMETIC_SUBTRACTION_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"

/**********************************************/

// In place subtraction in Z/(R^8 + 1)Z (xs = xs - ys)
__device__ void
device_p3_bigSub (usfixn64 *xm, usfixn64 *ym, usfixn64 *um)
{
//	usfixn64 um[8];
	short i, pos;
	unsigned short c = 0;
	usfixn64 num1, num2;
	num1 = 0;
	num2 = R - 1;

	for (i = 0; i <= 7; i++)
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
		for (i = 0; i < 8; i++)
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
			for (i = 1; i < 8; i++)
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
__device__ void
device_p3_bigPrimeSub_permutated_2 (usfixn64 *xm, usfixn64 * ym,
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
		for (i = 0; i < 8; i++)
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
			for (i = 1; i < 8; i++)
			{
//				um[i] = 0;
				xm[offset] = 0;
				offset += permutationStride;
			}
		}
	}
}

/**********************************************/
__device__ void
device_p3_bigSub_plain (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym)
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

/**********************************************/
__device__ void
device_p3_bigSub_uVector_plain_2 (uvector8 & __restrict__ xm,
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

/**********************************************/
__device__ __inline__ void
device_p3_bigPrimeSub_plain_inPlace (usfixn64 * __restrict__ xm,
																	usfixn64 * __restrict__ ym)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1;

//	num1 = 0;
//	num2 = R - 1;
	short i;
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
		}
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_bigPrimeSub_plain_ptx_v0 (usfixn64 * __restrict__ xm,
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
		for (i = 1; i < 8; i++)
			xm[i] = 0;
	}
}

/**********************************************/
__device__ void
device_p3_bigSubZero_3 (usfixn64 &l, usfixn64 &h, usfixn64 &carry)
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
	{
		c = 1;
		s = R - num1;
	}
	l = s;

	num1 = h + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
	{
		c = 1;
		s = R - num1;
	}
	h = s;

	num1 = carry + c;
//		if (xm[i] < num1) //there is not enough to do subtraction
	{
		c = 1;
		s = R - num1;
	}
	carry = s;
}

/**********************************************/
__global__ void
kernel_p3_subtraction_plain (usfixn64 * xs, usfixn64 * ys, usfixn64 *parameters)
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

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * 8);
	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * 8);
//	ud = (usfixn64*) ((char*) us + tid * sizeof(usfixn64) * 8);

//	device_bigSub_plain(xd, yd);
//	device_bigPrimeSub_plain_inPlace(xd,yd);
//	device_bigPrimeSub_permutated_2(xd,yd, 1);
	device_p3_bigPrimeSub_plain_ptx_v0 (xd, yd);
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
kernel_p3_subtraction_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
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
		device_p3_bigPrimeSub_permutated_2 (&xs[idx], &ys[idx], permutationStride);
}

/**********************************************/
#endif
