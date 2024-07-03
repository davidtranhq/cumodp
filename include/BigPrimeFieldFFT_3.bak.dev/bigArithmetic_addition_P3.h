#ifndef BIG_ARITHMETIC_ADDITION_H_P3
#define BIG_ARITHMETIC_ADDITION_H_P3

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"
#include "BigPrimeFieldFFT_3/bigArithmetic_cyclicShift_P3.h"

/**********************************************/

__device__ void
device_p3_shiftRight_uConstArray8 (uConstArray8_align8 & x, usfixn64 & tmp)
{
//	tmp = x.i[7];
//	x.i[7] = x.i[6];
//	x.i[6] = x.i[5];
//	x.i[5] = x.i[4];
//	x.i[4] = x.i[3];
//	x.i[3] = x.i[2];
//	x.i[2] = x.i[1];
//	x.i[1] = x.i[0];
//	x.i[0] = tmp;
	//device_negate_plain(tmp);
	asm("{\n\t"
			".reg .u64 t0 ;\n\t"
			"mov.u64 t0,%7;\n\t"
			"mov.u64 %7,%6;\n\t"
			"mov.u64 %6,%5;\n\t"
			"mov.u64 %5,%4;\n\t"
			"mov.u64 %4,%3;\n\t"
			"mov.u64 %3,%2;\n\t"
			"mov.u64 %2,%1;\n\t"
			"mov.u64 %1,%0;\n\t"
			"mov.u64 %0,t0;\n\t"
			"}"
			:"+l"(x.i[0]),"+l"(x.i[1]),"+l"(x.i[2]),"+l"(x.i[3]),"+l"(x.i[4]),"+l"(x.i[5]),"+l"(x.i[6]),"+l"(x.i[7])
			:);
}

/**********************************************/
__device__ void
device_p3_bigPrimeAdd_correct (usfixn64 *xm, usfixn64 *ym, usfixn64 *um)
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

/**********************************************/
__device__ __inline__ void
device_p3_bigPrimeAdd_plain (usfixn64 *__restrict__ xm, usfixn64 *__restrict__ ym)
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

/**********************************************/
__device__ __inline__ void
device_p3_bigPrimeAdd_plain_inPlace (usfixn64 * __restrict__ xm,
																	usfixn64 * __restrict__ ym)
{
	unsigned short c = 0;

	short pos1;
	usfixn64 num1, num2;

//	num1 = 0;
//	num2 = R - 1;
	short i;
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
device_p3_bigPrimeAdd_plain_ptx_v0 (usfixn64 * __restrict__ xm,
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
		}
		else
		{
			xm[pos1]--;
			for (i = 0; i < pos1; i++)
				xm[i] = R_MINUS_ONE;
		}
	}
}

/**********************************************/
__device__ void
device_p3_bigPrimeAdd_permutated (const usfixn64 *xm, const usfixn64 * ym, usfixn64 *um,
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

/**********************************************/
__device__ void
device_p3_bigPrimeAdd_permutated_ptx_v0 (usfixn64 * xm, usfixn64 *ym,
																			const usfixn64 permutationStride)
{
	unsigned short c = 0;
	usfixn64 num1 = 0;
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
	for (i = 0; i < 8; i++)
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
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;
			offset += permutationStride;
			xm[offset] = 0;

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

/**********************************************/
//xm = xm + ym
//ym = xm - ym
//__device__ inline void fft_base2_permutated(usfixn64 * __restrict__ xm,
//		usfixn64 * __restrict__ ym, const usfixn64 xIdx,const usfixn64 yIdx,
//		const usfixn64 permutationStride)
__device__ __inline__ void
device_p3_fft_base2_permutated (usfixn64 * xm, usfixn64 * ym,
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
		for (i = 0; i < 8; i++)
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
			for (i = 1; i < 8; i++)
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
		for (i = 0; i < 8; i++)
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
			for (i = 1; i < 8; i++)
			{
				//				um[i] = 0;
				ym[offset] = 0;
				offset += permutationStride;
				__syncthreads ();
			}
		}
	}
}

/**********************************************/
//using vectorized data structure + offset is computed more efficiently
__device__ inline void
device_p3_fft_base2_permutated_4 (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
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

/**********************************************/
__device__ inline void
device_p3_fft_base2_permutated_7 (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
												const short shNo, const usfixn64 permutationStride)
{
	short c = 0;
	usfixn64 num1;
//	num2;
	unsigned short c2 = 0;
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

	device_p3_cyclicShift_permutated_5 (ym, shNo, permutationStride);
	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		num1 = xm[idx + offset] + ym[idx + offset] + c;
//		num3 = ym[idx + offset] + c2;
//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;
//		xConst = xm[offset];
//		yConst = ym[offset];

//		num1 = xConst + yConst + c;
//		num3 = yConst + c2;

//		num1 = xm[offset] + ym[offset] + c;
//		num3 = ym[offset] + c2;

		num1 = xm[offset];		// +ym[offset] + c;
		num3 = ym[offset]; //+ c2;
		num1 = num1 + c + num3;
		num3 += c2;
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
		xm[offset] = num1;
		posAdd = (posAdd == -1 && num1 > 0) * i;
		posSub = (posSub == -1 && num3 > 0) * i;
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
			xm[offset] = R_MINUS_ONE;
			offset += permutationStride;
		}
		xm[offset]--;
	}

	if (c > 0 && posAdd == -1)
	{
		xm[offset] = ULMAX;
		offset += permutationStride;

//#pragma unroll (COEFFICIENT_SIZE-1)
		for (i = 1; i < COEFFICIENT_SIZE; i++)
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
	if (c2 > 0 && posSub == -1)
	{
		//			um[0] = ULMAX;
//			ym[idx] = ULMAX;
//			yConst.i[0] = ULMAX;
//		ym[offset] = ULMAX;
//		offset += permutationStride;
		ym[0] = ULMAX;
		offset = permutationStride;
//			offset = 0;
//#pragma unroll (COEFFICIENT_SIZE-1)
		for (i = 1; i < COEFFICIENT_SIZE; i++)
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

/**********************************************/
//using vectorized data structure + offset is computed more efficiently
__device__ void
device_p3_fft_base2_permutated_8 (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
												const short & shNo, const usfixn64 & permutationStride)
{
//	device_cyclicShift_permutated_5(ym, shNo, permutationStride);
//	device_cyclicShift_permutated_6(ym, shNo, permutationStride);
	device_p3_cyclicShift_permutated_7 (ym, shNo, permutationStride);
	device_p3_fft_base2_permutated (xm, ym, permutationStride);
//	fft_base2_vector8_permutated(xm, ym, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_addition_plain (usfixn64 * xs, usfixn64 * ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
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

//first solution, COEFFICIENT_SIZE iterations over all digits
//	device_bigPrimeAdd_plain(xd, yd);

	device_p3_bigPrimeAdd_plain_ptx_v0 (xd, yd);

	/*
	 * second solution, going on for every 2 coefficients,
	 * checking the carry at the end in bigPrimeAdd_check,
	 * should pass carry from each 2-step function, initially c=0
	 */
//	i = 0;
//	device_smallAdd2_plain(&xd[i],&yd[i],c);
//	i += 2;
//	device_smallAdd2_plain(&xd[i],&yd[i],c);
//	i += 2;
//	device_smallAdd2_plain(&xd[i],&yd[i],c);
//	i += 2;
//	device_smallAdd2_plain(&xd[i],&yd[i],c);
//	bigPrimeAdd_check(xd,c);
}

/**********************************************/
__global__ void
kernel_p3_addition_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
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
//	if (padding == 0)
	device_p3_bigPrimeAdd_permutated_ptx_v0 (&xs[idx], &ys[idx], permutationStride);

//	usfixn64 offset=tid;
//	uConstArray8_align8 x,y;
//	for(short i=0;i<8;i++)
//	{
//		x.i[i]=xs[offset];
//		y.i[i]=ys[offset];
//		offset+=permutationStride;
//	}
//	device_bigPrimeAdd_plain_inPlace_2_bits(x.i,y.i);
//
//	offset=tid;
//	for(short i=0; i<8;i++)
//	{
//		xs[offset]=x.i[i];
//		offset+=permutationStride;
//	}
}

/**********************************************/
#endif /* BIG_ARITHMETIC_ADDITION_H_P3 */
