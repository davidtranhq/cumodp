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
#ifndef BIG_ARITHMETIC_MULTIPLICATION_CU_
#define BIG_ARITHMETIC_MULTIPLICATION_CU_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"

#include "bigArithmetic_multiplication_P4_device.h"
//#include "BigPrimeFieldFFT_4/bigArithmetic_addition_P4.h"
//#include "BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.h"

/**********************************************/
__global__ void
kernel_p4_mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
												 const usfixn64* __restrict__ ys, usfixn64* parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 usfixn32* __restrict__ signVector)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0, c = 0;
	short i = 0, j = 0;
	usfixn64 sign = 0;
//	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
//			cVector[COEFFICIENT_SIZE];
//	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
//			cVectorSub[COEFFICIENT_SIZE];
//	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 y[COEFFICIENT_SIZE];
//	usfixn64 x[COEFFICIENT_SIZE];

//	usfixn64 lAdd, hAdd, cAdd;
//	usfixn64 lSub, hSub, cSub;
//	usfixn64 m0Array[COEFFICIENT_SIZE], m1Array[COEFFICIENT_SIZE];

//	__shared__ usfixn64 xshared[COEFFICIENT_SIZE*BLOCK_SIZE];
//	__shared__ usfixn64 xshared[COEFFICIENT_SIZE][BLOCK_SIZE];
//	usfixn64 shared_idx = threadIdx.x;

	usfixn64 gridStride = gridDim.x * blockDim.x;
//	usfixn64 gridStride = 8192;

	usfixn64 h0, h1, h2;
	usfixn64 u1;
	usfixn64 xOffset, xOffset_init = tid
			+ (COEFFICIENT_SIZE - 1) * permutationStride;
	for (tid; tid < n; tid += gridStride)
	{
		offset = tid;
//		shared_idx = threadIdx.x;
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
		for (i = COEFFICIENT_SIZE - 1; i >= 0; i--)
		{
//			x[i] = xs[offset];

//			xshared[shared_idx]=xs[tid+i*permutationStride];
//			xshared[threadIdx.x+i*BLOCK_SIZE]=xs[tid+i*permutationStride];
//			shared_idx+=BLOCK_SIZE;
//			xshared[i][threadIdx.x]=xs[offset];
//		}
//		for(i=0;i<COEFFICIENT_SIZE;i++)
//		{
//			y[(COEFFICIENT_SIZE - 1) - i] = ys[offset];
			y[i] = ys[offset];
			offset += permutationStride;
		}
//		for (i = 0; i < 8; i++)
//		{
//			lVectorSub[i] = 0;
//			hVectorSub[i] = 0;
//			cVectorSub[i] = 0;
//			lVector[i] = 0;
//			hVector[i] = 0;
//			cVector[i] = 0;
//		}

//		shared_idx = 15*BLOCK_SIZE + threadIdx.x;
		offset = tid + (permutationStride << 4);
//		offset = tid + COEFFICIENT_SIZE * (permutationStride);

		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		i=0;
		{
			offset -= permutationStride;
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
//				device_oneShiftRight (y);
				device_p4_oneShiftLeft (y);

			xOffset = xOffset_init;
//			shared_idx = 15*BLOCK_SIZE + threadIdx.x;
//			for (j = 8 - i; j < 8; j++)
//			for (j = 7; j >= 8 - i; j--)
			for (j = COEFFICIENT_SIZE - 1; j >= 0; j--)
			{
//				m=xs[j]*y[8-j];
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[7-j])
//				);

				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(xs[xOffset]),"l"(y[j])
				);

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[j])
//				);

//				asm("{\n\t"
//										"mul.lo.u64 %0,%2,%3;\n\t"
//										"mul.hi.u64 %1,%2,%3;\n\t"
//										"}"
//										:"=l"(m0),"=l"(m1)
//										:"l"(xshared[j][threadIdx.x]),"l"(y[j])
//								);

//
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0Array[j]),"=l"(m1Array[j])
//						:"l"(xs[xOffset]),"l"(y[j])
//				);

				xOffset -= permutationStride;
//				shared_idx -=BLOCK_SIZE;

//				if (j > 7 - i)
//					add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
//				else
//					add128 (u[0], u[1], m0, m1, u[2], u[0], u[1], u[2]);

//			}
//
//			for (j = 7; j >= 0; j--)
//			{
//				h0 = u[0];
//				h1 = u[1];
//				h2 = u[2];
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					h0 = uSub[0];
//					h1 = uSub[1];
//					h2 = uSub[2];
//				}
				h0 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[0] : uSub[0];
				h1 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[1] : uSub[1];
				h2 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[2] : uSub[2];

				/*****************************/
//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(h0),"=l"(c):"l"(h0),"l"(m0));

//				h0 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%1,%2;\n\t"
						"addc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(u1),"+l"(c):"l"(m1));
//				u1 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,%1,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
//				h1 = s;
//				h2 += c;
				asm("{\n\t"
						"add.cc.u64 %0,%0,%3;\n\t"
						"addc.cc.u64 %2,%2,%1;\n\t"
						"}"
						:"+l"(h1),"+l"(c),"+l"(h2):"l"(u1));
//								h1 = s;
//								h2 += c;

				/*****************************/
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					uSub[0] = h0;
//					uSub[1] = h1;
//					uSub[2] = h2;
//				}
//				else
//				{
//					u[0] = h0;
//					u[1] = h1;
//					u[2] = h2;
//				}
				s = (j > COEFFICIENT_SIZE - 1 - i);
				uSub[0] = (s) ? h0 : uSub[0];
				uSub[1] = (s) ? h1 : uSub[1];
				uSub[2] = (s) ? h2 : uSub[2];

				u[0] = (s) ? u[0] : h0;
				u[1] = (s) ? u[1] : h1;
				u[2] = (s) ? u[2] : h2;

				/*****************************/
			}

//			//			for (j = 0; j < 8 - i; j++)
//			for (j = 7 - i; j >= 0; j--)
//			{
//				//				m=xs[j]*y[8-j];
////				asm("{\n\t"
////						"mul.lo.u64 %0,%2,%3;\n\t"
////						"mul.hi.u64 %1,%2,%3;\n\t"
////						"}"
////						:"=l"(m0),"=l"(m1)
////						:"l"(x[j]),"l"(y[7-j])
////				);
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);
//				xOffset -= permutationStride;
//				c = u[2];
//				u[2] = 0;
//				add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
//				// 				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
//			}

//			device_p4_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
//						 cVectorSub[7 - i]);

//			device_p4_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
//			device_p4_toLHC (uSub[0], uSub[1], uSub[2], lSub, hSub, cSub);
//			device_p4_toLHC (u[0], u[1], u[2], lAdd, hAdd, cAdd);
//			offset = tid + (7 - i) * permutationStride;
			sign = 0;
			device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);

			lVector[offset] = u[0];
			hVector[offset] = u[1];
			cVector[offset] = u[2];
			signVector[offset] = sign;
//			hVectorSub[offset] = uSub[1];
//			cVectorSub[offset] = uSub[2];
//			offset=offset-permutationStride;
		}

//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);
	}
}
/**********************************************/
__global__ void
kernel_p4_mult_revised_8lhc_step1_v0 (const usfixn64 * __restrict__ xs,
														const usfixn64* __restrict__ ys,
														usfixn64* parameters,
														usfixn64* __restrict__ lVector,
														usfixn64* __restrict__ hVector,
														usfixn64* __restrict__ cVector,
														usfixn32* __restrict__ signVector)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0, c = 0;
	short i = 0, j = 0;
	usfixn64 sign = 0;
	usfixn64 lVectorLocal[COEFFICIENT_SIZE], hVectorLocal[COEFFICIENT_SIZE],
			cVectorLocal[COEFFICIENT_SIZE], signVectorLocal[COEFFICIENT_SIZE];
//	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
//			cVectorSub[COEFFICIENT_SIZE];
//	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 y[COEFFICIENT_SIZE];
//	usfixn64 x[COEFFICIENT_SIZE];

//	usfixn64 lAdd, hAdd, cAdd;
//	usfixn64 lSub, hSub, cSub;
//	usfixn64 m0Array[COEFFICIENT_SIZE], m1Array[COEFFICIENT_SIZE];

//	__shared__ usfixn64 xshared[COEFFICIENT_SIZE*BLOCK_SIZE];
//	__shared__ usfixn64 xshared[COEFFICIENT_SIZE][BLOCK_SIZE];
//		__shared__ usfixn64 yshared[COEFFICIENT_SIZE][BLOCK_SIZE];
//	usfixn64 shared_idx = threadIdx.x;

	usfixn64 gridStride = gridDim.x * blockDim.x;
//	usfixn64 gridStride = 8192;

	usfixn64 h0, h1, h2;
	usfixn64 u1;
	usfixn64 xOffset, xOffset_init = tid
			+ (COEFFICIENT_SIZE - 1) * permutationStride;
	for (tid; tid < n; tid += gridStride)
	{
		offset = tid;
//		shared_idx = threadIdx.x;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
//			x[i] = xs[offset];

//			xshared[shared_idx]=xs[tid+i*permutationStride];
//			xshared[threadIdx.x+i*BLOCK_SIZE]=xs[tid+i*permutationStride];
//			shared_idx+=BLOCK_SIZE;
//			xshared[i][threadIdx.x]=xs[offset];
//		}
//		for(i=0;i<COEFFICIENT_SIZE;i++)
//		{
			y[(COEFFICIENT_SIZE - 1) - i] = ys[offset];
//			yshared[(COEFFICIENT_SIZE - 1) - i][threadIdx.x]=ys[offset];
			offset += permutationStride;
		}
//		for (i = 0; i < 8; i++)
//		{
//			lVectorSub[i] = 0;
//			hVectorSub[i] = 0;
//			cVectorSub[i] = 0;
//			lVector[i] = 0;
//			hVector[i] = 0;
//			cVector[i] = 0;
//		}

//		shared_idx = 15*BLOCK_SIZE + threadIdx.x;
		offset = tid + (permutationStride << 4);
//		offset = tid + COEFFICIENT_SIZE * (permutationStride);

		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		i=0;
		{
			offset -= permutationStride;
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
//				device_oneShiftRight (y);
				device_p4_oneShiftLeft (y);

			xOffset = xOffset_init;
//			shared_idx = 15*BLOCK_SIZE + threadIdx.x;
//			for (j = 8 - i; j < 8; j++)
//			for (j = 7; j >= 8 - i; j--)
			for (j = COEFFICIENT_SIZE - 1; j >= 0; j--)
			{
//				m=xs[j]*y[8-j];
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[7-j])
//				);
//
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(xs[xOffset]),"l"(y[j])
				);

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(ys[((15-j+i)&0xF)*permutationStride+tid])
//				);

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[j])
//				);

//				asm("{\n\t"
//										"mul.lo.u64 %0,%2,%3;\n\t"
//										"mul.hi.u64 %1,%2,%3;\n\t"
//										"}"
//										:"=l"(m0),"=l"(m1)
//										:"l"(xshared[j][threadIdx.x]),"l"(y[j])
//								);

//
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0Array[j]),"=l"(m1Array[j])
//						:"l"(xs[xOffset]),"l"(y[j])
//				);

				xOffset -= permutationStride;
//				shared_idx -=BLOCK_SIZE;

//				if (j > 7 - i)
//					add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
//				else
//					add128 (u[0], u[1], m0, m1, u[2], u[0], u[1], u[2]);

//			}
//
//			for (j = 7; j >= 0; j--)
//			{
//				h0 = u[0];
//				h1 = u[1];
//				h2 = u[2];
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					h0 = uSub[0];
//					h1 = uSub[1];
//					h2 = uSub[2];
//				}
				h0 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[0] : uSub[0];
				h1 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[1] : uSub[1];
				h2 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[2] : uSub[2];

				/*****************************/
//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(h0),"=l"(c):"l"(h0),"l"(m0));

//				h0 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%1,%2;\n\t"
						"addc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(u1),"+l"(c):"l"(m1));
//				u1 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,%1,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
//				h1 = s;
//				h2 += c;
				asm("{\n\t"
						"add.cc.u64 %0,%0,%3;\n\t"
						"addc.cc.u64 %2,%2,%1;\n\t"
						"}"
						:"+l"(h1),"+l"(c),"+l"(h2):"l"(u1));
//								h1 = s;
//								h2 += c;

				/*****************************/
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					uSub[0] = h0;
//					uSub[1] = h1;
//					uSub[2] = h2;
//				}
//				else
//				{
//					u[0] = h0;
//					u[1] = h1;
//					u[2] = h2;
//				}
				s = (j > COEFFICIENT_SIZE - 1 - i);
				uSub[0] = (s) ? h0 : uSub[0];
				uSub[1] = (s) ? h1 : uSub[1];
				uSub[2] = (s) ? h2 : uSub[2];

				u[0] = (s) ? u[0] : h0;
				u[1] = (s) ? u[1] : h1;
				u[2] = (s) ? u[2] : h2;

				/*****************************/
			}

//			//			for (j = 0; j < 8 - i; j++)
//			for (j = 7 - i; j >= 0; j--)
//			{
//				//				m=xs[j]*y[8-j];
////				asm("{\n\t"
////						"mul.lo.u64 %0,%2,%3;\n\t"
////						"mul.hi.u64 %1,%2,%3;\n\t"
////						"}"
////						:"=l"(m0),"=l"(m1)
////						:"l"(x[j]),"l"(y[7-j])
////				);
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);
//				xOffset -= permutationStride;
//				c = u[2];
//				u[2] = 0;
//				add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
//				// 				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
//			}

//			device_p4_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
//						 cVectorSub[7 - i]);

//			device_p4_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
//			device_p4_toLHC (uSub[0], uSub[1], uSub[2], lSub, hSub, cSub);
//			device_p4_toLHC (u[0], u[1], u[2], lAdd, hAdd, cAdd);
//			offset = tid + (7 - i) * permutationStride;
			sign = 0;
			device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);

//			lVector[offset] = u[0];
//			hVector[offset] = u[1];
//			cVector[offset] = u[2];
//			signVector[offset] = sign;

			lVectorLocal[COEFFICIENT_SIZE - 1 - i] = u[0];
			hVectorLocal[COEFFICIENT_SIZE - 1 - i] = u[1];
			cVectorLocal[COEFFICIENT_SIZE - 1 - i] = u[2];
			signVectorLocal[COEFFICIENT_SIZE - 1 - i] = sign;
//			hVectorSub[offset] = uSub[1];
//			cVectorSub[offset] = uSub[2];
//			offset=offset-permutationStride;
		}

		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			lVector[offset] = lVectorLocal[i];
			hVector[offset] = hVectorLocal[i];
			cVector[offset] = cVectorLocal[i];
			signVector[offset] = signVectorLocal[i];
			offset += permutationStride;
		}
//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);
	}
}

/**********************************************/
__global__ void
kernel_p4_mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 const usfixn32* __restrict__ signVector)
{

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];

//	usfixn64 s0, s1, s2;
//	usfixn64 sign;
	usfixn64 l, h, c;
//	usfixn64 v0[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v1[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v0[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};
//		usfixn64 v1[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
		offset = tid;
		short i = 0;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
//		memset (v0, 0x00, 3 * sizeof(usfixn64));
//		memset (v1, 0x00, 3* sizeof(usfixn64));
//		s0 = lVector[offset];
//		s1 = hVector[offset];
//		s2 = cVector[offset];
//		sign = lVectorSub[offset];
			l = 0;
			h = 0;
			c = 0;
//		device_p4_toLHC (s0, s1, s2, l, h, c);
			device_p4_toLHC (lVector[offset], hVector[offset], cVector[offset], l, h, c);

//		v0[0] = l;
//		v0[1] = h;
//		v0[2] = c;
//		if (sign == 1)
			if (signVector[offset] == 1)
			{
				device_p4_bigSubZero_3 (l, h, c);

//			device_p4_bigSub_plain (v1, v0);
//			l = v1[0];
//			h = v1[1];
//			c = v1[2];
			}
			lVector[offset] = l;
			hVector[offset] = h;
			cVector[offset] = c;
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p4_mult_revised_8lhc_step3_v0 (usfixn64* __restrict__ xs,
														const usfixn64* __restrict__ lVector,
														const usfixn64 * __restrict__ hVector,
														const usfixn64* __restrict__ cVector,
														const usfixn32* __restrict__ signVector,
														const usfixn64 * __restrict__ parameters)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;
	usfixn64 n = parameters[7];

	uConstArray16_align8 v0;
	uConstArray16_align8 v1;
//	usfixn64 complement[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE };
	usfixn64 complement[COEFFICIENT_SIZE] =
		{ 0, 0, 0, R_MINUS_ONE,
		R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE,
		R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE,
		R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE };

	short j = 0;
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
//		v0.i[i]=0;
			v1.i[i] = 0;
		}

//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

		offset = tid;
//all l values [0:7]
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			v0.i[i] = lVector[offset];
			offset += permutationStride;
		}

		complement[0] = 0;
		complement[1] = 0;
		complement[2] = 0;
		for (j = 3; j < COEFFICIENT_SIZE; j++)
			complement[j] = R_MINUS_ONE;
		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (i > 0)
//			device_cyclicShift_plain (complement, 1);
//				device_cyclicShift_permutated_12 (complement);
				device_p4_cyclicShift_lhc_complement (complement, i);
			if (signVector[offset] == 1)
			{
				device_bigPrimeAdd_plain_inPlace (v0.i, complement);
//			device_bigPrimeAdd_plain_ptx_v0(v0.i,complement);
			}
			offset += permutationStride;
		}

		offset = tid + (permutationStride << 4) - permutationStride;

//	v1.i[0] = hVector[tid + 7 * permutationStride];
		v1.i[0] = hVector[offset];
//	device_p4_bigSub_plain(v0.i, v1.i);
		device_p4_bigPrimeSub_plain_inPlace  (v0.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}
//	v1.i[1] = cVector[tid + 7 * permutationStride];
//	v1.i[0] = cVector[tid + 6 * permutationStride];
		v1.i[1] = cVector[offset];
		offset -= permutationStride;
		v1.i[0] = cVector[offset];

//	device_p4_bigSub_plain(v0.i, v1.i);
		device_p4_bigPrimeSub_plain_inPlace  (v0.i, v1.i);

		offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}

		v1.i[0] = 0;
//positive h's [1:7]
		offset = tid;
		for (i = 1; i < COEFFICIENT_SIZE; i++)
		{
//		v1.i[i] = hVector[tid + (i - 1) * permutationStride];
			v1.i[i] = hVector[offset];
			offset += permutationStride;
		}
//		device_bigPrimeAdd_plain (v0.i, v1.i);
		device_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

		offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//			{
//				v1.i[i]=0;
//			}
		v1.i[0] = 0;
		v1.i[1] = 0;
		offset = tid;
		for (i = 2; i < COEFFICIENT_SIZE; i++)
		{
//		v1.i[i] = cVector[tid + (i - 2) * permutationStride];
			v1.i[i] = cVector[offset];
			offset += permutationStride;
		}
//		device_bigPrimeAdd_plain (v0.i, v1.i);
		device_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

//#################################################
//	uConstArray8_align8 v0_sub;
//	uConstArray8_align8 v1_sub;
//	uConstArray8_align8 v2_sub;

//	memset(v0_sub.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for (i = 0; i < 8; i++)
//	{
////				v0_sub.i[i]=0;
//		v1.i[i] = 0;
//	}
//	memset(v2_sub.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

//	############################# writing back to g-memory
		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
			offset += permutationStride;
		}
	}
}

/****************************************************/
__global__ void
kernel_p4_bigMult_plain (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
{
	//
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
	usfixn64 *xd, *yd, *xm, *ym;
//	usfixn64 *xd, *xm;

//	num1 = 0;
//	num2 = R - 1;

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * COEFFICIENT_SIZE);
	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	ud = (usfixn64*) ((char*) us + tid * sizeof(usfixn64) * 8);

//	device_cyclicShift_plain(xd,shiftNo);
//	bigMult(xs, ys, us);
//	device_bigMult_plain (xd, yd); //should be re-written based on P4 arithmetic
//	xd[0]=tid;
}
/****************************************************/

#endif
