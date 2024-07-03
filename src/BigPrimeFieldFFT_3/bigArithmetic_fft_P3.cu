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


#ifndef BIG_ARITHMETIC_FFT_CU_
#define BIG_ARITHMETIC_FFT_CU_

#include "../../include/BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "bigArithmetic_fft_P3_device.h"

extern __shared__ usfixn64 sharedMem[];

/**********************************************/
//FFT base 16 to be called from host
__global__ void
kernel_p3_base_fft16 (usfixn64 *xm, usfixn64 * parameters)
{
	//short rounds = 4; //log(16,2)

	short stride = 16; //stride/2
	stride /= 2;
	usfixn64 permutationStride = parameters[5];
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	device_base_fft16 (xm, permutationStride, tid);
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r1 (usfixn64 * xm,
													 const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	device_base_fft16_7_r1 (xm, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r21 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);

//#if (!NEW_CYCLIC_SHIFT_SHMEM )
//		device_base_fft16_7_r21(xm, permutationStride);
//#else
////	__shared__ usfixn64 tsh[BLOCK_SIZE][COEFFICIENT_SIZE / 2];
	usfixn64 tsh[COEFFICIENT_SIZE / 2];
	device_base_fft16_7_r21_shmem (xm, tsh, permutationStride);
//#endif
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r22 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	device_base_fft16_7_r22 (xm, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r31 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);

//#if (!NEW_CYCLIC_SHIFT_SHMEM )
//		device_base_fft16_7_r31(xm, permutationStride);
//#else
	__shared__ usfixn64 tsh[BLOCK_SIZE][COEFFICIENT_SIZE / 2];
	device_base_fft16_7_r31_shmem (xm, tsh[threadIdx.x], permutationStride);
//#endif
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r32 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	device_base_fft16_7_r32 (xm, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r41 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
//#if (!NEW_CYCLIC_SHIFT_SHMEM )
//	device_base_fft16_7_r41(xm, permutationStride);
//#else
	__shared__ usfixn64 tsh[BLOCK_SIZE][COEFFICIENT_SIZE / 2];
	device_base_fft16_7_r41_shmem (xm, tsh[threadIdx.x], permutationStride);
//#endif
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r42 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	device_base_fft16_7_r42 (xm, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r5 (usfixn64 * xm,
													 const usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
//	device_base_fft16_7_r5(xm, permutationStride);
	device_base_fft16_7_r5_v1 (xm, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_fft2_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
{
	short operation = parameters[0];
	usfixn64 permutationStride = parameters[5];
//	short shuffle = parameters[6];
	short padding = 0;	//parameters[?]
	usfixn64 idx;
	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);

	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
	//following indexing is slightly faster than above indexing

	idx = tid;
	if (padding == 0)
		device_p3_fft_base2_permutated (&xs[idx], &ys[idx], permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step1 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters)

{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 permutationStride = parameters[5];
	usfixn64 i, j;
	usfixn64 h, w;
	usfixn64 stride;

//	stride = l / K;
	stride = l;
	i = (tid % (l * K)) % stride;
	j = (tid % (l * K)) / stride;
	if ((i * j) == 0)
		return;

	h = (i * j) / stride; // h <K
	w = (i * j) % stride; // w < N/K

	device_p3_cyclicShift_twice_permutated_1 (&xs[tid], h, permutationStride);
}

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step2 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 permutationStride = parameters[5];
	usfixn64 i, j;
	usfixn64 h, w;
	usfixn64 stride;

//	stride = l / K;
	stride = l;
	i = (tid % (l * K)) % stride;
	j = (tid % (l * K)) / stride;
	if ((i * j) == 0)
		return;

	h = (i * j) / stride; // h <K
	w = (i * j) % stride; // w < N/K

	usfixn64 x[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 offset = tid;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		x[i] = xs[offset];
		offset += permutationStride;
	}

//	w=1;
	if (w > 0)
	{
		offset = w - 1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			y[i] = powers_omega[offset];
//			if (h==0 && w==1)
//				printf("i = %d, y[i]=%lu \n",i, y[i]);
			offset += stride;
		}
		device_p3_bigMult_plain (x, y);
//		if (h == 0 && w == 1)
//		{
//			for (i = 0; i < COEFFICIENT_SIZE; i++)
//				printf ("x[i]=%lu \n", x[i]);
//		}
	}

	offset = tid;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		xs[offset] = x[i];
//				xs[offset] = y[i];
		offset += permutationStride;
	}
}

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step21 (usfixn64 * xs,
																						 usfixn64* powers_omega, usfixn64 K,
																						 usfixn64 l, usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 permutationStride = parameters[5];
	short i, j;
	usfixn64 h, w;
	usfixn64 i0, j0;
	usfixn64 stride;

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
//		usfixn64 offset = 0;
//		usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0, c = 0;
//		short i = 0, j = 0;
	usfixn64 sign = 0;

//	stride = l / K;
	stride = l;
	i0 = (tid % (l * K)) % stride;
	j0 = (tid % (l * K)) / stride;
	if ((i0 * j0) == 0)
		return;

	h = (i0 * j0) / stride; // h <K
	w = (i0 * j0) % stride; // w < N/K
	if (w == 0)
		return;

	usfixn64 x[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 offset = tid;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		x[i] = xs[offset];
		offset += permutationStride;
	}

	usfixn64 h0, h1, h2;
	usfixn64 xOffset, xOffset_init = tid + 7 * permutationStride;
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)

//	if (w > 0)
	if (w == 0)
		return;
	{
		offset = w - 1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			y[COEFFICIENT_SIZE - i - 1] = powers_omega[offset];
			offset += stride;
		}
//			device_p3_bigMult_plain (x, y);
	}
	{
//		offset = tid;
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			x[i] = xs[offset];
//			y[7 - i] = ys[offset];
//			offset += permutationStride;
//		}

		offset = tid + (permutationStride << 3);

		for (i = 0; i < COEFFICIENT_SIZE; i++)
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
				device_p3_oneShiftLeft (y);

			xOffset = xOffset_init;
			for (j = 7; j >= 0; j--)
			{
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[j])
				);

				xOffset -= permutationStride;
				h0 = u[0];
				h1 = u[1];
				h2 = u[2];
				if (j > 7 - i)
				{
					h0 = uSub[0];
					h1 = uSub[1];
					h2 = uSub[2];
				}

				/*****************************/
				//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
				//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.cc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(s),"=l"(c):"l"(h0),"l"(m0));

				h0 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%1,%2;\n\t"
						"addc.cc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(s),"+l"(c):"l"(m1));
				usfixn64 u1 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.cc.u64 %1,%1,0;\n\t"
						"}"
						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
				h1 = s;
				h2 += c;

				/*****************************/
				if (j > 7 - i)
				{
					uSub[0] = h0;
					uSub[1] = h1;
					uSub[2] = h2;
				}
				else
				{
					u[0] = h0;
					u[1] = h1;
					u[2] = h2;
				}
			}
			sign = 0;
			device_p3_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
			lVector[offset] = u[0];
			hVector[offset] = u[1];
			cVector[offset] = u[2];
			signVector[offset] = sign;
		}
	}

//	w=1;

//	offset = tid;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		xs[offset] = x[i];
////				xs[offset] = y[i];
//		offset += permutationStride;
//	}
}

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step22 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters)
{

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];

	usfixn64 h, w;
	usfixn64 i0, j0;
	usfixn64 stride;

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	//		usfixn64 offset = 0;
	//		usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0;
	//		short i = 0, j = 0;
	usfixn64 sign = 0;

	//	stride = l / K;
	stride = l_input;
	short i = 0;
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
		i0 = (tid % (l_input * K)) % stride;
		j0 = (tid % (l_input * K)) / stride;
		if ((i0 * j0) == 0)
			//return;
			continue;

		h = (i0 * j0) / stride; // h <K
		w = (i0 * j0) % stride; // w < N/K

		if (w == 0)
//		return;
			continue;
//	usfixn64 s0, s1, s2;
//	usfixn64 sign;
		usfixn64 l, c;
//	usfixn64 v0[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v1[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v0[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};
//		usfixn64 v1[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};

		offset = tid;

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
//		device_p3_toLHC (s0, s1, s2, l, h, c);
			device_p3_toLHC (lVector[offset], hVector[offset], cVector[offset], l, h,
											 c);

//		v0[0] = l;
//		v0[1] = h;
//		v0[2] = c;
//		if (sign == 1)
			if (signVector[offset] == 1)
			{
				device_p3_bigSubZero_3 (l, h, c);

//			device_p3_bigSub_plain (v1, v0);
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
kernel_p3_twiddle_general_permutated_step23 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;
	usfixn64 n = parameters[7];

	usfixn64 h, w;
	usfixn64 i0, j0;
	usfixn64 stride;

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	//		usfixn64 offset = 0;
	//		usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0, c = 0;
	//		short i = 0, j = 0;
	usfixn64 sign = 0;

	//	stride = l / K;
	stride = l_input;
	i0 = (tid % (l_input * K)) % stride;
	j0 = (tid % (l_input * K)) / stride;
	if ((i0 * j0) == 0)
		return;

	h = (i0 * j0) / stride; // h <K
	w = (i0 * j0) % stride; // w < N/K
	if (w == 0)
		return;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
	usfixn64 complement[COEFFICIENT_SIZE] =
		{ 0, 0, 0, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE };

	short j = 0;
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
		for (i = 0; i < 8; i++)
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
				device_p3_cyclicShift_plain (complement, 1);
//				device_cyclicShift_permutated_12 (complement);
			if (signVector[offset] == 1)
			{
				device_p3_bigPrimeAdd_plain_inPlace (v0.i, complement);
//			device_p3_bigPrimeAdd_plain_ptx_v0(v0.i,complement);
			}
			offset += permutationStride;
		}

		offset = tid + (permutationStride << 3) - permutationStride;

//	v1.i[0] = hVector[tid + 7 * permutationStride];
		v1.i[0] = hVector[offset];
//	device_p3_bigSub_plain(v0.i, v1.i);
		device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

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

//	device_p3_bigSub_plain(v0.i, v1.i);
		device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

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
//		device_p3_bigPrimeAdd_plain (v0.i, v1.i);
		device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

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
//		device_p3_bigPrimeAdd_plain (v0.i, v1.i);
		device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

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

//	offset = tid;
////all l values [0:7]
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
////		v0_sub.i[i] = lVectorSub[tid + (i) * permutationStride];
//		v0_sub.i[i] = lVectorSub[offset];
//		offset += permutationStride;
//	}
//	v0_sub.i[7] = 0;
//
////	v1.i[0] = hVectorSub[tid + 7 * permutationStride];
////	device_p3_bigSub_plain(v0_sub.i, v1.i);
////	device_p3_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);
//
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = cVectorSub[tid + 6 * permutationStride];
////	v1.i[1] = cVectorSub[tid + 7 * permutationStride];
//
////	device_p3_bigPrimeAdd_plain(v1.i, v2_sub.i);
////	device_p3_bigSub_plain(v0_sub.i, v1.i);
//	device_p3_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);
//
//	offset = tid;
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = 0;
////positive h's [1:7]
//	offset = tid;
//	for (i = 1; i < COEFFICIENT_SIZE; i++)
//	{
////		v1.i[i] = hVectorSub[tid + (i - 1) * permutationStride];
//		v1.i[i] = hVectorSub[offset];
////		v0_sub.i[i] = hVectorSub[tid + (8-i) * permutationStride];
////		v0_sub.i[i] = hVectorSub[tid + (6+i) * permutationStride];
//		offset += permutationStride;
//	}
//	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);
//
//	offset = tid;
////positive c's [2:7]
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = 0;
//	v1.i[1] = 0;
//	offset = tid;
//	for (i = 2; i < COEFFICIENT_SIZE; i++)
//	{
////		v1.i[i] = cVectorSub[tid + (i - 2) * permutationStride];
//		v1.i[i] = cVectorSub[offset];
//		offset += permutationStride;
//	}
//	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);
//
////#################################################
////	device_p3_bigSub_plain(v0.i, v0_sub.i);
//	device_p3_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
////#################################################

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

/**********************************************/
//performing fft-16 on plain data on a single thread
__global__ void
kernel_p3_fft16_plain_singleThread (usfixn64 *xs, usfixn64 *ys)
{
	//Now, we only use one thread to do FFT16 16.
	//FFT16 16 has four rounds. Each round can be done in parallel.
	//int tid = blockIdx.x*blockDim.x + threadIdx.x;
	short i, j;
	usfixn64 *xd, *yd;

	//point initialization
	//xd = (usfixn64*)((char*)xs + tid*sizeof(usfixn64)*8);
	xd = (usfixn64*) ((char*) xs);
	yd = (usfixn64*) ((char*) ys);

	//first round
	//0  8
	device_p3_bigPrimeAdd_correct (&xd[0], &xd[64], &yd[0]);
	device_p3_bigSub (&xd[0], &xd[64], &yd[64]);
	//4  12
	device_p3_bigPrimeAdd_correct (&xd[32], &xd[96], &yd[32]);
	device_p3_bigSub (&xd[32], &xd[96], &yd[96]);
	//2  10
	device_p3_bigPrimeAdd_correct (&xd[16], &xd[80], &yd[16]);
	device_p3_bigSub (&xd[16], &xd[80], &yd[80]);
	//6  14
	device_p3_bigPrimeAdd_correct (&xd[48], &xd[112], &yd[48]);
	device_p3_bigSub (&xd[48], &xd[112], &yd[112]);
	//1  9
	device_p3_bigPrimeAdd_correct (&xd[8], &xd[72], &yd[8]);
	device_p3_bigSub (&xd[8], &xd[72], &yd[72]);
	//5  13
	device_p3_bigPrimeAdd_correct (&xd[40], &xd[104], &yd[40]);
	device_p3_bigSub (&xd[40], &xd[104], &yd[104]);
	//3  11
	device_p3_bigPrimeAdd_correct (&xd[24], &xd[88], &yd[24]);
	device_p3_bigSub (&xd[24], &xd[88], &yd[88]);
	//7  15
	device_p3_bigPrimeAdd_correct (&xd[56], &xd[120], &yd[56]);
	device_p3_bigSub (&xd[56], &xd[120], &yd[120]);
	for (i = 0; i < 128; i++)
	{
		xd[i] = yd[i];
	}

	//second round
	//0 4
	device_p3_bigPrimeAdd_correct (&xd[0], &xd[32], &yd[0]);
	device_p3_bigSub (&xd[0], &xd[32], &yd[32]);
	//8 12
	device_p3_cyclicShift (&xd[96], 4);
	device_p3_bigPrimeAdd_correct (&xd[64], &xd[96], &yd[64]);
	device_p3_bigSub (&xd[64], &xd[96], &yd[96]);

	//2  6
	device_p3_bigPrimeAdd_correct (&xd[16], &xd[48], &yd[16]);
	device_p3_bigSub (&xd[16], &xd[48], &yd[48]);
	//10 14
	device_p3_cyclicShift (&xd[112], 4);
	device_p3_bigPrimeAdd_correct (&xd[80], &xd[112], &yd[80]);
	device_p3_bigSub (&xd[80], &xd[112], &yd[112]);

	//1 5
	device_p3_bigPrimeAdd_correct (&xd[8], &xd[40], &yd[8]);
	device_p3_bigSub (&xd[8], &xd[40], &yd[40]);
	//9 13
	device_p3_cyclicShift (&xd[104], 4);
	device_p3_bigPrimeAdd_correct (&xd[72], &xd[104], &yd[72]);
	device_p3_bigSub (&xd[72], &xd[104], &yd[104]);

	//3  7
	device_p3_bigPrimeAdd_correct (&xd[24], &xd[56], &yd[24]);
	device_p3_bigSub (&xd[24], &xd[56], &yd[56]);
	//11  15
	device_p3_cyclicShift (&xd[120], 4);
	device_p3_bigPrimeAdd_correct (&xd[88], &xd[120], &yd[88]);
	device_p3_bigSub (&xd[88], &xd[120], &yd[120]);
	for (i = 0; i < 128; i++)
	{
		xd[i] = yd[i];
	}

	//third round
	//0 2
	device_p3_bigPrimeAdd_correct (&xd[0], &xd[16], &yd[0]);
	device_p3_bigSub (&xd[0], &xd[16], &yd[16]);

	//8 10
	device_p3_cyclicShift (&xd[80], 2);
	device_p3_bigPrimeAdd_correct (&xd[64], &xd[80], &yd[64]);
	device_p3_bigSub (&xd[64], &xd[80], &yd[80]);

	//4 6
	device_p3_cyclicShift (&xd[48], 4);
	device_p3_bigPrimeAdd_correct (&xd[32], &xd[48], &yd[32]);
	device_p3_bigSub (&xd[32], &xd[48], &yd[48]);

	//12 14
	device_p3_cyclicShift (&xd[112], 6);
	device_p3_bigPrimeAdd_correct (&xd[96], &xd[112], &yd[96]);
	device_p3_bigSub (&xd[96], &xd[112], &yd[112]);

	//1 3
	device_p3_bigPrimeAdd_correct (&xd[8], &xd[24], &yd[8]);
	device_p3_bigSub (&xd[8], &xd[24], &yd[24]);

	//9 11
	device_p3_cyclicShift (&xd[88], 2);
	device_p3_bigPrimeAdd_correct (&xd[72], &xd[88], &yd[72]);
	device_p3_bigSub (&xd[72], &xd[88], &yd[88]);

	//5 7
	device_p3_cyclicShift (&xd[56], 4);
	device_p3_bigPrimeAdd_correct (&xd[40], &xd[56], &yd[40]);
	device_p3_bigSub (&xd[40], &xd[56], &yd[56]);

	//13 15
	device_p3_cyclicShift (&xd[120], 6);
	device_p3_bigPrimeAdd_correct (&xd[104], &xd[120], &yd[104]);
	device_p3_bigSub (&xd[104], &xd[120], &yd[120]);
	for (i = 0; i < 128; i++)
	{
		xd[i] = yd[i];
	}

	//fourth
	//0 1
	device_p3_bigPrimeAdd_correct (&xd[0], &xd[8], &yd[0]);
	device_p3_bigSub (&xd[0], &xd[8], &yd[8]);

	//8 9
	device_p3_cyclicShift (&xd[72], 1);
	device_p3_bigPrimeAdd_correct (&xd[64], &xd[72], &yd[64]);
	device_p3_bigSub (&xd[64], &xd[72], &yd[72]);

	//4 5
	device_p3_cyclicShift (&xd[40], 2);
	device_p3_bigPrimeAdd_correct (&xd[32], &xd[40], &yd[32]);
	device_p3_bigSub (&xd[32], &xd[40], &yd[40]);

	//12 13
	device_p3_cyclicShift (&xd[104], 3);
	device_p3_bigPrimeAdd_correct (&xd[96], &xd[104], &yd[96]);
	device_p3_bigSub (&xd[96], &xd[104], &yd[104]);

	//2 3
	device_p3_cyclicShift (&xd[24], 4);
	device_p3_bigPrimeAdd_correct (&xd[16], &xd[24], &yd[16]);
	device_p3_bigSub (&xd[16], &xd[24], &yd[24]);

	//10 11
	device_p3_cyclicShift (&xd[88], 5);
	device_p3_bigPrimeAdd_correct (&xd[80], &xd[88], &yd[80]);
	device_p3_bigSub (&xd[80], &xd[88], &yd[88]);

	//6 7
	device_p3_cyclicShift (&xd[56], 6);
	device_p3_bigPrimeAdd_correct (&xd[48], &xd[56], &yd[48]);
	device_p3_bigSub (&xd[48], &xd[56], &yd[56]);

	//14 15
	device_p3_cyclicShift (&xd[120], 7);
	device_p3_bigPrimeAdd_correct (&xd[112], &xd[120], &yd[112]);
	device_p3_bigSub (&xd[112], &xd[120], &yd[120]);

	//0 0
	j = 0;
	for (i = 0; i < 8; i++)
	{
		xd[i] = yd[j++];
	}
	//1 8
	j = 64;
	for (i = 8; i < 16; i++)
	{
		xd[i] = yd[j++];
	}
	//2 4
	j = 32;
	for (i = 16; i < 24; i++)
	{
		xd[i] = yd[j++];
	}
	//3 12
	j = 96;
	for (i = 24; i < 32; i++)
	{
		xd[i] = yd[j++];
	}
	//4 2
	j = 16;
	for (i = 32; i < 40; i++)
	{
		xd[i] = yd[j++];
	}
	//5 10
	j = 80;
	for (i = 40; i < 48; i++)
	{
		xd[i] = yd[j++];
	}
	//6 6 omit
	j = 48;
	for (i = 48; i < 56; i++)
	{
		xd[i] = yd[j++];
	}
	//7 14
	j = 112;
	for (i = 56; i < 64; i++)
	{
		xd[i] = yd[j++];
	}
	//8 1
	j = 8;
	for (i = 64; i < 72; i++)
	{
		xd[i] = yd[j++];
	}
	//9 9 omit
	j = 72;
	for (i = 72; i < 80; i++)
	{
		xd[i] = yd[j++];
	}
	//10 5
	j = 40;
	for (i = 80; i < 88; i++)
	{
		xd[i] = yd[j++];
	}
	//11 13
	j = 104;
	for (i = 88; i < 96; i++)
	{
		xd[i] = yd[j++];
	}
	//12 3
	j = 24;
	for (i = 96; i < 104; i++)
	{
		xd[i] = yd[j++];
	}
	//13 11
	j = 88;
	for (i = 104; i < 112; i++)
	{
		xd[i] = yd[j++];
	}
	//14 7
	j = 56;
	for (i = 112; i < 120; i++)
	{
		xd[i] = yd[j++];
	}

	//15 15
	j = 120;
	for (i = 120; i < 128; i++)
	{
		xd[i] = yd[j++];
	}

}

/**********************************************/
__global__ void
kernel_p3_fft256_plain_2 (usfixn64 *xs, usfixn64* parameters)
{
	//1 block 256 threads
	usfixn64 i, h, w, pos, pos2;
	usfixn64 tid = threadIdx.x + blockDim.x * blockIdx.x;
	usfixn64 *xd;
	__shared__ usfixn64 xsm[INC * 8];
	__shared__ usfixn64 ysm[INC * 8];
	usfixn64 nCols;

	//0. import data
	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * 8);
//	pos = tid << 3; //tid*8
	pos = threadIdx.x << 3; //tid*8
	for (i = 0; i < 8; i++)
	{
		xsm[pos++] = xd[i];
//		xsm[pos++] = xs[tid*8+i];
	}
	__syncthreads ();  //copy the input to shared memory.

	//1. transposition
	transpose (xsm, ysm);  //transposition xsm into ysm
	__syncthreads ();

	//2. fft16
//	device_p3_fft16_plain(ysm, xsm); //calculate FFT16 of ysm. The final result is also stored in ysm. xsm is a temp array.
	device_p3_fft16_plain_3 (ysm, xsm, parameters); //calculate FFT16 of ysm. The final result is also stored in ysm. xsm is a temp array.
	__syncthreads ();

	//3. multiple w^{jk} and transposition
	transpose2 (ysm, xsm); //calculate ysm*w^{jk} and store the result into xsm
//	transpose2(xsm, ysm); //calculate ysm*w^{jk} and store the result into xsm
	__syncthreads ();

	//4. fft16 again
//	device_p3_fft16_plain(xsm, ysm); //calculate FFT16 of xsm. The final result is also stored in xsm. ysm is a temp array.
//	device_p3_fft16_plain_2(xs, xsm, ysm, parameters); //calculate FFT16 of xsm. The final result is also stored in xsm. ysm is a temp array.
	device_p3_fft16_plain_3 (xsm, ysm, parameters); //calculate FFT16 of xsm. The final result is also stored in xsm. ysm is a temp array.
	__syncthreads ();

	transpose (xsm, ysm);  //transposition xsm into ysm
	__syncthreads ();

	//5. write back, from xsm to xs
//	h = tid >> 4; //h=tid/16
//	w = tid - (h << 4); //w=h mod 16
	nCols = parameters[7] >> 4; // parameters[7]/16;
	h = tid >> 4; //h=tid/16
	w = tid - (h << 4); //w=h mod 16
//	pos2 = ((w * nCols) + h) << 3;
	pos2 = ((w << 4) + h) << 3;
	pos2 = (w * 16 + h) * 8;
	pos = threadIdx.x << 3;

	for (i = 0; i < 8; i++)
	{
//		xs[pos2++] = pos2++;//xsm[pos++]; //main line

//		xs[pos++] = xsm[pos++];
//		xs[(COEFFICIENT_SIZE * tid) + i] = pos++;//xsm[pos++];
		xs[(COEFFICIENT_SIZE * tid) + i] = ysm[pos++];
//		xs[pos2++] = pos++;
	}
}

/**********************************************/
__global__ void
kernel_p3_fft16_plain_2 (usfixn64 *xs, usfixn64 *parameters)
{
	__shared__ usfixn64 xsm[INC << 3];
	__shared__ usfixn64 ysm[INC << 3];

	device_p3_fft16_plain_2 (xs, xsm, ysm, parameters);
}

/**********************************************/

#endif
