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


#ifndef BIG_ARITHMETIC_FFT_DEVICE_H_
#define BIG_ARITHMETIC_FFT_DEVICE_H_

#include "../../include/BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "bigArithmetic_addition_P3_device.h"
#include "bigArithmetic_multiplication_P3_device.h"
#include "bigArithmetic_subtraction_P3_device.h"
#include "bigArithmetic_cyclicShift_P3_device.h"

extern __shared__ usfixn64 sharedMem[];

/**********************************************/
__device__ __inline__ void
device_p4_sub_single (usfixn64 &x, const usfixn64 y, short & carry)
{
//	short i, pos;
	short c = 0;
	usfixn64 num1;
	usfixn64 s = 0;
	num1 = 0;

	//step0
	num1 = y + carry;
//		if (xm[i] < num1) //there is not enough to do subtraction

	if (x < num1)
	{
		c = 1;
		s = R - num1 + x;
	}
	else
	{
		c = 0;
		s = x - num1;
	}
	x = s;
	carry = c;
}
/**********************************************/
__device__ __inline__ void
device_p4_rotate_permutated_shmem_step1 (usfixn64* __restrict__ xs,
																				 usfixn64 * __restrict__ tsh, short sn,
																				 usfixn64 permutationStride)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	usfixn64 permutationStride = parameters[5];
	tid = 0;
	usfixn64 offset = tid;
	usfixn16 newIdx = sn;
	int i, j, k;
//	uvector8 t;
	uvector4 t;
	usfixn64 stride = gridDim.x * blockDim.x;
	if (sn <= 0)
	{
		return;
	}
	if (sn > COEFFICIENT_SIZE)
	{
		return;
	}
	{
		newIdx = sn;
		offset = tid;
		//#pragma unroll
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			//#pragma unroll
			//			for (j = 0; j < COEFFICIENT_SIZE; j++)
			j = newIdx;
			{
				//				if (j == newIdx && j < COEFFICIENT_SIZE / 2)
				if (j < COEFFICIENT_SIZE / 2)
					//					t[j] = xs[offset];
					//					setUvector16Element (t, j, xs[offset]);
//					setUvector8Element (t, j, xs[offset]);
					setUvector4Element (t, j, xs[offset]);
				//				if (j == newIdx && j >= COEFFICIENT_SIZE / 2)
				else
					tsh[j - COEFFICIENT_SIZE / 2] = xs[offset];
				//					setUvector8Element (tsh, j-COEFFICIENT_SIZE/2, xs[offset]);
			}
			xs[offset] = 0;
			newIdx++;
			//			newIdx = newIdx & 0x1F;
//			newIdx = newIdx & 0xF;
//			newIdx = newIdx & 0x7;
			newIdx %= COEFFICIENT_SIZE;
			offset += permutationStride;
		}

		offset = tid;
		/***************************/
		//		for (i = 0; i < sn; i++)
		//		{
		//			if (i < COEFFICIENT_SIZE / 2)
		////				ys[offset] = t[i];
		////				ys[offset] = getUvector8Element (t, i);
		//				xs[offset] = getUvector8Element (t, i);
		//			else
		////				ys[offset] = tsh[threadIdx.x][i - COEFFICIENT_SIZE / 2];
		//				xs[offset] = tsh[threadIdx.x][i - COEFFICIENT_SIZE / 2];
		////				ys[offset] = getUvector8Element (tsh, i-COEFFICIENT_SIZE/2);
		//			offset += permutationStride;
		//		}
		for (i = 0; i < sn; i++)
		{
			if (i < COEFFICIENT_SIZE / 2)
			{
				//				ys[offset] = t[i];
				//				ys[offset] = getUvector8Element (t, i);
//				xs[offset] = getUvector8Element (t, i);
				xs[offset] = getUvector4Element (t, i);
			}
			else
			{
				//				ys[offset] = tsh[threadIdx.x][i - COEFFICIENT_SIZE / 2];
				xs[offset] = tsh[i - COEFFICIENT_SIZE / 2];
				//				ys[offset] = getUvector8Element (tsh, i-COEFFICIENT_SIZE/2);
			}
			offset += permutationStride;
		}

		/***************************/

		for (i = sn; i < COEFFICIENT_SIZE; i++)
		{
			if (i < COEFFICIENT_SIZE / 2)
				//				xs[offset] = t[i];
				xs[offset] = getUvector4Element (t, i);
			else
				xs[offset] = tsh[i - COEFFICIENT_SIZE / 2];
			//				xs[offset] = getUvector8Element (tsh, i-COEFFICIENT_SIZE/2);
			offset += permutationStride;
		}
//		for (i = sn; i < COEFFICIENT_SIZE / 2; i++)
//		{
//			//				xs[offset] = t[i];
////			xs[offset] = getUvector8Element (t, i);
//			xs[offset] = getUvector4Element (t, i);
//			offset += permutationStride;
//		}
//		for (i = COEFFICIENT_SIZE / 2; i < COEFFICIENT_SIZE; i++)
//		{
//			xs[offset] = tsh[i - COEFFICIENT_SIZE / 2];
//			//				xs[offset] = getUvector8Element (tsh, i-COEFFICIENT_SIZE/2);
//			offset += permutationStride;
//		}
	}
}
/**********************************************/
__device__ __inline__ void
device_p4_rotate_permutated_shmem_step2 (usfixn64 * __restrict__ xs, short sn,
																				 const usfixn64 permutationStride)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	tid = 0;
	usfixn64 offset = tid;
//	usfixn64 n = parameters[7];
//	int j = threadIdx.x;
//	short sn;
	int i;
	short c = 0;
	usfixn64 tmp = 0;
	short pos;
	usfixn64 tmpX;

//	sn = 6;
	usfixn64 stride = blockDim.x * gridDim.x;

//	for (tid; tid < n; tid += stride)
	{
//		sn = tid & 0x7;
		if (sn <= 0)
		{
			return;
//			continue;
		}
		if (sn > COEFFICIENT_SIZE)
		{
			return;
//			continue;
		}

		c = 0;
		pos = -1;
		offset = tid;

//		for (i = 0; i < sn; i++)
//		{
//			tmp = 0;
////		device_p4_sub_single (tmp, ys[offset], c);
//			device_p4_sub_single (tmp, xs[offset], c);
//			xs[offset] = tmp;
//			offset += permutationStride;
//		}
//
//		for (i = sn; i < COEFFICIENT_SIZE; i++)
//		{
//			tmp = 0;
//			device_p4_sub_single (xs[offset], tmp, c);
//			offset += permutationStride;
//		}
//

		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			tmpX = xs[offset];
			if (i < sn)
			{
				tmp = 0;
//		device_p4_sub_single (tmp, ys[offset], c);
				device_p4_sub_single (tmp, tmpX, c);
				tmpX = tmp;

			}

			else
			{
				tmp = 0;
				device_p4_sub_single (tmpX, tmp, c);
			}
			xs[offset] = tmpX;
			offset += permutationStride;
		}

//	return;
		if (c == 0)
			return;
//			continue;
//		continue;
		pos = -1;
		{
			offset = tid;
			for (i = 0; i < COEFFICIENT_SIZE; i++)
			{
				if (xs[offset] <= R_MINUS_ONE)
				{
					pos = i;
					break;
				}
				offset += permutationStride;
			}

			offset = tid;
			if (pos >= 0)
			{

				for (i = 0; i < pos; i++)
				{
					xs[offset] = 0;
					offset += permutationStride;
				}
//				offset = pos * permutationStride;
				//			um[pos]++;
				xs[offset]++;
			}
			else
			{
				//			um[0] = ULMAX;
				xs[offset] = ULMAX;
				offset = permutationStride;
				for (i = 1; i < COEFFICIENT_SIZE; i++)
				{
					//				um[i] = 0;
					xs[offset] = 0;
					offset += permutationStride;
				}
			}
		}
	}
}
/**********************************************/
/**********************************************/
__device__ __inline__ void
mulOmegaPower (usfixn64 *xs, short po) //multiplication by powers of omega (w), where w^16=r
{
	//w^(N/2k)=r and k=8 -> w^(256/16)=w^(16)=r
	short rp, wp;
	if (po <= 0)
	{
		return;
	}
	rp = po >> 4; //w^16=r
	wp = po - (rp << 4);

	rp > 8 ? device_p3_cyclicShift (xs, 8) : device_p3_cyclicShift (xs, rp);
	rp > 8 ? device_p3_cyclicShift (xs, rp - 8) : device_p3_cyclicShift (xs, 0);

	if (wp > 0)
	{
		bigMult (xs, W_0[wp - 1], xs);
	}

}

/**********************************************/
__device__ __inline__ void
transpose (usfixn64 *idata, usfixn64 *odata) //stride permutation for stride=16, transposes (h,w) to (w,h)
{
	usfixn64 *ip, *op;
	short i, pos, tid, w, h, pos2;

	tid = threadIdx.x;
	h = tid >> 4; //h=tid/16
	w = tid - (h << 4); //w=h mod 16
	pos = ((w << 4) + h) << 3;
	pos2 = tid << 3;

	ip = (usfixn64*) ((char*) idata + pos * sizeof(usfixn64)); //iregular
	op = (usfixn64*) ((char*) odata + pos2 * sizeof(usfixn64)); //regular

	for (i = 0; i < 8; i++)
	{
		op[i] = ip[i];
	}
	__syncthreads ();
}

/**********************************************/
__device__ __inline__ void
transpose2 (usfixn64 *idata, usfixn64 *odata) //stride permutation for stride=16, transposes (h,w) to (w,h), also multiplies by twiddle factor
{
	usfixn64 *ip, *op;
	short i, pos, tid, w, h, pos2;

	tid = threadIdx.x;
	pos = tid << 3;
	h = tid >> 4; //h=tid/16
	w = tid - (h << 4); //w=h mod 16
	pos2 = ((w << 4) + h) << 3;

	ip = (usfixn64*) ((char*) idata + pos * sizeof(usfixn64)); //iregular
	mulOmegaPower (ip, h * w);

	__syncthreads ();
	op = (usfixn64*) ((char*) odata + pos2 * sizeof(usfixn64)); //regular
	for (i = 0; i < 8; i++)
	{
		op[i] = ip[i];
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_fft16_plain (usfixn64 *xsm, usfixn64 *ysm)
{
	//Now, we only use one block to do $16$ FFT 16.
	short tid = threadIdx.x;
	short wid = (tid >> 5);  //warp no. [0,1,...,7]
	short bwid = (tid >> 6); //big warp no.[0,1,...,3]
	short sf = wid - (bwid << 1); //sub flag for odd warp[0,1]
	short wpn = ((tid - (wid << 5)) >> 3); //(tid-(wid*32))/8  [0,1,2,3]
	short wpid = (tid - ((tid >> 3) << 3)); //in [0,1,...,7]
	short posid = 0; //in[0,1,...,15]
	short i, pos, pos1;
	usfixn64 *xf, *yf;

	xf =
			(usfixn64*) ((char*) xsm + (bwid * 64 + wpn * 16) * sizeof(usfixn64) * 8); //
	yf =
			(usfixn64*) ((char*) ysm + (bwid * 64 + wpn * 16) * sizeof(usfixn64) * 8);

	//first round
	if (sf)
	{
		//device_p3_bigSub
		device_p3_bigSub (&xf[wpid * 8], &xf[(wpid + 8) * 8], &yf[(wpid + 8) * 8]);
	}
	else
	{
		//bigAdd
		device_p3_bigPrimeAdd_correct (&xf[wpid * 8], &xf[(wpid + 8) * 8],
																	 &yf[wpid * 8]);
	}
	__syncthreads ();
	if (sf)
	{
		pos = (wpid + 8) * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	else
	{
		pos = wpid * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	__syncthreads ();

	//second round
	posid = (wpid >= 4 ? wpid + 4 : wpid);
	if (sf > 0)
	{
		device_p3_cyclicShift (&xf[(posid + 4) * 8], ind2[wpid]);
	}
	__syncthreads ();
	if (sf)
	{
		//device_p3_bigSub
		device_p3_bigSub (&xf[posid * 8], &xf[(posid + 4) * 8],
											&yf[(posid + 4) * 8]);
	}
	else
	{
		//bigAdd
		device_p3_bigPrimeAdd_correct (&xf[posid * 8], &xf[(posid + 4) * 8],
																	 &yf[posid * 8]);
	}
	__syncthreads ();
	if (sf)
	{
		pos = (wpid + 8) * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	else
	{
		pos = wpid * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	__syncthreads ();

	//third round
	posid = (wpid >= 4 ? (wpid - 4) * 4 + 1 : wpid * 4);
	if (sf > 0)
	{
		device_p3_cyclicShift (&xf[(posid + 2) * 8], ind3[wpid]);
	}

	__syncthreads ();
	if (sf > 0)
	{
		//device_p3_bigSub
		device_p3_bigSub (&xf[posid * 8], &xf[(posid + 2) * 8],
											&yf[(posid + 2) * 8]);
	}
	else
	{
		//bigAdd
		device_p3_bigPrimeAdd_correct (&xf[posid * 8], &xf[(posid + 2) * 8],
																	 &yf[posid * 8]);
	}
	__syncthreads ();
	if (sf)
	{
		pos = (wpid + 8) * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	else
	{
		pos = wpid * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos + i];
		}
	}
	__syncthreads ();

	//fourth
	posid = (wpid << 1);
	if (sf > 0)
	{
		device_p3_cyclicShift (&xf[(posid + 1) * 8], ind4[wpid]);
	}
	__syncthreads ();
	if (sf)
	{
		//device_p3_bigSub
		device_p3_bigSub (&xf[posid * 8], &xf[(posid + 1) * 8],
											&yf[(posid + 1) * 8]);
	}
	else
	{
		//bigAdd
		device_p3_bigPrimeAdd_correct (&xf[posid * 8], &xf[(posid + 1) * 8],
																	 &yf[posid * 8]);
	}
	__syncthreads ();

	if (sf)
	{
		posid = wpid + 8;
		pos = posid * 8;
		pos1 = ind5[posid] * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos1 + i];
		}
	}
	else
	{
		posid = wpid;
		pos = posid * 8;
		pos1 = ind5[posid] * 8;
		for (i = 0; i < 8; i++)
		{
			xf[pos + i] = yf[pos1 + i];
		}
	}
	__syncthreads ();
}

/**********************************************/
__device__ __inline__ void
device_p3_fft16_plain_2 (usfixn64* xs, usfixn64* xsm, usfixn64 *ysm,
												 usfixn64 *parameters)
{
	usfixn64 bid = blockIdx.x;
	short tid = threadIdx.x;
	//			+ blockIdx.x*blockDim.x;

	//	int needBlock = (N >> 8); //N/256;
	usfixn64 needBlock = (parameters[7] >> 8); //N/256;
	//	short bn = ((needBlock >= BLOCK_DIM) ? needBlock / BLOCK_DIM : 1); //the number of block
	usfixn64 bn = ((needBlock >= MAX_DIM) ? needBlock / MAX_DIM : 1); //the number of block
	//	short i, j;
	usfixn64 i, j;
	//usfixn64 *A, *B;
	usfixn64 *A;
//		__shared__ usfixn64 xsm[INC << 3];
//		__shared__ usfixn64 ysm[INC << 3];
	int pos, pos1, bno; //tid*8

	bno = bid;
	pos = ((bno << 8) + tid) << 3;
	for (i = 0; i < bn; i++)
	{
		if (bno >= needBlock)
		{
			break;
		}
		pos1 = tid << 3;
		//pos = (bno<<8)+tid;  //each block deal 256 big numbers.
		A = (usfixn64 *) ((char*) xs + pos * sizeof(usfixn64));

		for (j = 0; j < 8; j++)
		{
			xsm[pos1++] = A[j];
		}
		__syncthreads ();

		device_p3_fft16_plain (xsm, ysm); //ysm is a temp array
		__syncthreads ();

		pos1 = (tid << 3);
		for (j = 0; j < 8; j++)
		{
			A[j] = xsm[pos1++];
		}
		bno += MAX_DIM;
		//		pos += 1048576;  //512block*256bignumber*8numbers
		pos += (1 << 26);  //maxDim block*256bignumber*8numbers = 2^26
	}
	__syncthreads ();
}

/**********************************************/
__device__ __inline__ void
device_p3_fft16_plain_3 (usfixn64* xsm, usfixn64 *ysm, usfixn64 *parameters)
{
	usfixn64 bid = blockIdx.x;
	short tid = threadIdx.x;
	//			+ blockIdx.x*blockDim.x;

	//	int needBlock = (N >> 8); //N/256;
	usfixn64 needBlock = (parameters[7] >> 8); //N/256;
	//	short bn = ((needBlock >= BLOCK_DIM) ? needBlock / BLOCK_DIM : 1); //the number of block
	usfixn64 bn = ((needBlock >= MAX_DIM) ? needBlock / MAX_DIM : 1); //the number of block
	//	short i, j;
	usfixn64 i, j;
	//usfixn64 *A, *B;
//		usfixn64 *A;
//		__shared__ usfixn64 xsm[INC << 3];
//		__shared__ usfixn64 ysm[INC << 3];
	int pos, pos1, bno; //tid*8

	bno = bid;
	pos = ((bno << 8) + tid) << 3;
	for (i = 0; i < bn; i++)
	{
		if (bno >= needBlock)
		{
			break;
		}
//			pos1 = tid << 3;
		//pos = (bno<<8)+tid;  //each block deal 256 big numbers.
//			A = (usfixn64 *) ((char*) xs + pos * sizeof(usfixn64));

//			for (j = 0; j < 8; j++)
//			{
//				xsm[pos1++] = A[j];
//			}
//			__syncthreads();

		device_p3_fft16_plain (xsm, ysm); //ysm is a temp array
		__syncthreads ();

//			pos1 = (tid << 3);
//			for (j = 0; j < 8; j++)
//			{
//				A[j] = xsm[pos1++];
//			}
		bno += MAX_DIM;
		//		pos += 1048576;  //512block*256bignumber*8numbers
		pos += (1 << 26);  //maxDim block*256bignumber*8numbers = 2^26
	}
	__syncthreads ();
}

/**********************************************/
__device__ __inline__ void
device_p3_swap_permutated (usfixn64 * __restrict__ xm,
													 usfixn64 * __restrict__ ym,
													 const usfixn64 permutationStride,
													 const usfixn64 coefficientSize)
{
	short i = 0;
	usfixn64 tmp = 0;
	usfixn64 offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		tmp = xm[offset];
		xm[offset] = ym[offset];
		ym[offset] = tmp;
		offset += permutationStride;
	}
}

/**********************************************/
__device__ __inline__ void
device_p3_swap_permutated_2 (usfixn64 * __restrict__ xm,
														 usfixn64 * __restrict__ ym,
														 const usfixn64 & permutationStride)
{
	short i = 0;
	usfixn64 tmp = 0;
	usfixn64 offset = 0;

	usfixn64 x[8], y[8];
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		tmp = xm[offset];
//		xm[offset] = ym[offset];
//		ym[offset] = tmp;
//		offset += permutationStride;
//	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		x[i] = xm[offset];
		offset += permutationStride;
	}

	offset = 0;
//	#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		y[i] = ym[offset];
		offset += permutationStride;
	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		ym[offset] = x[i];
		offset += permutationStride;
	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		xm[offset] = y[i];
		offset += permutationStride;
	}
}

/**********************************************/
//new code for performing fft16, 8 threads for fft-16
__device__ __inline__ void
device_base_fft16 (usfixn64 *xm, const usfixn64 permutationStride,
									 const usfixn64 tid)
{
	//short rounds = 4; //log(16,2)

	short stride = 16; //stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	device_p3_cyclicShift (&xm[idx + stride], shNo_FFT_16_r2[shNo]);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	device_p3_cyclicShift (&xm[idx + stride], shNo_FFT_16_r3[shNo]);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
	device_p3_cyclicShift (&xm[idx + stride], shNo_FFT_16_r4[shNo]);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;
}

/**********************************************/
//new code for performing fft16, 8 threads for fft-16
__device__ inline void
device_base_fft16_2 (usfixn64 * __restrict__ xm,
										 const usfixn64 permutationStride, const usfixn64 tid)
{
	//short rounds = 4; //log(16,2)

	short stride = 8; //stride=16/2
	usfixn64 idx = 0;

	short lstride = 3; //log 8
//	short i = tid % stride;
//	short shNo = tid % stride;
	short i = (tid & 0x7);
	short shNo = i;
//	short i = tid & (1 << lstride);
//	short shNo = tid & (1 << lstride);

	//round 1
//	idx = (tid / stride) * (stride << 1) + i;
	idx = ((tid >> lstride) << (lstride + 1)) + i;
	//no cyclic shift in round 1
//	fft_base2_permutated_4(&xm[idx], &xm[idx + stride], permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	lstride--;
//	i = i - (i >= stride) * stride;
	i = i - ((i >= stride) << lstride);

	//round 2
//	idx = (tid / stride) * (stride << 1) + i;
	idx = ((tid >> lstride) << (lstride + 1)) + i;

//	device_cyclicShift(&xm[idx + stride], shNo_FFT_16_r2[shNo]);
	device_p3_cyclicShift_permutated_5 (&xm[idx + stride], shNo_FFT_16_r2[shNo],
																			permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	lstride--;
//	i = i - (i >= stride) * stride;
	i = i - ((i >= stride) << lstride);

	//round 3
//	idx = (tid / stride) * (stride << 1) + i;
	idx = ((tid >> lstride) << (lstride + 1)) + i;
//	device_p3_cyclicShift(&xm[idx + stride], shNo_FFT_16_r3[shNo]);
	device_p3_cyclicShift_permutated_5 (&xm[idx + stride], shNo_FFT_16_r3[shNo],
																			permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + stride], permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
	stride >>= 1; //stride /=2;
	lstride--;
//	i = i - (i >= stride) * stride;
	i = i - ((i >= stride) << lstride);

	//round 4
//	idx = (tid / stride) * (stride << 1) + i;
	idx = ((tid >> lstride) << (lstride + 1)) + i;
//	device_p3_cyclicShift(&xm[idx + stride], shNo_FFT_16_r4[shNo]);
	device_p3_cyclicShift_permutated_5 (&xm[idx + stride], shNo_FFT_16_r4[shNo],
																			permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + stride], permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + stride],
																	permutationStride);
//	stride >>= 1; //stride /=2;
//	lstride--; is already 0 here
//	i = i - (i >= stride) * stride;
//	xm[idx]=tid;
//	short r, c;
//	r = idx / (8);
//	c = idx % (2);

	//	write a permutated_swap for this last step
	//last permutation 4 for fft16
//	shNo = (c * 2) + r;
//	usfixn64 t;
//	usfixn64 offset = 0;

//	usfixn64 coeffSize=8;
//	shNo = tid % 8;
//	/********************************************************************/
//	shNo = (tid & 0x7);
//	shNo = shNo_FFT_16_r5[shNo];
////	offset=(idx/16)*16;
////	device_swap(&xm[idx], &xm[idx + shNo], permutationStride, coeffSize);
//	if (shNo > 0)
//		device_p3_swap_permutated(&xm[idx], &xm[idx + shNo], permutationStride,
//				COEFFICIENT_SIZE);
//
////	shNo = idx % 16 + 1;
////	shNo = shNo_FFT_16_r5[shNo];
//
////	shNo = idx % 16+1;
//	shNo = (tid & 0x7) + 8;
//	shNo = shNo_FFT_16_r5[shNo];
//	idx++;
//	if (shNo > 0)
//		device_p3_swap_permutated(&xm[idx], &xm[idx + shNo], permutationStride,
//				COEFFICIENT_SIZE);
//	/********************************************************************/

	/********************************************************************/
	shNo = (tid & 0x7);
	i = shNo + 8;
	shNo = shNo_FFT_16_r5[shNo];
	//	offset=(idx/16)*16;
	//	device_swap(&xm[idx], &xm[idx + shNo], permutationStride, coeffSize);
	if (shNo > 0)
		device_p3_swap_permutated (&xm[idx], &xm[idx + shNo], permutationStride,
		COEFFICIENT_SIZE);

	//	shNo = idx % 16 + 1;
	//	shNo = shNo_FFT_16_r5[shNo];

	//	shNo = idx % 16+1;
//		shNo = (tid & 0x7) + 8;
	shNo = shNo_FFT_16_r5[i];
	idx++;
	if (shNo > 0)
		device_p3_swap_permutated (&xm[idx], &xm[idx + shNo], permutationStride,
		COEFFICIENT_SIZE);
	/********************************************************************/
////	device_swap(&xm[idx], &xm[idx + shNo], permutationStride, coeffSize);
//	device_swap(&xm[idx], &xm[offset + shNo], permutationStride, coeffSize);
//	for (i = 0; i < 8; i++)
//	{
//		t = xm[idx + offset + shNo];
//		xm[idx + offset + shNo] = xm[idx + offset];
//		xm[idx + offset] = t;
//		offset += permutationStride;
//	}
//
//	shNo = idx % 16 + 1;
//	shNo = shNo_FFT_16_r5[shNo];
//	for (i = 0; i < 8; i++)
//	{
//		t = xm[idx + offset + shNo];
//		xm[idx + offset + shNo] = xm[idx + offset];
//		xm[idx + offset] = t;
//		offset += permutationStride;
//	}
}

/**********************************************/
//new code for performing fft16, 8 threads for fft-16
//like fft_16_4 but accesses shared memory for reading shNo
__device__ __inline__ void
device_base_fft16_6 (usfixn64 * __restrict__ xm,
										 const usfixn64 & permutationStride)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	const short * __restrict__ shNo_FFT_16_r2_shmem,
//			const short * __restrict__ shNo_FFT_16_r3_shmem,
//			const short * __restrict__ shNo_FFT_16_r4_shmem,
//			const short *__restrict__ shNo_FFT_16_r5_shmem,
	//short rounds = 4; //log(16,2)

//	uConstShortArray8_align8 shNo_FFT_16_r2_shmem;
	usfixn64 shiftVector8;
	usfixn64 shiftValue;
	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
//	shNoBigShort = 0x0000000;
//	shNoBigShort = ((4<<40)|(4<<48)|(4<<56));
//	shNoBigShort = 32;
//	shNoBigChar = ((1<<8)|(1<<10)|(1<<12)|(1<<14));//0,0,0,0,4,4,4,4
	short shNoShort = 0;

//	short stride = 8; //stride=16/2
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;

//	short lstride = 3; //log 8
//	short i = tid % stride;
//	short shNo = tid % stride;
//	short i = (tid & 0x7);
	short i = (threadIdx.x & 0x7);
	short j = 0;
	short shNo = i;
	shNo = shNo << 3;
//	short shNo = (tid & 0x7);
//	short i = tid & (1 << lstride);
//	short shNo = tid & (1 << lstride);

	//round 1
//	idx = (tid / stride) * (stride << 1) + i;
//	idx = ((tid >> 3) << (3 + 1)) + i;
//	idx = ((tid >> 3) << (4)) + i;
//	idx = ((tid & 0xffff8) << (4)) + i;
	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;
	//no cyclic shift in round 1
//	fft_base2_permutated_4(&xm[idx], &xm[idx + 8], permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx2], permutationStride);
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx2], permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx+8], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx2], shNoShort, permutationStride);

//#################################
//	cyclic shift round 1

	j = 0;
	uvector8 ts, ys;
//	usfixn64 offset = 0;
//	idx += 8;
//	ys.i0 = xm[idx];
//	idx += permutationStride;
//	ys.i1 = xm[idx];
//	idx += permutationStride;
//	ys.i2 = xm[idx];
//	idx += permutationStride;
//	ys.i3 = xm[idx];
//	idx += permutationStride;
//	ys.i4 = xm[idx];
//	idx += permutationStride;
//	ys.i5 = xm[idx];
//	idx += permutationStride;
//	ys.i6 = xm[idx];
//	idx += permutationStride;
//	ys.i7 = xm[idx];
//
//	idx -= (permutationStride * 7);
//
//	ts.i0 = 0;
//	ts.i1 = 0;
//	ts.i2 = 0;
//	ts.i3 = 0;
//	ts.i4 = 0;
//	ts.i5 = 0;
//	ts.i6 = 0;
//	ts.i7 = 0;
//
//	j = 8 - shNoShort;
//	for (i = 0; i < shNoShort; i++)
//	{
//		setUvector8Element(ts, i, getUvector8Element(ys, j));
//		j++;
//	}
//
//	for (i = 7 - shNoShort; i >= 0; i--)
//	{
//		j = i + shNoShort;
//		setUvector8Element(ys, j, getUvector8Element(ys, i));
//	}
//	for (i = 0; i < shNoShort; i++)
//	{
//		setUvector8Element(ys, i, 0);
//	}
//	device_p3_bigSub_uVector_plain_2(ys, ts);
//
//	{
//		xm[idx] = ys.i0;
//		idx += (permutationStride);
//		xm[idx] = ys.i1;
//		idx += (permutationStride);
//		xm[idx] = ys.i2;
//		idx += (permutationStride);
//		xm[idx] = ys.i3;
//		idx += (permutationStride);
//		xm[idx] = ys.i4;
//		idx += (permutationStride);
//		xm[idx] = ys.i5;
//		idx += (permutationStride);
//		xm[idx] = ys.i6;
//		idx += (permutationStride);
//		xm[idx] = ys.i7;
//	}
//	idx = idx - (permutationStride * 7) - 8;

//#################################
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx + 8], permutationStride);
//#################################################
//	fft-base2 round 1
	short c = 0, c2 = 0;
	usfixn64 num1, num3 = 0;

	num1 = 0;
//	short i;
	//	usfixn64 pos = threadIdx.x;

//	usfixn64 offset = 0;
//	short pos1;
	short pos = 0;
	//	usfixn64 tmp = 0;
//		offset = 0;

//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
//			num1 = xm[offset] + ym[offset] + c;
//			num3 = ym[offset] + c2;
		num1 = xm[idx] + xm[idx + 8] + c;
		num3 = xm[idx + 8] + c2;

		//addition part
		if (num1 < xm[idx] || num1 < xm[idx + 8])	//there is overflow/truncation
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
		if (xm[idx] < num3) //there is not enough to do subtraction
		{
			c2 = 1;
			num3 = R - num3 + xm[idx];
		}
		else
		{
			c2 = 0;
			num3 = xm[idx] - num3;
		}
		xm[idx] = num1;
		xm[idx + 8] = num3;
//			offset += permutationStride;
		idx += permutationStride;
	}
	idx -= (8 * permutationStride);

//		offset = 0;
	if (c > 0)
	{
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] != 0)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}
		if (pos >= 0)	// shouldn't it be >0?
		{
//				offset = 0;

			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = R_MINUS_ONE;
				idx += permutationStride;
			}
			xm[idx]--;
			idx -= pos * permutationStride;
//				idx-=pos1 * permutationStride;
//				offset = pos1 * permutationStride;
			//xm[pos1*permutationStride+idx]--;
		}
		else
		{
			idx -= permutationStride * 8;
			//			xm[0] = ULMAX;
//				xm[offset] = ULMAX;
			xm[idx] = ULMAX;
			idx += permutationStride;
//				offset = permutationStride;
			for (i = 1; i < 8; i++)
			{
				//				xm[i] = 0;
//					xm[offset] = 0;
				xm[idx] = 0;
				idx += permutationStride;
			}
			idx -= (permutationStride * 8);
		}
	}

	idx += 8;
	if (c2 > 0)
	{
//			offset = 0;
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] < R_MINUS_ONE)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}

		if (pos >= 0)
		{
//				offset = 0;
			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = 0;
				idx += permutationStride;
			}
			xm[idx]++;
			idx -= pos * permutationStride;
//				offset = pos * permutationStride;
			//			um[pos]++;
//				ym[offset]++;
		}
		else
		{
//				offset = 0;
			//			um[0] = ULMAX;
//				ym[offset] = ULMAX;
//				offset = permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
			for (i = 1; i < 8; i++)
			{
				//				um[i] = 0;
				xm[idx] = 0;
				idx -= permutationStride;
			}
			xm[idx] = ULMAX;
		}
	}
//#################################################
	i = (tid & (0x07));
	__syncthreads ();
//	stride >>= 1; //stride /=2;
//	lstride--;
//	i = i - (i >= stride) * stride;
//	i = i - ((i >= 4) << 2); //i=i%4;
	i = i & 0x3;

	//round 2
//	idx = (tid / stride) * (stride << 1) + i;
//	idx = ((tid >> 2) << (2 + 1)) + i;
//	idx = ((tid >> 2) << (3)) + i;
	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;

//	shNo_FFT_16_r2_shmem.i[0] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[1] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[2] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[3] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[4] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[5] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[6] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[7] = 4; //0, 0, 0, 0, 4, 4, 4, 4

//	device_cyclicShift(&xm[idx + stride], shNo_FFT_16_r2[shNo]);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 4], shNo_FFT_16_r2[shNo], permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + 4], permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx + 4], shNo_FFT_16_r2_shmem.i[shNo],
//			permutationStride);

//	if (shNo == 5)
	shiftValue = (0x07);
//	shiftValue <<= (shNo << 3);
	shiftValue <<= (shNo);
//	if(tid==5)
	{
//		shNoBigShort >>=shNo;
//		shNoBigShort =;

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
		shNoShort = (shiftVector8 & shiftValue) >> (shNo);
//				>>32;
//		printf("shNo= %llu \n", shNoShort);
//	shNoChar = shNo_FFT_16_r2_shmem.i[shNo];
//		fft_base2_permutated_8(&xm[idx], &xm[idx + 4], shNoShort,
//		permutationStride
//		);
//		device_cyclicShift_permutated_7(&xm[idx2], shNoShort, permutationStride);
		//#################################
		//	cyclic shift round 2

		j = 0;
//		offset = 0;
		idx += 4;
		ys.i0 = xm[idx];
		idx += permutationStride;
		ys.i1 = xm[idx];
		idx += permutationStride;
		ys.i2 = xm[idx];
		idx += permutationStride;
		ys.i3 = xm[idx];
		idx += permutationStride;
		ys.i4 = xm[idx];
		idx += permutationStride;
		ys.i5 = xm[idx];
		idx += permutationStride;
		ys.i6 = xm[idx];
		idx += permutationStride;
		ys.i7 = xm[idx];

		idx -= (permutationStride * 7);

		ts.i0 = 0;
		ts.i1 = 0;
		ts.i2 = 0;
		ts.i3 = 0;
		ts.i4 = 0;
		ts.i5 = 0;
		ts.i6 = 0;
		ts.i7 = 0;

		j = 8 - shNoShort;
		for (i = 0; i < shNoShort; i++)
		{
			setUvector8Element (ts, i, getUvector8Element (ys, j));
			j++;
		}

		for (i = 7 - shNoShort; i >= 0; i--)
		{
			j = i + shNoShort;
			setUvector8Element (ys, j, getUvector8Element (ys, i));
		}
		for (i = 0; i < shNoShort; i++)
		{
			setUvector8Element (ys, i, 0);
		}
		device_p3_bigSub_uVector_plain_2 (ys, ts);

		{
			xm[idx] = ys.i0;
			idx += (permutationStride);
			xm[idx] = ys.i1;
			idx += (permutationStride);
			xm[idx] = ys.i2;
			idx += (permutationStride);
			xm[idx] = ys.i3;
			idx += (permutationStride);
			xm[idx] = ys.i4;
			idx += (permutationStride);
			xm[idx] = ys.i5;
			idx += (permutationStride);
			xm[idx] = ys.i6;
			idx += (permutationStride);
			xm[idx] = ys.i7;
		}
		idx = idx - (permutationStride * 7) - 4;
		//#################################################
		//	fft-base2 round 2
		c = 0;
		c2 = 0;
		num1 = 0;
		num3 = 0;
		pos = 0;
		//	usfixn64 tmp = 0;
		//		offset = 0;

		//#pragma unroll COEFFICIENT_SIZE
		for (i = 0; i < 8; i++)
		{
			//			num1 = xm[offset] + ym[offset] + c;
			//			num3 = ym[offset] + c2;
			num1 = xm[idx] + xm[idx + 4] + c;
			num3 = xm[idx + 4] + c2;

			//addition part
			if (num1 < xm[idx] || num1 < xm[idx + 4])	//there is overflow/truncation
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
			if (xm[idx] < num3) //there is not enough to do subtraction
			{
				c2 = 1;
				num3 = R - num3 + xm[idx];
			}
			else
			{
				c2 = 0;
				num3 = xm[idx] - num3;
			}
			xm[idx] = num1;
			xm[idx + 4] = num3;
			//			offset += permutationStride;
			idx += permutationStride;
		}
		idx -= (8 * permutationStride);

		//		offset = 0;
		if (c > 0)
		{
			pos = -1;
			for (i = 0; i < 8; i++)
			{
				if (xm[idx] != 0)
				{
					pos = i;
					break;
				}
				idx += permutationStride;
			}
			if (pos >= 0)	// shouldn't it be >0?
			{
				//				offset = 0;

				idx -= pos * permutationStride;
				for (i = 0; i < pos; i++)
				{
					xm[idx] = R_MINUS_ONE;
					idx += permutationStride;
				}
				xm[idx]--;
				idx -= pos * permutationStride;
				//				idx-=pos1 * permutationStride;
				//				offset = pos1 * permutationStride;
				//xm[pos1*permutationStride+idx]--;
			}
			else
			{
				idx -= permutationStride * 8;
				//			xm[0] = ULMAX;
				//				xm[offset] = ULMAX;
				xm[idx] = ULMAX;
				idx += permutationStride;
				//				offset = permutationStride;
				for (i = 1; i < 8; i++)
				{
					//				xm[i] = 0;
					//					xm[offset] = 0;
					xm[idx] = 0;
					idx += permutationStride;
				}
				idx -= (permutationStride * 8);
			}
		}

		idx += 4;
		if (c2 > 0)
		{
			//			offset = 0;
			pos = -1;
			for (i = 0; i < 8; i++)
			{
				if (xm[idx] < R_MINUS_ONE)
				{
					pos = i;
					break;
				}
				idx += permutationStride;
			}

			if (pos >= 0)
			{
				//				offset = 0;
				idx -= pos * permutationStride;
				for (i = 0; i < pos; i++)
				{
					xm[idx] = 0;
					idx += permutationStride;
				}
				xm[idx]++;
				idx -= pos * permutationStride;
				//				offset = pos * permutationStride;
				//			um[pos]++;
				//				ym[offset]++;
			}
			else
			{
				//				offset = 0;
				//			um[0] = ULMAX;
				//				ym[offset] = ULMAX;
				//				offset = permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
				for (i = 1; i < 8; i++)
				{
					//				um[i] = 0;
					xm[idx] = 0;
					idx -= permutationStride;
				}
				xm[idx] = ULMAX;
			}
		}
		//#################################################
		i = (tid & (0x07));
		//#################################
//		device_p3_fft_base2_permutated(&xm[idx], &xm[idx + 4], permutationStride);
	}
	__syncthreads ();
//	else
//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNo_FFT_16_r2_shmem.i[shNo],
//				permutationStride);
//	stride >>= 1; //stride /=2;
//	lstride--;
//	i = i - (i >= stride) * stride;
//	i = i - ((i >= 2) << 1);
	i = i & 0x1;

	shiftVector8 = (0x0606020204040000);
//	shiftValue = (0x07);	//should be repeated because of 64-bit registers
//	shiftValue <<= (shNo << 3);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shNo_FFT_16_r2_shmem.i[0] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[1] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[2] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[3] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[4] = 2; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[5] = 2; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[6] = 6; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[7] = 6; //0, 0, 0, 0, 4, 4, 4, 4

	//round 3
//	idx = (tid / stride) * (stride << 1) + i;
//	idx = ((tid >> 1) << (1 + 1)) + i;
//	idx = ((tid >> 1) << (2)) + i;
	idx = ((tid & 0xFFFFFE) << 1) + i;
//	device_p3_cyclicShift(&xm[idx + stride], shNo_FFT_16_r3[shNo]);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 2], shNo_FFT_16_r3[shNo],
//			permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + 2], permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx + 2], shNo_FFT_16_r2_shmem.i[shNo],
//			permutationStride);
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNo_FFT_16_r2_shmem.i[shNo],
//			permutationStride);
//	fft_base2_permutated_8(&xm[idx], &xm[idx + 2], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx2], shNoShort, permutationStride);
	//#################################
	//	cyclic shift round 3

	j = 0;
	idx += 2;
	ys.i0 = xm[idx];
	idx += permutationStride;
	ys.i1 = xm[idx];
	idx += permutationStride;
	ys.i2 = xm[idx];
	idx += permutationStride;
	ys.i3 = xm[idx];
	idx += permutationStride;
	ys.i4 = xm[idx];
	idx += permutationStride;
	ys.i5 = xm[idx];
	idx += permutationStride;
	ys.i6 = xm[idx];
	idx += permutationStride;
	ys.i7 = xm[idx];

	idx -= (permutationStride * 7);

	ts.i0 = 0;
	ts.i1 = 0;
	ts.i2 = 0;
	ts.i3 = 0;
	ts.i4 = 0;
	ts.i5 = 0;
	ts.i6 = 0;
	ts.i7 = 0;

	j = 8 - shNoShort;
	for (i = 0; i < shNoShort; i++)
	{
		setUvector8Element (ts, i, getUvector8Element (ys, j));
		j++;
	}

	for (i = 7 - shNoShort; i >= 0; i--)
	{
		j = i + shNoShort;
		setUvector8Element (ys, j, getUvector8Element (ys, i));
	}
	for (i = 0; i < shNoShort; i++)
	{
		setUvector8Element (ys, i, 0);
	}
	device_p3_bigSub_uVector_plain_2 (ys, ts);

	{
		xm[idx] = ys.i0;
		idx += (permutationStride);
		xm[idx] = ys.i1;
		idx += (permutationStride);
		xm[idx] = ys.i2;
		idx += (permutationStride);
		xm[idx] = ys.i3;
		idx += (permutationStride);
		xm[idx] = ys.i4;
		idx += (permutationStride);
		xm[idx] = ys.i5;
		idx += (permutationStride);
		xm[idx] = ys.i6;
		idx += (permutationStride);
		xm[idx] = ys.i7;
	}
	idx = idx - (permutationStride * 7) - 2;
	//#################################################
	//	fft-base2 round 3
	c = 0;
	c2 = 0;
	num1 = 0;
	num3 = 0;
	pos = 0;
	//	usfixn64 tmp = 0;
	//		offset = 0;

	//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
		//			num1 = xm[offset] + ym[offset] + c;
		//			num3 = ym[offset] + c2;
		num1 = xm[idx] + xm[idx + 2] + c;
		num3 = xm[idx + 2] + c2;

		//addition part
		if (num1 < xm[idx] || num1 < xm[idx + 2])	//there is overflow/truncation
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
		if (xm[idx] < num3) //there is not enough to do subtraction
		{
			c2 = 1;
			num3 = R - num3 + xm[idx];
		}
		else
		{
			c2 = 0;
			num3 = xm[idx] - num3;
		}
		xm[idx] = num1;
		xm[idx + 2] = num3;
		//			offset += permutationStride;
		idx += permutationStride;
	}
	idx -= (8 * permutationStride);

	//		offset = 0;
	if (c > 0)
	{
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] != 0)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}
		if (pos >= 0)	// shouldn't it be >0?
		{
			//				offset = 0;

			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = R_MINUS_ONE;
				idx += permutationStride;
			}
			xm[idx]--;
			idx -= pos * permutationStride;
			//				idx-=pos1 * permutationStride;
			//				offset = pos1 * permutationStride;
			//xm[pos1*permutationStride+idx]--;
		}
		else
		{
			idx -= permutationStride * 8;
			//			xm[0] = ULMAX;
			//				xm[offset] = ULMAX;
			xm[idx] = ULMAX;
			idx += permutationStride;
			//				offset = permutationStride;
			for (i = 1; i < 8; i++)
			{
				//				xm[i] = 0;
				//					xm[offset] = 0;
				xm[idx] = 0;
				idx += permutationStride;
			}
			idx -= (permutationStride * 8);
		}
	}

	idx += 2;
	if (c2 > 0)
	{
		//			offset = 0;
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] < R_MINUS_ONE)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}

		if (pos >= 0)
		{
			//				offset = 0;
			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = 0;
				idx += permutationStride;
			}
			xm[idx]++;
			idx -= pos * permutationStride;
			//				offset = pos * permutationStride;
			//			um[pos]++;
			//				ym[offset]++;
		}
		else
		{
			//				offset = 0;
			//			um[0] = ULMAX;
			//				ym[offset] = ULMAX;
			//				offset = permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
			for (i = 1; i < 8; i++)
			{
				//				um[i] = 0;
				xm[idx] = 0;
				idx -= permutationStride;
			}
			xm[idx] = ULMAX;
		}
	}
	//#################################################
	i = (tid & (0x07));
	//#################################
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx + 2], permutationStride);
	__syncthreads ();
//	stride >>= 1; //stride /=2;
//	lstride--;
//	i = i - (i >= stride) * stride;
//	i = i - ((i >= 1)<<0);
//	i = i - ((i >= 1));
//	i = i & 0x1;

	shiftVector8 = (0x0703050106020400);
//	shiftValue = (0x07);	//should be repeated because of 64-bit registers
//	shiftValue <<= (shNo << 3);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	i = 0;
//	shNo_FFT_16_r2_shmem.i[0] = 0; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[1] = 4; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[2] = 2; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[3] = 6; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[4] = 1; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[5] = 5; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[6] = 3; //0, 0, 0, 0, 4, 4, 4, 4
//	shNo_FFT_16_r2_shmem.i[7] = 7; //0, 0, 0, 0, 4, 4, 4, 4
	//round 4
//	idx = (tid / stride) * (stride << 1) + i;
//	idx = ((tid >> 0) << (0 + 1)) + i;
//	idx = ((tid) << (1)) + i;
	idx = (tid << 1);
//	device_p3_cyclicShift(&xm[idx + stride], shNo_FFT_16_r4[shNo]);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 1], shNo_FFT_16_r4[shNo],
//			permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx + 1], permutationStride);
//	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx + 1], shNo_FFT_16_r2_shmem.i[shNo],
//		permutationStride);
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNo_FFT_16_r2_shmem.i[shNo],
//			permutationStride);
//	fft_base2_permutated_8(&xm[idx], &xm[idx + 1], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx2], shNoShort, permutationStride);
	//#################################
	//	cyclic shift round 4

	j = 0;
	idx += 1;
	ys.i0 = xm[idx];
	idx += permutationStride;
	ys.i1 = xm[idx];
	idx += permutationStride;
	ys.i2 = xm[idx];
	idx += permutationStride;
	ys.i3 = xm[idx];
	idx += permutationStride;
	ys.i4 = xm[idx];
	idx += permutationStride;
	ys.i5 = xm[idx];
	idx += permutationStride;
	ys.i6 = xm[idx];
	idx += permutationStride;
	ys.i7 = xm[idx];

	idx -= (permutationStride * 7);

	ts.i0 = 0;
	ts.i1 = 0;
	ts.i2 = 0;
	ts.i3 = 0;
	ts.i4 = 0;
	ts.i5 = 0;
	ts.i6 = 0;
	ts.i7 = 0;

	j = 8 - shNoShort;
	for (i = 0; i < shNoShort; i++)
	{
		setUvector8Element (ts, i, getUvector8Element (ys, j));
		j++;
	}

	for (i = 7 - shNoShort; i >= 0; i--)
	{
		j = i + shNoShort;
		setUvector8Element (ys, j, getUvector8Element (ys, i));
	}
	for (i = 0; i < shNoShort; i++)
	{
		setUvector8Element (ys, i, 0);
	}
	device_p3_bigSub_uVector_plain_2 (ys, ts);

	{
		xm[idx] = ys.i0;
		idx += (permutationStride);
		xm[idx] = ys.i1;
		idx += (permutationStride);
		xm[idx] = ys.i2;
		idx += (permutationStride);
		xm[idx] = ys.i3;
		idx += (permutationStride);
		xm[idx] = ys.i4;
		idx += (permutationStride);
		xm[idx] = ys.i5;
		idx += (permutationStride);
		xm[idx] = ys.i6;
		idx += (permutationStride);
		xm[idx] = ys.i7;
	}
	idx = idx - (permutationStride * 7) - 1;
	//#################################################
	//	fft-base2 round 4
	c = 0;
	c2 = 0;
	num1 = 0;
	num3 = 0;
	pos = 0;
	//	usfixn64 tmp = 0;
	//		offset = 0;

	//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < 8; i++)
	{
		//			num1 = xm[offset] + ym[offset] + c;
		//			num3 = ym[offset] + c2;
		num1 = xm[idx] + xm[idx + 1] + c;
		num3 = xm[idx + 1] + c2;

		//addition part
		if (num1 < xm[idx] || num1 < xm[idx + 1])	//there is overflow/truncation
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
		if (xm[idx] < num3) //there is not enough to do subtraction
		{
			c2 = 1;
			num3 = R - num3 + xm[idx];
		}
		else
		{
			c2 = 0;
			num3 = xm[idx] - num3;
		}
		xm[idx] = num1;
		xm[idx + 1] = num3;
		//			offset += permutationStride;
		idx += permutationStride;
	}
	idx -= (8 * permutationStride);

	//		offset = 0;
	if (c > 0)
	{
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] != 0)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}
		if (pos >= 0)	// shouldn't it be >0?
		{
			//				offset = 0;

			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = R_MINUS_ONE;
				idx += permutationStride;
			}
			xm[idx]--;
			idx -= pos * permutationStride;
			//				idx-=pos1 * permutationStride;
			//				offset = pos1 * permutationStride;
			//xm[pos1*permutationStride+idx]--;
		}
		else
		{
			idx -= permutationStride * 8;
			//			xm[0] = ULMAX;
			//				xm[offset] = ULMAX;
			xm[idx] = ULMAX;
			idx += permutationStride;
			//				offset = permutationStride;
			for (i = 1; i < 8; i++)
			{
				//				xm[i] = 0;
				//					xm[offset] = 0;
				xm[idx] = 0;
				idx += permutationStride;
			}
			idx -= (permutationStride * 8);
		}
	}

	idx += 1;
	if (c2 > 0)
	{
		//			offset = 0;
		pos = -1;
		for (i = 0; i < 8; i++)
		{
			if (xm[idx] < R_MINUS_ONE)
			{
				pos = i;
				break;
			}
			idx += permutationStride;
		}

		if (pos >= 0)
		{
			//				offset = 0;
			idx -= pos * permutationStride;
			for (i = 0; i < pos; i++)
			{
				xm[idx] = 0;
				idx += permutationStride;
			}
			xm[idx]++;
			idx -= pos * permutationStride;
			//				offset = pos * permutationStride;
			//			um[pos]++;
			//				ym[offset]++;
		}
		else
		{
			//				offset = 0;
			//			um[0] = ULMAX;
			//				ym[offset] = ULMAX;
			//				offset = permutationStride;
//#pragma unroll (COEFFICIENT_SIZE-1)
			for (i = 1; i < 8; i++)
			{
				//				um[i] = 0;
				xm[idx] = 0;
				idx -= permutationStride;
			}
			xm[idx] = ULMAX;
		}
	}
	//#################################################
	i = (tid & (0x07));
	idx = (tid << 1);
	//#################################
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx + 1], permutationStride);
	__syncthreads ();
	/********************************************************************/

	shiftVector8 = (0xF9F7FBF900FE0200);
//	shiftVector8 = (0xF900000000000000);
	shiftValue = (0xFF);	//should be repeated because of 64-bit registers
//	shiftValue <<= (shNo << 3);
	shiftValue <<= (shNo);
//	shNoShort = (shiftVector8 & shiftValue)>>(shNo<< 3);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here
	//>>(shNo<<3);

//	shNo_FFT_16_r2_shmem.i[0] = 0;			//{ 0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[1] = 2;			//{ 0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[2] = -2;			//{ 0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[3] = 0;			//{ 0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[4] = -7;			//{-7,-5,-9,-7}
//	shNo_FFT_16_r2_shmem.i[5] = -5;			//{-7,-5,-9,-7}
//	shNo_FFT_16_r2_shmem.i[6] = -9;			//{-7,-5,-9,-7}
//	shNo_FFT_16_r2_shmem.i[7] = -7;			//{-7,-5,-9,-7}
//	shNo_FFT_16_r2_shmem.i[8] = 7;			//{7,9,5,7}
//	shNo_FFT_16_r2_shmem.i[9] = 9;			//{7,9,5,7}
//	shNo_FFT_16_r2_shmem.i[10] = 5;			//{7,9,5,7}
//	shNo_FFT_16_r2_shmem.i[11] = 7;			//{7,9,5,7}
//	shNo_FFT_16_r2_shmem.i[12] = 0;			//{0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[13] = 2;			//{0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[14] = -2;			//{0,2,-2,0}
//	shNo_FFT_16_r2_shmem.i[15] = 0;			//{0,2,-2,0}

//	shNo = (tid & 0x7);
//	i = shNo + 8;

	shNo >>= 3;
	i = shNo;
//	shNo = shNo_FFT_16_r2_shmem.i[shNo];
	shNo = shNoShort;
//	idx2 = idx + shNo;
//	if (tid < 8)
//	{
//		printf(" shNo = %d \n",shNo);
//	}
	//	offset=(idx/16)*16;
	//	device_swap(&xm[idx], &xm[idx + shNo], permutationStride, coeffSize);
//	shNo
	usfixn64 tmp = 0;
//	short j = 0;
	j = 0;
	if (shNo > 0)
	{
//		device_p3_swap_permutated_2(&xm[idx], &xm[idx + shNo], &permutationStride);
//		device_p3_swap_permutated_2(&xm[idx], &xm[idx2], permutationStride);
		for (j = 0; j < 8; j++)
		{
			tmp = xm[idx];
			xm[idx] = xm[idx + shNo];
			xm[idx + shNo] = tmp;
			idx += permutationStride;
//			idx2+=permutationStride;
		}
		idx = (tid << 1);
//		idx2=idx+shNo;
	}
	__syncthreads ();
	//	shNo = idx % 16 + 1;
	//	shNo = shNo_FFT_16_r5[shNo];

	shNo = i;
	shiftVector8 = (0x00FE020007050907);
//	shiftVector8 = (0xF900000000000000);
//	shiftValue = (0xFF);	//should be repeated because of 64-bit registers
//	shiftValue <<= (shNo << 3);
//	shNoShort = (shiftVector8 & shiftValue)>>(shNo<< 3);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here
	idx++;
	shNo = shNoShort;
//	idx2 = idx + shNo;

//	shNo_FFT_16_r2_shmem.i[0] = 7;			//{7,9,5,7}
//		shNo_FFT_16_r2_shmem.i[1] = 9;			//{7,9,5,7}
//		shNo_FFT_16_r2_shmem.i[2] = 5;			//{7,9,5,7}
//		shNo_FFT_16_r2_shmem.i[3] = 7;			//{7,9,5,7}
//		shNo_FFT_16_r2_shmem.i[4] = 0;			//{0,2,-2,0}
//		shNo_FFT_16_r2_shmem.i[5] = 2;			//{0,2,-2,0}
//		shNo_FFT_16_r2_shmem.i[6] = -2;			//{0,2,-2,0}
//		shNo_FFT_16_r2_shmem.i[7] = 0;			//{0,2,-2,0}
	//	shNo = idx % 16+1;
//		shNo = (tid & 0x7) + 8;
//	shNo = shNo_FFT_16_r2_shmem.i[i];
//	idx++;
//	idx2 = idx + 1;
	if (shNo > 0)
	{
//		device_p3_swap_permutated_2(&xm[idx], &xm[idx2], permutationStride);
		for (j = 0; j < 8; j++)
		{
			tmp = xm[idx];
			xm[idx] = xm[idx + shNo];
			xm[idx + shNo] = tmp;
			idx += permutationStride;
//					idx2+=permutationStride;
		}
	}

//		device_p3_swap_permutated_2(&xm[idx2], &xm[idx2 + shNo], &permutationStride);
//		device_p3_swap_permutated_2(&xm[idx], &xm[idx + shNo], &permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r1 (usfixn64 * xm, const usfixn64 & permutationStride)
{
	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
//	usfixn64 tid;
//	get_tid(tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
//	short shNo = i;

//	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;

//	idx = ((tid & 0xFFFFF8) << 1) + i;
	idx = (tid / 8) * 16 + i;
	idx2 = idx + 8;
//	idx =
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx+8], permutationStride);

	usfixn64 x[COEFFICIENT_SIZE], y[COEFFICIENT_SIZE], t[COEFFICIENT_SIZE];
	usfixn64 offset = idx;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		x[i] = xm[offset];
		t[i] = xm[offset];
		y[i] = xm[offset + 8];
		offset += permutationStride;
	}

	device_p3_bigPrimeAdd_plain (x, y);
	device_p3_bigSub_plain (t, y);

	offset = idx;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		xm[offset] = x[i];
		xm[offset + 8] = t[i];
		offset += permutationStride;
	}
//	printf("tid =%d, idx=%d \n ",tid,idx);

}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r1_v1 (usfixn64 * __restrict__ xm,
													 const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
//	short shNo = i;

//	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;

	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + 8], permutationStride);
//	printf("tid =%d, idx=%d \n ",tid,idx);

}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r21 (usfixn64 * xm, const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
//	short shNo = i;
	short shNo = (i << 3);
//	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

	i = i & 0x3;

	//round 2
	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	shNoShort = (0x0404040400000000ULL & shiftValue) >> (shNo);

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx+4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx+4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_9(&xm[idx + 4], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 4], shNoShort, permutationStride);
	device_p3_cyclicShift_permutated_11 (&xm[idx + 4], shNoShort,
																			 permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r21_shmem (usfixn64 * xm, usfixn64 * tsh,
															 const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
//	short shNo = i;
	short shNo = (i << 3);
//	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

	i = i & 0x3;

	//round 2
	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	shNoShort = (0x0404040400000000ULL & shiftValue) >> (shNo);

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx+4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx+4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_9(&xm[idx + 4], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 4], shNoShort, permutationStride);

	device_p4_rotate_permutated_shmem_step1 (&xm[idx + 4], tsh, shNoShort,
																					 permutationStride);
	device_p4_rotate_permutated_shmem_step2 (&xm[idx + 4], shNoShort,
																					 permutationStride);

//		device_p3_fft_base2_permutated(&xm[idx], &xm[idx+4], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r22 (usfixn64 * xm, const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
//	short shNo = i;
	short shNo = (i << 3);
//	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

	i = i & 0x3;

	//round 2
	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	shNoShort = (0x0404040400000000ULL & shiftValue) >> (shNo);

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	fft_base2_permutated_8(&xm[idx], &xm[idx+4], shNoShort, permutationStride);
//	device_cyclicShift_permutated_7(&xm[idx+4], shNoShort, permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + 4], permutationStride);
//	fft_base2_permutated_v1(&xm[idx], &xm[idx + 4], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r31 (usfixn64 * xm, const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
//	shiftVector8 = (0x0606020204040000);
	shNoShort = (0x0606020204040000 & shiftValue) >> (shNo);

	//round 3
	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//		device_cyclicShift_permutated_7(&xm[idx+2], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 2], shNoShort, permutationStride);
//	device_cyclicShift_permutated_9(&xm[idx + 2], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 2], shNoShort, permutationStride);
	device_p3_cyclicShift_permutated_11 (&xm[idx + 2], shNoShort,
																			 permutationStride);

}
/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r31_shmem (usfixn64 * xm, usfixn64 *tsh,
															 const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
//	shiftVector8 = (0x0606020204040000);
	shNoShort = (0x0606020204040000 & shiftValue) >> (shNo);

	//round 3
	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//		device_cyclicShift_permutated_7(&xm[idx+2], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 2], shNoShort, permutationStride);
//	device_cyclicShift_permutated_9(&xm[idx + 2], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 2], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_11(&xm[idx + 2], shNoShort, permutationStride);
	device_p4_rotate_permutated_shmem_step1 (&xm[idx + 2], tsh, shNoShort,
																					 permutationStride);
	device_p4_rotate_permutated_shmem_step2 (&xm[idx + 2], shNoShort,
																					 permutationStride);

}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r32 (usfixn64 * xm, const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	shNoShort = (0x0606020204040000 & shiftValue) >> (shNo);

	//round 3
	idx = ((tid & 0xFFFFFE) << 1) + i;
	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx + 2], permutationStride);
//	fft_base2_permutated_v1(&xm[idx], &xm[idx + 2], permutationStride);

}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r41 (usfixn64 * xm, const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shiftVector8 = (0x0703050106020400);
	shNoShort = (0x0703050106020400 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	/********************************/
	/*the main problem with r41 and r41 is very gross strided access to global memory
	 * as profiling shows the gld throughput is 50
	 * solution is L(?,2) strided permutation, this should be applied to input before r41 and r42, and then, in doing the inverse
	 */
	/********************************/

	idx = (tid << 1);
//	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/
//	device_cyclicShift_permutated_7(&xm[idx+1], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 1], shNoShort, permutationStride);
//		device_cyclicShift_permutated_9(&xm[idx + 1], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 1], shNoShort, permutationStride);
	device_p3_cyclicShift_permutated_11 (&xm[idx + 1], shNoShort,
																			 permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r41_shmem (usfixn64 * xm, usfixn64 *tsh,
															 const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shiftVector8 = (0x0703050106020400);
	shNoShort = (0x0703050106020400 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	/********************************/
	/*the main problem with r41 and r41 is very gross strided access to global memory
	 * as profiling shows the gld throughput is 50
	 * solution is L(?,2) strided permutation, this should be applied to input before r41 and r42, and then, in doing the inverse
	 */
	/********************************/

	idx = (tid << 1);
//	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/
//	device_cyclicShift_permutated_7(&xm[idx+1], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 1], shNoShort, permutationStride);
//		device_cyclicShift_permutated_9(&xm[idx + 1], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx + 1], shNoShort, permutationStride);
	//	device_p3_cyclicShift_permutated_11(&xm[idx + 1], shNoShort, permutationStride);
	device_p4_rotate_permutated_shmem_step1 (&xm[idx + 1], tsh, shNoShort,
																					 permutationStride);
	device_p4_rotate_permutated_shmem_step2 (&xm[idx + 1], shNoShort,
																					 permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r42 (usfixn64 * xm, const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shiftVector8 = (0x0703050106020400);
	shNoShort = (0x0703050106020400 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	idx = (tid << 1);
	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx2], permutationStride);
//	fft_base2_plain_inPlace(&xm[idx], &xm[idx2], permutationStride);
//	device_p3_fft_base2_permutated(&xm[idx], &xm[idx2], permutationStride);
//	fft_base2_permutated_4(&xm[idx], &xm[idx2], permutationStride);
//	fft_base2_permutated_v1(&xm[idx], &xm[idx2], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r41_v1 (usfixn64 * __restrict__ xm,
														const usfixn64 & permutationStride)
{

//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
//	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shiftVector8 = (0x0703050106020400);
	shNoShort = (0x0703050106020400 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	/********************************/
	/*the main problem with r41 and r41 is very gross strided access to global memory
	 * as profiling shows the gld throughput is 50
	 * solution is L(?,2) strided permutation, this should be applied to input before r41 and r42, and then, in doing the inverse
	 */
	/********************************/

//	idx = (tid << 1);
//	idx2 = idx + 1;
//	idx = (tid/8)*16+(tid%8);
//	idx = ((tid>>3)<<4)+(tid&0x7);
////	idx = (tid&(0xFFFFFFFFFFFFFFF8))+(tid&0x7);
//	usfixn64 idx2 = idx+8; //+8
//	idx = ((tid>>3)<<4)+(tid&0x7)+8;
	//	idx = (tid&(0xFFFFFFFFFFFFFFF8))+(tid&0x7);
//		usfixn64 idx2 = idx+8; //+8
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/
//	device_cyclicShift_permutated_7(&xm[idx+1], shNoShort, permutationStride);
//	device_cyclicShift_permutated_8(&xm[idx + 1], shNoShort, permutationStride);
//		device_cyclicShift_permutated_9(&xm[idx], shNoShort, permutationStride);
//	device_p3_cyclicShift_permutated_5(&xm[idx], shNoShort, permutationStride);
	idx = ((tid >> 3) << 4) + (tid & 0x7);
	if ((tid & 0x1) == 0)
		device_p3_cyclicShift_permutated_11 (&xm[idx], shNoShort,
																				 permutationStride);

	shNoShort++;
	idx += 8;
	if ((tid & 0x1) == 0)
		device_p3_cyclicShift_permutated_11 (&xm[idx], shNoShort,
																				 permutationStride);
//		device_cyclicShift_permutated_9(&xm[idx2], shNoShort, permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r42_v1 (usfixn64 * __restrict__ xm,
														const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
//	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	shiftVector8 = (0x0703050106020400);
	shNoShort = (0x0703050106020400 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

//	idx = (tid << 1);
//	idx2 = idx + 1;
//	idx = (tid/8)*16+(tid%8);
	idx = ((tid >> 3) << 4) + (tid & 0x7);
//	idx = (tid&(0xFFFFFFFFFFFFFFF8))+(tid&0x7);
	idx2 = idx + 8;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/
	device_p3_fft_base2_permutated (&xm[idx], &xm[idx2], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r5 (usfixn64 * xm, const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	idx = (tid << 1);
	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/

	shiftVector8 = (0xF9F7FBF900FE0200);
	shiftValue = (0xFF);	//should be repeated because of 64-bit registers
	shiftValue <<= (shNo);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here

	shNo >>= 3;
	i = shNo;
	shNo = shNoShort;
	idx2 = idx + shNo;
	if (shNo > 0)
		device_p3_swap_permutated_2 (&xm[idx], &xm[idx2], permutationStride);

	shNo = i;
	shiftVector8 = (0x00FE020007050907);

	shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here
	idx++;
	shNo = shNoShort;
	idx2 = idx + shNo;
	if (shNo > 0)
		device_p3_swap_permutated_2 (&xm[idx], &xm[idx2], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_base_fft16_7_r5_v1 (usfixn64 * xm, const usfixn64 & permutationStride)
{
//	usfixn64 tid = threadIdx.x + (blockIdx.x * blockDim.x);
	usfixn64 tid;
	get_tid (tid);
	usfixn64 shiftVector8;
//	shiftVector8 = (0x0404040400000000ULL);	//0,0,0,0,4,4,4,4
	short shNoShort = 0;
	usfixn64 idx = 0;
	usfixn64 idx2 = 0;
	short i = (tid & 0x7);
	short shNo = i;

	shNo = shNo << 3;
	//round 1
//	idx = ((tid & 0xFFFFF8) << 1) + i;
//	idx2 = idx + 8;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	i = i & 0x3;

	//round 2
//	idx = ((tid & 0xFFFFFC) << 1) + i;
//	idx2 = idx + 4;
	usfixn64 shiftValue = (0x07);
	shiftValue <<= (shNo);

//		shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//		fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
//	i = i & 0x1;

//	shiftVector8 = (0x0606020204040000);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

	//round 3
//	idx = ((tid & 0xFFFFFE) << 1) + i;
//	idx2 = idx + 2;

//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);

//	shiftVector8 = (0x0703050106020400);
//	shNoShort = (shiftVector8 & shiftValue) >> (shNo);

//	i = 0;

	//round 4

	idx = (tid << 1);
	idx2 = idx + 1;
//	fft_base2_permutated_8(&xm[idx], &xm[idx2], shNoShort, permutationStride);
	/********************************************************************/

	shiftVector8 = (0xF9F7FBF900FE0200);
	shiftValue = (0xFF);	//should be repeated because of 64-bit registers
	shiftValue <<= (shNo);
	shNoShort = (shiftVector8 & shiftValue) >> (shNo);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here

	shNo >>= 3;
	i = shNo;
	shNo = shNoShort;
	idx2 = idx + shNo;
	if (shNo > 0)
		device_p3_swap_permutated_2 (&xm[idx], &xm[idx2], permutationStride);

	shNo = i;
	shiftVector8 = (0x00FE020007050907);

	shNoShort = (shiftVector8 & shiftValue) >> (shNo << 3);
	if (shNoShort > 127)
		shNoShort -= (0x100); //needs work here
	idx++;
	shNo = shNoShort;
	idx2 = idx + shNo;
	if (shNo > 0)
		device_p3_swap_permutated_2 (&xm[idx], &xm[idx2], permutationStride);
}

/**********************************************/
__device__ __inline__ void
device_strided_permutation_single (usfixn64 *xm, usfixn64 *xm_shared,
																	 const usfixn64 stride, const usfixn64 offset)
{
	//threadIdx.x  and tid = offset
	short i = 0;
	short row = threadIdx.x >> 4; //row No.
	short column = threadIdx.x - (row << 4); //col No.

}

/**********************************************/
/**********************************************/
//only works if N/2 threads are being used for N coefficients
//this is compatible with fft_base_256
__device__ __inline__ void
device_strided_permutation_multiple (usfixn64 *xm, usfixn64 *xm_shared,
																		 const usfixn64 stride,
																		 const usfixn64 permutationStride,
																		 const short coefficientSize,
																		 const usfixn64 tid)
{
	short i = 0;
	short r = threadIdx.x / stride;
	short c = threadIdx.x - (r * stride);
	short idx = (c * stride) + r; //(r,c) to (c,r)

	usfixn64 offset = tid;
//#pragma unroll COEFFICIENT_SIZE
	//number of rounds = inputVectorSize / number of blocks
	// 2 rounds = input vector size of N / N/2 blocks
	for (i = 0; i < coefficientSize; i++)
	{
		xm_shared[threadIdx.x] = xm[offset];
		__syncthreads ();
		xm[offset] = xm_shared[idx];
		offset += permutationStride / 2;

		xm_shared[threadIdx.x] = xm[offset];
		__syncthreads ();
		xm[offset] = xm_shared[idx];
		offset += permutationStride / 2;
	}
}

/**********************************************/

#endif
