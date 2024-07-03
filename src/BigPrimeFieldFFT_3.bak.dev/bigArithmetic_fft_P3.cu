#ifndef BIG_ARITHMETIC_FFT_H_
#define BIG_ARITHMETIC_FFT_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_addition_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_multiplication_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_cyclicShift_P3.h"

//extern __shared__ usfixn64 sharedMem[];

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
__device__ void
device_p4_rotate_permutated_shmem_step1 (usfixn64* __restrict__ xs,
																				 usfixn64 * __restrict__ tsh, short sn,
																				 usfixn64 permutationStride)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	usfixn64 permutationStride = parameters[5];
	tid=0;
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
__device__ void
device_p4_rotate_permutated_shmem_step2 (usfixn64 * __restrict__ xs, short sn,
																				 const usfixn64 permutationStride)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	tid=0;
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void inline
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
__device__ void __inline__
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
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
__device__ void
device_strided_permutation_single (usfixn64 *xm, usfixn64 *xm_shared,
																	 const usfixn64 stride, const usfixn64 offset)
{
	//threadIdx.x  and tid = offset
	short i = 0;
	short row = threadIdx.x >> 4; //row No.
	short column = threadIdx.x - (row << 4); //col No.

}

/**********************************************/
//only works if N/2 threads are being used for N coefficients
//this is compatible with fft_base_256
__device__ void
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
