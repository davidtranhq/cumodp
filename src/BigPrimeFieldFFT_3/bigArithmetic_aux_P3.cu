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


#ifndef BIG_ARITHMETIC_AUX_H_P3
#define BIG_ARITHMETIC_AUX_H_P3

#include "../../include/BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "bigArithmetic_aux_P3_device.h"
/****************************************************/
extern __shared__ usfixn64 sharedMem[];


/**********************************************/
//
//__global__ void
//kernel_p3_permutation_16_plain (usfixn64 * xs, usfixn64* parameters);
///**********************************************/
////
//__global__ void
//kernel_p3_permutation_64_plain (usfixn64 * xs, usfixn64* parameters);
//
///**********************************************/
////
//__global__ void
//kernel_p3_permutation_256_plain (usfixn64 * xs, usfixn64 * ys,
//																 usfixn64* parameters, int lp, int mp);
//
///**********************************************/
//__global__ void
//kernel_p3_permutation_16_permutated (usfixn64 * __restrict__ xs,
//																		 usfixn64* parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_permutation_16_permutated_v1 (
//		usfixn64 * __restrict__ xs, const usfixn64* __restrict__ parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_permutation_16_permutated_v2 (
//		usfixn64 * __restrict__ xs, const usfixn64* __restrict__ parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_twiddle_16_vectorized_permutated_r1 (usfixn64 * xs,
//																							 usfixn64* parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_twiddle_16_general_vectorized_permutated_r2_join1 (
//		usfixn64 * __restrict__ xs, usfixn32 height, usfixn64* __restrict__ lVector,
//		usfixn64* __restrict__ hVector, usfixn64* __restrict__ cVector,
//		usfixn64* __restrict__ parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_twiddle_16_general_vectorized_permutated_r2_join2 (
//		usfixn64 * __restrict__ xs, usfixn32 height,
//		usfixn64* __restrict__ lVectorSub, usfixn64* __restrict__ hVectorSub,
//		usfixn64* __restrict__ cVectorSub, usfixn64* __restrict__ parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_twiddle_16_general_vectorized_permutated_r2_v1 (
//		usfixn64 * xs, usfixn64 * __restrict__ pow_omega, usfixn32 height,
//		usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
//		usfixn64* lVectorSub, usfixn64* hVectorSub, usfixn64* cVectorSub,
//		usfixn64* __restrict__ parameters);
//
///**********************************************/
__global__ void
kernel_p3_permutationTwiddle_16_plain (usfixn64 * xs, usfixn64* parameters);

/****************************************************/
__global__ void
kernel_p3_permutation_16_plain (usfixn64 * xs, usfixn64* parameters)
{
//	usfixn64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	short stride = 16;
	usfixn64 tmp[256 * COEFFICIENT_SIZE];			//stride^2
	short i, j, k;
	usfixn64 offset = 0;
	usfixn64 n = parameters[7];
	usfixn64 nBlocks = (n) / (stride * stride);
	short preIdx, postIdx;
	short s = 0;
	for (k = 0; k < nBlocks; k++)
	{
		offset = k * stride * stride * COEFFICIENT_SIZE;
		for (i = 0; i < stride; i++)
		{
			for (j = 0; j < stride; j++)
			{
				for (s = 0; s < COEFFICIENT_SIZE; s++)
					tmp[(i * stride + j) * COEFFICIENT_SIZE + s] = xs[offset
							+ ((i * stride) + (j)) * COEFFICIENT_SIZE + s];
			}
		}

		for (i = 0; i < stride; i++)
		{
			for (j = 0; j < stride; j++)
			{
				for (s = 0; s < COEFFICIENT_SIZE; s++)
					xs[offset + ((j * stride) + (i)) * COEFFICIENT_SIZE + s] = tmp[(j
							+ i * stride) * COEFFICIENT_SIZE + s];
//					xs[offset+((i*stride)+(j))*COEFFICIENT_SIZE+s] = (i+j*stride)*COEFFICIENT_SIZE+s;
//				xs[offset+(i*stride)+(j)] = (i);
			}
		}
	}
}

/**********************************************/
//not working
__global__ void
kernel_p3_permutation_64_plain (usfixn64 * xs, usfixn64* parameters)
{
//	usfixn64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	short stride = 16;
	usfixn64 tmp[256 * COEFFICIENT_SIZE];			//stride^2
	short i, j, k;
	usfixn64 offset = 0;
	usfixn64 n = parameters[7];
	usfixn64 nBlocks = (n) / (stride * stride);
	short preIdx, postIdx;
	short s = 0;
	for (k = 0; k < nBlocks; k++)
	{
		offset = k * stride * stride * COEFFICIENT_SIZE;
		for (i = 0; i < stride; i++)
		{
			for (j = 0; j < stride; j++)
			{
				for (s = 0; s < COEFFICIENT_SIZE; s++)
					tmp[(i * stride + j) * COEFFICIENT_SIZE + s] = xs[offset
							+ ((i * stride) + (j)) * COEFFICIENT_SIZE + s];
			}
		}

		for (i = 0; i < stride; i++)
		{
			for (j = 0; j < stride; j++)
			{
				for (s = 0; s < COEFFICIENT_SIZE; s++)
					xs[offset + ((i * stride) + (j)) * COEFFICIENT_SIZE + s] = tmp[(i
							+ j * stride) * COEFFICIENT_SIZE + s];
//					xs[offset+((i*stride)+(j))*COEFFICIENT_SIZE+s] = (i+j*stride)*COEFFICIENT_SIZE+s;
//				xs[offset+(i*stride)+(j)] = (i);
			}
		}
	}
}
/**********************************************/
//not working
__global__ void
kernel_p3_permutation_256_plain (usfixn64 * xs, usfixn64 * ys,
															usfixn64* parameters, int lp, int mp)
{
//	usfixn64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	short stride = parameters[2];
	short m = parameters[1];
	stride = lp;
	m = mp;
	usfixn64 permutationStride = parameters[5];
//	usfixn64 tmp[256 * COEFFICIENT_SIZE];			//stride^2
	usfixn64 i, j, k;
	usfixn64 r;
	usfixn64 offset = 0;
	usfixn64 n = parameters[7];
	usfixn64 nBlocks = (n) / (stride * m);
	short c = 0;
	usfixn64 idx = 0, outIdx = 0;
//	stride=16'
	{
		for (r = 0; r < nBlocks; r++)
		{
			offset = r * (stride * m) * COEFFICIENT_SIZE;
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < stride; j++)
				{
					idx = (i * stride + j) * COEFFICIENT_SIZE;
					outIdx = (i + j * m) * COEFFICIENT_SIZE;
					for (k = 0; k < COEFFICIENT_SIZE; k++)
						ys[offset + outIdx + k] = xs[offset + idx + k];
				}
			}
		}
	}
}
/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated (usfixn64 * __restrict__ xs,
																	usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 offset;
//	offset=tid;
//	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
	usfixn64 * shmem = sharedMem;
	short baseStride = 16;
//	short shIdx = threadIdx.x%baseStride;
	short shIdx = threadIdx.x & (0xF);
//	short shmemOffset=(threadIdx.x/baseStride)*(baseStride*baseStride);
	short shmemOffset = ((threadIdx.x >> 4) << 8);
	short i = 0;
	short c = 0;
//	offset = (tid/16)*256;
	offset = ((tid >> 4) << 8);
//	offset = ((blockIdx.x)<<8);
//	offset = (((threadIdx.x + blockIdx.x * blockDim.x)>>4)<<8);
//#pragma unroll COEFFICIENT_SIZE
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
		for (i = 0; i < baseStride; i++)
		{
//		shmem[shmemOffset+i*baseStride+shIdx]=xs[offset+shIdx+i*baseStride];
			shmem[shmemOffset + (i << 4) + shIdx] = xs[offset + shIdx + (i << 4)];
//		__syncthreads();
		}

		for (i = 0; i < baseStride; i++)
		{
//		xs[offset+shIdx+i*baseStride]=shmem[shmemOffset+i+shIdx*baseStride];
			xs[offset + shIdx + (i << 4)] = shmem[shmemOffset + i + (shIdx << 4)];
//		__syncthreads();
		}
		offset += permutationStride;
	}
}

/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated_v1 (usfixn64 * __restrict__ xs,
																		 const usfixn64* __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 offset;

//	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	usfixn64 * shmem = sharedMem;
	__shared__    __volatile__ usfixn64 shmem[16][33];
	usfixn64 n = parameters[7];

	short i, j, c;
	i = (tid >> 4);
	j = tid & 0xF;

	for (tid; tid < n; tid += blockDim.x * gridDim.x)
	{
		offset = tid;
//#pragma unroll COEFFICIENT_SIZE
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
			shmem[i][j] = xs[offset];
			__syncthreads ();
			xs[offset] = shmem[j][i];
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated_v2 (usfixn64 * __restrict__ xs,
																		 const usfixn64* __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 offset;

//	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	usfixn64 * shmem = sharedMem;
	__shared__    __volatile__ usfixn64 shmem[16][33];
	usfixn64 n = parameters[7];

	short i, j, c;
	i = (tid >> 4);
	j = tid & 0xF;

	for (tid; tid < n; tid += blockDim.x * gridDim.x)
	{
		offset = tid;
//#pragma unroll COEFFICIENT_SIZE
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
			shmem[i][j] = xs[offset];
			__syncthreads ();
			xs[offset] = shmem[j][i];
			offset += permutationStride;
		}
	}
}

/**********************************************/
//__global__ void
//kernel_p3_twiddle_16_general_vectorized_permutated_r2_join1 (
//		usfixn64 * __restrict__ xs, usfixn32 height, usfixn64* __restrict__ lVector,
//		usfixn64* __restrict__ hVector, usfixn64* __restrict__ cVector,
//		usfixn64* __restrict__ parameters)
//{
//	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
////	usfixn64 permutationStride = parameters[5];
////	usfixn64 offset = tid;
////	uConstArray8_align8 tmp;
////	short i = 0;
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
////		tmp.i[i] = xs[offset + i * permutationStride];
////	}
////	short h = (tid % 256) / 16;
////	short w = tid % 16;
//
////	short h = (tid &0xFF) / 16;
////	short h = ((tid & 0xFF) >> 4);
////	short w = (tid & 0xF);
//
////	tid = tid % (height * 16);
////	short h = tid / 16;
////	short w = tid % 16;
//
//	tid = tid % (height << 4);
//	short h = tid >> 4;
//	short w = tid & (0xF);
//
////	mulOmegaPower(tmp.i, h * w);
////		mulOmegaPower_vectorzied(tmp.i, h );
////		mulOmegaPower_vectorized(tmp.i, h );
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
////		xs[offset + i * permutationStride] = tmp.i[i];
////	}
////	__shared__ usfixn64 shmemOmega[128];
////	usfixn64 *shmemOmega;
////	short j=0;
////	if(threadIdx.x==0)
////	{
////
////			for(j=0;j<8;j++)
////				for(i=1;i<16;i++)
////		{
//////			shmemOmega[(8*i)+(j)]=W_0[i-1][j];
////				shmemOmega[(i)+(16*j)]=W_0[i-1][j];
////		}
////	}
////	__syncthreads();
//	h = h * w;
//
////	mulOmegaPower_vectorized_permutated_r2(&xs[tid], h * w, permutationStride);
////	mulOmegaPower_vectorized_permutated_r2_v1(xs, lVector, hVector, cVector, lVectorSub, hVectorSub, cVectorSub, h , permutationStride);
//	device_newMult15_join_r1_mulOmega (xs, lVector, hVector, cVector, parameters,
//																		 h);
////	   device_newMult15_join_r2_mulOmega(xs, lVectorSub, hVectorSub, cVectorSub, parameters,h);
//
//}

/**********************************************/
//__global__ void
//kernel_p3_twiddle_16_general_vectorized_permutated_r2_join2 (
//		usfixn64 * __restrict__ xs, usfixn32 height,
//		usfixn64* __restrict__ lVectorSub, usfixn64* __restrict__ hVectorSub,
//		usfixn64* __restrict__ cVectorSub, usfixn64* __restrict__ parameters)
//{
//	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
////	usfixn64 permutationStride = parameters[5];
////	usfixn64 offset = tid;
////	uConstArray8_align8 tmp;
////	short i = 0;
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
////		tmp.i[i] = xs[offset + i * permutationStride];
////	}
////	short h = (tid % 256) / 16;
////	short w = tid % 16;
//
////	short h = (tid &0xFF) / 16;
////	short h = ((tid & 0xFF) >> 4);
////	short w = (tid & 0xF);
//
////	tid = tid % (height * 16);
////	short h = tid / 16;
////	short w = tid % 16;
//
//	tid = tid % (height << 4);
//	short h = tid >> 4;
//	short w = tid & (0xF);
//
////	mulOmegaPower(tmp.i, h * w);
////		mulOmegaPower_vectorzied(tmp.i, h );
////		mulOmegaPower_vectorized(tmp.i, h );
////	for (i = 0; i < COEFFICIENT_SIZE; i++)
////	{
////		xs[offset + i * permutationStride] = tmp.i[i];
////	}
////	__shared__ usfixn64 shmemOmega[128];
////	usfixn64 *shmemOmega;
////	short j=0;
////	if(threadIdx.x==0)
////	{
////
////			for(j=0;j<8;j++)
////				for(i=1;i<16;i++)
////		{
//////			shmemOmega[(8*i)+(j)]=W_0[i-1][j];
////				shmemOmega[(i)+(16*j)]=W_0[i-1][j];
////		}
////	}
////	__syncthreads();
//	h = h * w;
//
////	mulOmegaPower_vectorized_permutated_r2(&xs[tid], h * w, permutationStride);
////	mulOmegaPower_vectorized_permutated_r2_v1(xs, lVector, hVector, cVector, lVectorSub, hVectorSub, cVectorSub, h , permutationStride);
////	device_newMult15_join_r1_mulOmega(xs, lVector, hVector, cVector, parameters,h);
//	device_newMult15_join_r2_mulOmega (xs, lVectorSub, hVectorSub, cVectorSub,
//																		 parameters, h);
//
//}
//
/**********************************************/
//__global__ void
//kernel_p3_permutationTwiddle_16_plain (usfixn64 * xs, usfixn64* parameters)
//{
//
////	usfixn64 tid = threadIdx.x + blockIdx.x*blockDim.x;
//	short stride = 16;
//	usfixn64 tmp[256 * COEFFICIENT_SIZE];			//stride^2
//	short i, j, k;
//	usfixn64 offset = 0;
//	usfixn64 n = parameters[7];
//	usfixn64 nBlocks = (n) / (stride * stride);
//	short idx;
//	short idx2;
//	short s = 0;
//	for (k = 0; k < nBlocks; k++)
//	{
//		offset = k * stride * stride * COEFFICIENT_SIZE;
//		for (i = 0; i < stride; i++)
//		{
//			for (j = 0; j < stride; j++)
//			{
//				for (s = 0; s < COEFFICIENT_SIZE; s++)
//					tmp[(i * stride + j) * COEFFICIENT_SIZE + s] = xs[offset
//							+ ((i * stride) + (j)) * COEFFICIENT_SIZE + s];
//				idx = i * j;
//				mulOmegaPower (&tmp[(i * stride + j) * COEFFICIENT_SIZE], idx);
//			}
//		}
//
//		for (i = 0; i < stride; i++)
//		{
//			for (j = 0; j < stride; j++)
//			{
//				for (s = 0; s < COEFFICIENT_SIZE; s++)
//					xs[offset + ((i + j * stride) * COEFFICIENT_SIZE) + s] = tmp[((i
//							* stride + j) * COEFFICIENT_SIZE) + s];
//			}
//		}
//	}
//}
//
/**********************************************/
__global__ void
kernel_p3_permutation_256_permutated (usfixn64 * xs, usfixn64 * ys,
																	 usfixn64* parameters, int lp, int mp)
{
//	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	usfixn64 permutationStride = parameters[5];
//	usfixn64 offset;
////	offset=tid;
////	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
////	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	usfixn64 * shmem = sharedMem;
//	short baseStride = 256;
//	baseStride = parameters[2];
//	short m = parameters[1]; //L(k,m)=L(256,16) over matrices of size (256*16)=4096
////m=16;
////	short shIdx = threadIdx.x%baseStride;
////	short shIdx = threadIdx.x;// & (0xF);
////		short shIdx = (tid%256);// & (0xF);
//	short shIdx = threadIdx.x;
////	usfixn64 gmemIdx = (tid / baseStride) * (baseStride * baseStride) + (tid % baseStride);
//
//	usfixn64 gmemIdx = (tid / baseStride) * (m * baseStride) + (threadIdx.x);
//	gmemIdx = threadIdx.x;
//	short j = threadIdx.x;
////	short shmemOffset=(threadIdx.x/baseStride)*(baseStride*baseStride);
////	short shmemOffset = ((threadIdx.x >> 4) << 8);
//
//	short i = 0;
//	short c = 0;
////	offset = (tid/16)*256;
//
////	offset = ((blockIdx.x)<<8);
////	offset = (((threadIdx.x + blockIdx.x * blockDim.x)>>4)<<8);
//
//	short r = 0;
////	short shOffset = 0;
//	short dimX = blockDim.x;
//	short dimY = (32 * 32) / dimX;
//	dimY = dimX;
////	short nRepeats = (baseStride) / ((32 * 32) / blockDim.x));
//	short nRepeats = (m + (dimY - 1)) / (dimY);
////#pragma unroll COEFFICIENT_SIZE
////	if(nRepeats ==1 )
////		dimY =m;
//
//	short blockWidth = (baseStride / dimX);
//	short b = (blockIdx.x) % blockWidth; //8 = baseStride/b
////	b=blockIdx.x;
////	b= (blockIdx.x/8)*8+blockIdx.x%8;
//
//	if (m >= 32)
//		for (c = 0; c < COEFFICIENT_SIZE; c++)
//		{
////		offset = c * permutationStride;
//			for (r = 0; r < nRepeats; r++)
//			{
//				offset = c * permutationStride + r * baseStride * dimY + b * dimX
//						+ baseStride * m * (blockIdx.x / (blockWidth));
////			for (i = 0; i <(32*32)/blockDim.x; i++)
//				for (i = 0; i < dimY; i++)
//				{
//					shmem[i + shIdx * dimY] = xs[offset + i * baseStride + gmemIdx];
////				xs[offset + i * baseStride + j] = r;//baseStride*m*(blockIdx.x/(blockWidth))*(blockWidth);
////				shmem[i + shIdx * dimY] = offset + i * baseStride + gmemIdx;
////				shmem[i + shIdx * dimY] = offset + gmemIdx;
////				shmem[i + shIdx * dimY] = i*dimX+ shIdx ;
//				}
////		}
//				offset = c * permutationStride + b * m * dimY + r * dimX
//						+ baseStride * m * (blockIdx.x / (blockWidth));
////			offset = c * permutationStride + r*dimY*m;
////		for (r = 0; r < nRepeats; r++)
////		{
////			offset = c * permutationStride + r*baseStride*dimY;
////			for (i = 0; i <(32*32)/blockDim.x; i++)
//				for (i = 0; i < dimY; i++)
//				{
////				ys[offset + i * dimY + gmemIdx] = baseStride*m*(blockIdx.x/(blockWidth))*(blockWidth);//shmem[i*dimX +shIdx];
//					ys[offset + i * m + gmemIdx] = shmem[i * dimX + shIdx];
////				ys[offset + i * baseStride + gmemIdx] = shmem[i+ dimY*shIdx]/8;
//				}
//			}
//		}
//
//
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 offset;
//	offset=tid;
//	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
	usfixn64 * shmem = sharedMem;
	short baseStride = 256;
	baseStride = parameters[2];
	short m = parameters[1]; //L(k,m)=L(256,16) over matrices of size (256*16)=4096
	m = mp;
	baseStride = lp;
//m=16;
//	short shIdx = threadIdx.x%baseStride;
//	short shIdx = threadIdx.x;// & (0xF);
//		short shIdx = (tid%256);// & (0xF);
	short shIdx = threadIdx.x;
//	usfixn64 gmemIdx = (tid / baseStride) * (baseStride * baseStride) + (tid % baseStride);

	usfixn64 gmemIdx = (tid / baseStride) * (m * baseStride) + (threadIdx.x);
	gmemIdx = threadIdx.x;
	short j = threadIdx.x;
//	short shmemOffset=(threadIdx.x/baseStride)*(baseStride*baseStride);
//	short shmemOffset = ((threadIdx.x >> 4) << 8);

	short i = 0;
	short c = 0;
//	offset = (tid/16)*256;

//	offset = ((blockIdx.x)<<8);
//	offset = (((threadIdx.x + blockIdx.x * blockDim.x)>>4)<<8);

	short r = 0;
//	short shOffset = 0;
	short dimX = blockDim.x;
	short dimY = (32 * 32) / dimX;
	dimY = dimX;
//	short nRepeats = (baseStride) / ((32 * 32) / blockDim.x));
	short nRepeats = (m + (dimY - 1)) / (dimY);
//#pragma unroll COEFFICIENT_SIZE
//	if(nRepeats ==1 )
//		dimY =m;

	short blockWidth = (baseStride / dimX);
	short b;
	b = (blockIdx.x) % blockWidth; //8 = baseStride/b
//	b=blockIdx.x;
//	b= (blockIdx.x/8)*8+blockIdx.x%8;

	usfixn64 offsetIn = c * permutationStride + r * baseStride * dimY + b * dimX
			+ baseStride * m * (blockIdx.x / (blockWidth));
	offsetIn += dimX;	// at end
	usfixn64 offsetOut = c * permutationStride + b * m * dimY
			+ baseStride * m * (blockIdx.x / (blockWidth));
	offsetOut += dimX;	// at end

	if (m >= 32)
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
			//		offset = c * permutationStride;
			for (r = 0; r < nRepeats; r++)
			{
				offset = c * permutationStride + r * baseStride * dimY + b * dimX
						+ baseStride * m * (blockIdx.x / (blockWidth));
				//			for (i = 0; i <(32*32)/blockDim.x; i++)
				for (i = 0; i < dimY; i++)
				{
					shmem[i + shIdx * dimY] = xs[offset + i * baseStride + gmemIdx];
					//				xs[offset + i * baseStride + j] = r;//baseStride*m*(blockIdx.x/(blockWidth))*(blockWidth);
					//				shmem[i + shIdx * dimY] = offset + i * baseStride + gmemIdx;
					//				shmem[i + shIdx * dimY] = offset + gmemIdx;
					//				shmem[i + shIdx * dimY] = i*dimX+ shIdx ;
				}
				//		}
				offset = c * permutationStride + b * m * dimY + r * dimX
						+ baseStride * m * (blockIdx.x / (blockWidth));
				//			offset = c * permutationStride + r*dimY*m;
				//		for (r = 0; r < nRepeats; r++)
				//		{
				//			offset = c * permutationStride + r*baseStride*dimY;
				//			for (i = 0; i <(32*32)/blockDim.x; i++)
				for (i = 0; i < dimY; i++)
				{
					//				ys[offset + i * dimY + gmemIdx] = baseStride*m*(blockIdx.x/(blockWidth))*(blockWidth);//shmem[i*dimX +shIdx];
					ys[offset + i * m + gmemIdx] = shmem[i * dimX + shIdx];
					//				ys[offset + i * baseStride + gmemIdx] = shmem[i+ dimY*shIdx]/8;
				}
			}
		}
}

/**********************************************/
__global__ void
kernel_p3_permutation_256_permutated_2 (usfixn64 * xs, usfixn64 * ys,
																		 usfixn64* parameters)
{
	int lP = 32;
	int mP = 32;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 offset;
//	offset=tid;
//	__shared__ usfixn64 shmem[256];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
//	__shared__ usfixn64 shmem[512];	// (blockSize/16)*256, means block size is 16 here, for blockSize=32, then shmem[512]
	usfixn64 * shmem = sharedMem;
//	short baseStride = parameters[2];
//	short m = parameters[1]; //L(k,m)=L(256,16) over matrices of size (256*16)=4096
	int m = mP;
	int l = lP;
//m=16;
//	short shIdx = threadIdx.x%baseStride;
//	short shIdx = threadIdx.x;// & (0xF);
//		short shIdx = (tid%256);// & (0xF);
	short shIdx = threadIdx.x;
//	usfixn64 gmemIdx = (tid / baseStride) * (baseStride * baseStride) + (tid % baseStride);

//	usfixn64 gmemIdx = (tid / baseStride) * (m * baseStride) + (threadIdx.x);
//	gmemIdx = threadIdx.x;
	short j = threadIdx.x;
	short gmemIdx = threadIdx.x % l;
	short i = 0;
	short c = 0;
//	offset = (tid/16)*256;

	short r = 0;
//	short shOffset = 0;
	short dimX = blockDim.x;
	short dimY = (2 * 1024) / dimX;

//	dimY = dimX;
//	short nRepeats = (baseStride) / ((32 * 32) / blockDim.x));
	short nRepeats = (m + (dimY - 1)) / (dimY);
//#pragma unroll COEFFICIENT_SIZE
//	if(nRepeats ==1 )
//		dimY =m;

	short blockWidth = (dimX / l);
	short b;
	b = ((blockIdx.x) / blockWidth) * blockWidth + threadIdx.x / l; //8 = baseStride/b
//	b=blockIdx.x;
//	b= (blockIdx.x/8)*8+blockIdx.x%8;
	short bIdx = 0;
	usfixn64 offsetIn, offsetOut;
	if (b * l * m < parameters[0])
		if (m >= 32)
			for (c = 0; c < COEFFICIENT_SIZE; c++)
			{
				for (r = 0; r < nRepeats; r++)
				{
					offsetIn = c * permutationStride + b * l * m + r * dimY * l;
					//			for (i = 0; i <(32*32)/blockDim.x; i++)
					for (i = 0; i < dimY; i++)
					{
						shmem[i + j * dimY] = xs[offsetIn + i * l + gmemIdx];	// *blockWidth+ gmemIdx];
//					xs[offsetIn + i * l+gmemIdx] = blockIdx.x;//offsetIn;//i+r*dimY;
					}
					//		}
					offsetOut = c * permutationStride + b * l * m + r * dimY * l;
					//			offset = c * permutationStride + r*dimY*m;
					//		for (r = 0; r < nRepeats; r++)
					//		{
					//			offset = c * permutationStride + r*baseStride*dimY;
					//			for (i = 0; i <(32*32)/blockDim.x; i++)
					for (i = 0; i < dimY; i++)
					{
						//				ys[offset + i * dimY + gmemIdx] = baseStride*m*(blockIdx.x/(blockWidth))*(blockWidth);//shmem[i*dimX +shIdx];
						ys[offsetOut + i * l + gmemIdx] = shmem[i * dimX + j];
						//				ys[offset + i * baseStride + gmemIdx] = shmem[i+ dimY*shIdx]/8;
					}
				}
			}
}

/**********************************************/

//general 256_permutation, equivalent of permutation L(k^e, 256)
//ys = L(xs)
__global__ void
kernel_p3_permutation_256_general_permutated_v0 (usfixn64 n,
																							usfixn64 * __restrict__ xs,
																							usfixn64 * __restrict__ ys,
																							usfixn64 * parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	__shared__            __volatile__ usfixn64 shmem[256][17];
//	usfixn64 n = (1024); //height
	usfixn64 offsetDigit = 0;
	usfixn64 offsetBlock = 0;
	short r, i, c, l, j;
	l = 16;
	j = threadIdx.x;
	usfixn16 shmemoffset = 0;
	usfixn16 h, w, height;

//	h=j/16;
//	w=j%16;
//	height=256/l;
	h = j >> 4;
	w = j & 0xF;
	height = 16;
	usfixn64 offset;
	usfixn64 offsetGrid = 0;
//	usfixn64 a[16];

	offsetBlock = 0;
	c = 0;
//	usfixn64 gridStep = (n);
//	if(gridStep<=15)
//		gridStep=1;
//	else
//		gridStep = (1<<(n-15+1));

	n = (1 << n);
//	for (usfixn16 step=0;step<gridStep;step++)
	{
//		offsetGrid = step*(1<<15);
//		if(tid < offsetGrid)
//			break;
//#pragma unroll COEFFICIENT_SIZE
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
//		offsetDigit = c * permutationStride + blockIdx.x*256*n;
			offsetDigit = c * permutationStride + (blockIdx.x * n << 8) + offsetGrid;
//		for (r = 0; r < n / l; r++)
//		for (r = 0; r < n/l; r++)
			{
				r = blockIdx.y;
//			offsetBlock = r * l*256;
				offsetBlock = (r * l) << 8;

				offset = offsetBlock + offsetDigit + j;
//			for (i = 0; i < l; i++)
//#pragma unroll 16
				for (i = 0; i < 16; i++)
				{
//				shmem[j][i] = xs[offsetBlock + offsetDigit + i * 256 + j];
//				shmem[j][i] = i;//xs[offsetBlock + offsetDigit + (i<<8) + j];
					shmem[j][i] = xs[offset];
					offset += 256;
				}
				__syncthreads ();

				offsetBlock = r * l;
				for (i = 0; i < l; i++)
				{

//				shmemoffset  =(j/l)+i*(256/l);
//				shmemoffset = h+i*height;
					shmemoffset = h + (i << 4);
//				ys[offsetBlock + offsetDigit + (j/l)*n + i*(256/l)*n + (j%l)]=shmem[shmemoffset][j%l];//shmem[][j%l];
//				ys[offsetBlock + offsetDigit + h*n + i*(height)*n + (w)]=shmem[shmemoffset][w];//shmem[][j%l];
					ys[offsetBlock + offsetDigit + h * n + ((i * n) << 4) + (w)] =
							shmem[shmemoffset][w];					//shmem[][j%l];
				}
//			offsetBlock +=4096;
			}
		}
	}

}

/**********************************************/
__global__ void
kernel_p3_permutation_2_general_permutated_v0 (
		usfixn16 n, usfixn64 * __restrict__ xs, usfixn64 * ys,
		const usfixn64 * __restrict__ parameters)
{

	//doing permutation equivalent to L(2^n,2)

	__shared__ usfixn64 shmem[BLOCK_SIZE];
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn16 j = threadIdx.x;
	usfixn64 permutationStride = parameters[5];
	usfixn16 c = 0;
//	i = 0;

	usfixn64 offset = tid;

//	usfixn16 idx = (j/16)*8+(j%8);
//	usfixn16 idx = ((j >> 4) << 3) + (j & 0x7);
//	usfixn16 idx = ((j >> 1) &(0xFFF8)) + (j & 0x7);
//	usfixn16 idx = ((j >> 1) &(0x8)) + (j & 0x7);
//	i = j % 16;
//	i = j & (0xF);

//	offset=(tid/8)*16+(tid%8);

//	usfixn16 idx = (j/16)*16+(j%16)/2;
//	usfixn16 idx = (j & 0xFFF0) + ((j & 0xF) >> 1);
////	usfixn64 x;
//
//	if ((j & 0x1) == 1)
//	{
//		idx += 8;
////			shmem[j >> 1] = x;//xs[offset];
////			shmem[idx] = x;//xs[offset];
//	}

	usfixn16 idx = (j & 0xFFF0) + ((j & 0xF) >> 1) + ((j & 0x1) << 3);
	//	usfixn64 x;

//		if ((j & 0x1) == 1)
//		{
//			idx += 8;
//	//			shmem[j >> 1] = x;//xs[offset];
//	//			shmem[idx] = x;//xs[offset];
//		}

//		else// ((j & 0x1) == 1)
//		{
////			shmem[(j >> 1)+8] = x;//xs[offset];
//			shmem[idx+8] = x;//xs[offset];
//		}

	usfixn64 y[COEFFICIENT_SIZE];
//#pragma unroll COEFFICIENT_SIZE
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
//		if ((j % 2) == 0)
//		{
//			shmem[j/2] = xs[offset];
//		}
//		x=xs[offset];

//		if ((j & 0x1) == 0)
//		{
////			shmem[j >> 1] = x;//xs[offset];
//			shmem[idx] = x;//xs[offset];
//		}
//		shmem[idx]=x;
		shmem[idx] = xs[offset];
//
//		else// ((j & 0x1) == 1)
//		{
////			shmem[(j >> 1)+8] = x;//xs[offset];
//			shmem[idx+8] = x;//xs[offset];
//		}
//		__syncthreads();
//		if(j<blockDim.x/2)
//		if((j%16)<8)
//		{
////			printf("j=%d, tid=%lu,\n",j,tid);
//			ys[offset]=shmem[(j/16)*8+(j%8)];
////			ys[offset]=j;
//			//shmem[(j/16)*8+(j%8)];
//		}
//		if (i < 8)
		{
			//			printf("j=%d, tid=%lu,\n",j,tid);
//			ys[offset] = shmem[idx];

//			ys[offset] = shmem[j];

//			__syncthreads();
			y[c] = shmem[j];
			//			ys[offset]=j;
			//shmem[(j/16)*8+(j%8)];
		}
//		else
//		{
//			ys[offset]=shmem[idx+128];
//		}
		offset += permutationStride;
//		__syncthreads();
	}

	offset = tid;

//#pragma unroll COEFFICIENT_SIZE
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
		ys[offset] = y[c];
		offset += permutationStride;
	}

//	__syncthreads();
//
//	offset = tid;
//	for (c = 0; c < COEFFICIENT_SIZE; c++)
//	{
////			if (j % 2 == 1)
////			{
////				shmem[j/2] = xs[offset];
////			}
//		if ((j & 0x1) == 1)
//		{
//			shmem[j >> 1] = xs[offset];
//		}
////		__syncthreads();
//		//		if(j<blockDim.x/2)
////			if((j%16)>=8)
////			{
////	//			ys[offset]=shmem[(j/16)*8+(j%8)];
////				ys[offset]=shmem[(j/16)*8+(j%8)];
//////				ys[offset]=j;
//////				offset;
////			}
//		if (i >= 8)
//		{
//			//			ys[offset]=shmem[(j/16)*8+(j%8)];
//			ys[offset] = shmem[idx];
//			//				ys[offset]=j;
//			//				offset;
//		}
//		offset += permutationStride;
////		__syncthreads();
//	}
//	doing L(k,16) = 4 L(8k,2)
//	find a way for quick computation of L(8k,2)
}

/**********************************************/
__global__ void
kernel_p3_transposition (usfixn64* xs, usfixn64 * ys, usfixn64 nPermutationBlocks,
											usfixn64 blockStride, usfixn64* parameters)
{
	__shared__ usfixn64 shmem[16][16];

	int blockDim = 16;
	usfixn64 x = blockIdx.x * blockDim + threadIdx.x;
	usfixn64 y = blockIdx.y * blockDim + threadIdx.y;
	usfixn64 width = gridDim.x * blockDim;
	usfixn64 height = gridDim.y * blockDim;
//	if (x==0 && y==0)
//	{
//		printf("height=%lu , width = %lu \n", height, width);
//	}
	usfixn64 permutationStride = parameters[5];
	usfixn64 n = parameters[7];
	short c = 0, j = 0;
	usfixn64 blockOffset = 0;
	usfixn64 blockNo;
	usfixn64 offset = 0;
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
		blockOffset = 0;
		offset = c * permutationStride;
		for (blockNo = 0; blockNo < nPermutationBlocks; blockNo++)
		{
			shmem[threadIdx.y][threadIdx.x] = xs[(y) * width + x + offset
					+ blockOffset];
//			__syncthreads ();
//		x = blockIdx.y * blockDim + threadIdx.x;  // transpose block offset
//		y = blockIdx.x * blockDim + threadIdx.y;
			ys[(y) + x * height + offset + blockOffset] =
					shmem[threadIdx.y][threadIdx.x];
			blockOffset += blockStride;
		}

	}
}

/**********************************************/
__global__ void
kernel_p3_cpy_d2d (usfixn64* __restrict__ dest, const usfixn64 * __restrict__ src,
								const usfixn64 n, const usfixn64* __restrict__ parameters)
{
//	__shared__ usfixn64 shmem[32][32];

//	int blockDim = 32;
	usfixn64 offset = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	short c = 0, j = 0;
	usfixn64 gridStride = blockDim.x * gridDim.x;
	usfixn64 tidInit = tid;
#pragma unroll
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
//		offset=c*permutationStride;
		for (tid = tidInit; tid < n; tid += gridStride)
		{
			dest[tid + offset] = src[tid + offset];
		}
		offset += permutationStride;
	}
}

/**********************************************/
void
host_p3_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd)
{
	kernel_p3_cpy_d2d<<<32, 256>>> (dest, src, count,
			fd->device_parameters);
}

/**********************************************/
//compute L_{k,m} for any "input size", k, and m are multiples of 32 (TILE_DIMENSION)
//computes input/(k*m) block permutations
void
host_p3_transposition (data* fd, usfixn64 k, usfixn64 m)
{

	//checks
	if (k % TILE_DIMENSION != 0)
	{
		printf ("k in host_transposition is not a multiple of 32! return!;\n");
		return;
	}
	else if (m % TILE_DIMENSION != 0)
	{
		printf ("m in host_transposition is not a multiple of 32! return!;\n");
		return;
	}
	else if (fd->inputVectorSize % TILE_DIMENSION != 0)
	{
		printf (
				"inputVectorSize in host_transposition is not a multiple of 32! return!;\n");
		return;
	}

	dim3 blockDims (TILE_DIMENSION, TILE_DIMENSION, 1);
	dim3 gridDims (k / TILE_DIMENSION, m / TILE_DIMENSION, 1);

	/***************************************************/
	usfixn64 blockStride = k * m;
	usfixn64 nPermutationBlocks = fd->inputVectorSize / (blockStride);
//	printf ("npermutationBlocks = %lu, blockStride= %lu \n", nPermutationBlocks,
//					blockStride);
	/***************************************************/
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
//	printf ("k/32=%lu, m/32=%lu \n", gridDims.x, gridDims.y);
	kernel_p3_transposition<<<gridDims, blockDims>>> (fd->device_xVector,
			fd->device_yVector,
			nPermutationBlocks,
			blockStride,
			fd->device_parameters);
	host_p3_cpy_d2d (fd->device_xVector, fd->device_yVector, fd->inputVectorSize,
								fd);
//	cudaMemcpy (fd->device_xVector, fd->device_yVector,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE,
//							cudaMemcpyDeviceToDevice);
}

/***************************************************/
#endif /*	end of BIG_ARITHMETIC_AUX_H_P3*/
