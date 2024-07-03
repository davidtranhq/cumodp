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


#ifndef BIG_ARITHMETIC_AUX_CU_P4
#define BIG_ARITHMETIC_AUX_CU_P4

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"
#include "bigArithmetic_aux_P4_device.h"
/****************************************************/
//__global__ void
//kernel_p4_fft2_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters)
//{
//	short operation = parameters[0];
//	usfixn64 permutationStride = parameters[5];
////	short shuffle = parameters[6];
//	short padding = 0;	//parameters[?]
//	usfixn64 idx;
//	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);
//
//	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
//	//following indexing is slightly faster than above indexing
//
//	idx = tid;
//	if (padding == 0)
//		fft_base2_permutated (&xs[idx], &ys[idx], permutationStride);
//}

/****************************************************/

__global__ void
kernel_p4_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters)
{

//	implement one of newMult functions here
}

/************************************************/

__global__ void
kernel_p4_transposition_v0 (usfixn64* xs, usfixn64 * ys,
												 usfixn64 nPermutationBlocks, usfixn64 blockStride,
												 usfixn64* parameters)
{
	__shared__ usfixn64 shmem[32][32];

	int blockDim = 32;
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

	for (usfixn64 blockNo = 0; blockNo < nPermutationBlocks; blockNo++)
	{
//		blockOffset = blockStride * blockNo;
//		if (x==0 && y==0)
//			{
//			usfixn64 maxstride= (y + j) + x * height + c * permutationStride ;
////				printf("blockoffset = %lu \n", blockOffset);
//			printf("maxstride = %lu \n", maxstride);
//			}
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
//		for (int j = 0; j < blockDim; j += blockDim)
			shmem[threadIdx.y + j][threadIdx.x] = xs[(y + j) * width + x
					+ c * permutationStride + blockOffset];

//			__syncthreads ();

//		x = blockIdx.y * blockDim + threadIdx.x;  // transpose block offset
//		y = blockIdx.x * blockDim + threadIdx.y;

			ys[(y + j) + x * height + c * permutationStride + blockOffset] =
					shmem[threadIdx.y + j][threadIdx.x];
		}
		blockOffset += blockStride;

	}
}
/********************************************************/
__global__ void
kernel_p4_transposition (usfixn64* xs, usfixn64 * ys, usfixn64 nPermutationBlocks,
											usfixn64 blockStride, usfixn64* parameters)
{
	__shared__ usfixn64 shmem[32][32];

	int blockDim = 32;
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

/********************************************************/

__global__ void
kernel_p4_cpy_d2d_v0 (usfixn64* __restrict__ dest,
									 const usfixn64 * __restrict__ src, const usfixn64 n,
									 const usfixn64* __restrict__ parameters)
{
//	__shared__ usfixn64 shmem[32][32];

//	int blockDim = 32;
	usfixn64 offset = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	short c = 0, j = 0;
	for (tid; tid < n; tid += blockDim.x * gridDim.x)
	{
//		tid = threadIdx.x + blockIdx.x * blockDim.x;
		offset = tid;
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
			dest[offset] = src[offset];
			offset += permutationStride;
		}
	}
}
/********************************************************/

__global__ void
kernel_p4_cpy_d2d (usfixn64* __restrict__ dest, const usfixn64 * __restrict__ src,
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

/********************************************************/
void
host_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd)
{
	kernel_p4_cpy_d2d <<<TILE_DIMENSION, 256>>> (dest, src, count,
																						fd->device_parameters);
}

#endif
