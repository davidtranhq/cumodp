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



//#include "mixed_radix_reverse.h"
typedef unsigned int usfixn32;
typedef unsigned long long int usfixn64;

__global__ void
kernel_crt_multiplications_plain (usfixn32 * vs, usfixn32 * mult2D,
																	usfixn32 * result)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn32 s[16] =
		{ 0 };
	usfixn32 m[16];
	usfixn32 nPrimes = 16;
	usfixn32 i, j, c;
	usfixn32 k = 0;

	j = tid % nPrimes;
	usfixn32 x = vs[tid];
	usfixn64 mult;

	for (tid = 0; tid < 16; tid++)
	{
		if (j > 0)
			for (i = 0; i < nPrimes; i++)
			{
				m[i] = mult2D[(j - 1) * 16 + i];
			}

		if (j == 0)
		{
			s[0] = x;
		}

//	for(i=0;i<nPrimes;i++)
		if (j > 0)
		{
			for (i = 0; i < j; i++)
			{
				mult = x * m[i]; //v(i)*m2D(j-1,i)
				s[i] += (mult & (0xFFFFFFFF));
				mult >>= 32;
				s[i + 1] += (mult);
			}
		}
//	for(i=0;i<nPrimes;i++)
		{
		}
	}
}
/****************************************************/
__global__ void
kernel_crt_multiplications (usfixn32 * __restrict__ vs,
														const usfixn32 * __restrict__ mult2D,
														usfixn32 * result, usfixn32* parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn32 s[17] =
		{ 0 };
	memset (s, 0x00, 17 * sizeof(usfixn32));
	usfixn32 m[16] =
		{ 0 };
	usfixn32 nPrimes = 16;
	usfixn32 i, j, c;
	usfixn32 k = 0;

//	j = tid % nPrimes;
//	j = tid % 16;
	j = threadIdx.x & (0xF);
	usfixn64 x = vs[tid];
	usfixn64 mult;

	usfixn64 sum = 0;
	usfixn32 permutationStride = parameters[5];
	usfixn32 offset = 0;

	usfixn32 n = parameters[0];
	if (tid >= n)
		return;

//	for (tid = 0; tid < 16; tid++)
	{
		if (j > 0)
			for (i = 0; i < nPrimes; i++)
			{
				m[i] = mult2D[(j - 1) * 16 + i];
			}

		if (j == 0)
		{
			s[0] = x;
		}

		usfixn32 m0;
		usfixn64 carry = 0;
//	for(i=0;i<nPrimes;i++)
		if (j > 0)
		{
			for (i = 0; i < j + 1; i++)
			{
				mult = usfixn64 (x * m[i]); //v(i)*m2D(j-1,i)
				m0 = (mult & 0xFFFFFFFF);
				sum = s[i] + m0 + carry;
//				sum = s[i] + (mult);
//				if (sum < s[i] || sum < m0)
//						carry++;
				s[i] = (sum & 0xFFFFFFFF);
				mult >>= 32;
				sum >>= 32;
//				sum=(sum & 0xFFFFFFFF);
				carry = (mult) + sum;
//				printf("tid = %d, mult=%lu, carry=%lu\n",tid,mult,carry);
			}
		}

//		usfixn32 permutationStride = 1024; //inputVectorSize
		offset = 0;
		for (i = 0; i < nPrimes + 1; i++)
		{
			result[tid + offset] = s[i];
			offset += permutationStride;
		}
	}
}
/****************************************************/

__device__ void
device_sum_17_u32 (usfixn32 * s, usfixn32 *r, usfixn32 step)
{
	usfixn32 i;
	usfixn32 sum;
	usfixn32 carry = 0;

	for (i = 0; i < step; i++)
	{
		sum = s[i] + r[i];
		if (sum < s[i] || sum < r[i])
			r[i + 1]++;
		r[i] = sum;
		s[i] = 0;

//		r[i]=i;
	}
}

/****************************************************/
__device__ void
device_sum_33_u32 (usfixn32 * s, usfixn32 *r, usfixn32 step)
{
	usfixn32 i;
	usfixn32 sum;
	usfixn32 carry = 0;

	for (i = 0; i < step; i++)
	{
		sum = s[i] + r[i];
		if (sum < s[i] || sum < r[i])
			r[i + 1]++;
		r[i] = sum;
		s[i] = 0;

//		r[i]=i;
	}
}
/****************************************************/
__global__ void
kernel_crt_multiplications_v1 (usfixn32 * __restrict__ vs,
															 const usfixn32 * __restrict__ mult2D,
															 usfixn32 * result, usfixn32* parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn32 s[17] =
		{ 0 };
	memset (s, 0x00, 17 * sizeof(usfixn32));
	usfixn32 m[16] =
		{ 0 };
	usfixn32 nPrimes = 16;
	usfixn32 i, j, c;
	usfixn32 k = 0;

	usfixn32 r[17] =
		{ 0 };
//	j = tid % nPrimes;
//	j = tid % 16;
	j = threadIdx.x & (0xF);

	usfixn64 mult;

	usfixn64 sum = 0;
	usfixn32 permutationStride = parameters[5];
	usfixn32 offset = 0;

	usfixn32 n = parameters[0];
	if (tid >= n)
		return;

//	for (tid = 0; tid < 16; tid++)

	usfixn32 m0;
	usfixn64 carry = 0;

//	for(i=0;i<17;i++)
//		s[i]=16-threadIdx.x;
	if (j > 0)
		for (i = 0; i < nPrimes; i++)
		{
//			m[i] = mult2D[(j - 1) * 16 + i];
			m[i] = i;
		}

	usfixn64 x;
	for (j = 0; j < 16; j++)
	{

		x = vs[tid + j * permutationStride];

		if (j == 0)
		{
			s[0] = x;
//			continue;
		}
		carry = 0;
		m0 = 0;
//	for(i=0;i<nPrimes;i++)
		if (j > 0)
		{
			for (i = 0; i < j + 1; i++)
			{
				mult = usfixn64 (x * m[i]); //v(i)*m2D(j-1,i)
				m0 = (mult & 0xFFFFFFFF);
				sum = s[i] + m0 + carry;
//				sum = s[i] + (mult);
//				if (sum < s[i] || sum < m0)
//					s[i+1]++;
//						carry++;
//				s[i] += (sum & 0xFFFFFFFF);
				s[i] = (sum & 0xFFFFFFFF);
				mult >>= 32;
				sum >>= 32;
//				sum=(sum & 0xFFFFFFFF);
				carry = (mult) + sum;
//				printf("tid = %d, mult=%lu, carry=%lu\n",tid,mult,carry);
			}
		}
//		__syncthreads();
		device_sum_17_u32 (s, r, j);
//		__syncthreads();
//		usfixn32 permutationStride = 1024; //inputVectorSize
//		offset = 0;
//		offset = j;
//#pragma unroll 17
//		for (i = 0; i < nPrimes + 1; i++)
	}

	offset = 0;
//#pragma unroll 17
	for (i = 0; i < nPrimes; i++)
	{
		result[tid + offset] = r[i];
		offset += permutationStride;
	}
}

/****************************************************/
__global__ void
kernel_crt_multiplications_32words_v1 (usfixn32 * __restrict__ vs,
																			 const usfixn32 * __restrict__ mult2D,
																			 usfixn32 * result, usfixn32* parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn32 s[33] =
		{ 0 };
	memset (s, 0x00, 33 * sizeof(usfixn32));
	usfixn32 m[32] =
		{ 0 };
	usfixn32 nPrimes = 32;
	usfixn32 i, j, c;
	usfixn32 k = 0;

	usfixn32 r[33] =
		{ 0 };
	j = tid % nPrimes;
//	j = tid % 16;
//	j = threadIdx.x & (0xF);

	usfixn64 mult;

	usfixn64 sum = 0;
	usfixn32 permutationStride = parameters[5];
	usfixn32 offset = 0;

	usfixn32 n = parameters[0];
	if (tid >= n)
		return;

//	for (tid = 0; tid < 16; tid++)

	usfixn32 m0;
	usfixn64 carry = 0;

//	for(i=0;i<17;i++)
//		s[i]=16-threadIdx.x;
	if (j > 0)
		for (i = 0; i < nPrimes; i++)
		{
//			m[i] = mult2D[(j - 1) * 16 + i];
			m[i] = i;
		}

	usfixn64 x;
	for (j = 0; j < 32; j++)
	{

		x = vs[tid + j * permutationStride];

		if (j == 0)
		{
			s[0] = x;
//			continue;
		}
		carry = 0;
		m0 = 0;
//	for(i=0;i<nPrimes;i++)
		if (j > 0)
		{
			for (i = 0; i < j + 1; i++)
			{
				mult = usfixn64 (x * m[i]); //v(i)*m2D(j-1,i)
				m0 = (mult & 0xFFFFFFFF);
				sum = s[i] + m0 + carry;
//				sum = s[i] + (mult);
//				if (sum < s[i] || sum < m0)
//					s[i+1]++;
//						carry++;
//				s[i] += (sum & 0xFFFFFFFF);
				s[i] = (sum & 0xFFFFFFFF);
				mult >>= 32;
				sum >>= 32;
//				sum=(sum & 0xFFFFFFFF);
				carry = (mult) + sum;
//				printf("tid = %d, mult=%lu, carry=%lu\n",tid,mult,carry);
			}
		}
//		__syncthreads();
		// device_sum_17_u32(s,r,j);
		device_sum_33_u32 (s, r, j);
//		__syncthreads();
//		usfixn32 permutationStride = 1024; //inputVectorSize
//		offset = 0;
//		offset = j;
//#pragma unroll 17
//		for (i = 0; i < nPrimes + 1; i++)
	}

	offset = 0;
//#pragma unroll 17
	for (i = 0; i < nPrimes; i++)
	{
		result[tid + offset] = r[i];
		offset += permutationStride;
	}
}

/****************************************************/
__global__ void
kernel_crt_additions (usfixn32 * __restrict__ vs,
											const usfixn32 * __restrict__ mult2D, usfixn32 * result,
											usfixn32 * __restrict__ parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;

	__shared__ usfixn64 sharedSum[256];
//	__shared__ usfixn32 sharedCarry[256];
	usfixn32 i, j, c, k;

	usfixn32 s = 0;
	usfixn32 offset = 0;
	usfixn32 permutationStride = parameters[5];
	j = threadIdx.x;
	sharedSum[j] = 0;
//	usfixn64 sum = 0;

	usfixn64 sum64 = 0;

	if (tid >= parameters[0])
		return;

	usfixn32 carry = 0;
//	j = j % 16;
	j = j & (0xF);
	for (c = 0; c < 16; c++)
	{
		carry = 0;
		s = result[tid + offset];
//		sharedSum[threadIdx.x] = 1;
		sharedSum[threadIdx.x] = result[tid + offset];
		if (j == 0)
		{
			sharedSum[threadIdx.x] += sum64;
		}
		result[tid + offset] = 0;
		__syncthreads ();
//		result[tid+offset]=c;

		if (j < 8)
		{
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 8]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+8]=carry;
			sharedSum[threadIdx.x] = sum64;

//			result[tid+offset]=(sum&(0xFFFFFFFF));
		}
		__syncthreads ();

		if (j < 4)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 4];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 4]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));

			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 4]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+4]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 2)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 2];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 2]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 2]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+2]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 1)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 1];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 1]);
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 1]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);

//			sharedSum[threadIdx.x] = sum64;

//			sum = sharedSum[threadIdx.x];
			result[tid + offset] = (sum64 & (0xFFFFFFFF));
//			sum >>= 32;
//			sum = 0;
			sum64 >>= 32;
//			carry = (sum64 >> 32);
//			result[tid + offset + permutationStride] += (sum64);
		}
		__syncthreads ();

//		if(j==0)
//		{
//			sum64=0;
//			for(short k=0;k<16;k++)
//				sum64+=result[tid+k+offset];
//			result[tid+offset]=sum64&(0xFFFFFFFF);
//			sum64>>=32;
//			result[tid+offset+permutationStride]+=sum64;
//		}
		offset += permutationStride;
	}
}
/****************************************************/
__global__ void
kernel_crt_additions_v1 (usfixn32 * __restrict__ vs, usfixn32 * result,
												 usfixn32 * __restrict__ parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;

	__shared__ usfixn64 sharedSum[256];
//	__shared__ usfixn32 sharedCarry[256];
	usfixn32 i, j, c, k;

	usfixn32 s = 0;
	usfixn32 offset = 0;
	usfixn32 permutationStride = parameters[5];
	j = threadIdx.x;
	sharedSum[j] = 0;
//	usfixn64 sum = 0;

	usfixn64 sum64;

	if (tid >= parameters[0])
		return;

	usfixn32 carry = 0;
//	j = j % 16;
	j = j & (0xF);

	sum64 = 0;
#pragma unroll 16
	for (c = 0; c < 16; c++)
	{
		carry = 0;
//		s = result[tid + offset];
//		sharedSum[threadIdx.x] = 1;
		sharedSum[threadIdx.x] = result[tid + offset];
		if (j == 0)
		{
			sharedSum[threadIdx.x] += sum64;
		}
		result[tid + offset] = 0;
		__syncthreads ();
//		result[tid+offset]=c;

		if (j < 8)
		{
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 8]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+8]=carry;
			sharedSum[threadIdx.x] = sum64;

//			result[tid+offset]=(sum&(0xFFFFFFFF));
		}
		__syncthreads ();

		if (j < 4)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 4];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 4]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));

			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 4]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+4]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 2)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 2];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 2]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 2]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+2]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 1)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 1];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 1]);
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 1]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);

//			sharedSum[threadIdx.x] = sum64;

//			sum = sharedSum[threadIdx.x];
			result[tid + offset] = (sum64 & (0xFFFFFFFF));
//			sum >>= 32;
//			sum = 0;
			sum64 >>= 32;
//			carry = (sum64 >> 32);
//			result[tid + offset + permutationStride] += (sum64);
		}
		__syncthreads ();

//		if(j==0)
//		{
//			sum64=0;
//			for(short k=0;k<16;k++)
//				sum64+=result[tid+k+offset];
//			result[tid+offset]=sum64&(0xFFFFFFFF);
//			sum64>>=32;
//			result[tid+offset+permutationStride]+=sum64;
//		}
		offset += permutationStride;
	}
}

/****************************************************/
__global__ void
kernel_crt_additions_v2 (usfixn32 * __restrict__ vs, usfixn32 * result,
												 usfixn32 * __restrict__ parameters)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;

	__shared__ usfixn64 sharedSum[256];
//	__shared__ usfixn32 sharedCarry[256];
	usfixn32 i, j, c, k;

	usfixn32 s = 0;
	usfixn32 offset = 0;
	usfixn32 permutationStride = parameters[5];
	j = threadIdx.x;
	sharedSum[j] = 0;
//	usfixn64 sum = 0;

	usfixn64 sum64;

	if (tid >= parameters[0])
		return;

	usfixn32 carry = 0;
//	j = j % 16;
	j = j & (0xF);

	sum64 = 0;
#pragma unroll 16
	for (c = 0; c < 16; c++)
	{
		carry = 0;
//		s = result[tid + offset];
//		sharedSum[threadIdx.x] = 1;
		sharedSum[threadIdx.x] = result[tid + offset];
		if (j == 0)
		{
			sharedSum[threadIdx.x] += sum64;
		}
		result[tid + offset] = 0;
		__syncthreads ();
//		result[tid+offset]=c;

		if (j < 8)
		{
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 8]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+8]=carry;
			sharedSum[threadIdx.x] = sum64;

//			result[tid+offset]=(sum&(0xFFFFFFFF));
		}
		__syncthreads ();

		if (j < 4)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 4];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 4]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));

			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 4]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+4]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 2)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 2];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 2]);
//			result[tid+offset]=(sum&(0xFFFFFFFF));
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 2]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);
//			carry = (sum64 >> 32);
//			sharedSum[threadIdx.x+2]=carry;
			sharedSum[threadIdx.x] = sum64;
		}
		__syncthreads ();

		if (j < 1)
		{
//			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 1];
//			sharedSum[threadIdx.x] += (sharedSum[threadIdx.x + 1]);
			sum64 = sharedSum[threadIdx.x] + (sharedSum[threadIdx.x + 1]);
//			sharedSum[threadIdx.x] = sum64 & (0xFFFFFFFF);

//			sharedSum[threadIdx.x] = sum64;

//			sum = sharedSum[threadIdx.x];
			result[tid + offset] = (sum64 & (0xFFFFFFFF));
//			sum >>= 32;
//			sum = 0;
			sum64 >>= 32;
//			carry = (sum64 >> 32);
//			result[tid + offset + permutationStride] += (sum64);
		}
		__syncthreads ();

//		if(j==0)
//		{
//			sum64=0;
//			for(short k=0;k<16;k++)
//				sum64+=result[tid+k+offset];
//			result[tid+offset]=sum64&(0xFFFFFFFF);
//			sum64>>=32;
//			result[tid+offset+permutationStride]+=sum64;
//		}
		offset += permutationStride;
	}
}

/****************************************************/

__global__ void
kernel_crt_additions_plain (usfixn32 *vs, usfixn32 * mult2D, usfixn32 * result)
{
	usfixn32 tid = threadIdx.x + blockIdx.x * blockDim.x;

	__shared__ usfixn32 sharedSum[128];
	usfixn32 i, j, c, k;

	usfixn32 s = 0;
	usfixn32 offset = 0;
	usfixn32 permutationStride = 1024;
	j = threadIdx.x;
	sharedSum[j] = 0;
	usfixn64 sum = 0;

	j = threadIdx.x;
	j = j % 16;

	for (c = 0; c < 16; c++)
	{
		for (tid = 0; tid < 16; tid++)
		{
			s = result[tid + offset];
			sharedSum[tid] = 1;
		}

		if (j < 8)
		{
			sum = sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 4];
		}
		__syncthreads ();

		if (j < 4)
		{
			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 2];
		}
		__syncthreads ();

		if (j < 2)
		{
			sum += sharedSum[threadIdx.x] + sharedSum[threadIdx.x + 1];
		}
		__syncthreads ();

		if (j < 1)
		{
//			sum=offset;
//			sum=123;
			result[threadIdx.x + offset] = (sum & (0xFFFFFFFF));
			sum >>= 32;
			result[threadIdx.x + offset + 1] += (sum);
		}
		__syncthreads ();

		offset += permutationStride;
	}
}
