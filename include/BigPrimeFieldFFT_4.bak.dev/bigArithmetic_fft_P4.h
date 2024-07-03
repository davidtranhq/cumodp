// for 16-word prime
#ifndef BIG_ARITHMETIC_FFT_H_
#define BIG_ARITHMETIC_FFT_H_

#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_addition_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_multiplication_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_cyclicShift_P4.h"

/************************************************/

__device__ void inline
device_p4_swap_permutated (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym,
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
/************************************************/

__device__ void __inline__
device_p4_swap_permutated_2 (usfixn64 * __restrict__ xm,
	usfixn64 * __restrict__ ym,
	const usfixn64 & permutationStride)
{
	short i = 0;
	usfixn64 tmp = 0;
	usfixn64 offset = 0;

	usfixn64 y[COEFFICIENT_SIZE];
//	           , y[COEFFICIENT_SIZE];
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		tmp = xm[offset];
//		xm[offset] = ym[offset];
//		ym[offset] = tmp;
//		offset += permutationStride;
//	}
//
//	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		x[i] = xm[offset];
//		offset += permutationStride;
//	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		y[i] = ym[offset];
		offset += permutationStride;
	}

	offset = 0;
//#pragma unroll COEFFICIENT_SIZE
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		ym[offset] = xm[offset];  //x[i];
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
/************************************************/
/************************************************/

//new code for performing fft32, 16 threads for fft-32
__device__ void
device_base_fft32 (usfixn64 *xm, const usfixn64 permutationStride,
	const usfixn64 tid)
{
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32; //stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r4[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r4[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	//round 5
	idx = (tid / stride) * (stride << 1) + i;
//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r5[shNo]);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r5[shNo],
		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1; //stride /=2;
	i = i - (i >= stride) * stride;

	idx = (tid / 16) * 32;
	i = tid % (16);
	idx = idx + i;
	//round 6 permutation
	shNo = shNo_FFT_32_r6[i];
	if (shNo > 0)
		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
	shNo = shNo_FFT_32_r6[i + 16];
	idx += 16;
	if (shNo > 0)
		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
//	if(shNo<0)
//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx-shNo], permutationStride);
}
/************************************************/
/************************************************/
//FFT base 32 to be called from host
__global__ void
kernel_p4_base_fft32 (usfixn64 *xm, usfixn64 * parameters)
{
	//short rounds = 4; //log(16,2)

//	short stride = 32; //stride/2
//	stride /= 2;
	usfixn64 permutationStride = parameters[5];
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 idx = 0;
//	short i = tid % stride;
//	short shNo = tid % stride;

	device_base_fft32 (xm, permutationStride, tid);
}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r1_DFT (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters)
{
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r2_shift (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters,
	const short * __restrict__ device_shNo_FFT_32_r2)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//																		permutationStride);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride],
		device_shNo_FFT_32_r2[shNo],
		permutationStride);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
//		stride >>= 1; //stride /=2;
//		i = i - (i >= stride) * stride;
}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r2_DFT (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//			device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//																				permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;
}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r3_shift (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters,
	const short * __restrict__ device_shNo_FFT_32_r3)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//	permutationStride
//	);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
//	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
//																		permutationStride);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride],
		device_shNo_FFT_32_r3[shNo],
		permutationStride);
//	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
//	stride >>= 1;//stride /=2;
//	i = i - (i >= stride) * stride;

}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r3_DFT (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//	permutationStride
//	);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
//																			permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;
}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r4_shift (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters,
	const short * __restrict__ device_shNo_FFT_32_r4)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//	permutationStride
//	);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
//	permutationStride
//	);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r4[shNo]);
//	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r4[shNo],
//																		permutationStride);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride],
		device_shNo_FFT_32_r4[shNo],
		permutationStride);
}

/***********************************************/
__global__ void
kernel_p4_base_fft32_r4_DFT (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters)
{

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
	//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
	//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//	permutationStride
//	);
	//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
	//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
//	permutationStride
//	);
	//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r4[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r4[shNo],
//		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;
}
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_shift (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters,
	const short * __restrict__ device_shNo_FFT_32_r5)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	short i = tid % stride;
	short shNo = tid % stride;

	//round 1
	idx = (tid / stride) * (stride << 1) + i;
	//no cyclic shift in round 1
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 2
	idx = (tid / stride) * (stride << 1) + i;
	//	printf("%d \n",idx);
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r2[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r2[shNo],
//																			permutationStride);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 3
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r3[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r3[shNo],
//																			permutationStride);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 4
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r4[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r4[shNo],
//																			permutationStride);
//		device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
	stride >>= 1;	//stride /=2;
	i = i - (i >= stride) * stride;

	//round 5
	idx = (tid / stride) * (stride << 1) + i;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r5[shNo]);
//	device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r5[shNo],
//																		permutationStride);
	device_p4_cyclicShift_permutated_11 (&xm[idx + stride],
		device_shNo_FFT_32_r5[shNo],
		permutationStride);

}
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_DFT (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
	idx = (tid << 1);
	stride = 1;
	//	cyclicShift (&xm[idx + stride], shNo_FFT_32_r5[shNo]);
//		device_p4_cyclicShift_permutated_11 (&xm[idx + stride], shNo_FFT_32_r5[shNo],
//		permutationStride);
	device_p4_fft_base2_permutated (&xm[idx], &xm[idx + stride], permutationStride);
//		stride >>= 1;//stride /=2;
//		i = i - (i >= stride) * stride;
//
//		idx = (tid / 16) * 32;
//		i = tid % (16);
//		idx = idx + i;
//		//round 6 permutation
//		shNo = shNo_FFT_32_r6[i];
//		if (shNo > 0)
//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
//		shNo = shNo_FFT_32_r6[i + 16];
//		idx += 16;
//		if (shNo > 0)
//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx + shNo], permutationStride);
	//	if(shNo<0)
	//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx-shNo], permutationStride);
}
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_permutation (usfixn64 * xm,
	const usfixn64 * __restrict__ parameters,
	const short * __restrict__ shNoListGlobal)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	//short rounds = 4; //log(16,2)
	//short rounds = 5; //log(32,2)
	short stride = 32;	//stride/2
	stride /= 2;
	usfixn64 idx = 0;
//	short i = tid % stride;
//	short shNo = tid % stride;

//	idx = (tid / 16) * 32;
//	i = tid % (16);
	short i, shNo;
	idx = (tid >> 4) << 5;
//	idx = tid & 0xFFFFFFFFFFFFFFF0L;
	i = tid & 0xF;
	idx = idx + i;
	//round 6 permutation
//	shNo = shNo_FFT_32_r6[i];
	shNo = shNoListGlobal[i];

	if (shNo > 0)
		device_p4_swap_permutated (&xm[idx], &xm[idx + shNo], permutationStride,
			COEFFICIENT_SIZE);
//	shNo = shNo_FFT_32_r6[i + 16];
	shNo = shNoListGlobal[i + 16];
	idx += 16;
	if (shNo > 0)
		device_p4_swap_permutated (&xm[idx], &xm[idx + shNo], permutationStride,
			COEFFICIENT_SIZE);
	//	if(shNo<0)
	//		device_p4_swap_permutated_2 (&xm[idx], &xm[idx-shNo], permutationStride);
}
/***********************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step1 (usfixn64 * xs, usfixn64* powers_omega,
	usfixn64 K, usfixn64 l,
	usfixn64* parameters)
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

	if (h == 0)
		return;
//	printf("h=%lu, tid=%lu\n",h,tid);
	usfixn64 offset = tid;
	usfixn64 stat = 0;
//	for(short c=0;c<COEFFICIENT_SIZE;c++)
//		{
//		if(xs[offset]==0)
//			stat++;
//		offset+=permutationStride;
//		}
//	if(stat==COEFFICIENT_SIZE)
//		{
//		printf("tid = %lu, i=%lu, j=%lu, h=%lu, w=%lu \n", tid, i,j, h,w);
//		return;
//		}
//	device_cyclicShift_twice_permutated_1 (&xs[tid], h, permutationStride);
	usfixn64 zero[COEFFICIENT_SIZE] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 x[COEFFICIENT_SIZE];
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		x[i] = xs[tid + i * permutationStride];
	}

	if (h >= COEFFICIENT_SIZE)
	{
//		device_p4_cyclicShift_permutated_11 (&xs[tid], COEFFICIENT_SIZE,
//																			permutationStride);
		device_p4_bigPrimeSub_permutated_2 (zero, x, 1);
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			xs[tid + i * permutationStride] = zero[i];
		}
		h -= COEFFICIENT_SIZE;
	}
	device_p4_cyclicShift_permutated_11 (&xs[tid], h, permutationStride);
}
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step2 (usfixn64 * xs, usfixn64* powers_omega,
	usfixn64 K, usfixn64 l,
	usfixn64* parameters)
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
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
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
//		device_bigMult_plain (x, y);
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

/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21_v0 (
	const usfixn64 * __restrict__ xs, const usfixn64* __restrict__ powers_omega,
	const usfixn64 K, const usfixn64 l, usfixn64* __restrict__ lVector,
	usfixn64* __restrict__ hVector, usfixn64* __restrict__ cVector,
	usfixn32* __restrict__ signVector, const usfixn64* __restrict__ parameters)
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
	usfixn64 u1;
//		short i = 0, j = 0;
	usfixn64 sign = 0;

//	stride = l / K;
	stride = l;
	i0 = (tid % (l * K)) % stride;
	j0 = (tid % (l * K)) / stride;
	if ((i0 * j0) == 0)
		return;

//	h = (i0 * j0) / stride; // h <K
	w = (i0 * j0) % stride; // w < N/K

//	usfixn64 x[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 offset = tid;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		x[i] = xs[offset];
//		offset += permutationStride;
//	}

	usfixn64 h0, h1, h2;
	usfixn64 xOffset;
//	xOffset_init = tid
//			+ (COEFFICIENT_SIZE - 1) * permutationStride;
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)

	if (w == 0)
		return;

//	if (w > 0)
	{
		offset = w - 1;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			y[COEFFICIENT_SIZE - i - 1] = powers_omega[offset];
//			if (tid == 620)
//						{
//							printf ("powers_omega[i]=%lu \n", powers_omega[offset]);
//						}
			offset += stride;
		}

//			device_bigMult_plain (x, y);
	}
//	else
//	{
//		y[COEFFICIENT_SIZE - 1] = 1;
//		for (i = 1; i < COEFFICIENT_SIZE; i++)
//				{
//					y[COEFFICIENT_SIZE - i - 1] = 0;
////					offset += stride;
//				}
//	}
	{
//		offset = tid;
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			x[i] = xs[offset];
//			y[7 - i] = ys[offset];
//			offset += permutationStride;
//		}

		//offset = tid + (permutationStride *COEFFICIENT_SIZE);
		offset = tid + (permutationStride << 4);
//
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			offset -= permutationStride;
//			u[0] = 0;
//			u[1] = 0;
//			u[2] = 0;
//			uSub[0] = 0;
//			uSub[1] = 0;
//			uSub[2] = 0;
//			if (i > 0)
//				//				device_oneShiftRight (y);
//				device_p4_oneShiftLeft (y);
//
//			xOffset = xOffset_init;
//			for (j = COEFFICIENT_SIZE - 1; j >= 0; j--)
//			{
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[j])
//				);
//
//				xOffset -= permutationStride;
//				h0 = u[0];
//				h1 = u[1];
//				h2 = u[2];
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					h0 = uSub[0];
//					h1 = uSub[1];
//					h2 = uSub[2];
//				}
//
//				/*****************************/
//				//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"=l"(c):"l"(h0),"l"(m0));
//
//				h0 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%1,%2;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(m1));
//				usfixn64 u1 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,%1,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
//				h1 = s;
//				h2 += c;
//
//				/*****************************/
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
//			}
//			sign = 0;
//			device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
//			lVector[offset] = u[0];
//			hVector[offset] = u[1];
//			cVector[offset] = u[2];
//			signVector[offset] = sign;
//		}

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

//			xOffset = xOffset_init;
			xOffset = tid + (COEFFICIENT_SIZE - 1) * permutationStride;
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

/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21 (const usfixn64 * xs,
	const usfixn64* powers_omega,
	const usfixn64 K, const usfixn64 l,
	usfixn64* lVector, usfixn64* hVector,
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
	usfixn64 u1;
//		short i = 0, j = 0;
	usfixn64 sign = 0;

//	stride = l / K;
//	usfixn64 x[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	usfixn64 lVectorLocal[COEFFICIENT_SIZE];
	usfixn64 hVectorLocal[COEFFICIENT_SIZE];
	usfixn64 cVectorLocal[COEFFICIENT_SIZE];
	usfixn64 signVectorLocal[COEFFICIENT_SIZE];
	usfixn64 offset = tid;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		x[i] = xs[offset];
//		offset += permutationStride;
//	}

	usfixn64 h0, h1, h2;
	usfixn64 xOffset;
//	xOffset_init = tid
//			+ (COEFFICIENT_SIZE - 1) * permutationStride;
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	usfixn64 gridStride = blockDim.x * gridDim.x;

	stride = l;
	usfixn64 lK = l * K;
	for (tid; tid < n; tid += gridStride)
	{
//		i0 = (tid % (l * K)) % stride;
//		j0 = (tid % (l * K)) / stride;
		i0 = (tid % (lK)) % stride;
		j0 = (tid % (lK)) / stride;
		if ((i0 * j0) == 0)
		{
//			return;
			continue;
		}

		//	h = (i0 * j0) / stride; // h <K
		w = (i0 * j0) % stride; // w < N/K

		if (w == 0)
		{
//			return;
			continue;
		}

//	if (w > 0)
		{
			offset = w - 1;
			for (i = 0; i < COEFFICIENT_SIZE; i++)
			{
				y[COEFFICIENT_SIZE - i - 1] = powers_omega[offset];
//			if (tid == 620)
//						{
//							printf ("powers_omega[i]=%lu \n", powers_omega[offset]);
//						}
				offset += stride;
			}

//			device_bigMult_plain (x, y);
		}
//	else
//	{
//		y[COEFFICIENT_SIZE - 1] = 1;
//		for (i = 1; i < COEFFICIENT_SIZE; i++)
//				{
//					y[COEFFICIENT_SIZE - i - 1] = 0;
////					offset += stride;
//				}
//	}
		{
//		offset = tid;
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			x[i] = xs[offset];
//			y[7 - i] = ys[offset];
//			offset += permutationStride;
//		}

			//offset = tid + (permutationStride *COEFFICIENT_SIZE);
//			offset = tid + (permutationStride << 4);
//
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			offset -= permutationStride;
//			u[0] = 0;
//			u[1] = 0;
//			u[2] = 0;
//			uSub[0] = 0;
//			uSub[1] = 0;
//			uSub[2] = 0;
//			if (i > 0)
//				//				device_oneShiftRight (y);
//				device_p4_oneShiftLeft (y);
//
//			xOffset = xOffset_init;
//			for (j = COEFFICIENT_SIZE - 1; j >= 0; j--)
//			{
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[j])
//				);
//
//				xOffset -= permutationStride;
//				h0 = u[0];
//				h1 = u[1];
//				h2 = u[2];
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					h0 = uSub[0];
//					h1 = uSub[1];
//					h2 = uSub[2];
//				}
//
//				/*****************************/
//				//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"=l"(c):"l"(h0),"l"(m0));
//
//				h0 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%1,%2;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(m1));
//				usfixn64 u1 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,%1,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
//				h1 = s;
//				h2 += c;
//
//				/*****************************/
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
//			}
//			sign = 0;
//			device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
//			lVector[offset] = u[0];
//			hVector[offset] = u[1];
//			cVector[offset] = u[2];
//			signVector[offset] = sign;
//		}

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

//			xOffset = xOffset_init;
				xOffset = tid + (COEFFICIENT_SIZE - 1) * permutationStride;
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

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
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
//			lVectorLocal[i]=u[0];
//			hVectorLocal[i]=u[1];
//			cVectorLocal[i]=u[2];
//			signVectorLocal[i]=sign;
				//			hVectorSub[offset] = uSub[1];
				//			cVectorSub[offset] = uSub[2];
				//			offset=offset-permutationStride;
			}
		}
	}

//	w=1;

//	offset = tid+(COEFFICIENT_SIZE-1)*permutationStride;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		lVector[offset] = lVectorLocal[i];
////		hVector[offset] = hVectorLocal[i];
////		cVector[offset] = cVectorLocal[i];
////		signVector[offset] = signVectorLocal[i];
////		xs[offset] = x[i];
////				xs[offset] = y[i];
//		offset -= permutationStride;
//	}
}
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21_lessReg (const usfixn64 * xs,
	const usfixn64* powers_omega,
	const usfixn64 K,
	const usfixn64 l,
	usfixn64* lVector,
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
//		usfixn64 offset = 0;
//		usfixn64 permutationStride = parameters[5];
//	usfixn64 m0 = 0, m1 = 0;
//	usfixn64 s = 0, c = 0;

//		short i = 0, j = 0;
//	usfixn64 sign = 0;

//	stride = l / K;
//	usfixn64 x[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 y[COEFFICIENT_SIZE] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//	usfixn64 lVectorLocal[COEFFICIENT_SIZE];
//	usfixn64 hVectorLocal[COEFFICIENT_SIZE];
//	usfixn64 cVectorLocal[COEFFICIENT_SIZE];
//	usfixn64 signVectorLocal[COEFFICIENT_SIZE];
	usfixn64 offset = tid;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		x[i] = xs[offset];
//		offset += permutationStride;
//	}

//	usfixn64 xOffset;
//	xOffset_init = tid
//			+ (COEFFICIENT_SIZE - 1) * permutationStride;
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	usfixn64 gridStride = blockDim.x * gridDim.x;

	stride = l;
	usfixn64 lK = l * K;
	for (tid; tid < n; tid += gridStride)
	{
//		i0 = (tid % (l * K)) % stride;
//		j0 = (tid % (l * K)) / stride;
		i0 = (tid % (lK)) % stride;
		j0 = (tid % (lK)) / stride;
		if ((i0 * j0) == 0)
		{
//			return;
			continue;
		}

		//	h = (i0 * j0) / stride; // h <K
		w = (i0 * j0) % stride; // w < N/K

		if (w == 0)
		{
//			return;
			continue;
		}

//	if (w > 0)
		{
			offset = w - 1;
			for (i = 0; i < COEFFICIENT_SIZE; i++)
			{
				y[COEFFICIENT_SIZE - i - 1] = powers_omega[offset];
//			if (tid == 620)
//						{
//							printf ("powers_omega[i]=%lu \n", powers_omega[offset]);
//						}
				offset += stride;
			}

//			device_bigMult_plain (x, y);
		}
//	else
//	{
//		y[COEFFICIENT_SIZE - 1] = 1;
//		for (i = 1; i < COEFFICIENT_SIZE; i++)
//				{
//					y[COEFFICIENT_SIZE - i - 1] = 0;
////					offset += stride;
//				}
//	}
		{
//		offset = tid;
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			x[i] = xs[offset];
//			y[7 - i] = ys[offset];
//			offset += permutationStride;
//		}

			//offset = tid + (permutationStride *COEFFICIENT_SIZE);
//			offset = tid + (permutationStride << 4);
//
//		for (i = 0; i < COEFFICIENT_SIZE; i++)
//		{
//			offset -= permutationStride;
//			u[0] = 0;
//			u[1] = 0;
//			u[2] = 0;
//			uSub[0] = 0;
//			uSub[1] = 0;
//			uSub[2] = 0;
//			if (i > 0)
//				//				device_oneShiftRight (y);
//				device_p4_oneShiftLeft (y);
//
//			xOffset = xOffset_init;
//			for (j = COEFFICIENT_SIZE - 1; j >= 0; j--)
//			{
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[j])
//				);
//
//				xOffset -= permutationStride;
//				h0 = u[0];
//				h1 = u[1];
//				h2 = u[2];
//				if (j > COEFFICIENT_SIZE - 1 - i)
//				{
//					h0 = uSub[0];
//					h1 = uSub[1];
//					h2 = uSub[2];
//				}
//
//				/*****************************/
//				//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"=l"(c):"l"(h0),"l"(m0));
//
//				h0 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%1,%2;\n\t"
//						"addc.cc.u64 %1,0,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(m1));
//				usfixn64 u1 = s;
//				asm("{\n\t"
//						"add.cc.u64 %0,%2,%3;\n\t"
//						"addc.cc.u64 %1,%1,0;\n\t"
//						"}"
//						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
//				h1 = s;
//				h2 += c;
//
//				/*****************************/
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
//			}
//			sign = 0;
//			device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
//			lVector[offset] = u[0];
//			hVector[offset] = u[1];
//			cVector[offset] = u[2];
//			signVector[offset] = sign;
//		}

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
				usfixn64 u[3] =
				{ 0, 0, 0 };
				usfixn64 uSub[3] =
				{ 0, 0, 0 };
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

//			xOffset = xOffset_init;
				usfixn64 xOffset;
				xOffset = tid + (COEFFICIENT_SIZE - 1) * permutationStride;
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

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);
//

					usfixn64 m0 = 0, m1 = 0;

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
					usfixn64 h0, h1, h2;
					h0 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[0] : uSub[0];
					h1 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[1] : uSub[1];
					h2 = (j <= COEFFICIENT_SIZE - 1 - i) ? u[2] : uSub[2];

					/*****************************/
					//				add128 (h0, h1, m0, m1, h2, h0, h1, h2);
					//				add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
					//	usfixn64 s = 0, c = 0;
					usfixn64 s = 0, c = 0;
					usfixn64 u1;
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
				usfixn64 sign = 0;
				sign = 0;
				device_p4_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);

				lVector[offset] = u[0];
				hVector[offset] = u[1];
				cVector[offset] = u[2];
				signVector[offset] = sign;
//			lVectorLocal[i]=u[0];
//			hVectorLocal[i]=u[1];
//			cVectorLocal[i]=u[2];
//			signVectorLocal[i]=sign;
				//			hVectorSub[offset] = uSub[1];
				//			cVectorSub[offset] = uSub[2];
				//			offset=offset-permutationStride;
			}
		}
	}

//	w=1;

//	offset = tid+(COEFFICIENT_SIZE-1)*permutationStride;
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
//		lVector[offset] = lVectorLocal[i];
////		hVector[offset] = hVectorLocal[i];
////		cVector[offset] = cVectorLocal[i];
////		signVector[offset] = signVectorLocal[i];
////		xs[offset] = x[i];
////				xs[offset] = y[i];
//		offset -= permutationStride;
//	}
}
/**********************************************************/
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step22 (usfixn64 * xs, usfixn64 K,
	usfixn64 l_input, usfixn64* lVector,
	usfixn64* hVector, usfixn64* cVector,
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

//	usfixn64 u[3] =
//		{ 0, 0, 0 };
//	usfixn64 uSub[3] =
//		{ 0, 0, 0 };
	//		usfixn64 offset = 0;
	//		usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0;
	//		short i = 0, j = 0;
	usfixn64 sign = 0;

	//	stride = l / K;
	stride = l_input;
	i0 = (tid % (l_input * K)) % stride;
	j0 = (tid % (l_input * K)) / stride;
	if ((i0 * j0) == 0)
		return;

//	h = (i0 * j0) / stride; // h <K
	w = (i0 * j0) % stride; // w < N/K
	if (w == 0)
		return;

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
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)
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
//				if(tid==620)
//										{printf("before sub3: l, h , c = %lu, %lu, %lu \n",l,h,c);
//										}
				device_p4_bigSubZero_3 (l, h, c);
//			device_p4_bigSub_plain (v1, v0);
//			l = v1[0];
//			h = v1[1];
//			c = v1[2];
//				if(tid==620)
//							{
//								printf("l, h , c = %lu, %lu, %lu \n",l,h,c);
//								printf("================================\n");
//							}
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
kernel_p4_twiddle_general_permutated_step23 (usfixn64 * xs, usfixn64 K,
	usfixn64 l_input, usfixn64* lVector,
	usfixn64* hVector, usfixn64* cVector,
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

//	usfixn64 u[3] =
//		{ 0, 0, 0 };
//	usfixn64 uSub[3] =
//		{ 0, 0, 0 };
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

//	h = (i0 * j0) / stride; // h <K
	w = (i0 * j0) % stride; // w < N/K
	if (w == 0)
		return;

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
//	for (tid; tid < n; tid += gridDim.x * blockDim.x)
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
//			if (tid == 620)
//			{
//				printf ("lvector[i]=%lu, hvector[i]=%lu, cvector[i]=%lu \n", lVector[offset], hVector[offset], cVector[offset]);
//			}
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
					device_p4_cyclicShift_permutated_12 (complement);
//				device_cyclicShift_lhc_complement (complement, i);
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
			device_p4_bigPrimeSub_plain_inPlace (v0.i, v1.i);

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
			device_p4_bigPrimeSub_plain_inPlace (v0.i, v1.i);

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
////	device_p4_bigSub_plain(v0_sub.i, v1.i);
////	device_p4_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);
//
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = cVectorSub[tid + 6 * permutationStride];
////	v1.i[1] = cVectorSub[tid + 7 * permutationStride];
//
////	device_bigPrimeAdd_plain(v1.i, v2_sub.i);
////	device_p4_bigSub_plain(v0_sub.i, v1.i);
//	device_p4_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);
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
//	device_bigPrimeAdd_plain (v0_sub.i, v1.i);
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
//	device_bigPrimeAdd_plain (v0_sub.i, v1.i);
//
////#################################################
////	device_p4_bigSub_plain(v0.i, v0_sub.i);
//	device_p4_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
////#################################################

//	############################# writing back to g-memory
			offset = tid;
			for (i = 0; i < COEFFICIENT_SIZE; i++)
			{
//			if (v0.i[i] >= R && i < COEFFICIENT_SIZE - 1)
//			{
//				v0.i[i] -= R;
//				v0.i[i + 1]++;
//			}
				xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
				offset += permutationStride;
			}
		}
	}

#endif
