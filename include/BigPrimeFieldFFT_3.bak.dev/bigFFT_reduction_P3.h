// for 8-word prime
#ifndef BIG_FFT_REDUCTION_H_
#define BIG_FFT_REDUCTION_H_

/***************************************************/
__global__ void
kernel_p3_reduction_8_v0 (usfixn64 * xs, usfixn64 * ys, usfixn64 prime,
												usfixn64 * precomputed_r_array, short index,
												usfixn64 * parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	short c = 0;
	usfixn64 reduced = 0;
	usfixn64 tmp = 0;
	for (c = 0; c < COEFFICIENT_SIZE; c++)
	{
		tmp = xs[offset] * precomputed_r_array[c];
		tmp = tmp % prime;
		reduced += tmp;
		offset += permutationStride;
	}
	reduced = reduced % prime;
	ys[tid + index * permutationStride] = reduced;
}

/***************************************************/
__global__ void
kernel_p3_reduction_8 (const usfixn64 * __restrict__ xs, usfixn64 * __restrict__ ys,
										 const usfixn64 * __restrict__ primeList, const usfixn64 * __restrict__ precomputed_r_array,
										 const usfixn64 * __restrict__ parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 permutationStride = parameters[5];
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	short c = 0;
	usfixn64 reduced = 0;
	usfixn64 tmp = 0;
	short i = 0;
//	usfixn64 reduced[32];

	usfixn64 x[COEFFICIENT_SIZE];
	usfixn64 currentPrime;
	usfixn64 t;

//	i=tid %16;
	i = tid & 0xF;
	currentPrime = primeList[i];

	usfixn64 gridStride = blockDim.x * gridDim.x;
	n = n << 4;
	for (tid; tid < n; tid += gridStride)
	{
		offset = tid / 16;
		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
			x[c] = xs[offset];
			offset += permutationStride;
		}

//	for (i = 0; i < 16; i++)
//		i = tid % 32;
//		i = tid & 0x1F; // it remains the same for all strides over a grid!
//		tid /= 32;
//		t=tid/32;
		reduced = 0;

		for (c = 0; c < COEFFICIENT_SIZE; c++)
		{
//			tmp = xs[offset] * precomputed_r_array[c];
			tmp = x[c] * c; //precomputed_r_array[c];
			tmp = tmp % currentPrime;
			reduced += tmp;
//			offset += permutationStride;
		}
		reduced = reduced % currentPrime;

//			tid = threadIdx.x + blockIdx.x * blockDim.x;
//		if (i < 16)
//			ys1[tid + i * permutationStride] = reduced;
//		else
//			ys2[tid + i * permutationStride] = reduced;
//		if (i < 16)
		ys[tid] = reduced;
//		else
//			ys2[tid] = reduced;
	}
}
/***************************************************/
/***************************************************/

#endif /* end of BIG_FFT_REDUCTION_H_ */
