#ifndef BIG_ARITHMETIC_AUX_H_P3
#define BIG_ARITHMETIC_AUX_H_P3

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
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

///**********************************************/
//__global__ void
//kernel_p3_permutation_256_permutated (usfixn64 * xs, usfixn64 * ys,
//																			usfixn64* parameters, int lp, int mp);
//
///**********************************************/
//__global__ void
//kernel_p3_permutation_256_permutated_2 (usfixn64 * xs, usfixn64 * ys,
//																				usfixn64* parameters);
//
///**********************************************/
////general 256_permutation, equivalent of permutation L(k^e, 256)
////ys = L(xs)
//__global__ void
//kernel_p3_permutation_256_general_permutated_v0 (usfixn64 n,
//																								 usfixn64 * __restrict__ xs,
//																								 usfixn64 * __restrict__ ys,
//																								 usfixn64 * parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_permutation_2_general_permutated_v0 (
//		usfixn16 n, usfixn64 * __restrict__ xs, usfixn64 * ys,
//		const usfixn64 * __restrict__ parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_transposition (usfixn64* xs, usfixn64 * ys,
//												 usfixn64 nPermutationBlocks, usfixn64 blockStride,
//												 usfixn64* parameters);
//
///**********************************************/
//__global__ void
//kernel_p3_cpy_d2d (usfixn64* __restrict__ dest,
//									 const usfixn64 * __restrict__ src, const usfixn64 n,
//									 const usfixn64* __restrict__ parameters);
//
///**********************************************/
//void
//host_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd);
//
///**********************************************/
////compute L_{k,m} for any "input size", k, and m are multiples of 32 (TILE_DIMENSION)
////computes input/(k*m) block permutations
//void
//host_transposition (data* fd, usfixn64 k, usfixn64 m);

/****************************************************/
__global__ void
kernel_p3_permutation_16_plain (usfixn64 * xs, usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_permutation_64_plain (usfixn64 * xs, usfixn64* parameters);

__global__ void
kernel_p3_permutation_256_plain (usfixn64 * xs, usfixn64 * ys,
																 usfixn64* parameters, int lp, int mp);

/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated (usfixn64 * __restrict__ xs,
																		 usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated_v1 (
		usfixn64 * __restrict__ xs, const usfixn64* __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_permutation_16_permutated_v2 (usfixn64 * __restrict__ xs,
																				const usfixn64* __restrict__ parameters);
/**********************************************/
__global__ void
kernel_p3_permutation_256_permutated (usfixn64 * xs, usfixn64 * ys,
																			usfixn64* parameters, int lp, int mp);
/**********************************************/
__global__ void
kernel_p3_permutation_256_permutated_2 (usfixn64 * xs, usfixn64 * ys,
																				usfixn64* parameters);
/**********************************************/
//general 256_permutation, equivalent of permutation L(k^e, 256)
//ys = L(xs)
__global__ void
kernel_p3_permutation_256_general_permutated_v0 (usfixn64 n,
																								 usfixn64 * __restrict__ xs,
																								 usfixn64 * __restrict__ ys,
																								 usfixn64 * parameters);
/**********************************************/
//__global__ void
//kernel_p3_permutation_256_general_permutated_v0 (usfixn64 n,
//																								 usfixn64 * __restrict__ xs,
//																								 usfixn64 * __restrict__ ys,
//																								 usfixn64 * parameters)
/**********************************************/
__global__ void
kernel_p3_permutation_2_general_permutated_v0 (
		usfixn16 n, usfixn64 * __restrict__ xs, usfixn64 * ys,
		const usfixn64 * __restrict__ parameters);
/**********************************************/
__global__ void
kernel_p3_transposition (usfixn64* xs, usfixn64 * ys,
												 usfixn64 nPermutationBlocks, usfixn64 blockStride,
												 usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_cpy_d2d (usfixn64* __restrict__ dest,
									 const usfixn64 * __restrict__ src, const usfixn64 n,
									 const usfixn64* __restrict__ parameters);
/**********************************************/
void
host_p3_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd);
/**********************************************/
//compute L_{k,m} for any "input size", k, and m are multiples of 32 (TILE_DIMENSION)
//computes input/(k*m) block permutations
void
host_p3_transposition (data* fd, usfixn64 k, usfixn64 m);
/**********************************************/
#endif /*	end of BIG_ARITHMETIC_AUX_H_P3*/
