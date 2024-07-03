#ifndef BIG_ARITHMETIC_ADDITION_H_P3
#define BIG_ARITHMETIC_ADDITION_H_P3

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_cyclicShift_P3.h"

/**********************************************/
///**********************************************/
//__device__ __inline__ void
//device_p3_bigPrimeSub_plain_ptx_v0_local_add (usfixn64 * __restrict__ xm,
//																							usfixn64 * __restrict__ ym);
//
///**********************************************/
///**********************************************/
//__device__ void
//device_p3_bigSub_uVector_plain_2_local_add (uvector8 & __restrict__ xm,
//																						const uvector8 & __restrict__
//																						ym);
///**********************************************/
//
////working and faster than permutated_3 and permutated_4
//__device__ void inline
//device_p3_cyclicShift_permutated_5_local_add (usfixn64 * __restrict__ xs,
//																							const short sn,
//																							const usfixn64 permutationStride);
//
///**********************************************/
//
//__device__ void inline
//device_p3_cyclicShift_permutated_5_local_add (usfixn64 * __restrict__ xs,
//																							const short sn,
//																							const usfixn64 permutationStride);
//
///**********************************************/
////computing cyclic shift with vectors of size 8 on register
//__device__ __inline__ void
//device_p3_cyclicShift_permutated_7_local_add (
//		usfixn64 * __restrict__ xs, const short & sn,
//		const usfixn64 & permutationStride);
//
///**********************************************/
//
//__device__ void
//device_p3_shiftRight_uConstArray8 (uConstArray8_align8 & x, usfixn64 & tmp);
//
///**********************************************/
//__device__ void
//device_p3_bigPrimeAdd_correct (usfixn64 *xm, usfixn64 *ym, usfixn64 *um);
//
///**********************************************/
//__device__ __inline__ void
//device_p3_bigPrimeAdd_plain (usfixn64 *__restrict__ xm,
//														 usfixn64 *__restrict__ ym);
//
///**********************************************/
//__device__ __inline__ void
//device_p3_bigPrimeAdd_plain_inPlace (usfixn64 * __restrict__ xm,
//																		 usfixn64 * __restrict__ ym);
//
///**********************************************/
//__device__ __inline__ void
//device_p3_bigPrimeAdd_plain_ptx_v0 (usfixn64 * __restrict__ xm,
//																		usfixn64 * __restrict__ ym);
//
///**********************************************/
//__device__ void
//device_p3_bigPrimeAdd_permutated (const usfixn64 *xm, const usfixn64 * ym,
//																	usfixn64 *um, const usfixn64 idx,
//																	const short permutationStride);
//
///**********************************************/
//__device__ void
//device_p3_bigPrimeAdd_permutated_ptx_v0 (usfixn64 * xm, usfixn64 *ym,
//																				 const usfixn64 permutationStride);
///**********************************************/
//
////xm = xm + ym
////ym = xm - ym
////__device__ inline void fft_base2_permutated(usfixn64 * __restrict__ xm,
////		usfixn64 * __restrict__ ym, const usfixn64 xIdx,const usfixn64 yIdx,
////		const usfixn64 permutationStride)
//__device__ __inline__ void
//device_p3_fft_base2_permutated (usfixn64 * xm, usfixn64 * ym,
//																const usfixn64 & permutationStride);
//
///**********************************************/
////using vectorized data structure + offset is computed more efficiently
//__device__ inline void
//device_p3_fft_base2_permutated_4 (usfixn64 * __restrict__ xm,
//																	usfixn64 * __restrict__ ym,
//																	const usfixn64 & permutationStride);
///**********************************************/
//__device__ inline void
//device_p3_fft_base2_permutated_7 (usfixn64 * __restrict__ xm,
//																	usfixn64 * __restrict__ ym, const short shNo,
//																	const usfixn64 permutationStride);
///**********************************************/
////using vectorized data structure + offset is computed more efficiently
//__device__ void
//device_p3_fft_base2_permutated_8 (usfixn64 * __restrict__ xm,
//																	usfixn64 * __restrict__ ym,
//																	const short & shNo,
//																	const usfixn64 & permutationStride);

/**********************************************/
__global__ void
kernel_p3_addition_plain (usfixn64 * xs, usfixn64 * ys, usfixn64 *parameters);

/**********************************************/
__global__ void
kernel_p3_addition_permutated (usfixn64 *xs, usfixn64 *ys,
															 usfixn64 *parameters);

/**********************************************/
#endif /* BIG_ARITHMETIC_ADDITION_H_P3 */
