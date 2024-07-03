#ifndef BIG_ARITHMETIC_MULTIPLICATION_H_
#define BIG_ARITHMETIC_MULTIPLICATION_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"

//#include "BigPrimeFieldFFT_3/bigArithmetic_addition_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"
/**********************************************/

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread (usfixn64 * xs, usfixn64* ys,
																		 usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread_8lhc (usfixn64 * xs, usfixn64* ys,
																					usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_mult_revised_v0 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_mult_revised_step1 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub);

/**********************************************/
__global__ void
kernel_p3_mult_revised_step2 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub);
/**********************************************/
__global__ void
kernel_p3_mult_revised_step3 (usfixn64* __restrict__ xs,
															usfixn64* __restrict__ lVector,
															usfixn64 * __restrict__ hVector,
															usfixn64* __restrict__ cVector,
															usfixn64* __restrict__ lVectorSub,
															usfixn64 * __restrict__ hVectorSub,
															usfixn64* __restrict__ cVectorSub,
															usfixn64 * __restrict__ parameters);

/************************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
																	 const usfixn64* __restrict__ ys,
																	 usfixn64* parameters,
																	 usfixn64* __restrict__ lVector,
																	 usfixn64* __restrict__ hVector,
																	 usfixn64* __restrict__ cVector,
																	 usfixn32* __restrict__ signVector);

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
																	 usfixn64* __restrict__ lVector,
																	 usfixn64* __restrict__ hVector,
																	 usfixn64* __restrict__ cVector,
																	 const usfixn32* __restrict__ signVector);

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step3 (usfixn64* __restrict__ xs,
																	 const usfixn64* __restrict__ lVector,
																	 const usfixn64 * __restrict__ hVector,
																	 const usfixn64* __restrict__ cVector,
																	 const usfixn32* __restrict__ signVector,
																	 const usfixn64 * __restrict__ parameters);

/************************************************/
__global__ void
kernel_p3_lhc (usfixn64 * xs, usfixn64* ys, usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_bigMult_plain (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters);

/**********************************************/
__global__ void
kernel_p3_bigMult_plain_2 (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters);

/**********************************************/
__global__ void
kernel_p3_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters);

/**********************************************/

//xs = xs*ys
//most efficient big multiplication
//not verified
//uses both shared memory, lmem, register, global mem
//try to fix negate and bigAddPlain to have static indexing
__global__ void
kernel_p3_newMult14 (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
										 usfixn64 * __restrict__ parameters, usfixn64* lVector,
										 usfixn64 *hVector, usfixn64* cVector);

/**********************************************/
__global__ void
kernel_p3_newMult14_incStep (usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_newMult14_join2 (usfixn64* __restrict__ xs,
													 usfixn64 * __restrict__ parameters,
													 usfixn64* __restrict__ lVector,
													 usfixn64 * __restrict__ hVector,
													 usfixn64* __restrict__ cVector);

/**********************************************/
__global__ void
kernel_p3_mulLong_revised (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_mulLong_plain (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters);

/**********************************************/

__global__ void
kernel_p3_mult_revised_singleThread (usfixn64 * xs, usfixn64* ys,
																		 usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread_8lhc (usfixn64 * xs, usfixn64* ys,
																					usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_mult_revised_v0 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_mult_revised_step1 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub);
/**********************************************/
__global__ void
kernel_p3_mult_revised_step2 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub);

/**********************************************/
__global__ void
kernel_p3_mult_revised_step3 (usfixn64* __restrict__ xs,
															usfixn64* __restrict__ lVector,
															usfixn64 * __restrict__ hVector,
															usfixn64* __restrict__ cVector,
															usfixn64* __restrict__ lVectorSub,
															usfixn64 * __restrict__ hVectorSub,
															usfixn64* __restrict__ cVectorSub,
															usfixn64 * __restrict__ parameters);

/************************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
																	 const usfixn64* __restrict__ ys,
																	 usfixn64* parameters,
																	 usfixn64* __restrict__ lVector,
																	 usfixn64* __restrict__ hVector,
																	 usfixn64* __restrict__ cVector,
																	 usfixn32* __restrict__ signVector);

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
																	 usfixn64* __restrict__ lVector,
																	 usfixn64* __restrict__ hVector,
																	 usfixn64* __restrict__ cVector,
																	 const usfixn32* __restrict__ signVector);

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step3 (usfixn64* __restrict__ xs,
																	 const usfixn64* __restrict__ lVector,
																	 const usfixn64 * __restrict__ hVector,
																	 const usfixn64* __restrict__ cVector,
																	 const usfixn32* __restrict__ signVector,
																	 const usfixn64 * __restrict__ parameters);

/************************************************/
__global__ void
kernel_p3_lhc (usfixn64 * xs, usfixn64* ys, usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_bigMult_plain (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters);

/**********************************************/
__global__ void
kernel_p3_bigMult_plain_2 (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters);
/**********************************************/
__global__ void
kernel_p3_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters);

/**********************************************/

/**********************************************/
//xs = xs*ys
//most efficient big multiplication
//not verified
//uses both shared memory, lmem, register, global mem
//try to fix negate and bigAddPlain to have static indexing
__global__ void
kernel_p3_newMult14 (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
										 usfixn64 * __restrict__ parameters, usfixn64* lVector,
										 usfixn64 *hVector, usfixn64* cVector);
/**********************************************/
__global__ void
kernel_p3_newMult14_incStep (usfixn64 * __restrict__ parameters);
/**********************************************/
__global__ void
kernel_p3_newMult14_join2 (usfixn64* __restrict__ xs,
													 usfixn64 * __restrict__ parameters,
													 usfixn64* __restrict__ lVector,
													 usfixn64 * __restrict__ hVector,
													 usfixn64* __restrict__ cVector);
/**********************************************/
__global__ void
kernel_p3_mulLong_revised (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters);
/**********************************************/
__global__ void
kernel_p3_mulLong_plain (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters);
/**********************************************/
#endif
