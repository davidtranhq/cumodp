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


#ifndef BIG_ARITHMETIC_MULTIPLICATION_H_
#define BIG_ARITHMETIC_MULTIPLICATION_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"

//#include "BigPrimeFieldFFT_3/bigArithmetic_addition_P3.h"
//#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread_8lhc (usfixn64 * xs, usfixn64* ys,
																usfixn64* parameters)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_v0 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_step1 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
										usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
										usfixn64* lVectorSub, usfixn64* hVectorSub,
										usfixn64* cVectorSub)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_step2 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
										usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
										usfixn64* lVectorSub, usfixn64* hVectorSub,
										usfixn64* cVectorSub)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_step3 (usfixn64* __restrict__ xs, usfixn64* __restrict__ lVector,
										usfixn64 * __restrict__ hVector,
										usfixn64* __restrict__ cVector,
										usfixn64* __restrict__ lVectorSub,
										usfixn64 * __restrict__ hVectorSub,
										usfixn64* __restrict__ cVectorSub,
										usfixn64 * __restrict__ parameters)
;

/************************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
												 const usfixn64* __restrict__ ys, usfixn64* parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 usfixn32* __restrict__ signVector)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 const usfixn32* __restrict__ signVector)
;

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step3 (usfixn64* __restrict__ xs,
												 const usfixn64* __restrict__ lVector,
												 const usfixn64 * __restrict__ hVector,
												 const usfixn64* __restrict__ cVector,
												 const usfixn32* __restrict__ signVector,
												 const usfixn64 * __restrict__ parameters)
;

/************************************************/
__global__ void
kernel_p3_lhc (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
;

/**********************************************/
__global__ void
kernel_p3_bigMult_plain (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
;

/**********************************************/
__global__ void
kernel_p3_bigMult_plain_2 (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
;

/**********************************************/
__global__ void
kernel_p3_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters)
;


/**********************************************/
//xs = xs*ys
//most efficient big multiplication
//not verified
//uses both shared memory, lmem, register, global mem
//try to fix negate and bigAddPlain to have static indexing
__global__ void
kernel_p3_newMult14 (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
					 usfixn64 * __restrict__ parameters, usfixn64* lVector,
					 usfixn64 *hVector, usfixn64* cVector)
;

/**********************************************/
__global__ void
kernel_p3_newMult14_incStep (usfixn64 * __restrict__ parameters)
;

/**********************************************/
__global__ void
kernel_p3_newMult14_join2 (usfixn64* __restrict__ xs, usfixn64 * __restrict__ parameters,
								 usfixn64* __restrict__ lVector,
								 usfixn64 * __restrict__ hVector,
								 usfixn64* __restrict__ cVector)
;

/**********************************************/
__global__ void
kernel_p3_mulLong_revised (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters)
;

/**********************************************/
__global__ void
kernel_p3_mulLong_plain (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters);
/**********************************************/
#endif
