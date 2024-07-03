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


// for 16-word prime
#ifndef BIG_ARITHMETIC_FFT_H_
#define BIG_ARITHMETIC_FFT_H_

#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"

/************************************************/
//FFT base 32 to be called from host
__global__ void
kernel_p4_base_fft32 (usfixn64 *xm, usfixn64 * parameters);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r1_DFT (usfixn64 * xm,
														 const usfixn64 * __restrict__ parameters);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r2_shift (
		usfixn64 * xm, const usfixn64 * __restrict__ parameters,
		const short * __restrict__ device_shNo_FFT_32_r2);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r2_DFT (usfixn64 * xm,
														 const usfixn64 * __restrict__ parameters);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r3_shift (
		usfixn64 * xm, const usfixn64 * __restrict__ parameters,
		const short * __restrict__ device_shNo_FFT_32_r3);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r3_DFT (usfixn64 * xm,
														 const usfixn64 * __restrict__ parameters);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r4_shift (
		usfixn64 * xm, const usfixn64 * __restrict__ parameters,
		const short * __restrict__ device_shNo_FFT_32_r4);

/***********************************************/
__global__ void
kernel_p4_base_fft32_r4_DFT (usfixn64 * xm,
														 const usfixn64 * __restrict__ parameters);
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_shift (
		usfixn64 * xm, const usfixn64 * __restrict__ parameters,
		const short * __restrict__ device_shNo_FFT_32_r5);
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_DFT (usfixn64 * xm,
														 const usfixn64 * __restrict__ parameters);
/***********************************************/
__global__ void
kernel_p4_base_fft32_r5_permutation (usfixn64 * xm,
																		 const usfixn64 * __restrict__ parameters,
																		 const short * __restrict__ shNoListGlobal);
/***********************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step1 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters);
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step2 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters);

/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21_v0 (
		const usfixn64 * __restrict__ xs, const usfixn64* __restrict__ powers_omega,
		const usfixn64 K, const usfixn64 l, usfixn64* __restrict__ lVector,
		usfixn64* __restrict__ hVector, usfixn64* __restrict__ cVector,
		usfixn32* __restrict__ signVector, const usfixn64* __restrict__ parameters);

/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21 (const usfixn64 * xs,
																						 const usfixn64* powers_omega,
																						 const usfixn64 K, const usfixn64 l,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step21_lessReg (
		const usfixn64 * xs, const usfixn64* powers_omega, const usfixn64 K,
		const usfixn64 l, usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
		usfixn32* signVector, usfixn64* parameters);
/**********************************************************/
/**********************************************************/
__global__ void
kernel_p4_twiddle_general_permutated_step22 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);
/**********************************************/

__global__ void
kernel_p4_twiddle_general_permutated_step23 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);

#endif
