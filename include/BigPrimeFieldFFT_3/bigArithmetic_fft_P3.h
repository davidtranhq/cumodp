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


#ifndef BIG_ARITHMETIC_FFT_H_
#define BIG_ARITHMETIC_FFT_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"

extern __shared__ usfixn64 sharedMem[];

/**********************************************/
//FFT base 16 to be called from host
__global__ void
kernel_p3_base_fft16 (usfixn64 *xm, usfixn64 * parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r1 (usfixn64 * xm,
													 const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r21 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r22 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r31 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r32 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r41 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r42 (usfixn64 * xm,
														const usfixn64 * __restrict__ parameters);

/**********************************************/
__global__ void
kernel_p3_base_fft16_3_r5 (usfixn64 * xm,
													 const usfixn64 * __restrict__ parameters);

/**********************************************/

/**********************************************/
//only works if N/2 threads are being used for N coefficients
//this is compatible with fft_base_256

/**********************************************/
__global__ void
kernel_p3_fft2_permutated (usfixn64 *xs, usfixn64 *ys, usfixn64 *parameters);

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step1 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters)

																						;

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step2 (usfixn64 * xs,
																						usfixn64* powers_omega, usfixn64 K,
																						usfixn64 l, usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step21 (usfixn64 * xs,
																						 usfixn64* powers_omega, usfixn64 K,
																						 usfixn64 l, usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step22 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_twiddle_general_permutated_step23 (usfixn64 * xs, usfixn64 K,
																						 usfixn64 l_input,
																						 usfixn64* lVector,
																						 usfixn64* hVector,
																						 usfixn64* cVector,
																						 usfixn32* signVector,
																						 usfixn64* parameters);

/**********************************************/
//performing fft-16 on plain data on a single thread
__global__ void
kernel_p3_fft16_plain_singleThread (usfixn64 *xs, usfixn64 *ys);

/**********************************************/
__global__ void
kernel_p3_fft256_plain_2 (usfixn64 *xs, usfixn64* parameters);

/**********************************************/
__global__ void
kernel_p3_fft16_plain_2 (usfixn64 *xs, usfixn64 *parameters);

/**********************************************/

#endif
