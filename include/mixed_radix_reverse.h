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



typedef unsigned int usfixn32;
typedef unsigned long long int usfixn64;

__global__ void kernel_crt_multiplications_plain(usfixn32 *vs, usfixn32 * mult2D, usfixn32 * result);
/****************************************************/
__global__ void kernel_crt_multiplications(usfixn32 * __restrict__ vs, const usfixn32 * __restrict__ mult2D, usfixn32 * result, usfixn32* parameters);
/****************************************************/
__global__ void kernel_crt_additions(usfixn32 * __restrict__ vs, const usfixn32 * __restrict__ mult2D, usfixn32 * result, usfixn32 * __restrict__ parameters);
/****************************************************/
__global__ void kernel_crt_additions_plain(usfixn32 *vs, usfixn32 * mult2D, usfixn32 * result);
/****************************************************/
__global__ void kernel_crt_multiplications_v1(usfixn32 * __restrict__ vs, const usfixn32 * __restrict__ mult2D,
	usfixn32 * result, usfixn32* parameters);
/****************************************************/
__global__ void
kernel_crt_multiplications_32words_v1 (usfixn32 * __restrict__ vs,
	const usfixn32 * __restrict__ mult2D,
	usfixn32 * result, usfixn32* parameters);
/****************************************************/
