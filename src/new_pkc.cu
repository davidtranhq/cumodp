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



#include "new_pkc.h"

/** Divide by 2^k
 * @polynominal: Coefficients in CRT manner
 * @powers: Powers of 2^k
 * @N: Size of polynominal
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void DivPow(sfixn* polynominal, sfixn* powers, int N, const int* m, int pnum) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        for (int i = 0; i < pnum; i++) {
                if (idx < N-1) {
                        polynominal[i*N+idx+1] = quo_mod(polynominal[i*N+idx+1], powers[i*N+idx], m[i]);
                }
        }
}


/** Multiply by 2^{k*n}
 * @polynominal: Coefficients in CRT manner
 * @powers: Powers of 2^k
 * @N: Size of polynominal
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void MulPow(sfixn* polynominal, sfixn* powers, int N, const int* m, int pnum) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        for (int i = 0; i < pnum; i++) {
                 if (idx < N) {
                        polynominal[i*N+idx] = mul_mod(powers[i*N+N-2], polynominal[i*N+idx], m[i]);
                }
        }
}
