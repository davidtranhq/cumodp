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



#ifndef _MIXED_RADIX_REP_H_
#define _MIXED_RADIX_REP_H_

#include "new_desbound_conf.h"

/** Transpose polynominal from r*c to c*r
 * @transpose: Transposed polynominal
 * @polynominal: Original polynominal
 * @r: Original rows
 * @c: Original columns
 */
__global__ void transpose(sfixn* transpose, sfixn* polynominal, int r, int c);


/** Compute elements needed in matrices
 ** m[i,j], where m^{-1}[i] = m[i,j] mod m[j] with i < j
 ** n[i,j] = m[j] - m[i,j]
 * @matrices: Array storing m[i,j] and n[i,j] continuously
 * @pnum: Number of primes
 * @m: Primes
 */
__global__ void initMRR(sfixn* matrices, int pnum, const int* m);


/** Convert coefficients into Mixed Radix Representaion
 ** Compute ( .. ((X A_{0}) A_{1}) .. A_{s-1})
 * @mrr: Coefficients in MRR manner
 * @N: Size of original polynominal
 * @matrices: Array storing m[i,j] and n[i,j] continuously
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void MRR(sfixn* mrr, int N, sfixn* matrices, const int* m, int pnum);


/** Compute the sum of signs changed within each thread block
 * @A: Array storing the sign changed
 * @N: Size of A
 */
__global__ void signsSum(int* A, int N, int TBS);


/** Check the sign of a MRR coefficient
 ** If a = 0, return 0
 ** If a < 0, return -1
 ** Otherwise, return 1
 */
__device__ __host__ __inline__
int checkSign(sfixn* a, const int* m, int pnum);


/** Cout number of sign changed
 * @signs: If a coefficient's sign is changed, store 1, otherwise store 0;
 *	Then perform prefix sum on it, and the last element is the number of signs changed
 * @A: Polynominal
 * @N: Size of polynominal
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void signChange(int* signs, sfixn* A, int N, const int* m, int pnum);


/** Compare a and b
 ** If a > b, return 1
 ** If a < b, return -1
 ** Otherwise, return 0 
 */
__device__ int compareMRR(sfixn* a, sfixn* b, const int* m, int pnum);


/** Compute the largest and smallest number among the coefficients within a thread block
 * @index: [0], store the index with smallest coefficient;
	 [N-1], store the index with largest coefficient
 * @A: Coefficients
 * @N: Size of A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void abslargestMRR(int* index, sfixn* A, int N, int TBS, const int* m, int pnum);

#endif // _MIXED_RADIX_REP_H_
