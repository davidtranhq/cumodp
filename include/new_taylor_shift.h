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



#ifndef _TAYLOR_SHIFT_H_
#define _TAYLOR_SHIFT_H_

#include "new_desbound_conf.h"
#include "new_prefix_mul.h"

/** Compute factorial sequence (x+1)^n, with each coefficient n!/[k! * (n-k)!]
 * @A: Array storing (x+c)^{2^i}, 0 \leq i < \log[2](N) continuously
 * @facseq: Factorial sequence storing 1!, ... (\frac{N}{2})!
 * @powersc: Powers of c
 * @N: Size of original polynominal
 * @m: Primes
 * @pnum: Number of primes
 ** Regard to (x+c)^{i}, the leading coefficient is always 1, such that it is not stored
 */
__global__ void binomialCoef(sfixn* A, sfixn* facseq, sfixn* powersc, int N, const int* m, int pnum);


/** Compute each multiplication of polynomials in each thread block
 * @Mgpu1: Polynomials
 * @Mgpu2: Polynomials after multiplication
 * @length_poly: Size of each polynomial
 * @poly_on_layer: # of polynomials doing multiplication
 * @threadsForAmul: = # of threads per thread block / mulInThreadBlock
 * @mulInThreadBlock: # of multiplications computed in a thread block
 * @p: Prime number
 */
__global__ void modifiedListPlainMulGpu(sfixn *Mgpu1, sfixn *Mgpu2, int length_poly, int poly_on_layer, int threadsForAmul, int mulInThreadBlock, int p);


/** Compuse input polynomial for multiplication in the sub-tree of Taylor Shift
 * @mul: Output polynomial
 * @poly: Original polynomial
 * @N: Size of polynomial
 * @facseq: Factorial sequence
 * @base: Size of part of the polynomial that muplies (c+x)^base
 */
__global__ void composeMulPoly(sfixn* mul, sfixn* poly, int N, sfixn* facseq, int base);


/** Addition of two polynomials
 * @A: Polynomial
 * @B: Polynomial
 * @N: Size of polynomial
 */
__global__ void addPoly(sfixn* A, sfixn* B, int N);


/** Shift polynomial by c
 * @poly: Input polynomial
 * @N: Size of polynomial
 * @c: Shift by c
 * @m: Primes
 * @pnum: Number of primes
 */
void TayShi(sfixn* poly, int N, sfixn c, const int* m, int pnum);

#endif // _TAYLOR_SHIFT_H_
