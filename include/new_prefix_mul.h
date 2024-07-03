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



#ifndef _PREFIX_MUL_H_
#define _PREFIX_MUL_H_

#include "new_desbound_conf.h"

/** Initialize an array for factorial sequence
 * @A: Arrary storing 1 2 .. N .. 1 .. N
 * @N: Size of n!
 * @pnum: Number of primes
 */
__global__ void initArray(sfixn* A, int N, const int* m, int pnum);

/** Initialize an array for powers, with same value **/
__global__ void initArray(sfixn* A, int N, const int* m, int pnum, sfixn base);


/** Compute prefix multiplication in each thread block
 * @A: Array {a_1, a_2, ..., a_B, a_{B+1}, ..., a_n}, and store the result {a_1, a_1*a_2, .., a_1*...a_B, a_{B+1}, ..., a_{x*B+1}*..*a_n}
 * @N: Size of array A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Prime
 */
__global__ void inclusivePreMul(sfixn* A, int N, int TBS, const int m);


/** Update array
 * @A: Updated array
 * @B: Array {b_1, ..., b_{N/2*TBS}}, which is the last element of each thread block in A from the previous step
 * @N: Size of array A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Prime
 */
__global__ void exclusivePreMul(sfixn* A, sfixn* B, int N, int TBS, const int m);


/** Colloect the last element of each thread block of previous step
 * @lasts: Last elements 
 * @A: Original polynomial
 * @nb: # of thread blocks in previous step
 * @TBS: Number of coefficients dealt by a thread block during previous step
 */
__global__ void collectLast(sfixn* lasts, sfixn* A, int nb, int TBS);


/** Compute the prefix multiplication
 * @A: Array {a_1, a_2, ..., a_n}, and store the result {a_1, a_1*a_2, .., a_1*..*a_n}
 * @N: Size of array A
 * @m: Primes
 * @pnum: Number of primes
 */
void prefixMul(sfixn* A, int N, const int* m, int pnum);

#endif // _PREFIX_MUL_H_
