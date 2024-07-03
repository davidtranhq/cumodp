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



#ifndef _DESBOUND_H_
#define _DESBOUND_H_

#include "new_desbound_conf.h"
#include "new_prefix_mul.h"
#include "new_taylor_shift.h"
#include "new_pkc.h"
#include "new_mixed_radix_rep.h"


/** Compute P_{k, c} = 2^{k*n} P(\frac{x+c}{2^{k}})
 * @h_A: Input polynomial storing as ( f_0 .. f_d ) mod prime[0] .. ( f_0 .. f_d ) mod prime[np-1]
 * @h_B: Ouput polynominal after computing PKC
 * @N: Number of coefficients of input polynominal
 * @pnum: Number of primes
 * @k: Input k
 * @c: Input c
 */
void pkc_cpu(sfixn* h_A, sfixn* h_B, int N, const int* m, int pnum, int k, sfixn c);


/** Represent coefficients in MRR manner, compute the number of signs changed, and the largest positive number and the smallest negative number
 * @h_A: Input polynominal, storing as ( f_0 .. f_d ) mod prime[0] .. ( f_0 .. f_d ) mod prime[np-1]
 * @h_B: Ouput polynominal after computing PKC
 * @N: Number of coefficients of input polynominal
 * @pnum: Number of primes
 * @sc: Number of signs changed
 * @mrrPosNorm: Largest positive number represented in MRR manner
 * @mrrNegNorm: Smallest negative number represented in MRR manner
 */
void mrr_cpu(sfixn* h_A, sfixn* h_B, int N, const int* m, int pnum, int* sc, sfixn* mrrPosNorm, sfixn* mrrNegNorm);


#endif // _DESBOUND_H_
