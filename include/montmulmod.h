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



///////////////////////////////////////////////////////////////////////////////
//   The montgomery reduction based modular multiplication
///////////////////////////////////////////////////////////////////////////////

#ifndef _MONT_MUL_MOD_
#define _MONT_MUL_MOD_

#include "inlines.h"
#include "defines.h"

/**
 * compute a = a * b quo 2^k and b = a * b mod 2^k 
 *
 * h, l stand for a, b
 * 
 * wsz = wordsize
 *
 * The requirements are 
 * (0) 0 <= a, b < 2^wsz
 * (1) 0 < k <= wsz
 * (2) a * b < 2^(wsz + k)
 *
 * The requirement (2) is needed, otherwise the quotient 
 * cannot fit into an sfixn.
 */
__device__ __host__ inline void mul_hi_lo(sfixn *h, sfixn *l, sfixn k) {
    // compue pd = [high32, low32] of a * b
    ulongfixnum pd = (ulongfixnum)(usfixn)(*h) * (ulongfixnum)(usfixn)(*l);
    // h gets pd quo 2^k
    *h = (sfixn)(pd >> k);
    // rem gets the last low32
    usfixn rem = (usfixn) pd;
    // reset the first 32 - k bits
    rem <<= (BASE - k);
    rem >>= (BASE - k);
    // convert back to a signed number 
    *l = (sfixn)rem;
}

/**
 * Given a and b such that 0 <= a, b < p, it computes
 *
 * a * b  
 * ------  mod p
 *   R 
 *
 * fourier_reduction(a, b) = fourier_reduction(b, a)
 *
 */
__device__ __host__  __inline__
sfixn fourier_reduction(sfixn a, sfixn b, const fprime_t * const fpp) {
    // q1 = a * b quo R --> a
    // r1 = a * b mod R --> b
    mul_hi_lo(&a, &b, fpp->rpow);
    
    // q2 = r1 * x quo R --> b
    // r2 = r1 * x mod R --> x
    sfixn p = fpp->p, x = p - 1;
    mul_hi_lo(&b, &x, fpp->rpow);
    
    // q1 - q2 --> a
    a -= b;

    // q3 = c * 2^npow * r2 quo R 
    // q3 = c * r2 / 2^(rpow - npow) --> x
    x *= fpp->c;
    x = (usfixn)x >> fpp->rnpow;

    // q1 - q2 + q3 --> a
    a += x;

    // a is in range (-p, 2*p)
    // correct the error as follows
    a += (a >> BASE_1) & p;
    a -= p;
    a += (a >> BASE_1) & p;
     
    return a;
}

/* convert a to a * R mod p */
/* f stands for Fourier or forward */
__device__ __host__ inline 
sfixn frep(sfixn a, const fprime_t * const fpp) {
    return mul_mod(a, fpp->r, fpp->p);
}

/**
 * the inverse of frep, a * R mod p --> a mod p;
 */
__device__ __host__ inline 
sfixn inv_frep(sfixn a, const fprime_t * const fpp) {
    return fourier_reduction(a, 1, fpp);
}

/**
 * compute a^e mod p.
 *
 * @a, integer in the Montgomery representation
 * @e, nonnegative integer
 *
 * The result is also in the Montgomery representation.
 */
__device__ __host__ inline 
sfixn fpow_mod(sfixn a, sfixn e, const fprime_t * const fpp) {
    sfixn x = fpp->r;
    if (e <= 0) return x;
    while (e) {
        if (e & 1) x = fourier_reduction(x, a, fpp);
        a = fourier_reduction(a, a, fpp);
        e >>= 1;
    }
    return x;
}

#endif // END OF FILE
