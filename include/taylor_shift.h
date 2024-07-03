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



# ifndef _TAYLOR_SHIFT_H_
# define _TAYLOR_SHIFT_H_
#include "taylor_shift_conf.h"

/* Important to notice :

  n  : number of coefficients of the polynomial considered
       (size of the polynomial)
 n-1 : degree of the polynomial considered
  p  : prime number, it must be greater than n
*/


/* ------------- Taylor_shift procedure modulo Zp ------------- */

/* This procedure is called s times by the CPU with s different prime
   numbers. It gives the CRT representation of the polynomial
   2^{k*n}*P( (x+a)/2^k ) modulo p (called PKC), combination of a
   Taylor_shift by a of a input polynomial P and multiplications for
   each coefficient by a power of 2^k. In each call, one column of X
   will be filled modulo p.

   PARAMETERS :
   ------------

     n : number of coefficients of the polynomial considered
         (size of the polynomial)
     e : e = log2(n) so n = 2^e
     p : prime number
     pinv : invert of p
     X : array of size n*s containing CRT and MRR representation of
         the output polynomial
     step : indicates what column of X we will fill
     s : number of prime numbers needed
     Polynomial_device : input polynomial
     a : integer for the Taylor_shift by a
     k : power of 2 considered
     Factorial_device : array containing the sequence (i!)_{i in [0;n]}
     Powers_a : array containing in this code :
                @ firstly : the sequence (a^i)_{i in [0;n]} mod p
                @ secondly : the sequence ((2^k)^i)_{i in [0;n]} mod p
     Monomial_shift_device : array containing the sequence of the coeffi-
                             cients of ((x+1)^{2^i})_{i in [0;e-1]} mod p
     Mgpu : array used in several intermediate computations mod p
     Polynomial_shift_device : array in 2 dimensions containing pairwise
                               intermediate computations to obtain the
                               Taylor_shift by 1 of a polynomial mod p
     fft_device : array used for Fast Fourier Transform mod p          */

void taylor_shift_Zp(sfixn n, sfixn e, sfixn p, double pinv, sfixn *X, sfixn step, sfixn s, sfixn *Polynomial, sfixn a, sfixn k, sfixn *Factorial_device, sfixn *Powers_a, sfixn *Monomial_shift_device, sfixn *Mgpu, sfixn **Polynomial_shift_device, sfixn *fft_device);

#endif // _TAYLOR_SHIFT_H_
