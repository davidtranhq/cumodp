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



#ifndef _TAYLOR_SHIFT_FFT_H_
#define _TAYLOR_SHIFT_FFT_H_

// store in Polynomial_shift the coeffs in fft we want to use
__global__ void mult_adjust_GPU(sfixn *Polynomial_shift, sfixn *fft, sfixn n, sfixn local_n, sfixn winv, sfixn p, double pinv);
  /*                   EXAMPLE


  --------------------------------------------------------

     ARRAY fft considered
       _______________ _______________ _______________ _______________
      |               |               |               |               |
      |  real coeffs  |    useless    |  real coeffs  |    useless    |
      |_______________|_______________|_______________|_______________|
           PART=0                       PART=2


     B = 2 * local_n =  size of a PART

  --------------------------------------------------------

     ARRAY Polynomial_shift_device considered
       _______________________________
      |               |               |
      |  real coeffs  |  real coeffs  |
      |_______________|_______________|
            part=0         part=1

     B = 2 * local_n = size of a part

                                                              */


// transfer at each step the polynomials which need to be multiplicated
__global__ void transfert_array_fft_GPU(sfixn *fft, sfixn *Mgpu, sfixn n, sfixn local_n);


#endif // _TAYLOR_SHIFT_FFT_H_
