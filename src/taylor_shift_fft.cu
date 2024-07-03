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



#include "taylor_shift_conf.h"
#include "taylor_shift_cpu.h"
#include "taylor_shift_kernel.h"
#include "taylor_shift.h"
#include "taylor_shift_fft.h"
#include "inlines.h"


// store in Polynomial_shift the coeffs in fft we want to use
__global__ void mult_adjust_GPU(sfixn *Polynomial_shift, sfixn *fft, sfixn n, sfixn local_n, sfixn winv, sfixn p, double pinv)
{

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

  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn B = 2*local_n;

  if (i < n)
    {
      sfixn part = i / B;
      sfixn pos  = i % B;
      sfixn bool1 = (sfixn) (pos != 0);

      Polynomial_shift[i] = bool1 * mul_mod(winv, fft[2*B*part + pos-1], p, pinv);
    }
}


// transfer at each step the polynomials which need to be multiplicated
__global__ void transfert_array_fft_GPU(sfixn *fft, sfixn *Mgpu, sfixn n, sfixn local_n)
{
  sfixn i = blockIdx.x * blockDim.x + threadIdx.x;
  sfixn part, pos, bool1, bool2;
  part = i / local_n;
  pos  = i % local_n;
  bool2 = part % 2;
  bool1 = 1 - bool2;

  if (i<2*n)
    fft[i] = bool1 * Mgpu[(part/2)*local_n + pos];
}
