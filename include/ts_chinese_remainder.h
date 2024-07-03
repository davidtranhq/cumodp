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



# ifndef _TS_CHINESE_REMAINDER_H_
# define _TS_CHINESE_REMAINDER_H_


// gives 'positions' in arrays D and L
// !!! not used any more in the code !!!
//__device__ void funcPOS(sfixn i, sfixn s, sfixn *out);


/* CONVERSION OF MODULAR NUMBERS WITH MIXED RADIX REPRESENTATION BY A MATRIX FORMULA

   We want to have b = (...((x A_1)A_2)...)A_{s-1} using A_k such that :

           [               |                                            ]
           [               |                                            ]
           [     I_{k-1}   |                   0                        ]
           [               |                                            ]
           [    ___________|_________________________________________   ]
           [               |                                            ]
           [               |  1   n_{k,k+1}  n_{k,k+2}  ...  n_{k,s}    ]
   A_k =   [               |  0   m_{k,k+1}      0      ...     0       ]
           [               |  0       0      m_{k,k+2}          :       ]
           [               |  :                  0   .          :       ]
           [       0       |  :                    .   .        :       ]
           [               |  :        ...           .   .      :       ]
           [               |  :                        .   .    :       ]
           [               |  0          ...             0   m_{k,s}    ]
           [               |                                            ]


    We don't really need to do matrix multiplication, but just to consider diagonal elements
    and the k-th lines of these matrices. We need to use vectors D and L.
*/


// procedure to compute the elements of vectors D and L
__global__ void createDL(sfixn *D, sfixn *L, sfixn* primes_GPU, sfixn s);


// procedure copying the elements of the Taylor shift by 1 in Z/pZ in X 
__global__ void copyX(sfixn *X, sfixn s, sfixn *Polynomial_shift, sfixn n, sfixn step);


// procedure computing (...(X A_{k0})...)A_{s-T}
__global__ void mixed_radix_phase1(sfixn* X, sfixn s, sfixn n, sfixn* primes_GPU, sfixn* D, sfixn* L);


// procedure computing (...(X A_{s-T+1})...)A_{s-1}
__global__ void mixed_radix_phase2(sfixn* X, sfixn s, sfixn n, sfixn *primes_GPU, sfixn* D, sfixn* L);


// procedure computing (...(X A_1)...)A_{s-1}
void recombine(sfixn* X, sfixn s, sfixn n, sfixn *primes_GPU);


// used for tests
void mixed_radix_phase_cpu(sfixn* X, sfixn s, sfixn n, sfixn* D, sfixn* L);


# endif // _TS_CHINESE_REMAINDER_H_
