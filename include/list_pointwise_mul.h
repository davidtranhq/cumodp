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



#ifndef _LIST_POINTWISE_MUL_H
#define _LIST_POINTWISE_MUL_H
#include <stdio.h>
#include <stdlib.h>


#include <stdio.h>
#include <stdlib.h>
#include "inlines.h"
#include "defines.h"
#include "types.h"
#include "printing.h"
#include "cudautils.h"
#include "rdr_poly.h"


#define N_TREAD 512 
#define PROCESS_PER_TRD 16
__global__ void list_pointwise_mul(sfixn *L, sfixn ln, sfixn p, double pinv, sfixn length_layer);
__global__ void list_expand_to_fft(sfixn *M, sfixn *L, sfixn length_poly, sfixn ln, sfixn start_offset, sfixn length_layer);
__global__ void list_truncate_for_invfft(sfixn *L, sfixn *K, sfixn ln);
__global__ void list_shrink_after_invfft(sfixn *L, sfixn *S, sfixn ln, sfixn length_result, sfixn length_layer);


__host__ void CPU_list_pointwise_mul(sfixn *L, sfixn ln, sfixn num_poly, sfixn p, double pinv);
__host__ void CPU_list_expand_to_fft(sfixn *M, sfixn *L, sfixn length_M, sfixn length_poly, sfixn ln, sfixn start_offset);
__host__ void CPU_list_truncate_for_invfft(sfixn *L, sfixn *K, sfixn ln, sfixn num_poly);
__host__ void CPU_list_shrink_after_invfft(sfixn *L, sfixn *S, sfixn ln, sfixn length_result, sfixn num_poly);
#endif

