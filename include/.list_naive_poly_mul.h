#ifndef _LIST_NAIVE_POLY_MUL_H_
#define _LIST_NAIVE_POLY_MUL_H_

#include <stdio.h>
#include <stdlib.h>

#include "inlines.h"
#include "defines.h"
#include "types.h"
#include "printing.h"
#include "cudautils.h"
#include "rdr_poly.h"


__device__ __global__ void list_poly_mul_ker(sfixn *M, sfixn length_poly, sfixn length_result, sfixn start_offset, sfixn result_offset, sfixn p, double pinv);
__device__ __global__ void list_poly_mul_ker_higher(sfixn *M, sfixn length_poly, sfixn length_result, sfixn start_offset, sfixn result_offset, sfixn p, double pinv);

#endif
