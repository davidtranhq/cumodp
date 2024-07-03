#ifndef _LIST_FAST_DIVISION_H_
#define _LIST_FAST_DIVISION_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "subproduct_tree.h"
#include "list_naive_poly_mul.h"
#include "list_pointwise_mul.h"
#include "list_inv_stockham.h"
#include "list_fft_poly_mul.h"
#include "power_inversion.h"

#include "list_poly_rev.h"
__device__ __global__ void poly_minus(sfixn *R, sfixn *G, sfixn length, sfixn p);

void list_fast_division(sfixn *A, sfixn length_poly_A, sfixn *B, sfixn length_poly_B, sfixn num_poly, sfixn p);
#endif
