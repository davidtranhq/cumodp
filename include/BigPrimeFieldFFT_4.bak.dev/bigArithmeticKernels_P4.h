// for 16-word prime

//kernel_operation_
#ifndef BIG_ARITHMETIC_KERNELS_H_P4
#define BIG_ARITHMETIC_KERNELS_H_P4


#include "BigPrimeFieldFFT_4/bigPrimeField_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_aux_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_addition_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_multiplication_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_cyclicShift_P4.h"
#include "BigPrimeFieldFFT_4/bigArithmetic_fft_P4.h"
#include "BigPrimeFieldFFT_4/bigFFT_reduction_P4.h"
/**********************************************/
//
//__global__ void
//kernel_transposition (usfixn64* xs, usfixn64 * ys, usfixn64 nPermutationBlocks,
//											usfixn64 blockStride, usfixn64* parameters);
///**********************************************/
//
//void
//host_cpy_d2d (usfixn64 * dest, usfixn64 * src, usfixn64 count, data* fd);

/*************************************************************
 * parameters[0] = operation;
 * parameters[1] = iterations;
 * parameters[2] = paddingMethod;
 * parameters[3] = dynamicMemSize;
 * parameters[4] = strideSize;
 * ?
 *************************************************************/
#endif /* BIG_ARITHMETIC_GLOBAL_H_ */
