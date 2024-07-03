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


#ifndef _FASTPOLYEVALUATION_H
#define _FASTPOLYEVALUATION_H

#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"

#include "cudautils.h"
#include "types.h"
#include "list_inv_stockham.h"
#include "subproduct_tree.h"
#include "list_stockham.h"

#define MAX_LEVEL 27

/*
 As our implementation is based on an adaptive algorithm
 described in the PhD thesis of Sardar Haque (Computer Science, UWO, 2013),
 we need to have different number of threads in a thread block for different purposesl.
 All Tmax, Tinv, Tdiv are the number of threads in a thread block for different scenario.
 For exmaple Tmax is the choice by default and it is used for all oridinary purposes.
 Tinv is used when we are computing the inverse from scratch using Newton iteration.
 */
const sfixn Tmax = 512;
const sfixn Tinv = 16;
const sfixn Tdiv = 128;
const sfixn Tmul = 512;

/*
 In some coefficient arithmetic like list of plain division,
 a thread is responsible for W of coefficient operations.
 */
const sfixn W = 4;

/*
 plainMulLimit is the level above which we do fast arithmetic or FFT based arithmetic.
 According to the Chapter 10 of the PhD thesis of Sardar Haque (Computer Science, UWO, 2013),
 it is H.
 */
const sfixn plainMulLimit = 8;

struct PolyEvalSteps
{
	sfixn *Ml[MAX_LEVEL], *Mr[MAX_LEVEL];
	sfixn *InvMl[MAX_LEVEL], *InvMr[MAX_LEVEL];
	sfixn *fL[MAX_LEVEL], *fR[MAX_LEVEL];
	float timeSubTree, timeSinvTree, timeEva;
};

void
subProductTree (sfixn k, sfixn p);
void
subInvTree (sfixn k, sfixn p);
void
fastEvaluation (sfixn *F, sfixn k, sfixn p, sfixn flag);
void
fastEvaluationLow (sfixn *F, sfixn k, sfixn p, sfixn flag);

struct PolyEvalSteps
fastEvaluation (sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag);
struct PolyEvalSteps
onlySubproductTree (sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag);
struct PolyEvalSteps
onlyTrees (sfixn *points, sfixn k, sfixn p, sfixn flag);

__global__ void
copyMgpu (sfixn *dest, sfixn *source, sfixn length_poly);

__global__ void
allZero (sfixn *X, sfixn n);

__global__ void
pointAdd2 (sfixn *dest, sfixn *source, sfixn l, sfixn n, sfixn p);

__global__ void
zeroInbetween (sfixn *X, sfixn *Y, sfixn n, sfixn l);

__global__ void
pointMul (sfixn *dest, sfixn *source, sfixn n, sfixn p);

__global__ void
scalarMul (sfixn *A, sfixn ninv, sfixn L, sfixn p);

__global__ void
pointAdd (sfixn *dest, sfixn *source, sfixn n, sfixn p);

__global__ void
listPlainMulGpu (sfixn *Mgpu1, sfixn *Mgpu2, sfixn length_poly,
								 sfixn poly_on_layer, sfixn threadsForAmul,
								 sfixn mulInThreadBlock, sfixn p);

__global__ void
listPolyinv (sfixn *Mgpu, sfixn *invMgpu, sfixn poly_on_layer, sfixn prime);

__global__ void
listReversePoly (sfixn *revMgpu, sfixn *Mgpu, sfixn length_poly,
								 sfixn poly_on_layer);

__global__ void
listCpLdZeroPoly (sfixn *B, sfixn *A, sfixn length_poly, sfixn poly_on_layer);

__global__ void
allNeg (sfixn *X, sfixn n, sfixn p);

__global__ void
listPolyDegInc (sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly,
								sfixn poly_on_layer, sfixn newLength);

__global__ void
listCpUpperCuda (sfixn *dest, sfixn *source, sfixn n, sfixn l);

__global__ void
listCpLowerCuda (sfixn *dest, sfixn *source, sfixn n, sfixn l);

__global__ void
list2wayCp (sfixn *dest, sfixn *source, sfixn l, sfixn totalLength, sfixn p);

__global__ void
leavesSubproductTree (sfixn *M1gpu, sfixn *Mgpu, sfixn numPoints,
											sfixn rightSubtree);
__global__ void
EvalistReversePolyDec (sfixn *Mgpu, sfixn *revMgpu, sfixn length_poly,
											 sfixn poly_on_layer);
__global__ void
EvalistPolyDegInc (sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly,
									 sfixn poly_on_layer, sfixn newLength);

__global__ void
EvalistZeroBetweenRev (sfixn *revMzero, sfixn *M, sfixn length_poly,
											 sfixn poly_on_layer);

__global__ void
subtractLowerSingle (sfixn *U, sfixn *F, sfixn *G, sfixn n, sfixn p);

__global__ void
pointwiseMulHalf (sfixn *A, sfixn *B, sfixn length_poly, sfixn poly_on_layer,
									sfixn p);
__global__ void
subtractLower (sfixn *R, sfixn *A, sfixn *BQ, sfixn length_poly,
							 sfixn poly_on_layer, sfixn p);

__global__ void
listPlainCudaDiv (sfixn *M, sfixn *F, sfixn start, sfixn length,
									sfixn threadsPerDiv, sfixn DivPerBlock, sfixn polyNum,
									sfixn P);
__global__ void
PlainCudaDiv (sfixn *M, sfixn *F, sfixn m, sfixn f, sfixn p);

#endif

