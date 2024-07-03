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



#ifndef _FASTPOLYINTERPOLATION_H
#define _FASTPOLYINTERPOLATION_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "fastPolyEvaluation.h"

struct PolyInterpolateSteps
{
	sfixn *lefts[MAX_LEVEL], *rights[MAX_LEVEL], *c, *I;
};

void linearCombination(sfixn k, sfixn p);
struct PolyInterpolateSteps interpolate(sfixn *c, sfixn k, sfixn p, sfixn flag);
struct PolyInterpolateSteps fastInterpolate(sfixn *v, sfixn *u, sfixn k, sfixn p, sfixn flag);
void computeI(sfixn *I, sfixn k, sfixn p, sfixn flag);
void computeIForSmallK(sfixn *I, sfixn k, sfixn p);
void computeIForBigK(sfixn *I, sfixn k, sfixn p);
void computeCs(sfixn *c, sfixn *v, sfixn k, sfixn p);
void generateM(sfixn *m, sfixn k, sfixn p);
void generateMForSmallK(sfixn *m, sfixn k, sfixn p);
void generateMForBigK(sfixn *m, sfixn k, sfixn p);
void fastEvaluate(sfixn *F, sfixn k, sfixn p, sfixn flag);
void fastEvaluationWithExistingTreesLow(sfixn *F, sfixn k, sfixn p, sfixn flag);
void fastEvaluationWithExistingTrees(sfixn *F, sfixn k, sfixn p, sfixn flag);

__global__ void addToLeading(sfixn *R, sfixn *portion1, sfixn *portion2, sfixn n, sfixn p);
__global__ void derivate(sfixn *f, sfixn n, sfixn p);
__global__ void computeC(sfixn *c, sfixn* v, sfixn *sinv, sfixn n, sfixn p);
__global__ void plainCrossMultiplications(sfixn *R, sfixn *M, sfixn *I, sfixn polyCount, sfixn threadsForAnOperation, sfixn operationsInABlock, sfixn p);
__global__ void change_elem(int *arr, int idx, int val);

#endif
