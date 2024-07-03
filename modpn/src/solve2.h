/* This file is part of the MODPN library

    MODPN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MODPN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Copyright, Marc Moreno Maza <moreno@csd.uwo.ca>
*/


#ifndef _SOLVE_2_H_
#define _SOLVE_2_H_



#include "IteratedResultant.h"
#include "RegularGcd.h"
#include "Types.h"
#include "Factorization.h"
#include "GCD.h"
#include "HGCD.h"
#include "MPMMTS.h"
#include "FMUL.h"


#include <time.h>



	//
/* CUDA related */
#if defined(_cuda_modp_)
#include "cumodp.h"
#endif

	//



sfixn testP(preFFTRep *f1, preFFTRep *f2, sfixn p);

regular_chain2 *EX_RegularChain2_Init(preFFTRep *f1, preFFTRep *f2);

regular_chain2 *EX_RegularChain2_Copy_Init(preFFTRep *f1, preFFTRep *f2);

void EX_RegularChain2_Free(void *element);

void EX_RegularChain2_Print(void *element);
void FILE_RegularChain2_Print(void *element, FILE *file);

LinkedQueue *modular_generic_solve2(preFFTRep *F1, preFFTRep *F2, 
    preFFTRep *g, preFFTRep *h, SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *modular_solve2_select(sfixn method, preFFTRep *F1, preFFTRep *F2,
    preFFTRep *g, MONTP_OPT2_AS_GENE *pPtr);

void irregularRCinsertion(computingChoince method, preFFTRep *F1, preFFTRep *F2, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc);

void ModularSolve2(computingChoince method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc);

void ModularGenericSolve2(computingChoince method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, polyBothInHostDev *h, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc);

void ModularGenericSolve2ZeroResultant(computingChoince method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, polyBothInHostDev *h, polyBothInHostDev *S0, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc);


void plainMultiplication( sfixn *a, sfixn db, sfixn *b, sfixn dc, sfixn *c, sfixn p);


sfixn* gcdCUDA(sfixn *dG, sfixn *A, sfixn dA, sfixn *B, sfixn dB, MONTP_OPT2_AS_GENE *pPtr);
sfixn* divCUDA(sfixn *dG, sfixn dA, sfixn *A, sfixn dB, sfixn *B, MONTP_OPT2_AS_GENE *pPtr);

sfixn* gcdCUDAStraight(sfixn *dG, sfixn *A, sfixn dA, sfixn *B, sfixn dB, MONTP_OPT2_AS_GENE *pPtr);
sfixn* divCUDAStraight(sfixn *dG, sfixn dA, sfixn *A, sfixn dB, sfixn *B, MONTP_OPT2_AS_GENE *pPtr);



sfixn is_divided_byCUDA(sfixn dF, sfixn *F, sfixn dT, sfixn *T, MONTP_OPT2_AS_GENE *pPtr);

sfixn isDividedWithCudaOption(computingChoince method, sfixn dF, sfixn *F, sfixn dT, sfixn *T, MONTP_OPT2_AS_GENE *pPtr);


void contentbivariate(computingChoince method,  polyBothInHostDev *G,  polyBothInHostDev *F1,  MONTP_OPT2_AS_GENE *pPtr);



sfixn KroneckerSubs( sfixn *A, polyBothInHostDev *B, sfixn Inc);


preFFTRep *invKroneck( sfixn *A, sfixn dA, sfixn e);






void bivarDivComponentwise(computingChoince method, polyBothInHostDev *A,  polyBothInHostDev *B, polyBothInHostDev *C, MONTP_OPT2_AS_GENE *pPtr);



void bivarMulComponentwise( polyBothInHostDev *A,  polyBothInHostDev *B, polyBothInHostDev *C, MONTP_OPT2_AS_GENE *pPtr);





LinkedQueue *modular_solve2(preFFTRep *F1, preFFTRep *F2, preFFTRep *g, 
    MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *EX_ModularSolve2(preFFTRep *F1, preFFTRep *F2, sfixn p, computingChoince method);

sfixn totalSizeRC(LinkedQueue *queue);
void makecompactRC(sfixn totalSize, sfixn *compactRC, LinkedQueue *queue);


#endif /* solve2.h */
