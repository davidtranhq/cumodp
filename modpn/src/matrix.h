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


/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __matrix_h
#define __matrix_h 

#include "Types.h"
#include "UniHensel.h"

#include <stdlib.h>            /* rand() and srand() functions               */
#include <assert.h>            /* assert()                                   */



void multi_mat_transpose (sfixn N, sfixn n, sfixn dm, sfixn * dims, sfixn * data);

void mat_transpose (sfixn, sfixn, sfixn * );

POLYVECTOR_SLG *randomPolyVec_SLG(int M, int GN, sfixn p, int N);

void freeVec_SLG(POLYVECTOR_SLG *vec_slg);

void  printJMatrix(POLYMATRIX * mat);

void freeJMatrix(POLYMATRIX * mat);

POLYMATRIX *createJMatrix(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);


POLYMATRIX *initJMatrix(sfixn N, int m, int n, TriSet *ts);


POLYMATRIX *mulJMatrix(sfixn N, POLYMATRIX * mat1, POLYMATRIX * mat2, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
randomJMatrix(POLYMATRIX * mat, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
INVJMatrix(sfixn N, POLYMATRIX * mat, TriSet * ts, TriRevInvSet * tris,  MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
createJMatrix_ForLifting(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
SLG2PolyVecotr(POLYVECTOR_SLG *polyVec_SLG, TriSet *ts, TriRevInvSet * tris, sfixn N, MONTP_OPT2_AS_GENE * pPtr);



POLYMATRIX *
createJMatrix_PolyForlifting(TriSet *ts, TriRevInvSet * tris,  sfixn N, MONTP_OPT2_AS_GENE * pPtr);


POLYMATRIX *
scalarMulJMatrix_1(sfixn r, POLYMATRIX * mat, MONTP_OPT2_AS_GENE * pPtr);

POLYMATRIX *
subJMatrix_1(sfixn N, POLYMATRIX * mat1, POLYMATRIX * mat2, MONTP_OPT2_AS_GENE * pPtr) ;


TriSet *
AddHs2TriSet(sfixn N, TriSet *ts, POLYMATRIX *Hs, MONTP_OPT2_AS_GENE * pPtr);



POLYMATRIX *
increaseMatrix_ForLifting(sfixn N, POLYMATRIX *mat, TriSet *ts);


void *BlockDecompose(void * PTR);



POLYVECTOR_SLG * example_1_PolyVec_SLG();


#endif

