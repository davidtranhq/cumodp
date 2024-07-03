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
#ifndef __FINTERP_h
#define __FINTERP_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FMUL.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include <math.h>

/*
 * Subproduct tree related functions
 */

subProdTree* subProdTreeCre( sfixn itemNo, sfixn itemSz, sfixn *items, sfixn p);

subProdTree** createArrOfSPTrees(sfixn m, sfixn *bounds, sfixn *pts, sfixn p);

void subProdTreeFree(subProdTree *tree);

void printSubProdTree(subProdTree *tree);

void mulNodes(sfixn dr, sfixn *srcAddr, sfixn ds, sfixn *Addr1, sfixn *Addr2,
     sfixn nodeSz, MONTP_OPT2_AS_GENE *pPtr);

void freeArrOfSPTTrees(sfixn m, subProdTree **trees);

void printArrOfSPTTrees(sfixn m, subProdTree **trees);

void FastEvaluation_1(sfixn resultsz, sfixn *result, sfixn degf, sfixn *fPtr,
     subProdTree *tree, sfixn p);

void freePTS_TREE(PTS_TREE *pt);

void fastEvalMulti_test(sfixn M, sfixn *dims, sfixn *E, preFFTRep* poly,
     subProdTree**trees, MONTP_OPT2_AS_GENE *pPtr);

sfixn* FastEvaluation(sfixn n, sfixn degf, sfixn *fPtr, subProdTree *tree, sfixn p);

sfixn* SlowEvaluation(sfixn degf, sfixn *fPtr, sfixn nopts, sfixn *pts, sfixn p);

sfixn* linearCombineModulus(sfixn *Cs, subProdTree *tree, sfixn p);

sfixn* fastInterp(sfixn* polyDg, sfixn n, sfixn *Us,  subProdTree *tree,
       sfixn *Vs, sfixn p);

sfixn* direvative(sfixn deg, sfixn *coef, sfixn p);

sfixn* createPts(sfixn start, sfixn m, sfixn *bounds, sfixn p);

preFFTRep* fastInterpMulti_test(sfixn N, sfixn M, sfixn *dims, sfixn *EEE,
           sfixn **UsPtr, subProdTree **trees, MONTP_OPT2_AS_GENE *pPtr);

sfixn** convertpts2ptsPtr(sfixn m, sfixn *bounds, sfixn *pts);

PTS_TREE* createGoodPtsForf1f2(sfixn N, sfixn d1, preFFTRep *f1, sfixn d2,
          preFFTRep *f2, sfixn m, sfixn *bounds, sfixn *dims1, sfixn *dims2,
          MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthSlice(sfixn ith, sfixn N, sfixn m, sfixn w, sfixn slicesz,
           sfixn *slicedims, sfixn Ssz, sfixn *S, PTS_TREE *pts_tree,
           MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthDthSlice(sfixn ith, sfixn dth, sfixn N, sfixn m, sfixn w,
           sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S,
           PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLC(int *nextiAddr, int start, sfixn N,
           sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims,
           sfixn Ssz, sfixn *S, PTS_TREE *pts_tree, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLT(int *nextiAddr, int start, sfixn N,
           sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims,
           sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE *pPtr);

/*
 * FFT based methods
 */

/*
 * Create a set of primitive roots of unity such that the leading coefficients
 * of f1 and f2 will not vanish after evaluation. New version.
 */
int findGoodRoots(sfixn N, preFFTRep *f1, preFFTRep *f2, sfixn *es, sfixn *dims,
    sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

/**
 * Old version, 
 */
int createGoodRootsForf1f2(sfixn N, sfixn M, sfixn d1, preFFTRep *f1, sfixn d2,
    preFFTRep *f2, sfixn *es, sfixn *dims, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

/*
 * Fill S with the evaluations on a grid. Creating a so-called SCube.
 * */
void getSubResultantChains(sfixn N, sfixn w, sfixn Ssz, sfixn *S, sfixn *Edims1,
     sfixn *E1, sfixn *Edims2, sfixn *E2, MONTP_OPT2_AS_GENE *pPtr);

void printAllSRS(sfixn no, sfixn w, sfixn *AllSRS);

/*
 * Returning the evaluation vector of the subresultant with index i.
 * S0 is the subresultant.
 * Also works for a TFT SCube.
 */
sfixn* get_ithSlice_fromSubResultantChains(sfixn ith, sfixn N, sfixn w,
       sfixn slicesz, sfixn *slicedims, sfixn Ssz, sfixn *S);
/*
 * Returning the evaluation vector of Coeff(Si, y^d), for i from 0 to w-1.
 * S0 is the resultant.
 * Also works for a TFT SCube.
 */
sfixn* get_ithDthSubSlice_fromSubResultantChains(sfixn ith, sfixn dth, sfixn N_1,
       sfixn w, sfixn subslicesz, sfixn *subslicedims, sfixn Ssz, sfixn *S);

int tracingNextCandidateSlice(int start, sfixn w, sfixn Ssz, sfixn *S);

void set_ithSliceZero_fromSubResultantChains(sfixn ith, sfixn N, sfixn w,
     sfixn slicesz, sfixn *slicedims, sfixn Ssz, sfixn *S);

void fastDftMulti_test(sfixn M, sfixn *es, sfixn *dims, sfixn *E, preFFTRep* poly,
     sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* fastInvDftMulti_test(sfixn N, sfixn M, sfixn *es, sfixn *dims,
           sfixn *EEE, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthSliceDFT(sfixn ith, sfixn N, sfixn m, sfixn w, sfixn slicesz,
           sfixn *slicees, sfixn *slicedims, sfixn Ssz, sfixn *S,
           sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthDthSliceDFT(sfixn ith, sfixn dth, sfixn N, sfixn m, sfixn w,
           sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims, sfixn Ssz,
           sfixn *S, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLCDFT(int *nextiAddr, int start, sfixn N,
           sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees,
           sfixn *subslicedims, sfixn Ssz, sfixn *S, sfixn *rootsPtr,
           MONTP_OPT2_AS_GENE *pPtr);

sfixn SlowEvaluation1pt(sfixn degf, sfixn *fPtr, sfixn pt, sfixn p);

void getQuotients(sfixn N, sfixn dd, sfixn Qsz, sfixn *Q, sfixn* Edims1,
     sfixn *E1, sfixn *Edims2, sfixn *E2, MONTP_OPT2_AS_GENE *pPtr, int opt);

void permuteSlice1toN(sfixn N, sfixn slicesz, sfixn *slicedims, sfixn *slice);

preFFTRep* EX_QuoMulti(preFFTRep *f1, preFFTRep *f2, sfixn N,
           MONTP_OPT2_AS_GENE *pPtr, int opt);

preFFTRep* interpNextCandidateSliceLCDefective(int *nextiAddr, int start,
           sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicedims,
           sfixn Ssz, sfixn *S,  PTS_TREE* pts_tree, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLCDFTDefective(int *nextiAddr, int start,
           sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees,
           sfixn *subslicedims, sfixn Ssz, sfixn *S, sfixn *rootsPtr,
           MONTP_OPT2_AS_GENE *pPtr);

/*
 *  TFT based methods
 * */

preFFTRep* TFTInterpMultiD(sfixn N, sfixn M, sfixn *es, sfixn *dims, sfixn *E,
           sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

void TFTEvalMultiD(sfixn M, sfixn *es, sfixn *dims, sfixn *E, preFFTRep *poly,
     sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

void getSubResultantChainsTFT(sfixn N, sfixn M, sfixn w, sfixn Ssz, sfixn *S,
     sfixn *Edimsq, sfixn *E1, sfixn *Edims2, sfixn*E2,
     MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthSliceTFT(sfixn ith, sfixn N, sfixn m, sfixn w,
           sfixn slicesz, sfixn *slicees, sfixn *slicedims, sfixn Ssz,
           sfixn *S, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpIthDthSliceTFT(sfixn ith, sfixn dth, sfixn N, sfixn m,
           sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims,
           sfixn Ssz, sfixn *S,  sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLCTFT(int *nextiAddr, int start, sfixn N,
           sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees, sfixn *subslicedims,
           sfixn Ssz, sfixn *S, sfixn *rootsPtr, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep* interpNextCandidateSliceLCTFTDefective(int *nextiAddr, int start,
           sfixn N, sfixn m, sfixn w, sfixn subslicesz, sfixn *subslicees,
           sfixn *subslicedims, sfixn Ssz, sfixn *S, sfixn *rootsPtr,
           MONTP_OPT2_AS_GENE * pPtr);
#endif
