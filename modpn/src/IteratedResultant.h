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
#ifndef __Iter_h
#define __Iter_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include <math.h>
#include <string.h>

/* CUDA related */
//#if defined(_cuda_modp_)
//#include "cumodp.h"
//#endif
#ifdef _cuda_modp_
#include "cumodp.h"
#endif

void EX_SCUBE_Free(SCUBE *scube);

void EX_SCUBE_Print(SCUBE *scube);

SCUBE *EX_SCUBE_Init(ChainType_t ct, void *chain);

SCUBE *EX_SubResultantChain(preFFTRep *f1, preFFTRep *f2, sfixn N, 
    MONTP_OPT2_AS_GENE *pPtr);

SCUBE* EX_SubResultantChainSelect(computingChoince method, preFFTRep *f1, preFFTRep *f2,
    sfixn N, MONTP_OPT2_AS_GENE *pPtr);

SCUBE* EX_SubResultantChainOpt(preFFTRep *f1, preFFTRep *f2, sfixn N,
            MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *EX_ResultantFromChain(SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *EX_Resultant_Multi(preFFTRep *f1, preFFTRep *f2, sfixn N, 
    MONTP_OPT2_AS_GENE *pPtr);

sfixn iteratedResultant_zerodim(preFFTRep *poly, TriSet *ts, 
    MONTP_OPT2_AS_GENE *pPtr);

sfixn *iteratedResultant_onedim(sfixn *resDgAddr, preFFTRep *poly, 
    TriSet *ts, sfixn bound, sfixn freeVarNo, MONTP_OPT2_AS_GENE *pPtr);

sfixn EX_WidthInSCUBE(SCUBE *scube);

preFFTRep *EX_IthDthCoeff(sfixn i, sfixn d, SCUBE *scube, 
    MONTP_OPT2_AS_GENE *pPtr);

preFFTRep *EX_IthSubres(sfixn i, SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr);

#endif
