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
#ifndef __Factorization_h
#define __Factorization_h 

#include "Types.h"
#include "generalFuncs.h"
#include "MPMMTS.h"
#include "FINTERP.h"
#include "HGCD.h"

sfixn *LcmPolyPair(sfixn *dgLcmAddr, sfixn d1, sfixn *f1, sfixn d2, sfixn *f2, 
    MONTP_OPT2_AS_GENE * pPtr );

sfixn *SquareFreeFact(sfixn *degR, sfixn degF, sfixn *FPtr, sfixn p);

sfixn *SquarefreePart(sfixn *degQAddr, sfixn degF, sfixn *FPtr, 
    MONTP_OPT2_AS_GENE *pPtr);

#endif
