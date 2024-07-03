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
#ifndef __IsInvertible
#define __IsInvertible

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include "LinkedList.h"
#include "IteratedResultant.h"
#include "RegularGcd.h"
#include <math.h>

LinkedQueue *isInvertible_zeroDim(preFFTRep *poly, TriSet *ts, 
    MONTP_OPT2_AS_GENE *pPtr);

/* 

Though the following functions were not declared in this header file,
these are defined and used in the c file. 
Sardar Haque on April 5, 2012 add these declarations 
in the header file.
*/

LinkedQueue* isInvertible_zeroDim_Inner(preFFTRep *poly, TriSet *ts, 
    MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *EX_RegularizeInitial(preFFTRep *InPoly, TriSet *InTs, 
    MONTP_OPT2_AS_GENE *pPtr);



#endif
