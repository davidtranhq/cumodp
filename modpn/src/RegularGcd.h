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
#ifndef __RegularGcd
#define __RegularGcd

#include "Types.h"
#include "generalFuncs.h"
#include "FINTERP.h"
#include "MPMMTS.h"
#include "SubResultantSeq.h"
#include "LinkedList.h"
#include "IteratedResultant.h"
#include "IsInvertible.h"
#include <math.h>

LinkedQueue *
EX_RegularizeList_1(LinkedQueue *RegQueue, 
                    LinkedQueue *ToCheckQueue, 
                    TriSet *ts, 
                    MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcd(preFFTRep *f1, 
              preFFTRep *f2, 
              TriSet *ts,  
              SCUBE *scube, 
              MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcdNew(preFFTRep *f1, 
                 preFFTRep *f2, 
                 TriSet *ts,
                 SCUBE *scube, 
                 MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcdImp(preFFTRep *f1, 
                 preFFTRep *f2, 
                 TriSet *ts,
                 SCUBE *scube, 
                 MONTP_OPT2_AS_GENE *pPtr);

LinkedQueue *
EX_RegularGcd_Wrapped(preFFTRep *f1, 
                      preFFTRep *f2, 
                      TriSet *ts,  
                      sfixn M, 
                      MONTP_OPT2_AS_GENE *pPtr);

#endif
