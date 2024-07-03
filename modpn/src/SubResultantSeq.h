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
#ifndef __SubResultantSeq_H
# define __SubResultantSeq_H

#include "Types.h"
#include "generalFuncs.h"
#include "MPMMTS.h"
#include "FDIV.h"

sfixn *
SubResultantSeq(sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);

sfixn *
SubResultantSeq_1(sfixn w, sfixn Ssz, sfixn *S, sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);


sfixn *
SubResultantSeq_1_new(sfixn w, sfixn Ssz, sfixn *S, sfixn dg1, sfixn *f1, sfixn dg2, sfixn *f2, MONTP_OPT2_AS_GENE *pPtr);

void
printSRS(sfixn w, sfixn *SRS);

sfixn EX_Resultant_Uni(sfixn dgP, sfixn *P, sfixn dgQ, sfixn *Q, MONTP_OPT2_AS_GENE *pPtr);


#endif 
