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
#ifndef __GCD_h
#define __GCD_h 

#include "Types.h"
#include "generalFuncs.h"
#include "FMUL.h"
#include "FDIV.h"




sfixn
gcd_Uni_1(sfixn * change, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr,MONTP_OPT2_AS_GENE * pPtr, sfixn );

sfixn
ExGcd_Uni_1(sfixn * change,sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr,MONTP_OPT2_AS_GENE * pPtr, sfixn * degS, sfixn * degT, sfixn * S1Ptr, sfixn * S2Ptr,sfixn * T1Ptr, sfixn * T2Ptr, sfixn);

void 
plainDiv(sfixn * RPtr, sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr );

void
ExGcd_Uni(sfixn * uPtr, sfixn *dC, sfixn * vPtr, sfixn *dD, sfixn * gcdPtr, sfixn *dG, 
	  sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);

void
Gcd_Uni(sfixn * gcdPtr, sfixn *dG, 
	sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);


sfixn *
EX_GCD_UNI(sfixn *dGAddr, 
	   sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);


int
ExGcd_Uni_RFR(sfixn d, sfixn * vPtr, sfixn *dD, sfixn * gcdPtr, sfixn *dG,  sfixn * fPtr, sfixn dA, sfixn * gPtr, sfixn dB, MONTP_OPT2_AS_GENE * pPtr);

void
normalize_1(sfixn deg, sfixn * cof, MONTP_OPT2_AS_GENE * pPtr);


#endif
