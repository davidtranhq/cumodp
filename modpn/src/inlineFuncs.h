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
#ifndef __inlineFuncs_h
#define __lnlineFuncs_h 

#include "CONSTANTS.h"



extern sfixn BASE;
extern sfixn BASE_1;
extern sfixn BASEHALF;

/**
 * AddMod:
 * @a: a int.
 * @b: b int.
 * @p: Prime number.
 * 
 * Return value: (a+b) mod p. 
 **/
static inline
sfixn AddMod(sfixn a, sfixn b, sfixn p){
    sfixn r = a + b;
    r -= p;
    r += (r >> BASE_1) & p;
    return r;
}

/**
 * SubMod:
 * @a: a int.
 * @b: b int.
 * @p: Prime number.
 * 
 * Return value: (a-b) mod p. 
 **/
static inline
sfixn SubMod(sfixn a, sfixn b, sfixn p){
    sfixn r = a - b;
    r += (r >> BASE_1) & p;
    return r;}

/**
 * NegMod:
 * @a: a int.
 * @p: Prime number.
 * 
 * Return value: (-a) mod p. 
 **/
static inline
sfixn NegMod(sfixn a,  sfixn p){
    sfixn r = - a;
    r += (r >> BASE_1) & p;
    return r;}

// This may make things slower!
/**
 * MulMod:
 * @a: An int.
 * @b: An int.
 * @n: An int.
 * 
 * Return value: (a*b) mod n. 
 **/
static  inline
sfixn MulMod(sfixn a, sfixn b, sfixn n)
{
    sfixn q, res;

    double ninv=1/(double)n;

    q  = (sfixn) ((((double) a) * ((double) b)) * ninv);
    res = a*b - q*n;
    res += (res >> BASE_1) & n;
    res -= n;
    res += (res >> BASE_1) & n;
    return res;
}




static  inline unsigned int
partialBitRev(register unsigned int x, int n)
{
    register unsigned int y = 0x55555555;
    x = (((x >> 1) & y) | ((x & y) << 1));
    y = 0x33333333;
    x = (((x >> 2) & y) | ((x & y) << 2));
    y = 0x0f0f0f0f;
    x = (((x >> 4) & y) | ((x & y) << 4));
    y = 0x00ff00ff;
    x = (((x >> 8) & y) | ((x & y) << 8));
    y =((x >> 16) | (x << 16));
    return y>>(32-n);
}





//=========================================================
// A technical step before using Montgomery trick.
//=========================================================
static inline
sfixn convertToMondMulR(sfixn r,  MONTP_OPT2_AS_GENE * pPtr){
    sfixn p=pPtr->P, R=(1L<<pPtr->Rpow)%p;
    R=(MulMod(r, R, p))<<pPtr->Base_Rpow;
    return R;
}

static inline
void nextMCoefData(preFFTRep * Ptr, sfixn N, sfixn M){
    DAT(Ptr)+=M*(CUMI(Ptr, N));
}





static inline
void setData(preFFTRep * Ptr, sfixn * dataPtr){ 
    DAT(Ptr)=dataPtr;
}

static inline
void backupData(preFFTRep * Ptr){
    TMPDAT(Ptr)=DAT(Ptr);
}

static inline
void restoreData(preFFTRep * Ptr){
    DAT(Ptr)=TMPDAT(Ptr);
}

static inline
void resumeData(preFFTRep * Ptr){ 
    DAT(Ptr)=DEFDAT(Ptr);
}

static inline
void decreaseOneDim(preFFTRep * Ptr){
    SIZ(Ptr)=CUMI(Ptr, N(Ptr));
    (N(Ptr))--;
}

static inline
void resumeDim(preFFTRep * Ptr){
    N(Ptr)=DFN(Ptr);
    SIZ(Ptr)=DFSIZ(Ptr);
    DAT(Ptr)=DEFDAT(Ptr);
}
static inline
void nextCoefData(preFFTRep * Ptr, sfixn N){
    DAT(Ptr)+=CUMI(Ptr, N);
}


#endif

