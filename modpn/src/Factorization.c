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
#include "Factorization.h"

#define printFact 1

extern int Interrupted;

// input: f1 !=0, f2 !=0;
// output: (*dgLcmAddr, Lcm)
/**
 * LcmPolyPair:
 * @dgLcmAddr: .
 * @d1: degree of f1.
 * @f1: the coefficient vector of f1.
 * @d2: degree of f2.
 * @f2: the coefficient vector of f2.
 * @pPtr: info about the prime number 'p'.
 *   assume f1, f2 are not zeros.
 * Return value: the Least Common Multiple of unvariate polynomials f1, f2.
 **/

sfixn *LcmPolyPair(sfixn *dgLcmAddr, sfixn d1, sfixn *f1, sfixn d2, sfixn *f2, 
        MONTP_OPT2_AS_GENE * pPtr )
{
    sfixn degR, *degRAddr=&degR, *RPtr, *Lcm, *Prod, dgProd, 
          *dgProdAddr=&dgProd, *Gcd, dGcd, *dGcdAddr=&dGcd;   

    if(d1<d2) return LcmPolyPair(dgLcmAddr, d2, f2, d1, f1, pPtr);

    degR=d2-1;
    RPtr=(sfixn *)my_calloc(degR+1, sizeof(sfixn));
    UniRem(degRAddr, RPtr, d1, f1, d2, f2, pPtr);
    if(isZeroVec(degR+1, RPtr)){
        my_free(RPtr);
        *dgLcmAddr=d1;
        return EX_copyVec_0_to_d(d1, f1);
    }
    my_free(RPtr);
    Prod=EX_Mont_Mul(dgProdAddr, d1, f1, d2, f2, pPtr);
    Gcd = EX_GCD_UNI(dGcdAddr, f1, d1, f2, d2, pPtr);
    Lcm=EX_UniQuo(dgLcmAddr, dgProd, Prod, dGcd, Gcd, pPtr);
    my_free(Prod);
    my_free(Gcd);
    return Lcm;
}

// input poly (degF, FPtr)
// output square free part ((*degQAddr), sqrFreePtr)
/**
 * SquareFreeFact:
 * @degQAddr: (output) a pointer to the degree of the output square-free.
 * @degF: The degree of unvariate polynomal F.
 * @FPtr: The coefficient vector of F.
 * @p: Prime. 
 *
 * Return value: The coefficient vector of square-free part. 
 **/
sfixn *SquareFreeFact(sfixn *degQAddr, sfixn degF, sfixn *FPtr, sfixn p)
{
    sfixn *gcd, *sqrFreePtr;
    sfixn dG=0, *dGAddr=&dG;
    sfixn *FPrime;
    MONTP_OPT2_AS_GENE  prime;
    EX_MontP_Init_OPT2_AS_GENE(&prime, p);
    FPrime=direvative(degF, FPtr, p);
    gcd=EX_GCD_UNI(dGAddr, FPtr, degF, FPrime, degF-1, &prime);
    my_free(FPrime);
    sqrFreePtr=EX_UniQuo(degQAddr, degF, FPtr, dG, gcd, &prime);
    return sqrFreePtr;
}

sfixn *SquarefreePart(sfixn *degQAddr, sfixn degF, sfixn *FPtr, 
    MONTP_OPT2_AS_GENE *pPtr)
{
    sfixn *gcd, *sqrFreePtr;
    sfixn dG = 0, *dGAddr = &dG;
    sfixn *FPrime;

    FPrime = direvative(degF, FPtr, pPtr->P);

    //degF = shrinkDegUni(degF, FPrime);
    
	
    // gcd = EX_GCD_UNI(dGAddr, FPtr, degF, FPrime, degF-1, pPtr);
    gcd = EX_GCD_Uni_Wrapper(dGAddr, FPtr, degF, FPrime, degF-1, pPtr);
    my_free(FPrime);
	
    // the following if condition is added by Sardar Haque on 20th November 2012
    if(dG == 0 && gcd[0] == 0)
    {
     	 sqrFreePtr =  (sfixn *) my_malloc((degF+1)*sizeof(sfixn));	
	 sfixn i;
	 for(i = 0; i <= degF; ++i) sqrFreePtr[i] = FPtr[i];
	 *degQAddr = degF;

    }
    else
    {
	sqrFreePtr = EX_UniQuo(degQAddr, degF, FPtr, dG, gcd, pPtr);
    }			
    return sqrFreePtr;
}

