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


#ifndef __utilityFuncs_h
#define __utilityFuncs_h 

#include "CONSTANTS.h"
#include "AS.h"
#include "inlineFuncs.h"



extern sfixn BASE;
extern sfixn BASE_1;
extern sfixn BASEHALF;


/**
 * PolyCleanData:
 * @prePtr:  A C-Cube polynomial.
 * 
 * make prePtr a zero polynomial.
 * Return value: 
 **/
//#ifndef _mcompile_
void PolyCleanData(preFFTRep * prePtr);


/**
 * egcd:
 * @x: An int.
 * @y: An int.
 * @ao: A pointer of int.
 * @bo: A pointer of int.
 * Extended Euclidean Gcd of machine integers.
 * Return value: (*bo)x + (*ao)y = *vo. 
 **/
extern void egcd (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo);




/**
 * inverseMod:
 * @n: A sfixn in Z/pZ.
 * @p: A prime.
 * 
 * computer the invers of n modulo p, where n is reduced with respect to p.
 * Return value: (1/n) mod p 
 **/

static inline sfixn inverseMod(sfixn n, sfixn p){
    sfixn a, b, v;
    egcd(n, p, &a, &b, &v);
    if (b < 0)
        b += p;
    return b % p;
}


/**
 * DivMod:
 * @a: A sfixn number in Z/pZ.
 * @b: A sfixn number in Z/pZ.
 * @p: A prime number.
 * Compute the faction a/b modulo p.
 * Return value: (a/b) mod p. 
 **/

static inline sfixn DivMod(sfixn a, sfixn b, sfixn p){
    return MulMod(a, inverseMod(b, p), p);
}

/**
 * QuoMod:
 * @a: A sfixn number in Z/pZ.
 * @b: A sfixn number in Z/pZ.
 * @p: A prime number.
 * Compute the faction a/b modulo p.
 * Return value: (a/b) mod p. 
 **/

static  inline sfixn QuoMod(sfixn a, sfixn b, sfixn p){
    return MulMod(a, inverseMod(b, p), p);
}

/**
 * setPolyOne:
 * @Ptr: A C-Cube polynomial.
 * 
 * Set 'Ptr' to a constant polynomial has value 1.
 * Return value: 
 **/

static inline void setPolyOne(preFFTRep *Ptr){
    PolyCleanData(Ptr);
    DATI(Ptr, 0)=1;
}

/**
 * PowerMod:
 * @a: A sfixn number in Z/nZ.
 * @ee: A sfixn number.
 * @n: A moduli.
 * Compute the power a^ee modulo n.
 * Return value: (a^ee) mod n. 
 **/
sfixn PowerMod(sfixn a, sfixn ee, sfixn n);


/**
 * MontMulMod_OPT2_AS_GENE:
 * @a: A fixnum.
 * @b: A fixnum.
 * @pPtr: Information for the prime number p.
 * 
 * A improved Montgamoney Trick only works for Fourier Primes.
 * Please see the JSC paper for details.
 * Return value: (a*b)/R mod p, where R is the next power of 2 of p.
 **/
sfixn MontMulMod_OPT2_AS_GENE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr);

sfixn MontMulMod_OPT2_AS_Double_GENE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr);

/**
 * MontMulMod_OPT2_AS_GENE_SPE:
 * @a: A fixnum.
 * @b: A fixnum.
 * @pPtr: Information for the prime number p.
 * 
 * A improved Montgamoney Trick only works for Fourier Primes.
 * Please see the JSC paper for details.
 * This rountine is for special Fourier Primes who are in the shape of N+1,
 * where N is a power of 2.
 * Return value: (a*b)/R mod p, where R is the next power of 2 of p.
 **/

sfixn MontMulMod_OPT2_AS_GENE_SPE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr);

sfixn MontMulMod_OPT2_AS_Double_GENE_SPE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr);

/**
 * MultiNumbPolyMul_1:
 * @r: r is a scalar.
 * @f: f is a C-Cube polynomial.
 * @pPtr: Information for the prime p.
 * 
 * Compute the product of r and f modulo p.
 * Return value: r*f mod p.
 **/

void MultiNumbPolyMul_1(sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr);

//=====================================================
//  copying data from one dense multivariate polynomial
//  to the one.
//=====================================================

void fromtofftRep(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs);

/**
 * zeroCoefp:
 * @coefPtr: A vector.
 * @coefSiz: Size of 'coefPtr'.
 * 
 * Test if vector 'coefPtr' of size 'coefSiz' is a zero vector.
 * Return value: 1 for true, 0 for false. 
 **/
int zeroCoefp(sfixn * coefPtr, sfixn coefSiz);


/**
 * shrinkDeg:
 * @deg: The original degree.
 * @poly: The coefficient vector.
 * @coefSiz: The size of a coefficient.
 * 
 * View the 'poly' as an univariate polynomial has 'deg'+1 coefficients.
 * Each coefficient has size of 'coefSiz'.
 * By ignoring the leading zeros, we compute the actual degree.
 * Return value: The actual degree.
 */
sfixn shrinkDeg(sfixn deg, sfixn * poly, sfixn coefSiz);

/*
Another implementation of shrinkdegree for bivariate polynomial
sfixn newShrinkDeg(sfixn deg1, sfixn deg2, sfixn * poly)
{
	register int i, j;
	i = (deg1+1)*(deg2+1) - 1;
	j = i;
	for(; i >= 0; --i)
	{
		if(poly[i] != 0) break;			
	}
	i = j - i;
	
	return (deg2 - (i/(deg1+1)));
	 	
}
*/

/**
 * shrinkDegUni:
 * @deg: The original degree.
 * @coef: The coefficient vector.
 * 
 * By ignoring the leading zeros, we compute the actual degree.
 * Return value: The actual degree.
 **/
sfixn shrinkDegUni(sfixn deg, sfixn * cof);


/**
 * setLeadingCoefMultiOne:
 * @N: X_N is the main variable of 'f'.
 * @f: a C-Cube polynomial.
 * 
 * Set f's leading coefficient to be 1.
 * Return value: f.
 **/

preFFTRep *setLeadingCoefMultiOne(sfixn N, preFFTRep* f);

void MultiNumbPolyMulMonicize_1(sfixn N, sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr);
void subEqDgPoly_inner_1(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector);

void subEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int selector);

// return 1 means It IS a Zero.
/**
 * zeroPolyp:
 * @polyPtr: A C-Cube polynomial.
 * 
 * To check if 'polyPtr' is a zero polynomial
 * Return value: 1 for true. 0 for false.
 **/
sfixn zeroPolyp(preFFTRep * polyPtr);


// return 1 means It IS a Zero.
/**
 * constantPolyp:
 * @polyPtr: A C-Cube polynomial.
 * 
 * To check if 'polyPtr' is a constant polynomial
 * Return value: 1 for true. 0 for false.
 **/
sfixn constantPolyp(preFFTRep * polyPtr);


/**
 * negatePoly_1:
 * @prePtr:  A C-Cube polynomial.
 * @p: A prime number.
 *
 * Negate all coefficients of 'prePtr'.
 * Return value: 
 **/
//#ifndef _mcompile_
void negatePoly_1(preFFTRep * prePtr, sfixn p);

void addEqDgPoly_inner(sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p);

/**
 * addEqDgPoly_1:
 * @N: Number of variables in 'Ptr1' and 'Ptr2'.
 * @Ptr1: A C-Cube polynomial.
 * @Ptr2: A C-Cube polynomial.
 * @p: A prime number.
 * 
 * Suppose 'Ptr1' and 'Ptr2' has the same dimension and size.
 * Compute the sum of them.
 * Return value: Ptr1 = Ptr1 + Ptr2;
 **/

void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);


// we will use the smaller dgs which are 2nd. but accums should used prepectively.
//#ifndef _mcompile_
void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p);

// this in-place fun suppose Ptr1 is the larger buffer on all dimensions.

/**
 * subPoly_1:
 * @N: Number of variables in 'Ptr1' and 'Ptr2'.
 * @Ptr1: A C-Cube polynomial.
 * @Ptr2: A C-Cube polynomial.
 * @p: A prime number.
 * 
 * Suppose 'Ptr1' has the same dimension and larger(or equal) size on each dimension.
 * Compute the difference of them.
 * Return value: Ptr1 = Ptr1 - Ptr2;
 **/
//#ifndef _mcompile_
void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p);

void increaseOneDim(preFFTRep * Ptr);

/**
 * MontMulMod:
 * @a: A sfixn.
 * @b: A sfixn.
 * @pPtr: Information for prime p.
 * 
 * Compute 'a'*'b' modulo p Based on improved Montgomery trick.
 * Return value a*b mod p: 
 **/

sfixn MontMulMod(sfixn a, sfixn b, MONTP_OPT2_AS_GENE *pPtr);

/**
 * coMulVec_1:
 * @co: A coefficient.
 * @deg: Degree of univeriate polynomial f.
 * @vec: Coefficient vector of the univeriate polynomial f.
 * @pPtr: The Info for a prime p.
 * 
 * compute the product of 'co' and 'f' modulo p.
 * Return value: co*f mod p. 
 **/
//#ifndef _mcompile_
void coMulVec_1(sfixn co, sfixn deg, sfixn *vec, MONTP_OPT2_AS_GENE *pPtr);

/**
 * coMulAddVec:
 * @co: A coefficient.
 * @deg: Degree of univeriate polynomial f1 and f2.
 * @vec1: Coefficient vector of the univeriate polynomial f1.
 * @vec2: Coefficient vector of the univeriate polynomial f2.
 * @pPtr: The Info for a prime p.
 * 
 * compute the product of 'co' and 'f' modulo p.
 * Return value: co*f mod p. 
 **/
//#ifndef _mcompile_
void coMulAddVec(sfixn co, sfixn deg, sfixn *vec1, sfixn *vec2, MONTP_OPT2_AS_GENE *pPtr);

/**
 * getLeadingCoefMulti:
 * @N: X_N is the main variable of f.
 * @co: co is an empty buffer.
 * @f: A C-Cube polynomial.
 * Copy the leading coefficient of 'f' into 'co'.
 * Return value: co -- the leading coefficient of 'f'.
 **/
preFFTRep *getLeadingCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f);

/**
 * setCoefMulti:
 * @N: X_N is the main variable of f.
 * @f: A C-Cube polynomial.
 * @co: A coefficient whoes main variable is X_{N-1}.
 * @j: a index number. 
 * Set f's j-th coefficient to be 'co'.
 * Return value: co.
 **/
//#ifndef _mcompile_
preFFTRep *setCoefMulti(sfixn N, preFFTRep* f, preFFTRep* co, sfixn j);

#endif
