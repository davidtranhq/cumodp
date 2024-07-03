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


#include "solve2.h"


/*
If computingChoince == GPUsmart
 these 4 values are threshold for GPU computing.
 If a polynomial's degree is less than threshold_gcdLow (threshold_divLow) or
 bigger than threshold_gcdHigh (threshold_divHigh) their
 gcd (division) will be done CPU.
*/
#define threshold_gcdLow  1000 // it should be at least 2
#define threshold_gcdHigh  30000
#define threshold_divLow  2000 // it should be at least 2
#define threshold_divHigh  40000

#define threshold_divDiff  200 // the difference between the degrees of the univariate polynomials
#define threshold_div  2000

#define gcdLowLevel1 2000
#define gcdHighLevel1 18000

#define gcdLowLevel2 5000
#define gcdHighLevel2 20000

#define gcdLowLevel3 10000
#define gcdHighLevel3 30000




#define PrintFlag 0
#define PrintSpecial 0


extern sfixn CUDA_TAG;

/*
input:
@rc: regular chain.

Returns a number which is the sum of all
coefficients over the regular chains and other 
numbers such as partial degrees total number of pairs 
in regular chain.

*/
sfixn totalSizeRC(LinkedQueue *queue)
{
	sfixn total = 1; 
	LinearNode *current;
	regular_chain2 *T;
        preFFTRep *Ptr;

        if(! EX_LinkedQueue_IsEmpty(queue))
	{
	        current = queue->front;
	        while( current != NULL)
		{
		    total += 4;
		    T = (regular_chain2*) (current->element);		

		    Ptr = (preFFTRep *) T->poly0;
		    if(Ptr != NULL)
		    {   
			if( (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1) > 1)	
				total += (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1);
		    	
		    }
			//else total += 1;	
       		        
	     	    Ptr = (preFFTRep *) T->poly1;
		    if(Ptr != NULL) 	
		    {
			if( (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1)  > 1)
				total += (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1);
		    }
		    //else total += 1;	

		    current=current->next;
        	}
    	}   
        
	return total;
}


void makecompactRC(sfixn totalSize, sfixn *compactRC, LinkedQueue *queue)
{
	sfixn i = 1, j =0;
        compactRC[0] = 0;

	LinearNode *current;
	regular_chain2 *T;
        preFFTRep *Ptr;

        if(! EX_LinkedQueue_IsEmpty(queue))
	{
	        current = queue->front;
	        while( current != NULL)
		{
		    ++compactRC[0];

		    T = (regular_chain2*) (current->element);		

		    Ptr = (preFFTRep *) T->poly0;
		    if(Ptr != NULL)   	
		    {
			compactRC[i++] = BUSZSI(Ptr, 1);
			compactRC[i++] = BUSZSI(Ptr, 2);
			if( (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1) > 1)
				for(j=0; j < (BUSZSI(Ptr, 1) + 1)*(BUSZSI(Ptr, 2) + 1); ++j)
					compactRC[i++] = DATI(Ptr, j);
		    }		
		    else
		    {
			compactRC[i++] = 0;
			compactRC[i++] = 0;
		    }
       		        
	     	    Ptr = (preFFTRep *) T->poly1;
		    if(Ptr != NULL)
		    {
			compactRC[i++] = BUSZSI(Ptr, 1);
			compactRC[i++] = BUSZSI(Ptr, 2);
			if( (BUSZSI(Ptr, 1)+1)*(BUSZSI(Ptr, 2)+1) > 1)
				for(j=0; j < (BUSZSI(Ptr, 1) + 1)*(BUSZSI(Ptr, 2) + 1); ++j)
					compactRC[i++] = DATI(Ptr, j);

		    }
		    else
		    {
			compactRC[i++] = 0;
			compactRC[i++] = 0;
		    }	
		    current=current->next;
        	}
    	}   
}


/**
 * Solving bivariate polynomial systems.
 *
 * WP, Sat Jan 1 21:15:07 EST 2011
 */

/**
 * next_regular_subres_index
 * 
 * @scube, the subresultant chain 
 * @i, the current index
 * 
 * Compute the smallest index j >= i, such that S[j] is a regular 
 * subresultant, w if no such a regular subresultant exists. 
 *
 * Assumption:
 *
 * (1) i >= 1 and i <= w
 */
sfixn next_regular_subres_index(sfixn i, SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr)
{
    preFFTRep *init=NULL;
    sfixn w, j;

    w = EX_WidthInSCUBE(scube);
    // the input polynomial with smaller degree is always regarded 
    // as a regular subresultant
    if (i == w) return w;
    for (j = i; j < w; ++j) {
        init = EX_IthDthCoeff(j, j, scube, pPtr);
        if (!zeroPolyp(init)) { return j; } 
    }
    return j;
}

/**
 * Specialized regular chain structure for the bivariate solver.
 */
regular_chain2 *EX_RegularChain2_Init(preFFTRep *f1, preFFTRep *f2) {
    regular_chain2 *T;
    T = (regular_chain2 *) my_malloc(sizeof(regular_chain2));
    T->poly0 = f1;
    T->poly1 = f2;
	
    return T;
    	
}

regular_chain2 *EX_RegularChain2_Copy_Init(preFFTRep *f1, preFFTRep *f2) {
    regular_chain2 *T;
    T = (regular_chain2 *) my_malloc(sizeof(regular_chain2));
    T->poly0 = (f1 != NULL) ? EX_CopyOnePoly(f1) : NULL; 
    T->poly1 = (f2 != NULL) ? EX_CopyOnePoly(f2) : NULL; 
    return T;
}

void EX_RegularChain2_Free(void *element) {
    regular_chain2 *T = (regular_chain2 *)element;
    if (T != NULL) {
        EX_FreeOnePoly(T->poly0);
        EX_FreeOnePoly(T->poly1);
        my_free(T);
    }
}

void FILE_RegularChain2_Print(void *element, FILE *file) {
    regular_chain2 *T = (regular_chain2*) element;
    if (T == NULL) fprintf(file, "[]");
    else{		
        fprintf(file,"[");
        fprintPoly(file, T->poly0);
        fprintf(file,", ");
        fprintPoly(file, T->poly1);
        fprintf(file,"]");
    }
}



void EX_RegularChain2_Print(void *element) {
    regular_chain2 *T = (regular_chain2*) element;
    if (T == NULL) printf("NULL regular_chain2\n");
    printf("\n{");
    printPoly(T->poly0);
    printf(", ");
    printPoly(T->poly1);
    printf("}\n");
}

/**
 * Check if F is a multiple of T.
 *
 * @F, univariate polynomial of degree df
 * @T, univariate polynomial of degree dT
 * 
 * returns 0 if not a multiple of T, 1 otherwise.
 *
 * Assumptions:
 *
 * (1) df is the true degree of F
 * (2) dT is the true degree of T
 * (3) df, dT > 0, nonconstant
 *
 */
sfixn is_divided_by(sfixn dF, sfixn *F, sfixn dT, sfixn *T, MONTP_OPT2_AS_GENE *pPtr) 
{
    sfixn *R, dR, result;
    if (dF < dT) return 0;
    if (DEBUG) { assert(dF >= 0 && dT > 0 && F[dF] && T[dT]); }

    R = EX_UniRem(&dR, dF, F, dT, T, pPtr);
    if(PrintFlag == 1)
    {
	printf("\n");
	printf("\n degree of the remainder: %d\n",dR);
	printf("\n");
	EX_printPolyUni(dR, R, 'a');

    }	
 
    // if dR = 0, either R is zero or R is a constant.
    if (dR > 0) {
        result = 0;
    } else {
        result = (R[0] == 0) ? 1 : 0;
    }

    my_free(R);
    return result;
}



sfixn is_divided_byCUDA(sfixn dF, sfixn *F, sfixn dT, sfixn *T, MONTP_OPT2_AS_GENE *pPtr)
{

		if(PrintFlag == 1)
		{
			printf("\nWorking isdivide() in CUDA?\n");
		}

	sfixn *R, *Q;
	Q = (sfixn *) my_malloc((dF - dT + 1)*sizeof(sfixn));
	R = (sfixn *) my_malloc(dT*sizeof(sfixn));

	float time = divPrimeField(F, T, R, Q, dF+1, dT+1, pPtr->P);
	//the following line is given for avoiding compiler warning
        time += 1.0;
	   

		if(PrintFlag == 1)
		{
			printf("Quo from CUDA: ");
			EX_printPolyUni(dF-dT+1, R, 'a');
			printf("\n Remainder from CUDA: \n");	
			EX_printPolyUni(dT, R, 'a');
		}

	my_free(Q);	
	
	sfixn i;
	for(i = dT-1; i >= 0; --i)
	{
		if(R[i] != 0)
			break;
	}
	my_free(R); 
	if(i == -1) return 1;
	  
	return 0;


}



/*
* Check if F is a multiple of T.

@method: whether it will do on CPU or GPU
@dF: degree of F
@F: univariate polynomial of degree df
@dT: degree of T
@T: univariate polynomial of degree dT
@pPtr: prime for field
 
returns 0 if not a multiple of T, 1 otherwise.

 * Assumptions:
 *
 * (1) df is the true degree of F
 * (2) dT is the true degree of T
 * (3) df, dT > 0, nonconstant
 *
*/
sfixn isDividedWithCudaOption(computingChoince method, sfixn dF, sfixn *F, sfixn dT, sfixn *T, MONTP_OPT2_AS_GENE *pPtr) 
{
	int dFF = shrinkDegUni(dF, F);
	int dTT = shrinkDegUni(dT, T);		
	dF = dFF; 
	dT = dTT;
	if(PrintFlag == 1)
	{
		printf("\nWorking isdivide() ?%d %d\n",dF,dT);
	}	
	if (dF < dT) return 0;
	if(dF == 0) return 1;

	if(CUDA_TAG == 0 || method == CPU )
		return is_divided_by(dF, F, dT, T, pPtr);
	if( dT <= 100 ||  dF <= 100  )
		return is_divided_by(dF, F, dT, T, pPtr);

	if(method == GPU)
		return is_divided_byCUDA(dF, F, dT, T, pPtr);

	if( dT < threshold_div && dF < threshold_div)
		return is_divided_by(dF, F, dT, T, pPtr);
   
	if(dT > threshold_divHigh || dF > threshold_divHigh)
		return is_divided_by(dF, F, dT, T, pPtr);

	return is_divided_byCUDA(dF, F, dT, T, pPtr);	
}






/**
 * modular_generic_solve2
 *
 * Input: F1, F2 in k[x, y] and g, h in k[x].
 * Output: regular chains (A1, B1), ..., (Ae, Be) such that
 *         V(F1, F2, g) \ V(h) = V(A1, B1) union ... union V(Ae, Be). 
 *  
 * The initials of F1 and F2 MAY NOT BE regular modulo the resultant R.
 * But one of the initials is regular modulo R. 
 *
 * Assumption:
 *
 * - mvar(F1) = mvar(F2).
 * - mdeg(F1) >= mdeg(F2).
 * - F1 and F2 are shrinked. 
 * - the resultant R of F1 and F2 in y is not zero. 
 *
 */
LinkedQueue *modular_generic_solve2(preFFTRep *F1, preFFTRep *F2, 
    preFFTRep *g, preFFTRep *h, SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr) 
{
    preFFTRep *R, *init, *initF2, *poly;
    sfixn *Rs, dRs;
    sfixn *hs, dhs;
    sfixn  *Q, dQ;
    sfixn  *T, dT;
    sfixn  *G, dG;
    sfixn w, i, j, dinit;
    LinkedQueue *resQueue;

    ///////////////////////////////////////////////////////////////////////////
    // Initialization & Setup
    ///////////////////////////////////////////////////////////////////////////

    // R is the resultant of F1 and F2, a univariate polynomial in x
    R = EX_IthDthCoeff(0, 0, scube, pPtr);

    if (DEBUG) {
        assert(N(R) == 1);
        assert(N(F1) == 2);
        assert(N(F2) == 2);
        assert(N(h) == 1);
        assert(N(g) == 1);
        assert(!zeroPolyp(R));
        assert(BUSZSI(R, 1) == shrinkDeg(BUSZSI(R, 1), DAT(R), CUMI(R, 1)));
        assert(BUSZSI(h, 1) == shrinkDeg(BUSZSI(h, 1), DAT(h), CUMI(h, 1)) );
        assert(BUSZSI(g, 1) == shrinkDeg(BUSZSI(g, 1), DAT(g), CUMI(g, 1)) );
        assert(shrinkDeg(BUSZSI(F1, 2), DAT(F1), CUMI(F1, 2)) >=  
               shrinkDeg(BUSZSI(F2, 2), DAT(F2), CUMI(F2, 2)));
    }
    
    // Rs is the squarefree part of R as a plain coefficient array 
    // the degree is given by dRs
    //
    // hs is the squarefree part of h as a plain coefficient array 
    // the degree is given by dhs
    
    Rs = SquarefreePart(&dRs, BUSZSI(R, 1), DAT(R), pPtr);
    hs = SquarefreePart(&dhs, BUSZSI(h, 1), DAT(h), pPtr);

    // Q = Rs / hs, the degree is given by dQ
    // T = gcd(Q, g), the degree is given by dT
    Q = EX_UniQuo(&dQ, dRs, Rs, dhs, hs, pPtr);
    
    // handle the case that g is zero
    if (zeroPolyp(g)) {
        T = Q, dT = dQ;
    } else {
        T = EX_GCD_Uni_Wrapper(&dT, Q, dQ, DAT(g), BUSZSI(g, 1), pPtr);
        my_free(Q);
    } 

    my_free(Rs);
    my_free(hs);

    ////////////////////////////////////////////////////
    // T is either 1 or non constant
    // if T is 1 then, no regular GCD to be computed.
    ////////////////////////////////////////////////////
    // if (DEBUG) assert(dT > 0);

    ///////////////////////////////////////////////////////////////////////////
    // Compute the regular GCDs
    ///////////////////////////////////////////////////////////////////////////
    w = EX_WidthInSCUBE(scube); 
    // w is the number of polynomials in the sub. res. chain
    // w is also the smallest of the main degees between F1 and F2
    i = 1;
    
    resQueue = EX_LinkedQueue_Init();
    initF2 = EX_getInitial(F2);

    while (dT > 0) {

        while (i <= w) {
            j = next_regular_subres_index(i, scube, pPtr);
            // the initial of the j-th regular subresultant
            // it must not be zero;
            if  (j < w) {
                init = EX_IthDthCoeff(j, j, scube, pPtr);
            } else {
                init = initF2;
            }

            // check if init is a multiple of T
            dinit = shrinkDeg(BUSZSI(init, 1), DAT(init), CUMI(init, 1));
            if (!is_divided_by(dinit, DAT(init), dT, T, pPtr)) break;
            ++i;
        }
        
        if (i > w) {
            EX_LinkedQueue_Enqeue(resQueue, 
                EX_RegularChain2_Init(
                    CreateUniPoly(dT, T), EX_CopyOnePoly(F1)));
            my_free(T);
            EX_FreeOnePoly(initF2);
            return resQueue;
        }

        G = EX_GCD_Uni_Wrapper(&dG, T, dT, DAT(init), dinit, pPtr);

        if (j < w)  {
            poly = EX_CopyOnePoly(EX_IthSubres(j, scube, pPtr));
        } else {
            poly = EX_CopyOnePoly(F2);
        }

        if (dG == 0) { 
            EX_LinkedQueue_Enqeue(resQueue, 
                EX_RegularChain2_Init(CreateUniPoly(dT, T), poly));
            my_free(T);
            EX_FreeOnePoly(initF2);
            return resQueue;
        }

        // Q = R / G
        Q = EX_UniQuo(&dQ, dT, T, dG, G, pPtr);
        EX_LinkedQueue_Enqeue(resQueue, 
            EX_RegularChain2_Init(CreateUniPoly(dQ, Q), poly));

        my_free(Q);
        my_free(T);
        T = G, dT = dG, i = j + 1;
    }
    EX_FreeOnePoly(initF2);
    return resQueue;
}

/**
 * Inner function for solving bivariate polynomial systems.
 *
 * @F1, bivariate polynomial in x, y
 * @F2, bivariate polynomial in x, y
 * @g,  univariate polynomial in x 
 * @pPtr, Fourier prime structure
 *
 * Assumptions:
 *
 * (1) g is non-constant, univariate 
 * (2) F1 and F2 are bivariate with positive partial degrees in y
 */
LinkedQueue *modular_solve2_inner(preFFTRep *F1, preFFTRep *F2, preFFTRep *g, 
    SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr)
{
    LinkedQueue *resQueue, *resQueue2;
    preFFTRep *h1, *h2, *t1, *t2, *hpoly, *gcdpoly;
    sfixn *h, dh, dh1, dh2;
    sfixn *gcd, dgcd, dg;

    h1 = EX_getInitial(F1); 
    h2 = EX_getInitial(F2); 
    dh1 = shrinkDegUni(BUSZSI(h1, 1), DAT(h1));
    dh2 = shrinkDegUni(BUSZSI(h2, 1), DAT(h2));
    h = EX_GCD_Uni_Wrapper(&dh, DAT(h1), dh1, DAT(h2), dh2, pPtr);
    hpoly = CreateUniPoly(dh, h);
    resQueue = modular_generic_solve2(F1, F2, g, hpoly, scube, pPtr);
    
    if (dh == 0) { 
        EX_FreeOnePoly(h1);
        EX_FreeOnePoly(h2);
        EX_FreeOnePoly(hpoly);
        my_free(h);
        return resQueue; 
    }

    t1 = EX_GetPolyTail(F1);
    t2 = EX_GetPolyTail(F2);
    dg = shrinkDeg(BUSZSI(g, 1), DAT(g), CUMI(g, 1));
    gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(g), dg, h, dh, pPtr);
    gcdpoly = CreateUniPoly(dgcd, gcd);

    resQueue2 = modular_solve2(t1, t2, gcdpoly, pPtr);
    EX_LinkedQueue_Concat_1(resQueue, resQueue2);

    EX_FreeOnePoly(h1);
    EX_FreeOnePoly(h2);
    EX_FreeOnePoly(hpoly);
    EX_FreeOnePoly(t1);
    EX_FreeOnePoly(t2);
    EX_FreeOnePoly(gcdpoly);
    my_free(h);
    my_free(gcd);
    EX_LinkedQueue_Free(resQueue2, EX_RegularChain2_Free);

    return resQueue;
}

/**
 * With option to select scube type.
 */
LinkedQueue *modular_solve2_select_inner(sfixn method, preFFTRep *F1, 
    preFFTRep *F2, preFFTRep *g, SCUBE *scube, MONTP_OPT2_AS_GENE *pPtr)
{
    LinkedQueue *resQueue, *resQueue2;
    preFFTRep *h1, *h2, *t1, *t2, *hpoly, *gcdpoly;
    sfixn *h, dh, dh1, dh2;
    sfixn *gcd, dgcd, dg;

    h1 = EX_getInitial(F1); 
    h2 = EX_getInitial(F2); 
    dh1 = shrinkDegUni(BUSZSI(h1, 1), DAT(h1));
    dh2 = shrinkDegUni(BUSZSI(h2, 1), DAT(h2));
    h = EX_GCD_Uni_Wrapper(&dh, DAT(h1), dh1, DAT(h2), dh2, pPtr);
    hpoly = CreateUniPoly(dh, h);
    resQueue = modular_generic_solve2(F1, F2, g, hpoly, scube, pPtr);
    
    if (dh == 0) { 
        EX_FreeOnePoly(h1);
        EX_FreeOnePoly(h2);
        EX_FreeOnePoly(hpoly);
        my_free(h);
        return resQueue; 
    }

    t1 = EX_GetPolyTail(F1);
    t2 = EX_GetPolyTail(F2);
    dg = shrinkDeg(BUSZSI(g, 1), DAT(g), CUMI(g, 1));
    gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(g), dg, h, dh, pPtr);
    gcdpoly = CreateUniPoly(dgcd, gcd);

    resQueue2 = modular_solve2_select(method, t1, t2, gcdpoly, pPtr);
    EX_LinkedQueue_Concat_1(resQueue, resQueue2);

    EX_FreeOnePoly(h1);
    EX_FreeOnePoly(h2);
    EX_FreeOnePoly(hpoly);
    EX_FreeOnePoly(t1);
    EX_FreeOnePoly(t2);
    EX_FreeOnePoly(gcdpoly);
    my_free(h);
    my_free(gcd);
    EX_LinkedQueue_Free(resQueue2, EX_RegularChain2_Free);

    return resQueue;
}

/**
 * Add a univariate polynomial f1 to a bivariate polynomial f2,
 * resulting a bivariate polynomial.
 */
preFFTRep *add_poly_uni_bi(preFFTRep *f1, preFFTRep *f2, sfixn p) 
{
    sfixn dgs[3], i, j, dx1, dx2, dy2, sz;
    preFFTRep *result;

    if (DEBUG) { assert(N(f1) == 1 && N(f2) == 2); }
    dx1 = shrinkDeg(BUSZSI(f1, 1), DAT(f1), CUMI(f1, 1));
	// There was a bug here. It is corrected by Sardar Haque on 26 November 2012
    dx2 = BUSZSI(f2, 1);//shrinkDeg(BUSZSI(f2, 1), DAT(f2), CUMI(f2, 1)); // it was done by shrinkDeg function previously
    dy2 = shrinkDeg(BUSZSI(f2, 2), DAT(f2), CUMI(f2, 2));

    dgs[0] = 0;
    dgs[1] = (dx1 > dx2) ? dx1 : dx2;
    dgs[2] = dy2;
    
    if (DEBUG) assert(dy2 > 0);
		
    result = EX_InitOnePoly(2, dgs);
		//
		if(PrintFlag == 1)	
		{
			printf("univariate: ");
			EX_Poly_Print(f1);
			printf("Bivariate: ");
			EX_Poly_Print(f2);
			printf("A polynomial is created with %d %d\n",dgs[1],dgs[1]);
			EX_Poly_Print(result);
		}
		//


    sz = dgs[1] + 1;
    for (i = 0; i <= dgs[2]; ++i) {
        for (j = 0; j <= dx2; ++j) {
            (result->data)[i * sz + j] = (f2->data)[i * (dx2+1) + j]; 	// There was a bug here. It is corrected by Sardar Haque on 26 November 2012
        }								// It was not (dx2+1) before in the inner loop.
    }

		//
		if(PrintFlag == 1)	
		{
			printf("Copy of the bivariate poly with some scaling: \n");
			EX_Poly_Print(result);
		}
		//

    for (i = 0; i <= dx1; ++i) {
        (result->data)[i] = AddMod((result->data)[i], (f1->data)[i], p);
    }

    return EX_NormalizePoly(result);
}

/**
 * Handle corner cases for bivariate solving.
 * 
 * Assumptions:
 *
 * (1) F1 is a univariate or bivariate polynomial
 * (2) F2 is a univariate or bivariate polynomial
 * (3) g is a univariate polynomial 
 */
LinkedQueue *modular_solve2(preFFTRep *F1, preFFTRep *F2, preFFTRep *g, 
    MONTP_OPT2_AS_GENE *pPtr)
{
    sfixn d1, d2, d3, k;
    SCUBE *scube;
    preFFTRep *R, *gcdpoly, *F3, *zpoly, *Sk; 
    sfixn *gcd, *gcd2, dgcd, dgcd2;
    LinkedQueue *resQueue = NULL;

    // the case that g is nonzero constant
    if (!zeroPolyp(g) && constantPolyp(g)) {
        if (DEBUG) printf("g is a nonzero constant\n");
        return EX_LinkedQueue_Init();
    }
    
    // the case that both F1 and F2 are univariate polynomials in x
    if (N(F1) == N(F2) && N(F2) == 1) { 
        if (DEBUG) printf("both F1 and F2 are univariate polynomials\n");
        d1 = shrinkDeg(BUSZSI(F1, 1), DAT(F1), CUMI(F1, 1));
        d2 = shrinkDeg(BUSZSI(F2, 1), DAT(F2), CUMI(F2, 1));

        resQueue = EX_LinkedQueue_Init();
        gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(F1), d1, DAT(F2), d2, pPtr);
        if (dgcd == 0) {
            my_free(gcd);
            return resQueue;  
        } 

        if (zeroPolyp(g)) {
            gcdpoly = CreateUniPoly(dgcd, gcd);
            EX_LinkedQueue_Enqeue(resQueue, 
                EX_RegularChain2_Init(gcdpoly, NULL));
            my_free(gcd);
            return resQueue;
        }
        
        d3 = shrinkDeg(BUSZSI(g, 1), DAT(g), CUMI(g, 1));
        gcd2 = EX_GCD_Uni_Wrapper(&dgcd2, gcd, dgcd, DAT(g), d3, pPtr);
        
        if (dgcd2 == 0) {
            my_free(gcd);
            my_free(gcd2);
            return resQueue;  
        }

        gcdpoly = CreateUniPoly(dgcd2, gcd2);
        EX_LinkedQueue_Enqeue(resQueue, EX_RegularChain2_Init(gcdpoly, NULL));

        my_free(gcd);
        my_free(gcd2);
        return resQueue;
    }

    if (N(F1) == 1) { 
        if (DEBUG) printf("F1 is a univariate polynomial\n");
        F3 = add_poly_uni_bi(F1, F2, pPtr->P);
        resQueue = modular_solve2(F3, F2, g, pPtr);
        EX_FreeOnePoly(F3);
        return resQueue;
    }

    if (N(F2) == 1) { 
        if (DEBUG) printf("F2 is a univariate polynomial\n");
        F3 = add_poly_uni_bi(F2, F1, pPtr->P);
        resQueue = modular_solve2(F1, F3, g, pPtr);
        EX_FreeOnePoly(F3);
        return resQueue;
    }
    
    // The case both F1 and F2 are bivariate
    if (DEBUG) printf("F1 and F2 are bivariate polynomials\n");
    d1 = shrinkDeg(BUSZSI(F1, 2), DAT(F1), CUMI(F1, 2));
    d2 = shrinkDeg(BUSZSI(F2, 2), DAT(F2), CUMI(F2, 2));
    if (d1 < d2) { return modular_solve2(F2, F1, g, pPtr); }
    
    // printf("d1 = %d, d2 = %d\n", d1, d2);
    // EX_Poly_Print(F1);
    // EX_Poly_Print(F2);
    
    // build the scube
    scube = EX_SubResultantChainOpt(F1, F2, 2, pPtr);
    R = EX_IthDthCoeff(0, 0, scube, pPtr);
    if (zeroPolyp(R)) { 
        if (DEBUG) printf("The resultant is zero\n"); 
        // S_k is the last nonzero regular subresultant of F1 and F2;
        k = next_regular_subres_index(1, scube, pPtr);
        zpoly = CreateZeroPoly();

        if (k == EX_WidthInSCUBE(scube)) {
            resQueue = modular_solve2(F2, g, zpoly, pPtr);
        } else {
            Sk = EX_IthSubres(k, scube, pPtr);
            resQueue = modular_solve2(Sk, g, zpoly, pPtr);
        }
        EX_SCUBE_Free(scube);
        EX_FreeOnePoly(zpoly);
        return resQueue;
    }

    // generic case
    resQueue = modular_solve2_inner(F1, F2, g, scube, pPtr);
    //EX_SCUBE_Print(scube);
    EX_SCUBE_Free(scube);
    return resQueue;
}


/*
return Quo(A,B)
The computation is done on GPU no matter how small or big the polynomials are.
@dG: degree of Quo(A,B)
@A: univariate polynomial 
@dA: degree of A
@B: univariate polynomial 
@dB: degree of B
@pPtr : the prime for finite field.
*/


sfixn* divCUDAStraight(sfixn *dG, sfixn dA, sfixn *A, sfixn dB, sfixn *B, MONTP_OPT2_AS_GENE *pPtr)
{
	sfixn dAA = shrinkDegUni(dA, A);
	sfixn dBB = shrinkDegUni(dB, B);		
	dA = dAA; 
	dB = dBB;
	if(dB > dA || CUDA_TAG == 0 )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);
	if( (dA <= 100) || (dB <= 100) )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);

	sfixn *R, *Q;
	Q = (sfixn *) my_malloc((dA - dB + 1)*sizeof(sfixn));
	R = (sfixn *) my_malloc(dB*sizeof(sfixn));
		//
		if(PrintSpecial == 1)
		{
			printf("A: ");	
			EX_printPolyUni(dA, A, 'x');
			printf("B: ");
			EX_printPolyUni(dB, B, 'x');
		}
		//
	float time = divPrimeField(A, B, R, Q, dA+1, dB+1, pPtr->P);
	//The following line is for avoiding compiler warnings
	time += 1.0;
		//
		if(PrintSpecial == 1)
		{
			printf("Q: ");	
			EX_printPolyUni(dA-dB, Q, 'x');
			printf("R: ");
			EX_printPolyUni(dB-1, R, 'x');
		}
		//
		//
		if(PrintSpecial == 1){
		printf("Returning from CUDA based Quo..\n");	
		}
		//

			     
	sfixn i;
	for(i = dA-dB; i >= 0; --i)
	{
			//
			if(PrintSpecial == 1)
			{		
				printf(" %d =%d ",i,Q[i]);
			}
			//
		if(Q[i] != 0)
			break;
	}
			//
			if(PrintSpecial == 1)
			{		
				printf("\n");
			}
			//
	if(i == -1) (*dG) = 0;
	else (*dG) = i;
	my_free(R);
	return Q;

}


/*
return Quo(A,B)
The computation is done on GPU if the degree 
of A or B is in the range given by threshold_divLow and threshold_divHigh.
@dG: degree of Quo(A,B)
@A: univariate polynomial 
@dA: degree of A
@B: univariate polynomial 
@dB: degree of B
@pPtr : the prime for finite field.
*/


sfixn* divCUDA(sfixn *dG, sfixn dA, sfixn *A, sfixn dB, sfixn *B, MONTP_OPT2_AS_GENE *pPtr)
{
	sfixn dAA = shrinkDegUni(dA, A);
	sfixn dBB = shrinkDegUni(dB, B);		
	dA = dAA; 
	dB = dBB;
	if( CUDA_TAG == 0 || dA < dB || dB <= 100 ||  dA <= 100  )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);
	if( dB < threshold_div )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);

	if( (dA - dB) > threshold_divDiff )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);

	if( dB  > threshold_divHigh || dA > threshold_divHigh )
		return EX_UniQuo(dG, dA, A, dB, B, pPtr);


	
		//
		if(PrintSpecial == 1){
		printf("CUDA based Quo is calling..\n");	
		}
		//
	
	sfixn *R, *Q;
	Q = (sfixn *) my_malloc((dA - dB + 1)*sizeof(sfixn));
	R = (sfixn *) my_malloc(dB*sizeof(sfixn));
		//
		if(PrintSpecial == 1)
		{
			printf("A: ");	
			EX_printPolyUni(dA, A, 'x');
			printf("B: ");
			EX_printPolyUni(dB, B, 'x');
		}
		//
	float time = divPrimeField(A, B, R, Q, dA+1, dB+1, pPtr->P);
	//The following line is for avoiding compiler warnings
	time += 1.0;
		//
		if(PrintSpecial == 1)
		{
			printf("Q: ");	
			EX_printPolyUni(dA-dB, Q, 'x');
			printf("R: ");
			EX_printPolyUni(dB-1, R, 'x');
		}
		//
		//
		if(PrintSpecial == 1){
		printf("Returning from CUDA based Quo..\n");	
		}
		//

			     
	sfixn i;
	for(i = dA-dB; i >= 0; --i)
	{
			//
			if(PrintSpecial == 1)
			{		
				printf(" %d =%d ",i,Q[i]);
			}
			//
		if(Q[i] != 0)
			break;
	}
			//
			if(PrintSpecial == 1)
			{		
				printf("\n");
			}
			//
	if(i == -1) (*dG) = 0;
	else (*dG) = i;
	my_free(R);
	return Q;
}


/*
return gcd(A,B)
The computation is done on GPU no matter how small or big the polynomials are.
@dG: degree of gcd(A,B)
@A: univariate polynomial 
@dA: degree of A
@B: univariate polynomial 
@dB: degree of B
@pPtr : the prime for finite field.

*/
sfixn* gcdCUDAStraight(sfixn *dG, sfixn *A, sfixn dA, sfixn *B, sfixn dB, MONTP_OPT2_AS_GENE *pPtr)
{
	if( CUDA_TAG == 0 )
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);

	sfixn dAA = shrinkDegUni(dA, A);
	sfixn dBB = shrinkDegUni(dB, B);		
	dA = dAA; 
	dB = dBB;

	if( (dA <= 100) || (dB <= 100) )
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);


	sfixn *result = NULL;
	sfixn min = dA, max = dB, flag = 0;

	


	if(dA <= dB) {min = dA; max = dB; flag = 0;}
	else         {min = dB; max = dA; flag = 1;}
	if(min == 0)
	{	
		if(flag == 1 && B[0] == 0) result = (sfixn *) my_malloc((max+1)*sizeof(sfixn));
		if(flag == 0 && A[0] == 0) result = (sfixn *) my_malloc((max+1)*sizeof(sfixn));
	}
	else         result = (sfixn *) my_malloc((min+1)*sizeof(sfixn));
	sfixn i;
	//for(i = 0; i <= min; ++i) result[i] = 0;
	float time;
	//
		if(PrintSpecial == 1){
		printf("CUDA based GCD is calling..\n");
		printf("A: ");	
		EX_printPolyUni(dA, A, 'x');
		printf("B: ");
		EX_printPolyUni(dB, B, 'x');
	
		//printf("min: %d, initializing: ", min);
		//EX_printPolyUni(min, result, 'x');
		}
		//		


	if(min == dB)	time = gcdPrimeField(&i, A, B, result, dA+1, dB+1, pPtr->P, 1);
	else time = gcdPrimeField(&i, B, A, result, dB+1, dA+1, pPtr->P, 1);
	//The following line is for avoiding compiler warnings
	time += 1.0;
		//
		if(PrintSpecial == 1){
		printf("A: ");	
		EX_printPolyUni(dA, A, 'x');
		printf("B: ");
		EX_printPolyUni(dB, B, 'x');

		printf("returning GPU gcd\n");
		printf("GCD(A,B): ");
		EX_printPolyUni(i, result, 'x');
		}
		//
	/*
	for(i = min; i >= 0; --i)
	{
		if(result[i] != 0 )
			break;
	}
	if(i == -1) (*dG) = 0;
	else 	(*dG) = i;
	*/
	(*dG) = i;
	return result;
	
}


/*
return gcd(A,B)
The computation is done on GPU if the degree 
of A or B is is in the range given by threshold_gcdLow and threshold_gcdHigh.
@dG: degree of gcd(A,B)
@A: univariate polynomial 
@dA: degree of A
@B: univariate polynomial 
@dB: degree of B
@pPtr : the prime for finite field.

*/


sfixn* gcdCUDA(sfixn *dG, sfixn *A, sfixn dA, sfixn *B, sfixn dB, MONTP_OPT2_AS_GENE *pPtr)
{
	sfixn dAA = shrinkDegUni(dA, A);
	sfixn dBB = shrinkDegUni(dB, B);		
	dA = dAA; 
	dB = dBB;

	if( CUDA_TAG == 0 )
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);
	if( (dA <= 100)||(dB <= 100) )
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);
	if( (dA >= threshold_gcdHigh)||(dB >= threshold_gcdHigh) )
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);
	


	sfixn levelA = 0, levelB = 0;

	if(dA >= gcdLowLevel1 && dA <= gcdHighLevel1) levelA = 1;
	else if(dA >= gcdLowLevel2 && dA <= gcdHighLevel2) levelA = 2;
	else if(dA >= gcdLowLevel3 && dA <= gcdHighLevel3) levelA = 3;


	if(dB >= gcdLowLevel1 && dB <= gcdHighLevel1) levelB = 1;
	else if(dB >= gcdLowLevel2 && dB <= gcdHighLevel2) levelB = 2;
	else if(dB >= gcdLowLevel3 && dB <= gcdHighLevel3) levelB = 3;
	
	if(levelA != levelB || levelA == 0)
		return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);





        sfixn *result;
        sfixn min = dA;
        if(dA < dB) min = dA;
        else        min = dB;
        result = (sfixn *) my_malloc((min+1)*sizeof(sfixn));
        sfixn i;
        //for(i = 0; i <= min; ++i) result[i] = 0;
        float time;
        //	
		if(PrintSpecial == 1){
		printf("CUDA based GCD is calling..\n");
		printf("A: ");	
		EX_printPolyUni(dA, A, 'x');
		printf("B: ");
		EX_printPolyUni(dB, B, 'x');
	
		//printf("min: %d, initializing: ", min);
		//EX_printPolyUni(min, result, 'x');
		}
		//		


	if(min == dB)	time = gcdPrimeField(&i, A, B, result, dA+1, dB+1, pPtr->P, 1);
	else time = gcdPrimeField(&i, B, A, result, dB+1, dA+1, pPtr->P, 1);
	//The following line is for avoiding compiler warnings
	time += 1.0;
		//
		if(PrintSpecial == 1){
		printf("returning GPU gcd\n");
		printf("GCD(A,B): ");
		EX_printPolyUni(i, result, 'x');
		}
		//
	/*
	for(i = min; i >= 0; --i)
	{
		if(result[i] != 0 )
			break;
	}
	if(i == -1) (*dG) = 0;
	else 	(*dG) = i;
	*/
	(*dG) = i;
	return result;
	

	//return EX_GCD_Uni_Wrapper(dG, A, dA, B, dB, pPtr);
}


/*
This function computes the inverse Kronocker substitution
input: one univariate polynomial along with its degree and an integer
output: a bivariate polynomial 

@A: The output univariate polynomial.
@dA: degree of A
@e: for modulas opeartion in inverse kronocker substituition.
*/
preFFTRep *invKroneck( sfixn *A, sfixn dA, sfixn e)
{
	// The A should have inner degree e-1
	sfixn dgs[3] = {0, e-1, dA/e};
	preFFTRep *B = EX_InitOnePoly(2, dgs); // creating B
	sfixn i, j, k;
	for(i = 0; i < ((BUSZSI(B, 2)+1) * e); ++i) 
		DATI(B, i) = 0; //initializing B
	for(i = 0; i <= dA; ++i)
	{
		if(A[i] != 0)
		{
			j = i/e;
			k = i%e;
			
			DATI(B,e*j+k) = A[i];
		}
	}
		
	
	return B;
}


/*
This function compute the Kronocker substitution of bivariate
polynomial B by putting x1 = x and x2 = x^Inc.
Where x2 is the higher variable.
The output univariate polynomial is stored by pointer A.
The degree of this univariate polynomial is computed and return.
@A: The output univariate polynomial.
@B: the input bivariate polynomial.
@Inc: Substitution by x2 = x^Inc.
*/
sfixn KroneckerSubs(sfixn *A, polyBothInHostDev *B, sfixn Inc)
{
	sfixn length = Inc*(BUSZSI(B->polyHost, 2)+1); // the length of A
		//
		if(PrintFlag == 1){
			printf("degrees of B: %d %d\n",BUSZSI(B->polyHost, 1),BUSZSI(B->polyHost, 2));
			printf("increment to %d\n", Inc);
			printf("so the length is : %d\n", length);
		}
		//
	
	sfixn i, j;
	for(i = 0; i < length; ++i) A[i] = 0; //initialization of A
	for(i = 0; i <= BUSZSI(B->polyHost, 2); ++i)
	{
		for(j = 0; j <= BUSZSI(B->polyHost, 1); ++j)
		{
			A[(Inc*i) + j] = DATI(B->polyHost, ((BUSZSI(B->polyHost, 1)+1)*i)+j );
				//
				if(PrintFlag == 1){
					printf("  (%d  %d)   ",(Inc*i) + j,((BUSZSI(B->polyHost, 1)+1)*i)+j);
				}
				//
			
		}
		
	}
			//
			if(PrintFlag == 1){
				printf("\n");
			}
			//
			//
		

	return shrinkDegUni(length-1, A);// deleting leading 0s	
}

/*
Compute A = B/C, where B is bivariate and C is univariate in x.
and each coefficient of B is divisible by C.
@method: for future use if GPU based method is choosen.
@A: The output biivariate polynomial.
@B: the input bivariate polynomial.
@C: the input univariate polynomial.
@pPtr : the prime for finite field.
*/
void bivarDivComponentwise(computingChoince method, polyBothInHostDev *A,  polyBothInHostDev *B, polyBothInHostDev *C, MONTP_OPT2_AS_GENE *pPtr)
{
		//
		if(PrintFlag == 1){
			EX_Poly_Print(B->polyHost);
			EX_Poly_Print(C->polyHost);
		}
		//
	//// degree of A should be dgs[1] in xa and dgs[2] in y
	sfixn dgs[3] = {0, BUSZSI(B->polyHost, 1) - BUSZSI(C->polyHost, 1), BUSZSI(B->polyHost, 2)};
	A->polyHost = EX_InitOnePoly(2, dgs);
	sfixn i, j;
	for(i = 0; i < ((BUSZSI(A->polyHost, 2)+1) * (BUSZSI(A->polyHost, 1)+1)); ++i)
		DATI(A->polyHost, i) = 0;  //initializing A as zero polynomial
			/*
			printf("%d %d\n",BUSZSI(A->polyHost, 1), BUSZSI(A->polyHost, 2) );
		
			for(i = 0; i < ((BUSZSI(A->polyHost, 2)+1) * (BUSZSI(A->polyHost, 1)+1)); ++i)
				printf("%d ",DATI(A->polyHost, i)); //DATI(A->polyHost, i) = 0;
			printf("\n");
			*/
	
	sfixn *a, *b, *c, *d, da, db, dc;
	c = DAT(C->polyHost);
	dc = shrinkDegUni(BUSZSI(C->polyHost, 1), c);	
	db = BUSZSI(B->polyHost, 1);
	
	for(i = 0; i <= BUSZSI(A->polyHost, 2); ++i)
	{
		// in ith iteration the code divide  y^i coefficient of B with C and store it in appropriate place of A.
		b = DAT(B->polyHost) + ((BUSZSI(B->polyHost, 1)+1)*i);
		if(method == CPU)		
			a = EX_UniQuo(&da, db, b, dc, c, pPtr);	
		else if(method == GPU) 
			a = divCUDAStraight(&da, db, b, dc, c, pPtr);
		else
			a = divCUDA(&da, db, b, dc, c, pPtr);
			
			/*
			printf("i: %d\n",i);
			EX_printPolyUni(db, b, 'x');
			EX_printPolyUni(dc, c, 'x');
			EX_printPolyUni(da, a, 'x');
			*/
		d = DAT(A->polyHost) + ((BUSZSI(A->polyHost, 1)+1)*i);		

		for(j = 0; j <= da; ++j)
			d[j] = a[j];
		my_free(a);
		
	}
		
}
/*
implementation of plain univariate polynomial multiplication: a = b*c mod p

@a: univariate polynomial
@db: degree of b
@b: univariate polynomial
@dc: degree of c
@c: univariate polynomial
@p: prime
*/
void plainMultiplication( sfixn *a, sfixn db, sfixn *b, sfixn dc, sfixn *c, sfixn p)
{
	sfixn i, j;
	for(i = 0; i <= db; ++i)
	{
		for(j = 0; j <= dc; ++j)
		{
			a[i+j] = AddMod(a[i+j], MulMod(b[i], c[j], p), p);
		}
	}

}

/*
Compute A = B*C, where B are bivariate and C is univariate in x.

@A: The output biivariate polynomial.
@B: the input bivariate polynomial.
@C: the input univariate polynomial.
@pPtr : the prime for finite field.


*/
void bivarMulComponentwise(polyBothInHostDev *A,  polyBothInHostDev *B, polyBothInHostDev *C, MONTP_OPT2_AS_GENE *pPtr)
{
		//
		if(PrintFlag == 1){
			EX_Poly_Print(B->polyHost);
			EX_Poly_Print(C->polyHost);
		}
		//
					//
		if(PrintFlag == 1){
			//EX_Poly_Print(A->polyHost);			
		}
		//

	// degree of A should be dgs[1] in xa and dgs[2] in y
	sfixn dgs[3] = {0, BUSZSI(B->polyHost, 1) + BUSZSI(C->polyHost, 1), BUSZSI(B->polyHost, 2)};
	A->polyHost = EX_InitOnePoly(2, dgs);

	sfixn i;
	for(i = 0; i < ((BUSZSI(A->polyHost, 2)+1) * (BUSZSI(A->polyHost, 1)+1)); ++i) DATI(A->polyHost, i) = 0; //initializing A as zero polynomial
			//
			if(PrintFlag == 1){
				printf("initialize multiplication result: ");
				EX_Poly_Print(A->polyHost);

			}
			//

	
	sfixn *a, *b, *c, da, db, dc;
	c = DAT(C->polyHost);
	dc = BUSZSI(C->polyHost, 1);
	db = BUSZSI(B->polyHost, 1);
	da = BUSZSI(A->polyHost, 1);

	for(i = 0; i <= BUSZSI(A->polyHost, 2); ++i)
	{
		// in ith iteration the code multiply C with y^i coefficient of B and store it in appropriate place of A.		

		b = DAT(B->polyHost) + ((BUSZSI(B->polyHost, 1)+1)*i);		
		a = DAT(A->polyHost) + ((BUSZSI(A->polyHost, 1)+1)*i);
			//
			if(PrintFlag == 1){
			printf("Before multiplication: \n");
			printf("1st component: ");EX_printPolyUni(db, b, 'x');
			printf("2nd component: ");EX_printPolyUni(dc, c, 'x');				
			printf("resultant component: ");EX_printPolyUni(da, a, 'x');
			}
			//		
			plainMultiplication( a, db, b, dc, c, pPtr->P);
			//
			if(PrintFlag == 1){
			printf("After multiplication:\n ");
			printf("1st component: ");EX_printPolyUni(db, b, 'x');
			printf("2nd component: ");EX_printPolyUni(dc, c, 'x');				
			printf("resultant component: ");EX_printPolyUni(da, a, 'x');
			}
			//
	}
		//
		if(PrintFlag == 1){
			EX_Poly_Print(A->polyHost);			
		}
		//

}

/*
Computing the content of a bivariate polynomial F1.
The result is stored in G.
@method: for future use if GPU based method is choosen.
@F1: The input bivariate polynomial.
@G: the content of F1
@pPtr : the prime for finite field.
*/
void contentbivariate(computingChoince method,  polyBothInHostDev *G,  polyBothInHostDev *F1,  MONTP_OPT2_AS_GENE *pPtr)
{
	sfixn *a, *b, *c, da, db, dc = -1;
	preFFTRep *g;

	//a contains the coefficient of y^0 of F1. this redundant code is for avoid compiler warnings
	a = DAT(F1->polyHost);
	da = shrinkDegUni(BUSZSI(F1->polyHost, 1), a);	

	if(BUSZSI(F1->polyHost, 2) >= 0)
	{
		//a contains the coefficient of y^0 of F1.
		a = DAT(F1->polyHost);
		da = shrinkDegUni(BUSZSI(F1->polyHost, 1), a);	
		
			//
			if(PrintFlag == 1){
				printf("the first coefficient: \n");
				EX_printPolyUni(da, a, 'x');
			}
			//

	}
	if(BUSZSI(F1->polyHost, 2) >= 1)
	{
		//b contains the coefficient of y^1 of F1.	

		b = DAT(F1->polyHost) + BUSZSI(F1->polyHost, 1) + 1;
		db = shrinkDegUni(BUSZSI(F1->polyHost, 1), b);
		if(method == CPU)
			c = EX_GCD_Uni_Wrapper(&dc, a, da, b, db, pPtr);
		else if(method == GPU)
			c = gcdCUDAStraight(&dc, a, da, b, db, pPtr);				
		else
			c = gcdCUDA(&dc, a, da, b, db, pPtr);		
			//
			if(PrintFlag == 1){
			printf("the second coefficient: \n");
			EX_printPolyUni(db, b, 'x');
			printf("the first gcd: \n");
			EX_printPolyUni(dc, c, 'x');
			}
			//

	}
	else// some how if the degree of F1 in y is negative, it will return the coefficient of y^0 
	{
		g = CreateUniPoly(da, a);
		G->polyHost = g;
		return;		
	}


	//c contains the gcd of the coefficient of y^0 and y^1 in F1.
	sfixn i;
	for(i = 2; i <= BUSZSI(F1->polyHost, 2) && dc >= 0; ++i)
	{
		// in i-th iteration c = gcd(c, coefficient of y^i)
		a = c;
		da = dc;
		b = DAT(F1->polyHost) + (BUSZSI(F1->polyHost, 1) + 1)*i;
		db = shrinkDegUni(BUSZSI(F1->polyHost, 1), b);
		if(method == CPU)
			c = EX_GCD_Uni_Wrapper(&dc, a, da, b, db, pPtr);
		else if(method == GPU)
			c = gcdCUDAStraight(&dc, a, da, b, db, pPtr);
		else
			c = gcdCUDA(&dc, a, da, b, db, pPtr);

			//
			if(PrintFlag == 1){
			EX_printPolyUni(da, a, 'x');
			EX_printPolyUni(db, b, 'x');
			EX_printPolyUni(dc, c, 'x');
			}
			//
		//my_free(a);
	}
	if(dc >= 0)
	{
		g = CreateUniPoly(dc, c);
		G->polyHost = g;
		//my_free(c);
	}
}

/*
@method:=  flag for computing option between CPU and GPU. see Types.h for details.
@f1:= f1 is a univariate polynomial of x,  where coefficients are in a field.
@f2:= f2 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@pPtr:= it is a prime number that defines the coefficient field.
@rc:= The regular chain that will be returned by this function. 

This function checks whether f1 and f2 are regular or not.
If it is regular then it is added to rc.
Otherwise it ( ( f1/gcd(f1, lc(f2))), f2 ) is added to rc
and the function make a recursive call with gcd(f1, lc(f2)) and tail(f2)

*/

void irregularRCinsertion(computingChoince method, preFFTRep *f1, preFFTRep *f2, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc){
        sfixn realDeg2 = shrinkDeg(BUSZSI(f2, 2), DAT(f2), BUSZSI(f2, 1)+1);
	BUSZSI(f2, 2) = realDeg2;
	realDeg2 = shrinkDegUni(BUSZSI(f1, 1), DAT(f1));
	BUSZSI(f1, 1) = realDeg2;
	
	if(BUSZSI(f1, 1) == 0){
		if( DATI(f1,0) == 0) EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(EX_CopyOnePoly(f1), EX_CopyOnePoly(f2)));	
		return;
	}
		
	if(BUSZSI(f2, 2) == 0){ // if f2 is univariate in x	
		realDeg2 = shrinkDegUni(BUSZSI(f2, 1), DAT(f2));
		BUSZSI(f2, 1) = realDeg2;

			

		sfixn dgcd, *gcd;
		if(method == CPU)     gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(f2),  BUSZSI(f2, 1),  DAT(f1),  BUSZSI(f1, 1), pPtr);
		else if(method == GPU)gcd = gcdCUDAStraight(&dgcd, DAT(f2),  BUSZSI(f2, 1),  DAT(f1),  BUSZSI(f1, 1), pPtr);
		else	              gcd = gcdCUDA(&dgcd, DAT(f2),  BUSZSI(f2, 1),  DAT(f1),  BUSZSI(f1, 1), pPtr);
		realDeg2 = shrinkDegUni(dgcd, gcd);
		
		if(realDeg2 > 0) EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(NULL, CreateUniPoly(dgcd, gcd) ));		
		my_free(gcd);
		return;
	}	

	preFFTRep *lc2 = EX_getInitial(f2);
	realDeg2 = shrinkDegUni(BUSZSI(lc2, 1), DAT(lc2));
	BUSZSI(lc2, 1) = realDeg2;

	sfixn dg, *g;
	if(method == CPU)	g = EX_GCD_Uni_Wrapper(&dg, DAT(lc2), BUSZSI(lc2, 1), DAT(f1),  BUSZSI(f1, 1), pPtr);
	else if(method == GPU)	g = gcdCUDAStraight(&dg, DAT(lc2), BUSZSI(lc2, 1), DAT(f1),  BUSZSI(f1, 1), pPtr);
	else			g = gcdCUDA(&dg, DAT(lc2), BUSZSI(lc2, 1), DAT(f1),  BUSZSI(f1, 1), pPtr);
	realDeg2 = shrinkDegUni(dg, g); 
	dg = realDeg2;

	EX_FreeOnePoly(lc2);
	if(dg == 0){		
		EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(EX_CopyOnePoly(f1), EX_CopyOnePoly(f2)));
		my_free(g);
		return;
	}
	
	sfixn dd, *d;
	if(method == CPU)	d = EX_UniQuo(&dd, BUSZSI(f1, 1), DAT(f1),  dg, g, pPtr);
	else if(method == GPU)  d = divCUDAStraight(&dd, BUSZSI(f1, 1), DAT(f1),  dg, g, pPtr);
	else			d = divCUDA(&dd, BUSZSI(f1, 1), DAT(f1),  dg, g, pPtr);
	realDeg2 = shrinkDegUni(dd, d); 
	dd = realDeg2;

	if(dd > 0) EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init( CreateUniPoly(dd, d), EX_CopyOnePoly(f2)));			
		
	my_free(d);
	irregularRCinsertion(method, CreateUniPoly(dg, g), EX_GetPolyTail(f2) , pPtr, rc);
	my_free(g);	
}


/*
@method:=  flag for computing option between CPU and GPU. see Types.h for details.
@F1:= F1 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@F2:= F2 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@g:= g is a univariate polynomial in x, where coefficients are in a field.
@S0:= the first nonzero polynomial in the subresultant of F1 and F2. Note that, F1 and F2 have positive degree in y such that they have zero resultant in y.
@h := h is a univariate polynomial in x, where coefficients are in a field. h(x) = gcd(lc(F1), lc(F2))(x) != 0 
@pPtr:= it is a prime number that defines the coefficient field.
@rc:= The regular chain that will be returned by this function. 

Output:
On output this function add triaangular sets  {(A_1, B_1),...(A_e, B_e)} in rc
such that,
V(F1,F2,g) = V(A_1,B_1),..V(A_e,B_e).

For Details, please see 
Marc Moreno Maza and Wei Pan, Solving Bivariate Polynomial Systems on a GPU,
J. Phys.: Conf. Ser. 341 012022


*/

void ModularGenericSolve2ZeroResultant(computingChoince method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, polyBothInHostDev *h, polyBothInHostDev *S0, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc)
{		
	polyBothInHostDev *G1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G1->status = -1;
	contentbivariate( method, G1, F1, pPtr); // G1 = content of F1 is computed
			//
			if(PrintFlag == 1){
				printf("cont(F1): \n");
				EX_printPolyUni(BUSZSI(G1->polyHost, 1), DAT(G1->polyHost), 'x');			
			}
			//
	sfixn *a1, da1;
	//if(BUSZSI(G1->polyHost, 1) > 0)
	//{
	polyBothInHostDev *G2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G2->status = -1;
	contentbivariate( method, G2, F2, pPtr); // G2 = content of F2 is computed
			//
			if(PrintFlag == 1){
				printf("cont(F2): \n");
				EX_printPolyUni(BUSZSI(G2->polyHost, 1), DAT(G2->polyHost), 'x');			
			}
			//
		//a1 = gcd(G1,G2)
	if(method == CPU) 
		a1 = EX_GCD_Uni_Wrapper(&da1, DAT(G1->polyHost), BUSZSI(G1->polyHost, 1),DAT(G2->polyHost), BUSZSI(G2->polyHost, 1),pPtr);
	else if(method == GPU)
		a1 = gcdCUDAStraight(&da1, DAT(G1->polyHost), BUSZSI(G1->polyHost, 1),DAT(G2->polyHost), BUSZSI(G2->polyHost, 1),pPtr);	
	else
		a1 = gcdCUDA(&da1, DAT(G1->polyHost), BUSZSI(G1->polyHost, 1),DAT(G2->polyHost), BUSZSI(G2->polyHost, 1),pPtr);	
	EX_FreeOnePoly(G1->polyHost);
	EX_FreeOnePoly(G2->polyHost);  
	//}
	//else
	//{
	
	//	da1 = 0;
	//	a1 = (sfixn *) my_malloc(sizeof(sfixn));
	//	a1[0] = 1;
	//	EX_FreeOnePoly(G1->polyHost);
	//}
			//
			if(PrintFlag == 1){
				printf("gcd of cont(F1) and cont(F2):\n ");
				EX_printPolyUni(da1, a1, 'x');			
			}
			//
	//gcd of G1 and G2 is stored in a1. This is the gcd(F1,F2)
	preFFTRep *g1;
	g1 = CreateUniPoly(da1, a1);
	my_free(a1);
	polyBothInHostDev *G5 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G5->status = -1;
	G5->polyHost = g1;
	//G5 = a1
	polyBothInHostDev *G3 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G3->status = -1;
	contentbivariate( method, G3, S0, pPtr); // G3 = content of S0
			//
			if(PrintFlag == 1){
				printf("partial degrees of S0: %d %d\n",BUSZSI(S0->polyHost, 1),BUSZSI(S0->polyHost, 2));
				printf("cont(S) : \n");
				EX_printPolyUni(BUSZSI(G3->polyHost, 1), DAT(G3->polyHost), 'x');			
			}
			//

	polyBothInHostDev *G4;
	G4 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G4->status = -1;
			//
			if(PrintFlag == 1){
				printf("componentwise division:\n");
			}
			//
	bivarDivComponentwise( method, G4,  S0, G3, pPtr); // G4 = S0 is divided by G3 (content of S0) component wise.
	EX_FreeOnePoly(G3->polyHost);
			//
			if(PrintFlag == 1){
				EX_Poly_Print(G4->polyHost);
			}
			//
			

	polyBothInHostDev *G6 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G6->status = -1;
			//
			if(PrintFlag == 1){
				printf("componentwise multiplication\n");
			}
			//
	bivarMulComponentwise(G6,  G4, G5, pPtr); //G6 = componentwise multiplication of G4 by G5
	EX_FreeOnePoly(G4->polyHost);
	EX_FreeOnePoly(G5->polyHost);
	//Finally, G6 contains the gcd(F1, F2)
	// G6 is computed as follows: component_wise_mul( gcd( content(F1), content(F2) ),   componentwise_division(S0,content(S0)) )
			//
			if(PrintFlag == 1){
				printf("Adding g and G into regular chain (paper)... \n");
				EX_Poly_Print(g->polyHost);			
				EX_Poly_Print(G6->polyHost);			
			}
			//
	//Adding g and G6 into regular chain
	irregularRCinsertion(method, g->polyHost, G6->polyHost, pPtr, rc);
			//
			if(PrintFlag == 1){
				printf("done \n");
			}
			//
	
	sfixn *a2, da2;
	sfixn length = (BUSZSI(F1->polyHost, 1)+1)*(BUSZSI(G6->polyHost, 2)+1);
		//
		if(PrintFlag == 1){
			printf("inner degree of F1: %d outer degree of G (paper): %d total length: %d\n",BUSZSI(F1->polyHost, 1), BUSZSI(G6->polyHost, 2), length );
		}
	
	a2 = (sfixn *)my_malloc(length*sizeof(sfixn));
	da2 = KroneckerSubs(a2, G6, BUSZSI(F1->polyHost, 1)+1); //a2= Kronecker substituition of G6 considering the lower degree of F1
			//
			if(PrintFlag == 1){
				printf("Kronocker substituition\n");
				EX_Poly_Print(F1->polyHost);
				EX_Poly_Print(G6->polyHost);						
				EX_printPolyUni(da2, a2, 'x');
			}
			//
	sfixn *a4, da4;	
	length = (BUSZSI(F1->polyHost, 1)+1)*(BUSZSI(F1->polyHost, 2)+1);
	// F1 is already represented as Kronecker substituition
	// so we can do a univariate division  F1 div a2
	if(method == CPU)
		a4 = EX_UniQuo(&da4,  shrinkDegUni(length-1, DAT(F1->polyHost)), DAT(F1->polyHost), da2, a2, pPtr);
	else if(method == GPU) 
		a4 = divCUDAStraight(&da4,  shrinkDegUni(length-1, DAT(F1->polyHost)), DAT(F1->polyHost), da2, a2, pPtr);
	else
		a4 = divCUDA(&da4,  shrinkDegUni(length-1, DAT(F1->polyHost)), DAT(F1->polyHost), da2, a2, pPtr);
	preFFTRep *A4 = invKroneck( a4, da4, BUSZSI(F1->polyHost, 1)+1); // inverse kronecker substituition
		//
			if(PrintFlag == 1){			
				EX_printPolyUni(da4, a4, 'x');
				EX_Poly_Print(A4);
			}
			//
	// so A4 = F1 div G6
	my_free(a2);
	my_free(a4); 

	sfixn *a3, da3;
	length = (BUSZSI(F2->polyHost, 1)+1)*(BUSZSI(G6->polyHost, 2)+1);
		//
		if(PrintFlag == 1){
			printf("inner degree of F2: %d outer degree of G (paper): %d total length: %d\n",BUSZSI(F2->polyHost, 1), BUSZSI(G6->polyHost, 2), length );
		}
	a3 = (sfixn *)my_malloc(length*sizeof(sfixn));
	da3 = KroneckerSubs( a3, G6, BUSZSI(F2->polyHost, 1)+1);//a3= Kronecker substituition of G6 considering the lower degree of F2
			//
			if(PrintFlag == 1){
				printf("Kronocker substituition\n");
				EX_Poly_Print(F2->polyHost);
				EX_Poly_Print(G6->polyHost);						
				EX_printPolyUni(da3, a3, 'x');
			}
			//
	sfixn *a5, da5;
	length = (BUSZSI(F2->polyHost, 1)+1)*(BUSZSI(F2->polyHost, 2)+1);
	// F2 is already represented as Kronecker substituition
	// so we can do a univariate division  F2 div a3
	if(method == CPU)
		a5 = EX_UniQuo(&da5,  shrinkDegUni(length-1, DAT(F2->polyHost)), DAT(F2->polyHost), da3, a3, pPtr);
	else if(method == GPU) 
		a5 = divCUDAStraight(&da5,  shrinkDegUni(length-1, DAT(F2->polyHost)), DAT(F2->polyHost), da3, a3, pPtr);
	else
		a5 = divCUDA(&da5,  shrinkDegUni(length-1, DAT(F2->polyHost)), DAT(F2->polyHost), da3, a3, pPtr);
	preFFTRep *A5 = invKroneck( a5, da5, BUSZSI(F2->polyHost, 1)+1);
			//
			if(PrintFlag == 1){			
				EX_printPolyUni(da5, a5, 'x');
				EX_Poly_Print(A5);
			}
			//
	// so A5 = F2 div G6
	my_free(a3);
	my_free(a5);

	polyBothInHostDev *G7 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G7->status = -1;
	
	polyBothInHostDev *G8 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	G8->status = -1;
	



	// A4 A5
	sfixn *x1, *x2, *x3, *x4;
	sfixn dx1, dx2, dx3, dx4;
	preFFTRep *X4; 

	if(BUSZSI(A4, 2) == 0)
	{
			//
			if(PrintFlag == 1){
				printf("F1 became univariate\n");
			}			
			//
		x1 = SquarefreePart(&dx1, shrinkDegUni(BUSZSI(A4, 1), DAT(A4)), DAT(A4), pPtr); // if A4 is univariate in x, x1 = sqrfreepart(A4)
			//
			if(PrintFlag == 1){
				printf("F1:= squareFree(F1):\n");
				EX_printPolyUni(dx1, x1, 'x');
			}			
			//

		if(method == CPU)
			x2 = EX_GCD_Uni_Wrapper(&dx2, x1, dx1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);	//x2 = gcd(x1,h)
		else if(method == GPU)
			x2 = gcdCUDAStraight(&dx2, x1, dx1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);
		else
			x2 = gcdCUDA(&dx2, x1, dx1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);
			//
			if(PrintFlag == 1){
				printf("F1 quo h:\n");
				EX_printPolyUni(dx2, x2, 'x');
			}			
			//
		if(method == CPU)	
			x3 = EX_UniQuo(&dx3, dx1, x1, dx2, x2, pPtr);
		else if(method == GPU) 
			x3 = divCUDAStraight(&dx3, dx1, x1, dx2, x2, pPtr);
		else
			x3 = divCUDA(&dx3, dx1, x1, dx2, x2, pPtr);
		//x3 = x1 quo x2
		if(method == CPU)
			x4 = EX_GCD_Uni_Wrapper(&dx4, x3, dx3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);
		else if(method == GPU)
			x4 = gcdCUDAStraight(&dx4, x3, dx3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);	
		else
			x4 = gcdCUDA(&dx4, x3, dx3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);	
		//x4 = gcd(g,x3)	
		my_free(x1);	
		my_free(x2);
		my_free(x3);
		if(dx4 == 0)
		{
			my_free(x4); //if x4 is a constant or in other words if F1 becomes constant 
			return;
		}
		else
		{
			EX_FreeOnePoly(A4);
			X4 = CreateUniPoly(dx4, x4);
			my_free(x4);
			G7->polyHost = X4;
		}
	}
	else
	{
		G7->polyHost = A4;	
	}
	// So G7 contains either F1 quo gcd(F1,F2) or 
	//  the successive of following operations (if F1 is univariate in x): F1 = F1 quo gcd(F1,F2), F1=squarefreepart(F1), F1= F1 quo gcd(F1,h), F1 = gcd(F1,g)
	sfixn *y1, *y2, *y3, *y4;
	sfixn dy1, dy2, dy3, dy4;
	preFFTRep *Y4;
	
	if(BUSZSI(A5, 2) == 0)
	{
			//
			if(PrintFlag == 1){
				printf("F2 became univariate\n");
				
			}			
			//
		y1 = SquarefreePart(&dy1, shrinkDegUni(BUSZSI(A5, 1), DAT(A5)), DAT(A5), pPtr);// if A5 is univariate in x, y1 = sqrfreepart(A5)
			//
			if(PrintFlag == 1){
				printf("F2:= squareFree(F2):\n");
				EX_printPolyUni(dy1, y1, 'x');
			}			
			//

		if(method == CPU)
			y2 = EX_GCD_Uni_Wrapper(&dy2, y1, dy1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);//y2 = gcd(y1,h)
		else if(method == GPU)
			y2 = gcdCUDAStraight(&dy2, y1, dy1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);
		else
			y2 = gcdCUDA(&dy2, y1, dy1, DAT(h->polyHost), BUSZSI(h->polyHost, 1),pPtr);
			//
			if(PrintFlag == 1){
				printf("F2 quo h:\n");
				EX_printPolyUni(dy2, y2, 'x');
			}			
			//
		//y3 = y1 quo y2
		if(method == CPU)
			y3 = EX_UniQuo(&dy3, dy1, y1, dy2, y2, pPtr);
		else if(method == GPU) 
			y3 = divCUDAStraight(&dy3, dy1, y1, dy2, y2, pPtr);		
		else
			y3 = divCUDA(&dy3, dy1, y1, dy2, y2, pPtr);
		
		if(method == CPU)
			y4 = EX_GCD_Uni_Wrapper(&dy4, y3, dy3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);
		else if(method == GPU)
			y4 = gcdCUDAStraight(&dy4, y3, dy3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);
		else
			y4 = gcdCUDA(&dy4, y3, dy3, DAT(g->polyHost), BUSZSI(g->polyHost, 1),pPtr);
		//y4 = gcd(g,y3)
		my_free(y1);	
		my_free(y2);
		my_free(y3);
		if(dy4 == 0)
		{
			my_free(y4); //if y4 is a constant or in other words if F2 becomes constant 
			return;
		}
		else
		{
			EX_FreeOnePoly(A5);
			Y4 = CreateUniPoly(dy4, y4);
			my_free(y4);
			G8->polyHost = Y4;
		}			
	}
	else
	{
		G8->polyHost = A5;
	}
	// So G8 contains either F2 quo gcd(F1,F2) or 
	//  the successive of following operations (if F2 is univariate in x): F2 = F2 quo gcd(F1,F2), F2=squarefreepart(F2), F2= F2 quo gcd(F2,h), F2 = gcd(F2,g)
	if(BUSZSI(G7->polyHost, 2) == 0 && BUSZSI(G8->polyHost, 2) == 0) // if both F1 and F2 are univariate in x
	{
		if(PrintFlag == 1)
		{
 			printf("both F1 and F2 are univariate so returning.. x\n");
		} 
		return; 
	}
	
	if(BUSZSI(G7->polyHost, 2) == 0)
	{
		if(PrintFlag == 1)
		{
 			printf("F1 becomes univariate so adding (F1, F2) into rc.. \n");
		} 
		irregularRCinsertion(method, G7->polyHost, G8->polyHost, pPtr, rc); // if F1 becomes univariate in x
		return;
	}
	
	if(BUSZSI(G8->polyHost, 2) == 0)
	{
		if(PrintFlag == 1)
		{
 			printf("F2 becomes univariate so adding (F2,F1) into rc.. \n");
		} 
		irregularRCinsertion(method, G8->polyHost, G7->polyHost, pPtr, rc); // if F2 becomes univariate in x
		return;
	}
	
	EX_FreeOnePoly(G6->polyHost);

	sfixn rDeg7 = shrinkDeg(BUSZSI(G7->polyHost, 2), DAT(G7->polyHost), BUSZSI(G7->polyHost, 1)+1);
	BUSZSI(G7->polyHost, 2) = rDeg7;
	
	sfixn rDeg8 = shrinkDeg(BUSZSI(G8->polyHost, 2), DAT(G8->polyHost), BUSZSI(G8->polyHost, 1)+1);
	BUSZSI(G8->polyHost, 2) = rDeg8;

	if(rDeg7 >= rDeg8) ModularGenericSolve2(method, G7, G8, g, h, pPtr, rc);
	else ModularGenericSolve2(method, G8, G7, g, h, pPtr, rc);

	EX_FreeOnePoly(G7->polyHost);
	EX_FreeOnePoly(G8->polyHost);
}
/*
@method:=  flag for computing option between CPU and GPU. see Types.h for details.
@F1:= F1 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@F2:= F2 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@g:= g is a univariate polynomial in x, where coefficients are in a field.
@h := h is a univariate polynomial in x, where coefficients are in a field. h(x) = gcd(lc(F1), lc(F2))(x) != 0 
@pPtr:= it is a prime number that defines the coefficient field.
@rc:= The regular chain that will be returned by this function. 

Output:
On output this function add triaangular sets  {(A_1, B_1),...(A_e, B_e)} in rc
such that,
V(F1,F2,g) = V(A_1,B_1),..V(A_e,B_e).

For Details, please see 
Marc Moreno Maza and Wei Pan, Solving Bivariate Polynomial Systems on a GPU,
J. Phys.: Conf. Ser. 341 012022


*/

void ModularGenericSolve2(computingChoince method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, polyBothInHostDev *h, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc)
{
	sfixn *hs, dhs;
	hs = SquarefreePart(&dhs, BUSZSI(h->polyHost, 1), DAT(h->polyHost), pPtr); // hs= squarefree part of h
	sfixn *gs, dgs;
		//
		if(PrintFlag == 1)
		{
			printf("g: ");
			EX_Poly_Print(g->polyHost);
		}
		//
	if (!zeroPolyp(g->polyHost))
	{
		//
		if(PrintFlag == 1)
		{
			printf("g is not a zero poly\n");		
		}
		//
		gs = SquarefreePart(&dgs, BUSZSI(g->polyHost, 1), DAT(g->polyHost), pPtr); // gs = squarefree part of g provided g is not a zero polynomial.
	}
	else
	{ 
		//
		if(PrintFlag == 1)
		{
			printf("g is  a zero poly\n");		
		}
		//
		gs =  (sfixn *) my_malloc(sizeof(sfixn)); gs[0] = 0; dgs = 0;	 // if g is a zero polynomial then its square free part should be also a zero polynomial.
	}
	sfixn *gcd, dgcd;
	//computing gcd(gs,hs)
	if(method == CPU)
		gcd = EX_GCD_Uni_Wrapper(&dgcd, hs, dhs, gs, dgs, pPtr);
	else if(method == GPU)
		gcd = gcdCUDAStraight(&dgcd, hs, dhs, gs, dgs, pPtr);
	else
		gcd = gcdCUDA(&dgcd, hs, dhs, gs, dgs, pPtr);
		//
		if(PrintFlag == 1){
			printf("g: ");
			EX_Poly_Print(g->polyHost);
			printf("degree of gcd: %d and lsb of sqrfree(g): %d\n", dgs, gs[0]);
			printf("CONSIDER all x as x_1\n");
			printf("h = sqfree(h) : ");
			EX_printPolyUni(dhs, hs, 'x');
			printf("g = sqfree(g) : ");
			EX_printPolyUni(dgs, gs, 'x');
			printf("gcd(g,h): ");
			EX_printPolyUni(dgcd, gcd, 'x');
		}
		//
	sfixn *Q, dQ;
	preFFTRep *gNew;
	polyBothInHostDev *Gnew;
	// if gcd(gs,hs) is not a zero poly, we discard the leading 0s to be sure its true degree
	// finally, Gnew = gcd(gs,hs)
	if (!zeroPolyp(g->polyHost))
	{
			//
			if(PrintFlag == 1){
				printf("g is not zero\n");
			}
			//
		if(method == CPU)
			Q = EX_UniQuo(&dQ, dgs, gs, dgcd, gcd, pPtr);
		else if(method == GPU) 
			Q = divCUDAStraight(&dQ, dgs, gs, dgcd, gcd, pPtr);
		else
			Q = divCUDA(&dQ, dgs, gs, dgcd, gcd, pPtr);
		dQ = shrinkDegUni(dQ, Q);	
		gNew = CreateUniPoly(dQ, Q);
					
	}
	else	gNew = CreateZeroPoly();			
	Gnew = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	Gnew->status = -1;
	Gnew->polyHost = gNew;
		//
		if(PrintFlag == 1){
			printf("g = g quo gcd(g,h): ");
			EX_printPolyUni(BUSZSI(gNew, 1), DAT(gNew), 'x');
		}
		//
	my_free(gs);
	my_free(gcd);
	// check Gnew is zero or constant poly.
	if (BUSZSI(Gnew->polyHost, 1) == 0)// && DATI(Gnew->polyHost, 0) != 0 )
	{
				//
				if(PrintFlag == 1){
					printf("g is constant or zero\n");
				}
				//
		if(!zeroPolyp(Gnew->polyHost)) // if Gnew is a constant poly return
		{			
				//
				if(PrintFlag == 1){
					printf("g is constant returning from ModularGenericSolve2 \n");
				}
				//
			my_free(hs);
			EX_FreeOnePoly(gNew); 				    	
			return;
		}
	}
			/*
			if(PrintFlag == 1)
			{
				printf("Subresultant is going to be computed using %d method\n",method);
			}
			*/


	//Gnew is either zero or its degree is atleast 1.
        // sub-resultant chain is computed
	SCUBE *S = EX_SubResultantChainSelect(method, F2->polyHost, F1->polyHost, 2, pPtr);
			//
			if(PrintFlag == 1){
				printf("Subresultant is successfully computed\n");
			}	
			//	
	// S0 is the resultant
	preFFTRep *S0 = EX_IthDthCoeff(0, 0, S, pPtr);
				//
				if(PrintFlag == 1){
					printf("resultant of degree %d, R : ",BUSZSI(S0, 1));
					EX_printPolyUni(BUSZSI(S0, 1), DAT(S0), 'x');
				}
				//
	//We keep hs into H
	preFFTRep *h1 = CreateUniPoly(dhs, hs);
	my_free(hs);	    	
        polyBothInHostDev *H = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
        H->status = -1;
        H->polyHost = h1;	
	sfixn i, j;			
	polyBothInHostDev *Szero;
	preFFTRep *szero;
	//ModularGenericSolve2ZeroResultant() is called if S0 is a zero polynomial
	if (zeroPolyp(S0))
	{		
			//
			if(PrintFlag == 1){
				printf("Calling ModularGenericSolve2ZeroResultant as resultant is zero\n");
			}
			//		
		//the 1st nonzero polynomial of subresultant chain is the j-th polynomial in it from bottom.
		// We will store the j-th polynomial in Szero.
		j = next_regular_subres_index(1, S, pPtr);
		Szero = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
		Szero->status = -1;
		//If no such j-th polynomial from bottom exists in subresultant chain, we will copy F1 in Szero.
		if(j >= EX_WidthInSCUBE(S)) 
		{
				//
				if(PrintFlag == 1){
					printf("ModularGenericSolve2ZeroResultant will be called with F1 as as lowest nonzero polynomial in subresultant chain\n");
				}
				//
			szero = EX_CopyOnePoly(F1->polyHost);
		}	
		else
			szero = EX_IthSubres(j, S, pPtr);
		Szero->polyHost = szero;
			//
			if(PrintFlag == 1){
			printf("the 1st polynomial(%dth from bottom of subresultant chain) has degree: %d %d\n", j, BUSZSI(Szero->polyHost, 1), BUSZSI(Szero->polyHost, 2));
			}			
			//
		ModularGenericSolve2ZeroResultant(method, F1, F2, Gnew, H, Szero, pPtr, rc);
 		EX_FreeOnePoly(S0);
		EX_SCUBE_Free(S);
		EX_FreeOnePoly(h1);
		EX_FreeOnePoly(gNew); 
		return;
	}
	sfixn *r, dr;
	r = SquarefreePart(&dr, BUSZSI(S0, 1), DAT(S0), pPtr);// computing square free part of resultant and store it in r
	sfixn *Q1, dQ1;
	// Q1 = r quo H
	if(method == CPU)
		Q1 = EX_UniQuo(&dQ1, dr, r, BUSZSI(H->polyHost, 1), DAT(H->polyHost), pPtr); 
	else if(method == GPU) 
		Q1 = divCUDAStraight(&dQ1, dr, r, BUSZSI(H->polyHost, 1), DAT(H->polyHost), pPtr);
	else
		Q1 = divCUDA(&dQ1, dr, r, BUSZSI(H->polyHost, 1), DAT(H->polyHost), pPtr);
					//
					if(PrintFlag == 1){
					printf("R = sqfree(R) : ");
					EX_printPolyUni(dr, r, 'x');
					printf("R = R quo h : ");
					EX_printPolyUni(dQ1, Q1, 'x');
					//printf("R = gcd(R,g) : ");
					//EX_printPolyUni(dgcd, gcd1, 'x');
					}
					//
	sfixn *gcd1, dgcd1;
	polyBothInHostDev *R;
	preFFTRep *r2;
	//r2 = gcd(Q1,gnew)
	if(!zeroPolyp(Gnew->polyHost))
	{
		if(method == CPU)
			gcd1 = EX_GCD_Uni_Wrapper(&dgcd1, Q1, dQ1, DAT(Gnew->polyHost), BUSZSI(Gnew->polyHost, 1) , pPtr);
		else if(method == GPU)
			gcd1 = gcdCUDAStraight(&dgcd1, Q1, dQ1, DAT(Gnew->polyHost), BUSZSI(Gnew->polyHost, 1) , pPtr);
		else
			gcd1 = gcdCUDA(&dgcd1, Q1, dQ1, DAT(Gnew->polyHost), BUSZSI(Gnew->polyHost, 1) , pPtr);
					//
					if(PrintFlag == 1){
					printf("g is not a zero polynomial. R = gcd(R,g) : ");
					EX_printPolyUni(dgcd, gcd1, 'x');
					}
					//
		r2 = CreateUniPoly(dgcd1, gcd1);
	       	my_free(gcd1);
	}
	else
	{
		r2 = CreateUniPoly(dQ1, Q1);
					//
					if(PrintFlag == 1){
					printf("g is zero. so no more modification to R is required.\n");
					}
					//
	}

	R = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
	R->status = -1;
        R->polyHost = r2 ;
	dQ1 = shrinkDegUni(BUSZSI(R->polyHost, 1), DAT(R->polyHost));
	BUSZSI(R->polyHost, 1) = dQ1;	
					//
					if(PrintFlag == 1){
					printf("Final R is : ");
					EX_Poly_Print(R->polyHost);
					}
					//			
	my_free(Q1);
	my_free(r);
	EX_FreeOnePoly(gNew);
	EX_FreeOnePoly(h1);
	if(dQ1 == 0)
	{
				//
				if(PrintFlag == 1){
					printf("Back from ModularGenericSolve2 as R is constant. \n");		
				}
				//
		EX_FreeOnePoly(S0); 
		EX_SCUBE_Free(S);		
		return;
	}
	i = 1;
	sfixn a1;
	sfixn w = EX_WidthInSCUBE(S); //w = width of the subresultant chain
	sfixn *g4, dg4, *g5, dg5;
	preFFTRep *T1, *T2, *T3, *T10, *copyDummy1, *copyDummy2;
			//
			if(PrintFlag == 1){
			printf("The width of subresultant chain: %d. Degree(F2, y) %d \n starting loop,\n",w, BUSZSI(F2->polyHost, 2));
			}
			//	
	while( !(constantPolyp(R->polyHost)) )
	{	
				//
				if(PrintFlag == 1){
				printf("i= %d\n",i);
				}
				//
		// the following redundant line of code is written to avoid compiler warning
		j = next_regular_subres_index(0, S, pPtr);
		T1 = F1->polyHost;
		while(i <= BUSZSI(F2->polyHost, 2)) // R T1
		{
			j = next_regular_subres_index(i, S, pPtr);
				//
				if(PrintFlag == 1){
				printf("j= %d\n",j);				
				}
				//
                        if(j < w) // if j is a polynomial in subresultant chain, where i <= j < deg(F1, y) and j is minimum
			{			
				T1 = EX_IthDthCoeff(j, j, S, pPtr);
				a1 = shrinkDegUni(BUSZSI(T1, 1), DAT(T1));	
				BUSZSI(T1, 1) = a1;
					//
					if(PrintFlag == 1){
						printf("j<w\n");						
					}				
					//
			}
			else // if j is a polynomial in subresultant chain, where i <= j = deg(F2, y) 
			{
				T1 = EX_getInitial(F2->polyHost);
				a1 = shrinkDegUni(BUSZSI(T1, 1), DAT(T1));
				BUSZSI(T1, 1) = a1;
						//
					if(PrintFlag == 1){
						printf("j>=w\n");						
					}				
					//
			}
					//
					if(PrintFlag == 1){
					printf("The leading coefficient from subresultant chain: ");
					EX_printPolyUni(a1, DAT(T1), 'a');
					}
					//
			//if(is_divided_by(BUSZSI(T1, 1), DAT(T1), BUSZSI(R->polyHost, 1), DAT(R->polyHost), pPtr) == 0) 
			// if R  not divides T1, we got the correct j, break the last inner loop
			if(isDividedWithCudaOption(method,BUSZSI(T1, 1), DAT(T1), BUSZSI(R->polyHost, 1), DAT(R->polyHost), pPtr ) == 0 )
			{
					//
					if(PrintFlag == 1){
					printf(" not divisible..breaking the inner loop..\n");
					}
					//
					break;
			}
			else // R does  divides T1 so check for the next poly in subresultant chain by increasing i
			{
					//
					if(PrintFlag == 1){
					printf(" divisible..\n");
					}
					//
				 ++i;													
			}
		}
				//
				if(PrintFlag == 1){
				printf("i = %d\n", i);
				}
				//
		//adding current R and F1 to the regular chain.
		if(i > BUSZSI(F2->polyHost, 2))
		{
				//
				if(PrintFlag == 1){
				printf("adding current R and F1 to the regular chain. Back from ModularGenericSolve2\n");
				}
				//
			copyDummy1 = EX_CopyOnePoly(R->polyHost);
			copyDummy2 = EX_CopyOnePoly(F1->polyHost);
	        	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(copyDummy1, copyDummy2));
			EX_FreeOnePoly(S0); 
			EX_SCUBE_Free(S);		
			return;			
		}
		// g4 = gcd(R, T1) or gcd(R,lc(jth subresultant polynomial)
 		if(method == CPU)
			g4 = EX_GCD_Uni_Wrapper(&dg4, DAT(T1), BUSZSI(T1, 1), DAT(R->polyHost), BUSZSI(R->polyHost, 1) , pPtr);
		else if(method == GPU)
			g4 = gcdCUDAStraight(&dg4, DAT(T1), BUSZSI(T1, 1), DAT(R->polyHost), BUSZSI(R->polyHost, 1) , pPtr);
		else
			g4 = gcdCUDA(&dg4, DAT(T1), BUSZSI(T1, 1), DAT(R->polyHost), BUSZSI(R->polyHost, 1) , pPtr);
			//
			if(PrintFlag == 1)
			{
				printf("G = gcd(R,lc(jth subresultant polynomial)) : ");
				EX_printPolyUni(dg4, g4, 'x');
			}
			// T2 contains the j-th poly of subresultant chain or F1 or F2 (based on j and w).
		if(j < w)	T2 = EX_IthSubres(j, S, pPtr);     
		else 		T2 = EX_CopyOnePoly(F2->polyHost);
			//
			if(PrintFlag == 1)
			{
				printf("%d th polynomial of subresultant is : ",j);
				EX_Poly_Print(T2);
			}
			//
		if(dg4 == 0)//g4 is constant. Adding (R, jth subresultant polynomial) in the regular chain.Back from ModularGenericSolve2
		{
				//
				if(PrintFlag == 1){
				printf("G is constant. Adding (R, jth subresultant polynomial) in the regular chain.Back from ModularGenericSolve2 \n");
				}
				//
			copyDummy1 = EX_CopyOnePoly(R->polyHost);
			copyDummy2 = EX_CopyOnePoly(T2);
	        	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(copyDummy1, copyDummy2));
			EX_FreeOnePoly(R->polyHost);
			EX_FreeOnePoly(S0); 
			EX_SCUBE_Free(S);		
			return;		
		}
		if(method == CPU)
			g5 = EX_UniQuo(&dg5, BUSZSI(R->polyHost, 1), DAT(R->polyHost), dg4, g4, pPtr);
		else if(method == GPU) 
			g5 = divCUDAStraight(&dg5, BUSZSI(R->polyHost, 1), DAT(R->polyHost), dg4, g4, pPtr);
		else
			g5 = divCUDA(&dg5, BUSZSI(R->polyHost, 1), DAT(R->polyHost), dg4, g4, pPtr);
		T10 = CreateUniPoly(dg4, g4); //T10 = g4
		my_free(g4);
		T3 = CreateUniPoly(dg5, g5); // T3 = R quo T10
		my_free(g5); 
			//
			if(PrintFlag == 1){
			printf("G is not constant. R quo G: ");
			EX_Poly_Print(T3);
			printf(" Adding (R quo G, jth subresultant polynomial) in the regular chain. \n");
			}
			//
		//Adding (R quo G, lc(jth subresultant polynomial)) in the regular chain.
		copyDummy1 = EX_CopyOnePoly(T3);
		copyDummy2 = EX_CopyOnePoly(T2);				
        	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(copyDummy1, copyDummy2));
		R->polyHost = T10;
		i = j+1;
		dg5 = shrinkDegUni(BUSZSI(R->polyHost, 1), DAT(R->polyHost));
		//updating R
		BUSZSI(R->polyHost, 1) = dg5;
				//
				if(PrintFlag == 1){
				printf("new R : ");
				EX_Poly_Print(R->polyHost);
				}
				//
		if(dg5 == 0)		
		{
				//
				if(PrintFlag == 1){
					printf("Back from ModularGenericSolve2 as R is constant. \n");		
				}
				//
			break;
		}
	}
				//
				if(PrintFlag == 1){
					if(i > BUSZSI(F2->polyHost, 2)) printf("Back from ModularGenericSolve2 as i > deg(F2,y) \n");		
				}
				//	
	EX_FreeOnePoly(S0); 
	EX_SCUBE_Free(S);	
}


/*
Input:
@method:= flag for computing option between CPU and GPU. see Types.h for details.
@F1:= F1 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@F2:= F2 is a bivariate polynomial of x and y, where x < y, where coefficients are in a field.
@g:= g is a univariate polynomial in x, where coefficients are in a field. It is a zero polynomial when ModularSolve2() is called from top level function (not recursively).
@pPtr:= it is a prime number that defines the coefficient field.
@rc:= The regular chain that will be returned by this function. Intially it is an empty list.

Output:
On output this function add triaangular sets  {(A_1, B_1),...(A_e, B_e)} in rc
such that,
V(F1,F2,g) = V(A_1,B_1),..V(A_e,B_e).

For Details, please see 
Marc Moreno Maza and Wei Pan, Solving Bivariate Polynomial Systems on a GPU,
J. Phys.: Conf. Ser. 341 012022
*/

void ModularSolve2(computingChoince  method, polyBothInHostDev *F1, polyBothInHostDev *F2, polyBothInHostDev *g, MONTP_OPT2_AS_GENE *pPtr, LinkedQueue *rc)
{
		//
		if(PrintFlag == 1){
			printf("F1: ");
			EX_Poly_Print(F1->polyHost);
			printf("F2: ");
			EX_Poly_Print(F2->polyHost);
			printf("g: ");
			EX_Poly_Print(g->polyHost);
		}
		//
    // leading 0 coefficient with respect to variable y are discarded for both  F1 and F2.
    sfixn realDeg1 = shrinkDeg(BUSZSI(F1->polyHost, 2), DAT(F1->polyHost), BUSZSI(F1->polyHost, 1)+1); 
    sfixn realDeg2 = shrinkDeg(BUSZSI(F2->polyHost, 2), DAT(F2->polyHost), BUSZSI(F2->polyHost, 1)+1);
    BUSZSI(F1->polyHost, 2) = realDeg1;
    BUSZSI(F2->polyHost, 2) = realDeg2;	
	//
	if(PrintFlag == 1){
		printf("\n deg(F1, big var): %d deg(F2, big var): %d\n", realDeg1, realDeg2);
		printf("deg(F1, small var): %d deg(F2, small var): %d\n", BUSZSI(F1->polyHost, 1), BUSZSI(F2->polyHost, 1));
	}
	//

   sfixn realDegin;
   //If F1 is univariate in x then leading 0 coefficients with respect to x are discarded for F1
   if(BUSZSI(F1->polyHost, 2) == 0)
   {
	realDegin = shrinkDegUni(BUSZSI(F1->polyHost, 1), DAT(F1->polyHost));
	BUSZSI(F1->polyHost, 1) = realDegin;
   }	
   //If F2 is univariate in x then leading 0 coefficients with respect to x are discarded for F2
   if(BUSZSI(F2->polyHost, 2) == 0)	 	
   {
	realDegin = shrinkDegUni(BUSZSI(F2->polyHost, 1), DAT(F2->polyHost));
	BUSZSI(F2->polyHost, 1) = realDegin;
   }	
    //If all F1, F2 and g are zeros then return {(0,0)}
    if( (BUSZSI(F1->polyHost, 1) == 0) && (BUSZSI(F2->polyHost, 1) == 0) && (BUSZSI(F1->polyHost, 2) == 0) && (BUSZSI(F2->polyHost, 2) == 0) && zeroPolyp(g->polyHost) )
    {		
		//
		if(PrintFlag == 1){
			printf("\n All are zeros\n");
		}
		//
	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(NULL, NULL));
        return;
    }
    // If any of F1, F2, and g are constants then return
    if( constantPolyp(F1->polyHost) || constantPolyp(F2->polyHost) || constantPolyp(g->polyHost) )
    {
		//
		if(PrintFlag == 1){
			printf("\n Any of F1, or F2 or g is constant\n");
		}
		//
        return;	
    }
	
    preFFTRep *p3, *h1, *h2, *hpoly, *t1, *t2, *copyDummy;	
    polyBothInHostDev *F3,  *h, *T1, *T2, *P3;
    sfixn *gcd, *gcd1, dgcd, dgcd1, d1, d2, d3;	
    	    //
	    //sfixn i;
	    // 	   
    //if both F1 and F2 are univariate in x and if gcd(F1,F2,g) is not constant then add {(gcd(F1,F2,g),NULL)}  to the triangular set
    if( (BUSZSI(F1->polyHost, 2) <= 0) && (BUSZSI(F2->polyHost, 2) <= 0) )
    {
			//
			if(PrintFlag == 1){
				printf("Both of the polynomials are univariate in x\n");
			}
			//
	d1 = shrinkDegUni(BUSZSI(F1->polyHost, 1), DAT(F1->polyHost));
	d2 = shrinkDegUni(BUSZSI(F2->polyHost, 1), DAT(F2->polyHost));
			//
		if(PrintFlag == 1){
			printf("shrinking degree is done: %d %d\n", d1, d2);
		}
		//
	//compute gcd(F1,F2)
	if(method == CPU)
		gcd  = EX_GCD_Uni_Wrapper(&dgcd,  DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);
	else if(method == GPU)
		gcd = gcdCUDAStraight(&dgcd,  DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);
	else
		gcd  = gcdCUDA(&dgcd,  DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);
		//
		if(PrintFlag == 1){
			printf("gcd of F1 and F2 is  done:\n");

		}
		//
	d3 = shrinkDegUni(BUSZSI(g->polyHost, 1), DAT(g->polyHost));
	//compute gcd(gcd(F1,F2),g)
	if(method == CPU)		
		gcd1 = EX_GCD_Uni_Wrapper(&dgcd1, DAT(g->polyHost),  d3, gcd, dgcd, pPtr);		
        else if(method == GPU)
		gcd1 = gcdCUDAStraight(&dgcd1, DAT(g->polyHost),  d3, gcd, dgcd, pPtr);			
	else
		gcd1 = gcdCUDA(&dgcd1, DAT(g->polyHost),  d3, gcd, dgcd, pPtr);			
		//
		if(PrintFlag == 1){
					EX_printPolyUni(BUSZSI(g->polyHost, 1), DAT(g->polyHost), 'x');
					EX_printPolyUni(dgcd1, gcd1, 'x');
		}
		//		
	//if gcd(F1,F2,g) is not constant
	if(dgcd1 > 0)
	{
		p3 = CreateUniPoly(dgcd1, gcd1);
        	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(p3, NULL));
	}

        my_free(gcd);
	my_free(gcd1);	
	return;	
    }
    //if both of them are univariate in y and if gcd(F1,F2) is not constant then add {(g, gcd(F1,F2))} to the triangular set
    if( (BUSZSI(F1->polyHost, 1) <= 0) && (BUSZSI(F2->polyHost, 1) <= 0) )
    {
		//
		if(PrintFlag == 1){
			printf("Both of the polynomials are univariate in y\n");
		}
		//
	d1 = shrinkDegUni(BUSZSI(F1->polyHost, 2), DAT(F1->polyHost));
	d2 = shrinkDegUni(BUSZSI(F2->polyHost, 2), DAT(F2->polyHost));
		//
		if(PrintFlag == 1){
			printf("shrinking degree is done\n");
			printf("\n %d %d \n",d1, d2);
		}
		//
		
	if(method == CPU)
		gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);
        else if(method == GPU)
		gcd = gcdCUDAStraight(&dgcd, DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);			
	else
		gcd = gcdCUDA(&dgcd, DAT(F1->polyHost), d1, DAT(F2->polyHost), d2, pPtr);
	if(dgcd > 0)
	{
		p3 = CreateUniPolyY(dgcd, gcd);
		copyDummy = EX_CopyOnePoly(g->polyHost);
        	EX_LinkedQueue_Enqeue(rc, EX_RegularChain2_Init(copyDummy, p3));
	}		
	my_free(gcd);		
	return;	
    }
   //if F1 is univariate in y  then call ModularSolve2() recursively with F1+F2 and F2
    if( BUSZSI(F1->polyHost, 2) <= 0  )
    {
		//
		if(PrintFlag == 1){		
			printf("F1 is univariate in x. so we are adding F1 and F2.\n ");
		}
		//
	p3 = add_poly_uni_bi(F1->polyHost, F2->polyHost, pPtr->P);
        F3 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
        F3->status = -1;
        F3->polyHost = p3;
        ModularSolve2(method, F2, F3 , g, pPtr, rc);
	EX_FreeOnePoly(p3);
	return;	
    }
    //if F2 is univariate in y then call ModularSolve2() recursively with F1 and F1+F2
    if( BUSZSI(F2->polyHost, 2) <= 0 )
    {
		//
		if(PrintFlag == 1){		
			printf("F2 is univariate in x. so we are adding F1 and F2.\n ");
		}
		//

	p3 = add_poly_uni_bi(F2->polyHost, F1->polyHost,  pPtr->P);
        F3 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
        F3->status = -1;
        F3->polyHost = p3;
        ModularSolve2(method, F1, F3 , g, pPtr, rc);
	EX_FreeOnePoly(p3);
	return;	
    }
	    //
	    if(PrintFlag == 1){		
		    printf("both of the polynomials are bivariate\n");
	    }
	    //
    //deg(F1,y) should be less than or equal to deg(F2,y) otherwise swap(F1,F2) and call ModularSolve2() recursively 
    if( BUSZSI(F1->polyHost, 2) > BUSZSI(F2->polyHost, 2))
    {
	    if(PrintFlag == 1){		
	    printf("F1 and F2 are interchanged\n");
	    }	
            ModularSolve2(method, F2, F1 , g, pPtr, rc);
	    return;
    }
    //Both F1 and F2 are bivariate and deg(F1,y) <= deg(F2,y) 	
    if( BUSZSI(F1->polyHost, 2) <= BUSZSI(F2->polyHost, 2))	
    {
	    h1 = EX_getInitial(F1->polyHost); //h1 = lc(F1)
	    h2 = EX_getInitial(F2->polyHost); //h2 = lc(F2)
	    d1 = shrinkDegUni(BUSZSI(h1, 1), DAT(h1)); //cancelling leading 0s of h1
	    BUSZSI(h1, 1) = d1;
	    d2 = shrinkDegUni(BUSZSI(h2, 1), DAT(h2));//cancelling leading 0s of h2
	    BUSZSI(h2, 1) = d2;
	    //computing gcd(h1,h2)
	    if(method == CPU)
	    	gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(h1), d1, DAT(h2), d2, pPtr);
	    else if(method == GPU)
		gcd = gcdCUDAStraight(&dgcd, DAT(h1), d1, DAT(h2), d2, pPtr);
	    else	
		gcd = gcdCUDA(&dgcd, DAT(h1), d1, DAT(h2), d2, pPtr);

	    dgcd = shrinkDegUni(dgcd, gcd);		   
	    hpoly = CreateUniPoly(dgcd, gcd);
	    h = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
            h->status = -1;
            h->polyHost = hpoly; // h= gcd(lc(F1),lc(F2))
		//
		if(PrintFlag == 1){
			printf("leading coeff of F1 (x means x_1): ");		
			EX_printPolyUni(BUSZSI(h1, 1), DAT(h1), 'x');
			printf("leading coeff of F2 (x means x_1): ");
			EX_printPolyUni(BUSZSI(h2, 1), DAT(h2), 'x');	
			printf(" the degree of lc(F1), lc(F2) and deg(gcd(lc(F1),lc(F2))): %d %d %d\n",d1,d2, dgcd);
			printf("gcd of leading coefficient (x means x_1): ");
			EX_printPolyUni(dgcd, gcd, 'x');
			//EX_printPolyUni(BUSZSI(hpoly, 1), DAT(hpoly), 'x');		
			printf("calling ModularGenericSolve2\n");
		}
		//
	    ModularGenericSolve2(method, F2, F1, g, h, pPtr, rc);
		//
		if(PrintFlag == 1){		
			printf("back from ModularGenericSolve2. degree of the gcd of leading coeff: %d\n", dgcd);
		}
		//	
	    // h1 h2 gcd hpoly
	    EX_FreeOnePoly(h1);
	    EX_FreeOnePoly(h2); 
            my_free(gcd);
	    // if gcd(lc(F1),lc(F2)) is not constant
	    if(dgcd > 0)
	    {
		t1 = EX_GetPolyTail(F1->polyHost);// t1 = tail(F1) or Poly F1 without lc
        	T1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
                T1->status = -1;
                T1->polyHost = t1; 

	        t2 = EX_GetPolyTail(F2->polyHost); // t2 = tail(F2) or Poly F2 without lc
	        T2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
                T2->status = -1;
                T2->polyHost = t2;

		d1 = shrinkDegUni(BUSZSI(g->polyHost, 1), DAT(g->polyHost));//cancelling leading 0s of g
	        d2 = shrinkDegUni(BUSZSI(h->polyHost, 1), DAT(h->polyHost));//cancelling leading 0s of h
			//
			if(PrintFlag == 1){
				printf("tail of F1: ");	
				EX_Poly_Print(T1->polyHost);
				printf("tail of F2: ");	
				EX_Poly_Print(T2->polyHost);
				printf("g: ");	
				EX_Poly_Print(g->polyHost);
				printf("h: ");	
				EX_Poly_Print(h->polyHost);	
			}
			//
		//computing gcd(g,h)
	        if(method == CPU)
	        	gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(g->polyHost), d1, DAT(h->polyHost), d2, pPtr);
		else if(method == GPU)
			gcd = gcdCUDAStraight(&dgcd, DAT(g->polyHost), d1, DAT(h->polyHost), d2, pPtr);
		else
			gcd = gcdCUDA(&dgcd, DAT(g->polyHost), d1, DAT(h->polyHost), d2, pPtr);
		dgcd = shrinkDegUni(dgcd, gcd);
	        p3 = CreateUniPoly(dgcd, gcd);
		P3 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
                P3->status = -1;
                P3->polyHost = p3;
			//
			if(PrintFlag == 1){	
				printf("gcd of g and h: ");
				EX_Poly_Print(P3->polyHost);
			}
			//
	        ModularSolve2(method, T1, T2, P3, pPtr, rc);
		EX_FreeOnePoly(t1);
         	EX_FreeOnePoly(t2);  
		EX_FreeOnePoly(p3); 
                my_free(gcd);
	    }		     		    
	    EX_FreeOnePoly(hpoly); 
	    return;
    }	
}


/**
 * With option to select scube type
 */

LinkedQueue *modular_solve2_select(sfixn method, preFFTRep *F1, preFFTRep *F2,
    preFFTRep *g, MONTP_OPT2_AS_GENE *pPtr)
{
	
    sfixn d1, d2, d3, k;
    SCUBE *scube;
    preFFTRep *R, *gcdpoly, *F3, *zpoly, *Sk; 
    sfixn *gcd, *gcd2, dgcd, dgcd2;
    LinkedQueue *resQueue;

    // the case that g is nonzero constant
    if (!zeroPolyp(g) && constantPolyp(g)) {
        if (DEBUG) printf("g is a nonzero constant\n");
        return EX_LinkedQueue_Init();
    }
    
    // the case that both F1 and F2 are univariate polynomials in x
    if (N(F1) == N(F2) && N(F2) == 1) { 
        if (DEBUG) printf("both F1 and F2 are univariate polynomials\n");
        d1 = shrinkDeg(BUSZSI(F1, 1), DAT(F1), CUMI(F1, 1));
        d2 = shrinkDeg(BUSZSI(F2, 1), DAT(F2), CUMI(F2, 1));

        resQueue = EX_LinkedQueue_Init();
        gcd = EX_GCD_Uni_Wrapper(&dgcd, DAT(F1), d1, DAT(F2), d2, pPtr);
        if (dgcd == 0) {
            my_free(gcd);
            return resQueue;  
        } 

        if (zeroPolyp(g)) {
            gcdpoly = CreateUniPoly(dgcd, gcd);
            EX_LinkedQueue_Enqeue(resQueue, 
                EX_RegularChain2_Init(gcdpoly, NULL));
            my_free(gcd);
            return resQueue;
        }
        
        d3 = shrinkDeg(BUSZSI(g, 1), DAT(g), CUMI(g, 1));
        gcd2 = EX_GCD_Uni_Wrapper(&dgcd2, gcd, dgcd, DAT(g), d3, pPtr);
        
        if (dgcd2 == 0) {
            my_free(gcd);
            my_free(gcd2);
            return resQueue;  
        }

        gcdpoly = CreateUniPoly(dgcd2, gcd2);
        EX_LinkedQueue_Enqeue(resQueue, EX_RegularChain2_Init(gcdpoly, NULL));

        my_free(gcd);
        my_free(gcd2);
        return resQueue;
    }

    if (N(F1) == 1) { 
        if (DEBUG) printf("F1 is a univariate polynomial\n");
        F3 = add_poly_uni_bi(F1, F2, pPtr->P);
        resQueue = modular_solve2_select(method, F3, F2, g, pPtr);
        EX_FreeOnePoly(F3);
        return resQueue;
    }

    if (N(F2) == 1) { 
        if (DEBUG) printf("F2 is a univariate polynomial\n");
        F3 = add_poly_uni_bi(F2, F1, pPtr->P);
        resQueue = modular_solve2_select(method, F1, F3, g, pPtr);
        EX_FreeOnePoly(F3);
        return resQueue;
    }
    
    // the case both F1 and F2 are bivariate
    if (DEBUG) printf("F1 and F2 are bivariate polynomials\n");
    d1 = shrinkDeg(BUSZSI(F1, 2), DAT(F1), CUMI(F1, 2));
    d2 = shrinkDeg(BUSZSI(F2, 2), DAT(F2), CUMI(F2, 2));
    if (d1 < d2) { return modular_solve2_select(method, F2, F1, g, pPtr); }
    
    // build the scube with various types of code
    // method == 0 ==> gpu_fft + subprodtree
    // method == 1 ==> cpu_fft + subprodtree
    computingChoince myChoice;
   if(method ==0 )    myChoice = GPU;
   if(method == 1 )   myChoice = CPU;
    scube = EX_SubResultantChainSelect(myChoice, F1, F2, 2, pPtr);

    R = EX_IthDthCoeff(0, 0, scube, pPtr);
    if (zeroPolyp(R)) {
        if (DEBUG) printf("The resultant is zero\n"); 
        // S_k is the last nonzero regular subresultant of F1 and F2;
        k = next_regular_subres_index(1, scube, pPtr);
        zpoly = CreateZeroPoly();

        if (k == EX_WidthInSCUBE(scube)) {
            resQueue = modular_solve2_select(method, F2, g, zpoly, pPtr);
        } else {
            Sk = EX_IthSubres(k, scube, pPtr);
            resQueue = modular_solve2_select(method, Sk, g, zpoly, pPtr);
        }
        EX_FreeOnePoly(zpoly);
        EX_SCUBE_Free(scube);
        return resQueue;
    }

    resQueue = modular_solve2_select_inner(method, F1, F2, g, scube, pPtr);
    EX_SCUBE_Free(scube);
    return resQueue;
}


// This function will test whether the prime, p is good for this bivariate system (f1, f2)
// 1 -- for good.
// 0 -- for bad.

sfixn testP(preFFTRep *f1, preFFTRep *f2, sfixn p)
{    
	sfixn n, m, k;
	n = BUSZSI(f1, 1)+1;
	m = BUSZSI(f1, 2)+1;
	//printf("degree of first poly: %d %d\n",n-1,m-1);
	//printf("degree of second poly: %d %d\n",BUSZSI(f2, 1),BUSZSI(f2, 2));

	if(BUSZSI(f1, 1) < BUSZSI(f2, 1)) 	n = BUSZSI(f2, 1)+1;
	if(BUSZSI(f1, 2) < BUSZSI(f2, 2)) 	m = BUSZSI(f2, 2)+1;
	//printf("max n m : %d %d\n",n,m);
	n = 2*n*m;
	//printf("n: %d \n",n);
	p = p-1;
	//printf("p: %d \n",p);
	k = 1;
	while( ((p>>1)<<1) == p)
	{
		k = 2*k;
		p= p >>1;
		//printf("p k: %d %d \n",p, k);
		
	
	}
	if(k > n) return 1;
	return 0; 
}


/*
* NEW BIVARIATE SOLVER
* IMPLEMENTED BY SARDAR HAQUE ON 13TH NOVEMBER 2012
* Wrapper function, exported 
 * 
 * @F1: poly in Zp[x, y]
 * @F2: poly in Zp[x, y]
 * @p: fourier prime number
 * 
 * @return, a list of regular_chain2 objects stored in a linked-queue.
 *
 * A possible output is 
 *
 * [(A_1, B_1),  (A_2, NULL), (A_3, B_3), (NULL, B_4)]
 *
 * where A_i are non-constant univariate polynomials and
 * B_i are bivariate polynomials with positive degree in y.
*/

LinkedQueue *EX_ModularSolve2(preFFTRep *f1, preFFTRep *f2, sfixn p, computingChoince method) {

  
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    preFFTRep *g = CreateZeroPoly();
    //LinkedQueue *resQueue = NULL; this is a bug

    LinkedQueue *resQueue = EX_LinkedQueue_Init();
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	//
	//if(PrintFlag == 1){		
	//	EX_Poly_Print(f1);
	//	EX_Poly_Print(f2);
	//}

    
    polyBothInHostDev *F1, *F2, *Zerop;
    F1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
    F1->status = -1;
    F1->polyHost = f1; 

    F2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
    F2->status = -1;
    F2->polyHost = f2; 

    Zerop = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
    Zerop->status = -1;
    Zerop->polyHost = g;

    ModularSolve2(method, F1, F2, Zerop, pPtr, resQueue);


    //resQueue = modular_solve2(F1, F2, g, pPtr);
    EX_Poly_Free(g);
    return resQueue;
}



/**
 * OLD BIVARIATE SOLVER 
 * THIS IS IMPLEMENTED AND TESTED BY XIL LI AND WEI PAN
 * THIS IS MADE COMMENT BY SARDAR HAQUE ON 13TH NOVEMBER 2012
 
 * Wrapper function, exported 
 * 
 * @F1, poly in Zp[x, y]
 * @F2, poly in Zp[x, y]
 * @p, fourier prime number
 * 
 * @return, a list of regular_chain2 objects stored in a linked-queue.
 *
 * A possible output is 
 *
 * [(A_1, B_1),  (A_2, NULL), (A_3, B_3), (NULL, B_4)]
 *
 * where A_i are non-constant univariate polynomials and
 * B_i are bivariate polynomials with positive degree in y.
 */
/*
LinkedQueue *EX_ModularSolve2(preFFTRep *F1, preFFTRep *F2, sfixn p) {
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    preFFTRep *g = CreateZeroPoly();
    LinkedQueue *resQueue = NULL;
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    //EX_Poly_Print(F1);
    //EX_Poly_Print(F2);
    resQueue = modular_solve2(F1, F2, g, pPtr);
    EX_Poly_Free(g);
    return resQueue;
}
*/
///////////////////////////END OF FILE/////////////////////////////////////////
