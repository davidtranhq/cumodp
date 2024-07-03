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

#include "FDIV.h"

extern int Interrupted;
//===================================================================
// Newton Inverse.
//===================================================================

// type: exported.
// Input: degG=n-1, degF<l where x^l is the modulus.
//        degF< l implies F is reduced r.w.t x^l.
//        r=ceiling(log_2(l)), n=2^r.
// Output: G with degG=n-1. G is the inverse of F modulo x^r.
//         =>  G is also the inverse of F modulo x^l.
// note: good for all Fourier Primes.
/**
 * modularInvPM: Computes the inverse of F modulo x^r 
 *               This assumes F has a non-zero trailing coefficient.
 *               This works for all Fourier Primes.
 * @degG: degree of G (output) 
 * @GPtr: coefficient vector of G (output) 
 * @degF: degree of F
 * @FPtr: coefficient vector of F
 * @r: 
 * @n: equals 2^r
 * @pPtr: prime number structure
 * 
 * Using the Middle product trick.
 * So the running time is expected to be in 2 M(n) + o(M(n)) machine 
 * word operations.
 * Return value: G, the inverse of F modulo x^r 
 **/   
sfixn *modularInvPM(sfixn degG, sfixn * GPtr, sfixn degF, sfixn * FPtr, 
        sfixn r, sfixn n, MONTP_OPT2_AS_GENE * pPtr)
{
    int i,j;
    sfixn nn, halfnn;
    sfixn * rootsPtr=(sfixn *)my_calloc(n, sizeof(sfixn)), 
          * tmpGPtr=(sfixn *)my_calloc(n, sizeof(sfixn)), 
          * tmpFPtr=(sfixn *)my_calloc(n, sizeof(sfixn));


    GPtr[0]=1;

    for(i=1;i<=r;i++){
        nn=1<<i;
        halfnn=nn>>1;
        EX_Mont_GetNthRoots_OPT2_AS_GENE(i, nn, rootsPtr, pPtr);
        EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpGPtr, nn-1, GPtr, pPtr);
        if(degF>=(nn-1)) 
            EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpFPtr, nn-1, FPtr, pPtr);
        else    
            EX_Mont_DFT_OPT2_AS_GENE( nn, i, rootsPtr, tmpFPtr, degF, FPtr, pPtr);
        EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
        EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
        for(j=0;j<halfnn; j++) tmpFPtr[j]=tmpFPtr[j+halfnn];
        for(j=halfnn; j<nn; j++) tmpFPtr[j]=0;
        EX_Mont_DFT_OPT2_AS_GENE_1 ( nn, i, rootsPtr, tmpFPtr, pPtr );
        EX_Mont_PairwiseMul_OPT2_AS_R(nn, tmpFPtr, tmpGPtr, pPtr);
        EX_Mont_INVDFT_OPT2_AS_GENE_R_1 (nn, i, rootsPtr, tmpFPtr, pPtr);
        for(j=halfnn; j<nn; j++) GPtr[j]=SubMod(GPtr[j],tmpFPtr[j-halfnn],pPtr->P);
    }
    my_free (rootsPtr);
    my_free (tmpGPtr);
    my_free (tmpFPtr);
    return GPtr;
}

//===================================================================
// Fast Division.
//===================================================================
// type: Only for benchmark use. 
// Input: A, B. Suppose deg(A)>=deg(B)
//        B is monic.
// Ouput: (Q , R), A=QB+R. 
// note: good for all Fourier Primes.
/**
 * fastDiv_bench: Computes the division of A by B, returning the
 *                quotient Q and the coefficient vector of the remainder R.
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degQ: (output) degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: the quotient and the remainder of A by B
 **/   
double fastDiv_bench(sfixn * RPtr, sfixn degQ, sfixn * QPtr, sfixn degA, 
        sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    register sfixn j;
    sfixn * FPtr, * GPtr;
    sfixn degF, power1,power2,power3, n1,n2,n3,l1,l2,l3, tmp, sz1, sz2, sz3;
    sfixn dg, da;

    degF=degA-degB;
    if(degF<0) { //printf("Attension: degA<degB!\n");
        for(j=0;j<=degA;j++) RPtr[j]=APtr[j];
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return 0;
    }

    n1=1; power1=0; l1=degF+1; tmp=l1;
    while(tmp){tmp>>=1; n1<<=1; power1++;}

    n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
    while(tmp){tmp>>=1; n2<<=1; power2++;}

    n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
    while(tmp){tmp>>=1; n3<<=1; power3++;}

    dg=da=degF;
    sz1=sz2=l2;
    if(sz1<n1) sz1=n1;
    if(sz2<(degA+1)) sz2=degA+1;

    //degG=n1-1;
    sz3=degF+1;
    if(sz3<degB+1) sz3=degB+1;
    if(sz3<sz2) sz3=sz2;
    if(sz3<n2) sz3=n2;
    if(sz2<n3) sz3=n3;
    sz1=sz3;

    //  1. F=rev_m(B);
    FPtr=(sfixn * )my_calloc(sz3 ,sizeof(sfixn));
    FPtr=reverseUni(degB, FPtr, BPtr);

    //  2. get G such that GF=1 mod x^n;   my_free(rootsPtr);
    GPtr=(sfixn * )my_calloc(sz1, sizeof(sfixn));

    //===============================================================>
    GPtr=modularInvPM(n1-1, GPtr, degF, FPtr, power1, n1, pPtr);
    //===============================================================>

    //  3. rev_n(A)*G mod n; 
    FPtr=reverseUni(degA, FPtr, APtr);
    EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
    QPtr=reverseUni(degQ, QPtr, FPtr);
    cleanVec(n3-1, FPtr);
    EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
    for(j=0; j<l3; j++) RPtr[j]=SubMod(APtr[j], FPtr[j],pPtr->P);
    my_free(FPtr);
    my_free(GPtr); 

    return 0; 
} 

// type: exported 
// Input: A, B. Suppose deg(A)>=deg(B)
//        B is monic.
// Ouput: (Q , R), A=QB+R. 
// note: good for all Fourier Primes.
/**
 * fastDiv: Computes the division of A by B, returning the
 *                quotient Q and the coefficient vector of the remainder R.
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degQ: (output) degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: the quotient and the remainder of A by B
 **/   
void fastDiv(sfixn * RPtr, sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    register sfixn j;
    sfixn * FPtr, * GPtr;
    sfixn degF, power1,power2,power3, n1,n2,n3,l1,l2,l3, tmp, sz1, sz2, sz3;
    sfixn dg, da;

    degF=degA-degB;
    if(degF<0) { 
        //printf("Attension: degA<degB!\n");
        for(j=0;j<=degA;j++) RPtr[j]=APtr[j];
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }

    n1=1; power1=0; l1=degF+1; tmp=l1;
    while(tmp){tmp>>=1; n1<<=1; power1++;}

    n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
    while(tmp){tmp>>=1; n2<<=1; power2++;}

    n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
    while(tmp){tmp>>=1; n3<<=1; power3++;}

    dg=da=degF;
    sz1=sz2=l2;
    if(sz1<n1) sz1=n1;
    if(sz2<(degA+1)) sz2=degA+1;

    //degG=n1-1;
    sz3=degF+1;
    if(sz3<degB+1) sz3=degB+1;
    if(sz3<sz2) sz3=sz2;
    if(sz3<n2) sz3=n2;
    if(sz2<n3) sz3=n3;
    sz1=sz3;

    //  1. F=rev_m(B);
    FPtr=(sfixn * )my_calloc(sz3 ,sizeof(sfixn));
    FPtr=reverseUni(degB, FPtr, BPtr);

    //  2. get G such that GF=1 mod x^n;   my_free(rootsPtr);
    GPtr=(sfixn * )my_calloc(sz1, sizeof(sfixn));

    //===============================================================>
    GPtr=modularInvPM(n1-1, GPtr, degF, FPtr, power1, n1, pPtr);
    //===============================================================>

    //  3. rev_n(A)*G mod n; 
    FPtr=reverseUni(degA, FPtr, APtr);
    EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
    QPtr=reverseUni(degQ, QPtr, FPtr);
    cleanVec(n3-1, FPtr);
    EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
    for(j=0; j<l3; j++) RPtr[j]=SubMod(APtr[j], FPtr[j],pPtr->P);
    my_free(FPtr);
    my_free(GPtr); 
} 

/**
 * fastQuo: Computes the quotient  in the Euclidean
 *          (fast) division of A by B
 * @degQ: (output) degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: the quotient  in the Euclidean division of A by B.
 **/   
void fastQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    register sfixn j;
    sfixn *FPtr, *GPtr, *RPtr;
    sfixn degF, power1,power2,power3, n1,n2,n3,l1,l2,l3, tmp, sz1, sz2, sz3;
    sfixn dg, da;

    degF=degA-degB;

    if(degF<0) { 
        //printf("Attension: degA<degB!\n");
        //for(j=0;j<=degA;j++) RPtr[j]=APtr[j];
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }
    RPtr=(sfixn *)my_calloc(degB, sizeof(sfixn));
    n1=1; power1=0; l1=degF+1; tmp=l1;
    while(tmp){tmp>>=1; n1<<=1; power1++;}
    n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
    while(tmp){tmp>>=1; n2<<=1; power2++;}
    n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
    while(tmp){tmp>>=1; n3<<=1; power3++;}

    // sz1, sz2 are the size of vector for GPtr, ArevPtr. 
    // They are used twice in this function!
    dg=da=degF;
    sz1=sz2=l2;
    if(sz1<n1) sz1=n1;
    if(sz2<(degA+1)) sz2=degA+1;

    //degG=n1-1;
    sz3=degF+1;
    if(sz3<degB+1) sz3=degB+1;
    if(sz3<sz2) sz3=sz2;
    if(sz3<n2) sz3=n2;
    if(sz2<n3) sz3=n3;
    sz1=sz3;

    //  1. F=rev_m(B);
    FPtr=(sfixn * )my_calloc(sz3 ,sizeof(sfixn));
    FPtr=reverseUni(degB, FPtr, BPtr);

    //  2. get G such that GF=1 mod x^n;   my_free(rootsPtr);
    GPtr=(sfixn * )my_calloc(sz1, sizeof(sfixn));

    //===============================================================>
    GPtr=modularInvPM(n1-1, GPtr, degF, FPtr, power1, n1, pPtr);
    //===============================================================>

    //  3. rev_n(A)*G mod n; 
    FPtr=reverseUni(degA, FPtr, APtr);
    EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
    QPtr=reverseUni(degQ, QPtr, FPtr);

    my_free(FPtr);
    my_free(GPtr); 
    my_free(RPtr);
} 

/* The following code has bug for the following input:

p1:= 0*a^0+384276638*a^1+0*a^2+0*a^3+0*a^4+0*a^5+0*a^6+0*a^7+0*a^8+0*a^9+0*a^10+0*a^11+0*a^12+150                                                                                               
888267*a^13+0*a^14+0*a^15+0*a^16+0*a^17+0*a^18+0*a^19+0*a^20+0*a^21+0*a^22+0*a^23+0*a^24+2882222*                                                                                               
 a^25+0*a^26+0*a^27+0*a^28+0*a^29+0*a^30+0*a^31+0*a^32+0*a^33+0*a^34+0*a^35+0*a^36+378505063*a^37+                                                                                               
0*a^38+0*a^39+0*a^40+0*a^41+0*a^42+0*a^43+0*a^44+0*a^45+0*a^46+0*a^47+0*a^48+204579778*a^49+0*a^50
+0*a^51+0*a^52+0*a^53+0*a^54+0*a^55+0*a^56+0*a^57+0*a^58+0*a^59+0*a^60+314897937*a^61;

p2:= 0+0*a^1+0*a^2+0*a^3+0*a^4+0*a^5+0*a^6+0*a^7+0*a^8+0*a^9+0*a^10+0*a^11+0*a^12+0*a^13+0*a^14+0                                                                                                
 *a^15+0*a^16+0*a^17+0*a^18+0*a^19+0*a^20+0*a^21+0*a^22+0*a^23+0*a^24+0*a^25+0*a^26+0*a^27+0*a^28+\                                                                                                
 0*a^29+0*a^30+0*a^31+0*a^32+0*a^33+0*a^34+0*a^35+0*a^36+0*a^37+0*a^38+0*a^39+0*a^40+0*a^41+0*a^42\                                                                                                
 +0*a^43+0*a^44+0*a^45+0*a^46+0*a^47+0*a^48+0*a^49+0*a^50+0*a^51+0*a^52+0*a^53+0*a^54+0*a^55+0*a^5\                                                                                                
 6+0*a^57+0*a^58+0*a^59+0*a^60+0*a^61+0*a^62+0*a^63+0*a^64+0*a^65+0*a^66+0*a^67+0*a^68+0*a^69+0*a^\                                                                                                
 70+0*a^71+0*a^72+0*a^73+0*a^74+0*a^75+0*a^76+0*a^77+0*a^78+0*a^79+0*a^80+0*a^81+0*a^82+0*a^83+0*a\                                                                                                
 ^84+0*a^85+0*a^86+0*a^87+0*a^88+0*a^89+0*a^90+0*a^91+0*a^92+0*a^93+0*a^94+0*a^95+0*a^96+0*a^97+0*\                                                                                                
 a^98+0*a^99+0*a^100+0*a^101+0*a^102+0*a^103+0*a^104+0*a^105+0*a^106+0*a^107+0*a^108+0*a^109+0*a^1\                                                                                                
 10+0*a^111+0*a^112+0*a^113+0*a^114+0*a^115+0*a^116+0*a^117+0*a^118+0*a^119+0*a^120+0*a^121+0*a^12\                                                                                                
 2+0*a^123+0*a^124+0*a^125+0*a^126+0*a^127+0*a^128+0*a^129+0*a^130+0*a^131+0*a^132+0*a^133+0*a^134\                                                                                                
 +0*a^135+0*a^136+0*a^137+0*a^138+0*a^139+0*a^140+0*a^141+0*a^142+0*a^143+0*a^144+0*a^145+0*a^146+\                                                                                                
 0*a^147+0*a^148+0*a^149+0*a^150+0*a^151+0*a^152+0*a^153+0*a^154+0*a^155+0*a^156+0*a^157+0*a^158+0\                                                                                                
 *a^159+0*a^160+0*a^161+0*a^162+0*a^163+0*a^164+0*a^165+0*a^166+0*a^167+0*a^168+0*a^169+0*a^170+0*\                                                                                                
 a^171+0*a^172+0*a^173+0*a^174+0*a^175+0*a^176+0*a^177+0*a^178+0*a^179+0*a^180+0*a^181+0*a^182+0*a\                                                                                                
 ^183+0*a^184+0*a^185+0*a^186+0*a^187+0*a^188+0*a^189+0*a^190+0*a^191+0*a^192+0*a^193+0*a^194+0*a^\                                                                                                
 195+0*a^196+0*a^197+0*a^198+0*a^199+0*a^200+0*a^201+0*a^202+0*a^203+0*a^204+0*a^205+0*a^206+0*a^2\                                                                                                
 07+0*a^208+0*a^209+0*a^210+0*a^211+0*a^212+0*a^213+0*a^214+0*a^215+0*a^216+0*a^217+0*a^218+0*a^21\                                                                                                
 9+0*a^220+0*a^221+0*a^222+0*a^223+0*a^224+0*a^225+0*a^226+0*a^227+0*a^228+0*a^229+0*a^230+0*a^231\                                                                                                
 +0*a^232+0*a^233+0*a^234+0*a^235+0*a^236+0*a^237+0*a^238+0*a^239+0*a^240+0*a^241+0*a^242+0*a^243+\                                                                                                
 0*a^244+0*a^245+0*a^246+0*a^247+0*a^248+0*a^249+0*a^250+0*a^251+0*a^252+0*a^253+0*a^254+0*a^255+0\                                                                                                
 *a^256+0*a^257+0*a^258+0*a^259+0*a^260+0*a^261+0*a^262+0*a^263+0*a^264+0*a^265+0*a^266+0*a^267+0*\                                                                                                
 a^268+0*a^269+0*a^270+0*a^271+0*a^272+0*a^273+0*a^274+0*a^275+0*a^276+0*a^277+0*a^278+0*a^279+0*a\                                                                                                
 ^280+0*a^281+0*a^282+0*a^283+0*a^284+0*a^285+0*a^286+0*a^287+0*a^288+0*a^289+0*a^290+0*a^291+0*a^\                                                                                                
 292+0*a^293+0*a^294+0*a^295+0*a^296+0*a^297+0*a^298+0*a^299+0*a^300+0*a^301+0*a^302+0*a^303+0*a^3\                                                                                                
 04+0*a^305+0*a^306+0*a^307+0*a^308+0*a^309+0*a^310+0*a^311+0*a^312+0*a^313+0*a^314+0*a^315+0*a^31\                                                                                                
 6+0*a^317+0*a^318+0*a^319+0*a^320+0*a^321+0*a^322+0*a^323+0*a^324+0*a^325+0*a^326+0*a^327+0*a^328\                                                                                                
 +0*a^329+0*a^330+0*a^331+0*a^332+0*a^333+0*a^334+0*a^335+0*a^336+0*a^337+0*a^338+0*a^339+0*a^340+\                                                                                                
 0*a^341+0*a^342+0*a^343+0*a^344+0*a^345+0*a^346+0*a^347+0*a^348+0*a^349+0*a^350+0*a^351+0*a^352+0\                                                                                                
 *a^353+0*a^354+0*a^355+0*a^356+0*a^357+0*a^358+0*a^359+0*a^360+0*a^361+0*a^362+0*a^363+0*a^364+0*\                                                                                                
 a^365+0*a^366+0*a^367+0*a^368+0*a^369+0*a^370+0*a^371+0*a^372+0*a^373+0*a^374+0*a^375+0*a^376+0*a\                                                                                                
 ^377+0*a^378+0*a^379+0*a^380+0*a^381+0*a^382+0*a^383+0*a^384+0*a^385+0*a^386+0*a^387+0*a^388+0*a^\                                                                                                
 389+0*a^390+0*a^391+0*a^392+0*a^393+0*a^394+0*a^395+0*a^396+0*a^397+0*a^398+0*a^399+0*a^400+0*a^4\                                                                                                
 01+0*a^402+0*a^403+0*a^404+0*a^405+0*a^406+0*a^407+0*a^408+0*a^409+0*a^410+0*a^411+0*a^412+0*a^41\                                                                                                
 3+0*a^414+0*a^415+0*a^416+0*a^417+0*a^418+0*a^419+0*a^420+0*a^421+0*a^422+0*a^423+0*a^424+0*a^425\                                                                                                
 +0*a^426+0*a^427+0*a^428+0*a^429+0*a^430+0*a^431+0*a^432+0*a^433+0*a^434+0*a^435+0*a^436+0*a^437+\                                                                                                
 0*a^438+0*a^439+0*a^440+0*a^441+0*a^442+0*a^443+0*a^444+0*a^445+0*a^446+0*a^447+0*a^448+0*a^449+0\                                                                                                
 *a^450+0*a^451+0*a^452+0*a^453+0*a^454+0*a^455+0*a^456+0*a^457+0*a^458+0*a^459+0*a^460+0*a^461+11\                                                                                                
 5261378*a^462+0*a^463+0*a^464+0*a^465+0*a^466+0*a^467+0*a^468+0*a^469+0*a^470+0*a^471+0*a^472+0*a\                                                                                                
 ^473+91947322*a^474+0*a^475+0*a^476+0*a^477+0*a^478+0*a^479+0*a^480+0*a^481+0*a^482+0*a^483+0*a^4\                                                                                                
 84+0*a^485+89487732*a^486+0*a^487+0*a^488+0*a^489+0*a^490+0*a^491+0*a^492+0*a^493+0*a^494+0*a^495\                                                                                                
 +0*a^496+0*a^497+348167327*a^498+0*a^499+0*a^500+0*a^501+0*a^502+0*a^503+0*a^504+0*a^505+0*a^506+\                                                                                                
 0*a^507+0*a^508+0*a^509+446294540*a^510+0*a^511+0*a^512+0*a^513+0*a^514+0*a^515+0*a^516+0*a^517+0\                                                                                                
 *a^518+0*a^519+0*a^520+0*a^521+362066393*a^522+0*a^523+0*a^524+0*a^525+0*a^526+0*a^527+0*a^528+0*\                                                                                                
 a^529+0*a^530+0*a^531+0*a^532+0*a^533+449244192*a^534+0*a^535+0*a^536+0*a^537+0*a^538+0*a^539+0*a\                                                                                                
 ^540+0*a^541+0*a^542+0*a^543+0*a^544+0*a^545+384641460*a^546+0*a^547+0*a^548+0*a^549+0*a^550+0*a^\                                                                                                
 551+0*a^552+0*a^553+0*a^554+0*a^555+0*a^556+0*a^557+465066292*a^558+0*a^559+0*a^560+0*a^561+0*a^5\                                                                                                
 62+0*a^563+0*a^564+0*a^565+0*a^566+0*a^567+0*a^568+0*a^569+256009023*a^570+0*a^571+0*a^572+0*a^573+0*a^574+0*a^575+0*a^576+0*a^577+0*a^578+0*a^579+0*a^580+0*a^581+391664336*a^582;

bug found by Sardar Haque on April 26th 2013.

 */
// type: exported 
// Input: A, B. Suppose deg(A)>=deg(B)
//        B is monic.
// Ouput: (Q , R), A=QB+R. 
// note: good for all Fourier Primes.
/**
 * fastRem: Computes the remainder of A by B.
 * @degRAddr: (output) degree of the remainder
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value:  remainder of A by B
 **/   
void fastRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn j;
    sfixn * FPtr, * GPtr, *QPtr;
    sfixn degQ, degF,  power1,power2,power3, n1,n2,n3,l1,l2,l3, tmp, sz1, sz2, sz3;
    sfixn dg, da;

    degF=degA-degB;
    degQ=degF;


    if(degF<0) { 
        //printf("Attension: degA<degB!\n");
        for(j=0;j<=degA;j++) RPtr[j]=APtr[j];
        //for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }

    QPtr=(sfixn *)my_calloc(degQ+1, sizeof(sfixn));

    n1=1; power1=0; l1=degF+1; tmp=l1;
    while(tmp){tmp>>=1; n1<<=1; power1++;}

    n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
    while(tmp){tmp>>=1; n2<<=1; power2++;}

    n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
    while(tmp){tmp>>=1; n3<<=1; power3++;}

    // sz1, sz2 are the size of vector for GPtr, ArevPtr. 
    // They are used twice in this function!
    dg=da=degF;
    sz1=sz2=l2;
    if(sz1<n1) sz1=n1;
    if(sz2<(degA+1)) sz2=degA+1;

    //degG=n1-1;
    sz3=degF+1;
    if(sz3<degB+1) sz3=degB+1;
    if(sz3<sz2) sz3=sz2;
    if(sz3<n2) sz3=n2;
    if(sz2<n3) sz3=n3;
    sz1=sz3;
    //  1. F=rev_m(B);
    FPtr=(sfixn * )my_calloc(sz3 ,sizeof(sfixn));
    FPtr=reverseUni(degB, FPtr, BPtr);
    //  2. get G such that GF=1 mod x^n;   my_free(rootsPtr);
    GPtr=(sfixn * )my_calloc(sz1, sizeof(sfixn));
    //===============================================================>
    GPtr=modularInvPM(n1-1, GPtr, degF, FPtr, power1, n1, pPtr);
    //===============================================================>
    //  3. rev_n(A)*G mod n; 
    FPtr=reverseUni(degA, FPtr, APtr);
    EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
    QPtr=reverseUni(degQ, QPtr, FPtr);
    cleanVec(n3-1, FPtr);
    EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
    for(j=0; j<=(*degRAddr); j++) RPtr[j]=SubMod(APtr[j], FPtr[j],pPtr->P);
    while((! RPtr[*degRAddr]) && ((*degRAddr)>0)) (*degRAddr)--;
    my_free(QPtr);
    my_free(FPtr);
    my_free(GPtr); 
} 


/**
 * UniRem: Computes the remainder of A by B.
 * @degRAddr: (output) degree of the remainder
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 *               Using classical division for low quotient degree
 *               and fast division otherwise
 *               This works for all Fourier Primes.
 * 
 * Return value: remainder of A by B
 **/   
void UniRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    if (degB==0){ *degRAddr=0; return;}
    /* 	the following lines are made comment by Sardar Haque on April 26th 2013
        as we found bug in fastRem() code
    if(degA-degB>60)
        fastRem(degRAddr, RPtr, degA, APtr, degB, BPtr, pPtr);
    else
    */	plainRem(degRAddr, RPtr, degA, APtr, degB, BPtr, pPtr);
}

/**
 * UniQuo: Computes the quotient  in the Euclidean
 *          (fast) division of A by B
 * @degQ: (output) degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: the quotient  in the Euclidean division of A by B.
 **/   
void UniQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    int i;
    sfixn inv, *UB, p = pPtr->P;

    if (degB==0){
        assert(BPtr[0]!=0);
        inv=inverseMod(BPtr[0], pPtr->P);
        for(i=0; i<=degA; i++) QPtr[i]=MulMod(inv, APtr[i], pPtr->P);
        return;
    }

    if (degA - degB >= 128) {
        assert(degQ == degA - degB);
        if (BPtr[degB] != 1) {
            inv = inverseMod(BPtr[degB], p);
            UB = (sfixn *)my_malloc(sizeof(sfixn)*(degB + 1));
            for (i = 0; i <= degB; ++i) { UB[i] = MulMod(BPtr[i], inv, p); }
            fastQuo(degQ, QPtr, degA, APtr, degB, UB, pPtr);
            for (i = 0; i <= degQ; ++i) { QPtr[i] = MulMod(QPtr[i], inv, p); }
            my_free(UB);
        } else {
            fastQuo(degQ, QPtr, degA, APtr, degB, BPtr, pPtr);
        }
    } else{
        plainQuo(degQ, QPtr, degA, APtr, degB, BPtr, pPtr);
    }
}

/**
 * UniPseuQuo: Computes the pseudo-quotient  in the 
 *          pseudo-division of A by B
 * @degQ: (output) degree of the pseudo-quotient
 * @QPtr: (output) Coefficient vector  of the pseudo-quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: pseudo-quotient  in the pseudo-division of A by B.
 **/   
void UniPseuQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    int i;
    sfixn inv;

    if (degB==0){
        assert(BPtr[0]!=0);
        inv=inverseMod(BPtr[0], pPtr->P);
        for(i=0; i<=degA; i++) QPtr[i]=MulMod(inv, APtr[i], pPtr->P);
        return;
    }
    PlainPseudoQuotient(&degQ, QPtr, degA, APtr,degB, BPtr, pPtr);
}

/**
 * EX_UniRem: Computes the remainder of A by B (univariate monic division)
 * @degRAddr: (output) degree of the remainder
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 *               Using classical division for low quotient degree
 *               and fast division otherwise
 *               This works for all Fourier Primes.
 * 
 * Return value: the coefficient vector of the remainder of A by B
 **/ 
sfixn *EX_UniRem(sfixn *degRAddr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    int i;
    sfixn * tmpRPtr, *RPtr;
    *degRAddr=degB-1;
    tmpRPtr=(sfixn *)my_calloc((*degRAddr)+1, sizeof(sfixn));
    UniRem(degRAddr, tmpRPtr, degA, APtr, degB, BPtr, pPtr);
    RPtr=(sfixn *)my_calloc((*degRAddr)+1, sizeof(sfixn));
    for(i=0; i<=(*degRAddr); i++) RPtr[i]=tmpRPtr[i];
    my_free(tmpRPtr);
    return RPtr;
}

/**
 * EX_UniQuo: Computes the quotient  in the Euclidean
 *          (fast) univariate division of A by B
 * @degQAddr: (output) degree of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: the coefficient vector of the quotient  
 *                in the Euclidean division of A by B.
 **/   
sfixn * EX_UniQuo(sfixn *degQAddr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    int i;
    sfixn *tmpQPtr, *QPtr;

	
    //*degQAddr=degA;
    *degQAddr=degA-degB; // fixed by WP Mon Jan 10 23:21:28 EST 2011
    tmpQPtr=(sfixn *)my_calloc((*degQAddr)+1, sizeof(sfixn));
    UniQuo(*degQAddr, tmpQPtr, degA, APtr, degB, BPtr, pPtr);
    *degQAddr=shrinkDegUni(*degQAddr, tmpQPtr);
    QPtr=(sfixn *)my_calloc((*degQAddr)+1, sizeof(sfixn));
    for(i=0; i<=(*degQAddr); i++) QPtr[i]=tmpQPtr[i];
    my_free(tmpQPtr);
    return QPtr;
}

// type: exported, in-place.
// Input: A, B. Suppose deg(A)>=deg(B)
//        B is monic.
// Ouput: (Q , R), A=QB+R. 
// note: good for all Fourier Primes.
//  A is over-written by the remainder
// fast-division
void fastDiv_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, sfixn degB, 
        sfixn * BPtr, sfixn * BRevInvPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn j;
    sfixn * FPtr, * GPtr;
    sfixn degF,  power1,power2,power3, n1,n2,n3,l1,l2,l3, tmp, sz1, sz2, sz3;
    sfixn dg, da;

    degF=degA-degB;

    if(degF<0) {
        //printf("Attension: degA<degB!\n");
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }

    n1=1; power1=0; l1=degF+1; tmp=l1;
    while(tmp){tmp>>=1; n1<<=1; power1++;}

    n2=1; power2=0; l2=(l1<<1)-1;tmp=l2;
    while(tmp){tmp>>=1; n2<<=1; power2++;}

    n3=1; power3=0;  l3=degB+degQ+1, tmp=l3;
    while(tmp){tmp>>=1; n3<<=1; power3++;}

    // sz1, sz2 are the size of vector for GPtr, ArevPtr. 
    // They are used twice in this function!
    dg=da=degF;
    sz1=sz2=l2;
    if(sz1<n1) sz1=n1;
    if(sz2<(degA+1)) sz2=degA+1;

    //degG=n1-1;
    sz3=degF+1;
    if(sz3<degB+1) sz3=degB+1;
    if(sz3<sz2) sz3=sz2;
    if(sz3<n2) sz3=n2;
    if(sz2<n3) sz3=n3;
    sz1=sz3;

    //  1. F=rev_m(B);
    FPtr=(sfixn * )my_calloc(sz3 ,sizeof(sfixn));
    FPtr=reverseUni(degB, FPtr, BPtr);

    //  2. get G such that GF=1 mod x^n;   my_free(rootsPtr);
    GPtr=(sfixn * )my_calloc(sz1, sizeof(sfixn));
    for(j=0;j<=dg; j++) GPtr[j]=BRevInvPtr[j];

    //  3. rev_n(A)*G mod n; 
    FPtr=reverseUni(degA, FPtr, APtr);
    EX_Mont_FFTMul_OPT2_AS_GENE_1(n2, power2, da+dg, da, FPtr, dg, GPtr, pPtr);
    QPtr=reverseUni(degQ, QPtr, FPtr);
    cleanVec(n3-1, FPtr);
    EX_Mont_FFTMul_OPT2_AS_GENE(n3, power3, degQ+degB, FPtr, degQ, QPtr, degB, BPtr, pPtr);
    for(j=0; j<l3; j++) APtr[j]=SubMod(APtr[j], FPtr[j],pPtr->P);
    for(j=l3; j<=degA; j++) APtr[j]=0;

    my_free(FPtr);
    my_free(GPtr); 
} 

//===================================================================
// fmecg(res,e,r,p2) finds X :       res - r * X**e * p2
//===================================================================

// deg(res)>=deg(p2)+e
// note: good for all Fourier Primes.  
void fmedg_1(sfixn degRes, sfixn * resPtr, sfixn e, sfixn r, 
        sfixn degp2, sfixn * p2Ptr, MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn i, p=pPtr->P, R=(1L<<pPtr->Rpow)%p;
    R=(MulMod(r, R, p))<<pPtr->Base_Rpow;
    for(i=0; i<=degp2; i++) 
        resPtr[i+e]=SubMod(resPtr[i+e], MontMulMod_OPT2_AS_GENE(p2Ptr[i], R, pPtr), p);
}

//===================================================================
// Plain Division.
//===================================================================
// note: good for all Fourier Primes.
/**
 * plainDivMonic_1: Computes the monic division of A by B, producing the
 *                quotient Q and the coefficient vector of the remainder R.
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degQ: (output) degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 *               Use the classical / plain algorithm.
 *               In place: coeff vector of A is over-written 
 *               by that of the  remainder.
 *               The degree of the remainder is not returned.
 * 
 * Return value: the quotient and the remainder of A by (monic) B
 **/ 
void plainDivMonic_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn i,j;
    sfixn tmp;

    if((degA-degB)<0) {
        //printf("Attension: degA<degB!\n");
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }

    for(i=degA-degB; i>=0; i--){
        if (APtr[degB+i] != 0){
            QPtr[i]=APtr[degB+i];
            tmp=MontMulMod_OPT2_AS_GENE(QPtr[i], pPtr->R2BRsft, pPtr)<<(pPtr->Base_Rpow);
            for(j=0; j<=degB; j++) 
                APtr[i+j]=SubMod(APtr[i+j], MontMulMod_OPT2_AS_GENE(BPtr[j],tmp,pPtr), pPtr->P);
        } else{QPtr[i]=0;}
    }   
}

/**
 * plainRemMonic_1: Computes the remainder in the monic division of A by B.
 8                  A is overwritten.
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 *               Use the classical / plain algorithm.
 *               In place: coeff vector of A is over-written 
 *               by that of the  remainder.
 *               The degree of the remainder is not returned.
 * 
 * Return value: void
 **/ 
void plainRemMonic_1(sfixn degA, sfixn * APtr, sfixn degB, sfixn * BPtr, 
        MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn i,j;
    sfixn tmp, p, R, SFT;

    if((degA-degB)<0) { return;}

    p=pPtr->P; R=(1L<<pPtr->Rpow)%p; SFT=pPtr->Base_Rpow;
    R=MulMod(R,R,p)<<SFT;
    for(i=degA-degB; i>=0; i--){
        if (APtr[degB+i] != 0){
            tmp=MontMulMod_OPT2_AS_GENE(APtr[degB+i], R, pPtr)<<SFT;
            for(j=0; j<=degB; j++) 
                APtr[i+j]=SubMod(APtr[i+j], 
                        MontMulMod_OPT2_AS_GENE(BPtr[j],tmp,pPtr), pPtr->P);
        }
    } 
}

// note: good for all Fourier Primes.
/**
 * plainDiv_1: Computes the plain division of A by B, returning the
 *                quotient Q and the coefficient vector of the remainder R.
 * @degQ:          degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: (overwritten) Coefficient vector  of A (input) / R (output)
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 *               The coeff vector of A is overwritten by 
 *               that of the remainder. 
 * 
 * Return value: void
 **/   
void plainDiv_1(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    register sfixn i,j;


  /* The first line was here before April 1, 2013. Sardar Haque deleted degR as it is not in used.*/
    //sfixn degR, u, tmp, p, R, U, SFT;
      sfixn  u, tmp, p, R, U, SFT;

    if((degA-degB)<0) {
        for(j=0;j<=degQ;j++) QPtr[j]=0;
        return;
    }

 /* The first line was here before April 1, 2013. Sardar Haque deleted degR=degA as it is not in used.*/
    //degR=degA, u=inverseMod(BPtr[degB],pPtr->P);
	         u=inverseMod(BPtr[degB],pPtr->P);


    p=pPtr->P; R=(1L<<pPtr->Rpow)%p; SFT=pPtr->Base_Rpow;
    U=MulMod(u,R,p)<<SFT;
    R=MulMod(R, R, p)<<SFT;

    for(i=degA-degB; i>=0; i--){
        if (APtr[degB+i] != 0){
            QPtr[i]=MontMulMod_OPT2_AS_GENE(APtr[degB+i],U,pPtr);
            tmp=MontMulMod_OPT2_AS_GENE(QPtr[i], R, pPtr)<<SFT;
            for(j=0; j<=degB; j++) 
                APtr[i+j]=SubMod(APtr[i+j], 
                        MontMulMod_OPT2_AS_GENE(BPtr[j],tmp,pPtr), pPtr->P);
        }
        else{QPtr[i]=0;}
    } 
}


// note: good for all Fourier Primes.
/**
 * plainDiv: Computes the plain division of A by B, returning the
 *                quotient Q and the coefficient vector of the remainder R.
 * @RPtr: (output) Coefficient vector  of the remainder
 * @degQ:          degree of the quotient
 * @QPtr: (output) Coefficient vector  of the quotient
 * @degA: degree of A
 * @APtr: Coefficient vector  of A 
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 * 
 *               This works for all Fourier Primes.
 * 
 * Return value: void
 **/   
void plainDiv(sfixn * RPtr,sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    sfixn * tmpA=(sfixn *)my_calloc(degA+1, sizeof(sfixn));
    register sfixn i, degR;
    if(degA<degB) return;
    for(i=0; i<=degA; i++) tmpA[i]=APtr[i];
    plainDiv_1(degQ, QPtr, degA, tmpA, degB, BPtr, pPtr);
    degR=degB-1;
    while((! tmpA[degR]) && (degR>0)) degR--;
    for(i=0; i<=degR; i++) RPtr[i]=tmpA[i];
    my_free(tmpA);
}

void plainRem(sfixn *degRAddr, sfixn * RPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    sfixn * tmpA, *tmpQ, degQ;
    register sfixn i,j;
    if(degA<degB) {for(j=0;j<=degA;j++) RPtr[j]=APtr[j]; return;}
    tmpA=(sfixn *)my_calloc(degA+1, sizeof(sfixn));
    for(i=0; i<=degA; i++) tmpA[i]=APtr[i];
    degQ=degA-degB;
    tmpQ=(sfixn *)my_calloc(degQ+1, sizeof(sfixn));
    plainDiv_1(degQ, tmpQ, degA, tmpA, degB, BPtr, pPtr);
    while((! tmpA[* degRAddr]) && ((*degRAddr)>0)) (*degRAddr)--;
    for(i=0; i<=(*degRAddr); i++) RPtr[i]=tmpA[i];
    my_free(tmpQ);
    my_free(tmpA);
}

// note: good for all Fourier Primes.
void plainQuo(sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    sfixn * tmpA=(sfixn *)my_calloc(degA+1, sizeof(sfixn));
    register sfixn i;
    if(degA<degB) return;
    for(i=0; i<=degA; i++) tmpA[i]=APtr[i];
    plainDiv_1(degQ, QPtr, degA, tmpA, degB, BPtr, pPtr);
    my_free(tmpA);
}


void plainDivNew(sfixn * RPtr,sfixn degQ, sfixn * QPtr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr )
{
    sfixn * tmpA=(sfixn *)my_calloc(degA+1, sizeof(sfixn));
    register sfixn i, degR;
    if(degA<degB) { 
        //printf("Attension: degA<degB!\n");
        for(i=0;i<=degA;i++) RPtr[i]=APtr[i];
        for(i=0;i<=degQ;i++) QPtr[i]=0;
        return;
    }
    for(i=0; i<=degA; i++) tmpA[i]=APtr[i];
    plainDiv_1(degQ, QPtr, degA, tmpA, degB, BPtr, pPtr);
    degR=degB-1;
    while((! tmpA[degR]) && (degR>0)) degR--;
    for(i=0; i<=degR; i++) RPtr[i]=tmpA[i];
    my_free(tmpA);
}

/**
 * PlainPseudoRemainder: Computes the pseudo-remainder of A by B.
 * @degRAddr: (output) degree of the pseudo-remainder
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 *               Using classical pseudo-division.
 *               This works for all Fourier Primes.
 * 
 * Return value: void
 **/   
void PlainPseudoRemainder(sfixn *degRAddr, sfixn * RPtr, sfixn degA, 
        sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    int i;
    sfixn * monicBPtr, Blc, BlcInv, e, co;
    if(degA<degB) {
        //printf("In PlainPseudoRemainder(), degA is smaller than degB which is unexpected!\n");
        //fflush(stdout);
        Interrupted=1;
        return;
    } 
    monicBPtr=(sfixn *)my_calloc(degB+1, sizeof(sfixn));
    Blc=BPtr[degB];
    BlcInv=inverseMod(Blc, pPtr->P);
    for(i=0; i<degB; i++) monicBPtr[i]=MulMod(BlcInv, BPtr[i], pPtr->P);
    monicBPtr[degB]=1;
    plainRem(degRAddr, RPtr, degA, APtr, degB, monicBPtr, pPtr);
    e=degA - degB +1;
    co=PowerMod(Blc, e, pPtr->P);
    for(i=0; i<=(*degRAddr); i++) RPtr[i]=MulMod(co, RPtr[i], pPtr->P);
    my_free(monicBPtr);
}


/**
 * PlainPseudoQuotient: Computes the pseudo-quotient of A by B.
 * @degQAddr: (output) degree of the pseudo-remainder
 * @RPtr: (output) Coefficient vector for the reaminder (size degB)
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 *               Using classical pseudo-division.
 *               This works for all Fourier Primes.
 * 
 * Return value: void
 **/   
void PlainPseudoQuotient(sfixn *degQAddr, sfixn * QPtr, sfixn degA, 
        sfixn * APtr, sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    int i;
    sfixn * monicBPtr, Blc, BlcInv, e, co;
    if(degA<degB) {
        //printf("In PlainPseudoRemainder(), degA is smaller than degB which is unexpected!\n");
        //fflush(stdout);
        Interrupted=1;
        return;
    } 
    monicBPtr=(sfixn *)my_calloc(degB+1, sizeof(sfixn));
    Blc=BPtr[degB];
    BlcInv=inverseMod(Blc, pPtr->P);
    for(i=0; i<degB; i++) monicBPtr[i]=MulMod(BlcInv, BPtr[i], pPtr->P);
    monicBPtr[degB]=1;
    plainQuo(*degQAddr, QPtr, degA, APtr, degB, monicBPtr, pPtr);

    e=degA - degB +1;
    co=PowerMod(Blc, e, pPtr->P);
    for(i=0; i<=(*degQAddr); i++) QPtr[i]=MulMod(co, QPtr[i], pPtr->P);
    my_free(monicBPtr);
}

/**
 * EX_PQuo_Uni: Computes the pseudo-quotient of A by B.
 * @degQAddr: (output) degree of the pseudo-remainder
 * @degA: degree of A
 * @APtr: Coefficient vector  of A
 * @degB: degree of B
 * @BPtr: Coefficient vector  of B 
 * @pPtr: prime number structure
 *               Using classical pseudo-division.
 *               This works for all Fourier Primes.
 * 
 * Return value: Coefficient vector for the pseudo-reaminder (accurate size)
 **/  
sfixn * EX_PQuo_Uni(sfixn *degQAddr, sfixn degA, sfixn * APtr, 
        sfixn degB, sfixn * BPtr, MONTP_OPT2_AS_GENE * pPtr)
{
    int i;
    sfixn * tmpQPtr, *QPtr;
    *degQAddr=degA;
    tmpQPtr=(sfixn *)my_calloc((*degQAddr)+1, sizeof(sfixn));
    PlainPseudoQuotient(degQAddr, tmpQPtr, degA, APtr, degB, BPtr, pPtr);
    *degQAddr=shrinkDegUni(*degQAddr, tmpQPtr);
    QPtr=(sfixn *)my_calloc((*degQAddr)+1, sizeof(sfixn));
    for(i=0; i<=(*degQAddr); i++) QPtr[i]=tmpQPtr[i];
    my_free(tmpQPtr);
    return QPtr;
}

/////////////////////////////END OF FILE //////////////////////////////////////
