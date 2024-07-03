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


#include "utilityFunctions.h"




extern void egcd (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo)  {
    sfixn tmp;
    sfixn A,B,C,D,u,v,q;

    u = y; v = x; 
    A=1; B=0; 
    C=0; D=1;

    do {
        q = u / v;
        tmp = u;
        u = v;
        v = tmp - q*v;
        tmp = A;
        A = B;
        B = tmp - q*B;
        tmp = C;
        C = D;
        D = tmp - q*D;
    } while (v != 0);
    *ao=A; 
    *bo=C; 
    *vo=u;
}




sfixn PowerMod(sfixn a, sfixn ee, sfixn n)
{
    sfixn x, y;

    usfixn e;

    if (ee == 0) return 1;

    if (ee < 0)
        e = - ((usfixn) ee);
    else
        e = ee;

    x = 1;
    y = a;
    while (e) {
        if (e & 1) x = MulMod(x, y, n);
        y = MulMod(y, y, n);
        e = e >> 1;
    }

    if (ee < 0) x = inverseMod(x, n);

    return x;
}


sfixn
MontMulMod_OPT2_AS_GENE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr){
    sfixn q2=pPtr->c_sft, c;

    // b has been pre-leftShifted Base_Rpow bits.
    MulHiLoUnsigned(&a, &b);
    // now a keeps M1 quo R.  
    // q2=c<<{Npow}.
    MulHiLoUnsigned(&b, &q2);
    // now b keeps M3 quo R;
    a-=b;
    // load c;
    c=pPtr->c;
    // be careful the sign of q2.
    q2=(usfixn)q2>>pPtr->Base_Npow;
    q2*=c;
    // now q2 keeps M2 quo R;
    // load P.
    c=pPtr->P;
    a+=q2;
    // compute (M1+M2-M3) quo R.
    // the result is in range  -p<a<2p.
    // we need following 3 lines to reduce the error.
    a += (a >> BASE_1) & c;
    a -= c;
    a += (a >> BASE_1) & c;
    return a;
}


sfixn
MontMulMod_OPT2_AS_Double_GENE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr){
    sfixn q2=pPtr->c_sft, q22=q2;

    // b has been pre-leftShifted Base_Rpow bits.
    MulHiLoUnsigned(&a, &b);
    MulHiLoUnsigned(&b, &q2);
    a-=b;
    MulHiLoUnsigned(&x, &y);
    MulHiLoUnsigned(&y, &q22);
    x-=y;

    y=pPtr->Base_Npow;
    b=pPtr->c;

    // be careful the sign of q2.
    q2=(usfixn)q2>>y;
    q2*=b;
    a+=q2;

    q22=(usfixn)q22>>y;
    y=pPtr->P;
    q22*=b;
    x+=q22;
    // now q2 keeps M2 quo R;
    // load P.


    // compute (M1+M2-M3) quo R.
    // the result is in range  -p<a<2p.
    // we need following 3 lines to reduce the error.


    a += (a >> BASE_1) & y;
    a -= y;
    a += (a >> BASE_1) & y;


    x += (x >> BASE_1) & y;
    x -= y;
    x += (x >> BASE_1) & y;

    *r2nd=x;

    return a;
}


sfixn
MontMulMod_OPT2_AS_GENE_SPE(sfixn a, sfixn b, MONTP_OPT2_AS_GENE * pPtr){
    sfixn cpow=pPtr->c_pow, q3=pPtr->R_Npow, p=pPtr->Base_Rpow; 
    // b has been pre-leftShifted Base_Rpow bits.
    MulHiLoUnsigned(&a, &b);
    // now a keeps M1 quo R.

    b=((usfixn)b)>>p;
    b+=b<<cpow;
    p=pPtr->R_Nmask;
    q3=b>>q3;
    a-=q3;
    q3=pPtr->N2_Rpow;
    b&=p;
    b+=b<<cpow;
    b<<=q3;
    p=pPtr->P;
    a+=b;
    // the result is in range  -p<a<2p.
    // we need following 3 lines to reduce the error.
    a += (a >> BASE_1) & p;
    a -= p;
    a += (a >> BASE_1) & p;
    return a;
}


sfixn
MontMulMod_OPT2_AS_Double_GENE_SPE(sfixn * r2nd, sfixn a, sfixn b, sfixn x, sfixn y, MONTP_OPT2_AS_GENE * pPtr){

    sfixn cpow=pPtr->c_pow, q3=pPtr->R_Npow,  p=pPtr->Base_Rpow; 
    // b has been pre-leftShifted Base_Rpow bits.
    MulHiLoUnsigned(&a, &b);
    b=((usfixn)b)>>p;
    b+=b<<cpow;

    MulHiLoUnsigned(&x, &y);
    // now a keeps M1 quo R.
    y=((usfixn)y)>>p;
    y+=y<<cpow;

    p=pPtr->R_Nmask;


    a-=b>>q3;
    x-=y>>q3;

    q3=pPtr->N2_Rpow;


    b&=p;
    y&=p;

    b+=b<<cpow;
    b<<=q3;
    a+=b;

    y+=y<<cpow;
    y<<=q3;

    p=pPtr->P;

    x+=y;
    // the result is in range  -p<a<2p.
    // we need following 3 lines to reduce the error.
    a += (a >> BASE_1) & p;
    a -= p;
    a += (a >> BASE_1) & p;

    x += (x >> BASE_1) & p;
    x -= p;
    x += (x >> BASE_1) & p;

    *r2nd=x;

    return a;

}


void MultiNumbPolyMul_1(sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr){
    register sfixn i;
    sfixn R;
    R=convertToMondMulR(r, pPtr);
    for(i=0; i<SIZ(f); i++)  DATI(f, i)=MontMulMod_OPT2_AS_GENE( DATI(f, i), R, pPtr);
}

 
void fromtofftRep(sfixn N, sfixn * rccum, sfixn * res, sfixn * ccum,  sfixn * dgs, sfixn * coeffs){
    int i;
    sfixn d;
    int tmpRes=0, tmpCoeffs=0;
    //printf("N=%ld, Using fromtofftRep().\n", N);
    if(N==0){
        res[0]=coeffs[0];
        return;
    }
    d=dgs[N];
    if(N==1){
        for(i=0; i<=d; i++){ res[i]=coeffs[i];}
        return;
    }
    for(i=0; i<=d; i++){
        tmpCoeffs=i*ccum[N];
        tmpRes=i*rccum[N];
        fromtofftRep(N-1, rccum, res+tmpRes, ccum, dgs, coeffs+tmpCoeffs); 
    }
}


int zeroCoefp(sfixn * coefPtr, sfixn coefSiz){
    register int i;
    for(i=0;i<coefSiz;i++){ if (coefPtr[i]) return 0;}
    return 1;
}



sfixn shrinkDeg(sfixn deg, sfixn * poly, sfixn coefSiz){
    sfixn * tmpPtr=poly+deg*coefSiz;
    while((zeroCoefp(tmpPtr, coefSiz))&& deg>0) {deg--; tmpPtr-=coefSiz;}
    return deg;
}

 
sfixn shrinkDegUni(sfixn deg, sfixn * cof){
    while((deg>0) && (!cof[deg])) deg--;
    return deg;
}


preFFTRep *setLeadingCoefMultiOne(sfixn N, preFFTRep* f){
    sfixn d;
    register sfixn i;
    d=shrinkDeg(BUSZSI(f, N), DAT(f), CUMI(f, N));
    nextMCoefData(f,N,d);
    for(i=0; i<CUMI(f,N); i++) DATI(f, i)=0;
    DATI(f,0)=1;
    DAT(f)=DEFDAT(f);
    return f;
}

void MultiNumbPolyMulMonicize_1(sfixn N, sfixn r, preFFTRep * f,  MONTP_OPT2_AS_GENE * pPtr){
    MultiNumbPolyMul_1(r, f,  pPtr);
    setLeadingCoefMultiOne(N,f);
}

void subEqDgPoly_inner_1(sfixn N, sfixn * dgs, sfixn * accum, sfixn * data1, sfixn * data2, sfixn p, int selector)
{
    int i, offset=0;
    if(N==1){
        if(selector==1){
            for(i=0;i<=dgs[1];i++) data1[i]=SubMod(data1[i],data2[i],p);}
        else{
            for(i=0;i<=dgs[1];i++) data2[i]=SubMod(data1[i],data2[i],p);}
        return;}
        for(i=0; i<=dgs[N]; i++){
            offset=accum[N]*i;
            subEqDgPoly_inner_1(N-1, dgs, accum, data1+offset, data2+offset, p, selector);
        } 
}

void subEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p, int selector){
    subEqDgPoly_inner_1(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p, selector);
}


sfixn zeroPolyp(preFFTRep * polyPtr){
    register int i;
    if((!(N(polyPtr))) && (!(DATI(polyPtr,0)))) return 1; 
    for(i=0;i<(SIZ(polyPtr));i++){ if (DATI(polyPtr, i)) return 0;}
    return 1;
}




sfixn constantPolyp(preFFTRep * polyPtr){
    register int i;
    if((!(N(polyPtr))) && ((DATI(polyPtr, 0))==1)) return 1; 
    for(i=1;i<(SIZ(polyPtr));i++){ if (DATI(polyPtr, i)) return 0;}
    /* change by Anis on September 24th 2012
     At this point this function returns 1
     But it fails to differentiate between zero poly and constant poly.
     Example, an univariate zero poly fails the first if condition. And
     in the loop it is not executed anything. So return 1.
     it should be like the following:
    */
    if( DATI(polyPtr, 0) != 0 ) return 1; 
    return 0;
}


void PolyCleanData(preFFTRep * prePtr){
    register int i;
    for(i=0;i<SIZ(prePtr);i++){
        DATI(prePtr, i)=0;
    }
}

void negatePoly_1(preFFTRep * prePtr, sfixn p){
    register int i;
    for(i=0;i<SIZ(prePtr);i++){
        DATI(prePtr, i)=p-DATI(prePtr, i);
    }
}

void addEqDgPoly_inner(sfixn N, sfixn *dgs, sfixn *accum, sfixn *data1, sfixn *data2, sfixn p)
{
    int i, offset=0;
    if(N==1){
        for(i=0;i<=dgs[1];i++) data1[i]=AddMod(data1[i],data2[i],p);
        return;}
        for(i=0; i<=dgs[N]; i++){
            offset=accum[N]*i;
            addEqDgPoly_inner(N-1, dgs, accum, data1+offset, data2+offset, p);
        } 
}

void addEqDgPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
    addEqDgPoly_inner(N, BUSZS(Ptr1), CUM(Ptr1), DAT(Ptr1), DAT(Ptr2), p);
}

void subPoly_inner_1 (sfixn N, sfixn * accum1, sfixn * dgs2, sfixn * accum2, sfixn * data1, sfixn * data2, sfixn p)
{
    int i, offset1=0, offset2=0;
    if(N==1){
        for(i=0;i<=dgs2[1];i++) data1[i]=SubMod(data1[i],data2[i],p);
        return;}
        for(i=0; i<=dgs2[N]; i++){
            offset1=accum1[N]*i;
            offset2=accum2[N]*i;
            subPoly_inner_1(N-1, accum1, dgs2, accum2, data1+offset1, data2+offset2, p);
        } 

}

void subPoly_1(sfixn N, preFFTRep * Ptr1, preFFTRep * Ptr2, sfixn p){
    subPoly_inner_1(N, CUM(Ptr1), BUSZS(Ptr2), CUM(Ptr2), DAT(Ptr1), DAT(Ptr2), p);
}

void increaseOneDim(preFFTRep * Ptr){
    (N(Ptr))++;
    SIZ(Ptr)=(SIZ(Ptr))*( BUSZSI(Ptr, N(Ptr)) +1);
}

sfixn MontMulMod(sfixn a, sfixn b, MONTP_OPT2_AS_GENE *pPtr){
    sfixn brsft;
    brsft=MontMulMod_OPT2_AS_GENE(b, pPtr->R2BRsft, pPtr)<<(pPtr->Base_Rpow);
    return MontMulMod_OPT2_AS_GENE(a,brsft,pPtr);
}

void coMulVec_1(sfixn co, sfixn deg, sfixn *vec, MONTP_OPT2_AS_GENE *pPtr){
    int i;
    sfixn tmp;
    tmp=(MontMulMod_OPT2_AS_GENE(co,pPtr->R2BRsft,pPtr))<<pPtr->Base_Rpow;
    deg = shrinkDegUni(deg, vec);
    for(i=0; i<=deg; i++){
        vec[i] = MontMulMod_OPT2_AS_GENE(vec[i],tmp,pPtr);
    }

}

void coMulAddVec(sfixn co, sfixn deg, sfixn *vec1, sfixn *vec2, MONTP_OPT2_AS_GENE *pPtr){
    int i;
    sfixn tmp;
    tmp=(MontMulMod_OPT2_AS_GENE(co,pPtr->R2BRsft,pPtr))<<pPtr->Base_Rpow;
    deg = shrinkDegUni(deg, vec2);
    for(i=0; i<=deg; i++){
        vec1[i] = AddMod(vec1[i], MontMulMod_OPT2_AS_GENE(vec2[i], tmp, pPtr), pPtr->P);
    }

}

preFFTRep *getLeadingCoefMulti(sfixn N, preFFTRep* co, preFFTRep* f){
    sfixn d;
    d=shrinkDeg(BUSZSI(f, N), DAT(f), CUMI(f, N));
    backupData(f);
    decreaseOneDim(f); 
    nextMCoefData(f,N,d);
    fromtofftRep(N-1,  CUM(co), DAT(co), CUM(f), BUSZS(f), DAT(f));
    increaseOneDim(f);
    restoreData(f);
    return co;
}

preFFTRep *setCoefMulti(sfixn N, preFFTRep* f, preFFTRep* co, sfixn j){
    register sfixn i;
    nextMCoefData(f,N,j);

    for(i=0; i<SIZ(co); i++) DATI(f, i)=DATI(co, i);
    DAT(f)=DEFDAT(f);
    //printf("f=\n");
    //printPoly(f);
    //fflush(stdout);
    return co;
}
