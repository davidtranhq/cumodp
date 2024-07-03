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


#include "Types.h"
#include "generalFuncs.h"
#include "MultiDFFT.h"
#include "MPMMTS.h"
#include "FINTERP.h"
#include "Factorization.h"
#include "MapleCConverter.h"
#include "solve2.h"
#include <math.h>
#include "solve2.h"
#include <time.h>
#include "HGCD.h"

#include <stdio.h>
#include <sys/types.h>



extern sfixn FP[];
extern sfixn CUDA_TAG;


void resultant3_tst(int argc, char **argv) {
    sfixn p =  919601153;
    //sfixn p = 257;
    sfixn N = 3, d1 = 6, d2 = 6, d3 = 6;
    sfixn dgs1[] = {0, 6, 6, 6};
    sfixn dgs2[] = {0, 6, 6, 6};
    preFFTRep *F1;
    preFFTRep *F2;
    preFFTRep *result = NULL;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    double time0, time1;
    clock_t t0, t1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    if (argc == 2) {
        d3 = d2 = d1 = atoi(argv[1]);
        dgs1[1] = dgs2[1] = d1;
        dgs1[2] = dgs2[2] = d2;
        dgs1[3] = dgs2[3] = d3;
    } 

    F1 = EX_randomPoly(N, dgs1, p);
    F2 = EX_randomPoly(N, dgs2, p);
    
    CUDA_TAG = 1;
    t0 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t0 = clock() - t0;
    //EX_Poly_Print(F1);
    //EX_Poly_Print(F2);
    //EX_Poly_Print(result);
    EX_FreeOnePoly(result);
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;

    CUDA_TAG = 0;
    t1 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t1 = clock() - t1;
    EX_FreeOnePoly(result);
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    
    printf("d1 = %d, d2 = %d, d3 = %d, GPU = %.3f, CPU = %.3f\n", d1, d2, d3, time0, time1);
}

void resultant2_tst(int argc, char **argv) {
    sfixn p = 919601153;
    //sfixn p = 257;
    sfixn N = 2, d1 = 10, d2 = 10;
    sfixn dgs1[] = {0, 10, 10};
    sfixn dgs2[] = {0, 10, 10};
    preFFTRep *F1;
    preFFTRep *F2;
    preFFTRep *result = NULL;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    double time0, time1;
    clock_t t0, t1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    if (argc == 2) {
        d1 = atoi(argv[1]);
        d2 = d1; 
        dgs1[1] = dgs1[2] = d1;
        dgs2[1] = dgs2[2] = d2;
    }

    F1 = EX_randomPoly(N, dgs1, p);
    F2 = EX_randomPoly(N, dgs2, p);
    
    CUDA_TAG = 1;
    t0 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t0 = clock() - t0;
    //EX_Poly_Print(F1);
    //EX_Poly_Print(F2);
    //EX_Poly_Print(result);
    EX_FreeOnePoly(result);
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;

    CUDA_TAG = 0;
    t1 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t1 = clock() - t1;
    EX_FreeOnePoly(result);
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    
    printf("d1 = %d, d2 = %d, GPU = %.3f, CPU = %.3f\n", d1, d2, time0, time1);

}

void uniquo_tst(int argc, char ** argv){
    sfixn *A, dA = 9, *B, dB = 9, *Q, dQ;
    sfixn i;
    sfixn p = 919601153;
    //sfixn p = 257;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    clock_t t0, t1;
    double time0, time1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    if (argc == 2) {
        dB = atoi(argv[1]);
        dA = dB; 
    } else if (argc == 3) {
        dA = atoi(argv[1]);
        dB = atoi(argv[2]);
    }

    A = (sfixn *)my_malloc(sizeof(sfixn)*(dA + 1));
    B = (sfixn *)my_malloc(sizeof(sfixn)*(dB + 1));
    for (i = 0; i <= dA; ++i) A[i] = rand() % p;
    for (i = 0; i <= dB; ++i) B[i] = rand() % p;
    B[dB] = 1;
    
    dQ = dA - dB;
    Q = (sfixn *)my_malloc(sizeof(sfixn)*(dQ + 1));
    
    t0 = clock();
    fastQuo(dQ, Q, dA, A, dB, B, pPtr);
    t0 = clock() - t0;
    dQ = shrinkDegUni(dQ, Q);
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
    EX_printPolyUni(dQ, Q, 'x');

    dQ = dA - dB;
    t1 = clock(); 
    plainQuo(dQ, Q, dA, A, dB, B, pPtr);
    t1 = clock() - t1;
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    printf("dA = %d, dB = %d, fastQuo = %.3f, plainQuo = %.3f\n", dA, dB, time0, time1);
    EX_printPolyUni(dA, A, 'x');
    EX_printPolyUni(dB, B, 'x');
    EX_printPolyUni(dQ, Q, 'x');

    my_free(A);
    my_free(B);
    my_free(Q);
}

void HGCD_tst(int argc, char** argv) 
{
    //sfixn AA[] = { 22, 199, 234, 97, 205, 19, 200, 144, 13, 177};
    //sfixn BB[] = {112, 68, 245, 142, 35, 130, 124, 240, 174, 26};
    sfixn *A, dA = 9, *B, dB = 9, *G, dG;
    sfixn i;
    //sfixn p = 919601153;
    sfixn p = 257;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    clock_t t0, t1;
    double time0, time1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	/*
    if (argc == 2) {
        dA = dB = atoi(argv[1]);
    } else if (argc == 3) {
        dA = atoi(argv[1]);
        dB = atoi(argv[2]);
    }
	*/
	dA = 500;
	dB = 1000;
    A = (sfixn *)my_malloc(sizeof(sfixn)*(dA + 1));
    B = (sfixn *)my_malloc(sizeof(sfixn)*(dB + 1));
    for (i = 0; i <= dA; ++i) A[i] = rand() % p;
    for (i = 0; i <= dB; ++i) B[i] = rand() % p;
    //for (i = 0; i <= dA; ++i) A[i] = AA[i];
    //for (i = 0; i <= dB; ++i) B[i] = BB[i];
    //EX_printPolyUni(dA, A, 'x');
    //EX_printPolyUni(dB, B, 'x');

    t0 = clock();
    G = HalfGCD(&dG, A, dA, B, dB, pPtr);
    t0 = clock() - t0;
    EX_printPolyUni(dG, G, 'x');
    my_free(G);
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
      
    t1 = clock();
    G = EX_GCD_UNI(&dG, A, dA, B, dB, pPtr);
    t1 = clock() - t1;
    EX_printPolyUni(dG, G, 'x');
    my_free(G);
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    
    printf("dA = %d, dB = %d, HGCD = %.3f, PlainGCD = %.3f\n", dA, dB, time0, time1);

    my_free(A);
    my_free(B);
}






// n := 5;
// unknowns := [x, y];
// sys := [ x^(2*n) + a*y^n - y, y^(2*n) +b * x^n - x ];
void hass_solve2(int argc, char** argv) 
{
    sfixn i, j, a = 1, b = 1, n = 60, p = 919601153;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    sfixn dgs1[3] = {0, 120, 60};
    sfixn dgs2[3] = {0, 60, 120};
    LinkedQueue *resQueue = NULL;
    clock_t t0, t1;
    double time0, time1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    if (argc == 2) {
        n = atoi(argv[1]);
        dgs1[1] = dgs2[2] = 2 * n;
        dgs1[2] = dgs2[1] = n;
    } else if (argc == 4) {
        n = atoi(argv[1]);
        dgs1[1] = dgs2[2] = 2 * n;
        dgs1[2] = dgs2[1] = n;
        a = atoi(argv[2]);
        b = atoi(argv[3]);
    }

    preFFTRep *F1 = EX_InitOnePoly(2, dgs1);
    preFFTRep *F2 = EX_InitOnePoly(2, dgs2);

    (F1->data)[2*n] = 1;
    (F1->data)[2*n+1] = p - 1;
    (F1->data)[(2*n+1) * n] = a;

    (F2->data)[1] = p - 1;
    (F2->data)[n] = b;
    (F2->data)[(n + 1) * 2 * n] = 1;

    preFFTRep *zerop = CreateZeroPoly();

    t0 = clock();
    resQueue = modular_solve2_select(0, F2, F1, zerop, pPtr);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
    //EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f  ||  ", time0);

    //t1 = clock();
    //resQueue = modular_solve2(1, F2, F1, zerop, pPtr);
    //t1 = clock() - t1;
    //time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    ////EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    //EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    //printf(" total = %.3f\n", time1);

    EX_FreeOnePoly(F1);
    EX_FreeOnePoly(F2);
    EX_FreeOnePoly(zerop);
}


void newBinarySolve2(int argc, char** argv)
{
    sfixn p =  469762049;
    sfixn dgs1[3] = {0, 60, 60};
    sfixn dgs2[3] = {0, 60, 60};

	    CUDA_TAG = 1;
    sfixn whatUwant = 0, method = 1, i;
    FILE *a, *b, *c;

    if(argc > 1)
    {
	a = fopen("inputFile.txt","w");
	//fprintf(a,"/home/shaque/SysPool/Systems/%s/sys.txt ",argv[1]);
	//fprintf(a,"/home/shaque/SysPool/Systems/%s/sysinfo.txt",argv[1]);
	fprintf(a,"%s/sys.txt ",argv[1]);
	fprintf(a,"%s/sysinfo.txt",argv[1]);


	fclose(a);

	if(argc > 2) p = atoi(argv[2]);
	
	b = fopen("prime.txt","w");
	fprintf(b,"%d",p);
	fclose(b);

	i = system("maple -q binSolve.mpl");

	c = fopen("outputFile.txt","r");
	if(fscanf(c,"%d",&i) == 0) printf("invalid entries in outputFile.txt\n");
	fclose(c);
	//printf("%s\n",argv[1]);
	if(i == -1)
	{
		printf("check the polynomial system. It might not be a bivariate system..\n");
		return;
	}
	c = fopen("outputFile.txt","r");
	if(fscanf(c,"%d %d %d %d %d",&i, &dgs1[1], &dgs1[2], &dgs2[1], &dgs2[2]) == 0) printf("invalid entries in outputFile.txt\n");
	fclose(c);

    }		
    else
    {
	printf("no file is passed. returning...");
	return;
    }
	
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    preFFTRep *f1 = EX_randomPoly(2, dgs1, p);
    preFFTRep *f2 = EX_randomPoly(2, dgs2, p);
    preFFTRep *zerop = CreateZeroPoly();

    polyBothInHostDev *F1, *F2, *Zerop;
	a = fopen("Poly1.dat","r");
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f1, i)) == 0) printf("Polynomial is missing\n");
	}
	fclose(a);

	a = fopen("Poly2.dat","r");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f2, i))== 0) printf("Polynomial is missing\n");
	}
	fclose(a);

	//for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i) printf(" %d ",DATI(f1, i));
	//printf("\n");
	//for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i) printf(" %d ",DATI(f2, i));
	//printf("\n");
	//*/
   F1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F1->status = -1;
   F1->polyHost = f1; 

   F2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F2->status = -1;
   F2->polyHost = f2; 

   Zerop = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   Zerop->status = -1;
   Zerop->polyHost = zerop; 

   LinkedQueue *resQueue = EX_LinkedQueue_Init();

   clock_t t0, t1;
   double time0, time1;

    if(argc > 3) method = atoi(argv[3]);
	
    t0 = clock();
    ModularSolve2(method, F1, F2, Zerop, pPtr, resQueue);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
  
    //FILE *outFile = fopen("rc.txt","w");
    //fprintf(outFile,"p := %d;\n",p);
   // FILE_LinkedQueue_Print(resQueue, outFile);

    sfixn totalSize= totalSizeRC(resQueue);
    printf("\nTotal size: %d\n",totalSize);

    //fclose(outFile);
   //EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
   EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
   // printf(" total = %.3f  || \n ", time0);


    //i = system("maple -q verifyWithMaple.mpl");

    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);
    EX_FreeOnePoly(zerop);





}

void newSolve2(int argc, char** argv) 
{		
    sfixn p = 11;// 469762043;
    sfixn dgs1[3] = {0, 60, 60};
    sfixn dgs2[3] = {0, 60, 60};
    sfixn whatUwant = 0, method = 1, i;

    if(argc > 1) dgs1[1] = atoi(argv[1]);
    if(argc > 2) dgs1[2] = atoi(argv[2]);
    if(argc > 3) dgs2[1] = atoi(argv[3]);
    if(argc > 4) dgs2[2]= atoi(argv[4]);
    if(argc > 5) whatUwant = atoi(argv[5]);
    if(argc > 6) method = atoi(argv[6]);

    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    preFFTRep *f1 = EX_randomPoly(2, dgs1, p);
    preFFTRep *f2 = EX_randomPoly(2, dgs2, p);
    preFFTRep *zerop = CreateZeroPoly();

    polyBothInHostDev *F1, *F2, *Zerop;

    if(whatUwant == 0)
    {
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i)
		DATI(f1, i) = rand()% p;
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
		DATI(f2, i) = rand()% p;
    }
    else
    {
	FILE *a = fopen("PDs.dat","w");
	fprintf(a, "%d %d %d %d %d",p, dgs1[1], dgs1[2], dgs2[1], dgs2[2]);
	fclose(a);
	i = system("maple -q solve2Maple.mm");
	
	a = fopen("Poly1.dat","r");
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f1, i)) ==0) printf("ERROR in reading polynomial from file\n");
	}
	fclose(a);

	a = fopen("Poly2.dat","r");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f2, i)) == 0) printf("ERROR in reading polynomial from file\n");
	}
	fclose(a);
    }	
	///*
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i) printf(" %d ",DATI(f1, i));
	printf("\n");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i) printf(" %d ",DATI(f2, i));
	printf("\n");
	//*/
   F1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F1->status = -1;
   F1->polyHost = f1; 

   F2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F2->status = -1;
   F2->polyHost = f2; 

   Zerop = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   Zerop->status = -1;
   Zerop->polyHost = zerop; 

   LinkedQueue *resQueue = EX_LinkedQueue_Init();

   clock_t t0, t1;
   double time0, time1;

	
    t0 = clock();
    ModularSolve2(method, F1, F2, Zerop, pPtr, resQueue);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;

    EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f  || \n ", time0);
    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);
    EX_FreeOnePoly(zerop);


}



void random_dense_solve2(int argc, char** argv) 
{
    sfixn p = 469762049;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    sfixn dgs1[3] = {0, 60, 60};
    sfixn dgs2[3] = {0, 60, 60};
    LinkedQueue *resQueue = NULL;
    clock_t t0, t1;
    double time0, time1;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    if (argc == 2) {
        dgs1[1] = dgs1[2] = dgs2[1] = dgs2[2] = atoi(argv[1]);
    } else if (argc == 5) {
        dgs1[1] = atoi(argv[1]);
        dgs1[2] = atoi(argv[2]);
        dgs2[1] = atoi(argv[3]);
        dgs2[2] = atoi(argv[4]);
    }
    srand(time(NULL));
    preFFTRep *F1 = EX_randomPoly(2, dgs1, p);
    preFFTRep *F2 = EX_randomPoly(2, dgs2, p);
    preFFTRep *zerop = CreateZeroPoly();

    t0 = clock();
    resQueue = modular_solve2_select(0, F1, F2, zerop, pPtr);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
    //EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f  ||  ", time0);

    t1 = clock();
    resQueue = modular_solve2_select(1, F1, F2, zerop, pPtr);
    t1 = clock() - t1;
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    //EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f\n", time1);

    EX_FreeOnePoly(F1);
    EX_FreeOnePoly(F2);
    EX_FreeOnePoly(zerop);
}

void solve2_tst() {
    sfixn p = 257;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    sfixn dgs1[3] = {0, 3, 4};
    sfixn dgs2[3] = {0, 2, 3};
    LinkedQueue *resQueue = NULL;
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    preFFTRep *F1 = EX_randomPoly(2, dgs1, p);
    preFFTRep *F2 = EX_randomPoly(2, dgs2, p);
    preFFTRep *zerop = CreateZeroPoly();

    printf("Input poly 1 : ");
    EX_Poly_Print(F1);
    printf("\nInput poly 2 : ");
    EX_Poly_Print(F2);
    printf("\n");

    resQueue = modular_solve2_select(1, F1, F2, zerop, pPtr);

    EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);

    EX_FreeOnePoly(F1);
    EX_FreeOnePoly(F2);
    EX_FreeOnePoly(zerop);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
}

void poly_tst() {
    sfixn p = 257;
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    sfixn dgs[] = {0, 3, 3};
    preFFTRep *zerop, *init, *tail, *poly1;
    sfixn i, j;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    poly1 = EX_InitOnePoly(2, dgs);
    for (j = 0; j <= dgs[2]; ++j) {
        for (i = 0; i <= dgs[1]; ++i) {
            poly1->data[j * (dgs[1] + 1) + i] = i;
        }
    }

    init = EX_getInitial(poly1);
    tail = EX_GetPolyTail(poly1);
    zerop = CreateZeroPoly();
    
    EX_Poly_Print(poly1);
    EX_Poly_Print(init);
    EX_Poly_Print(tail);

    EX_FreeOnePoly(init);
    EX_FreeOnePoly(tail);
    EX_FreeOnePoly(poly1);
    EX_FreeOnePoly(zerop);
}


void regcd() {
    sfixn p = 40961;
    sfixn dgs1[3] = {0, 2, 3};
    sfixn dgs2[3] = {0, 2, 2};
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    TriSet *T;
    sfixn M = 2;
    sfixn N = 3;
    sfixn dgbound = 3;
    LinkedQueue *resQueue;

    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
    T = EX_initRandomTriSet(N, M, pPtr);

    preFFTRep *f1 = EX_randomPoly(M, dgs1, p);
    preFFTRep *f2 = EX_randomPoly(M, dgs2, p);
    preFFTRep *R = EX_Resultant_Multi(f1, f2, M, pPtr);
    
    EX_FreeOnePoly(ELEMI(T, 1));
    ELEMI(T, 1) = R;
    BDSI(T, 1) = BUSZSI(R, 1) - 1; // main degree of R is BUSZSI(R, 1) 
    

    resQueue = EX_RegularGcd_Wrapped(f1, f2, T, M, pPtr);

    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);
    EX_freeTriSet(T);
    EX_LinkedQueue_Free(resQueue, EX_RegularPair_Free);
}




// testing Univariate Multiplication (Plain, FFT, TFT).
void
UniMul(){
  
  sfixn p;
  int isOne1, isOne2;
  sfixn degbound, degA, *A, degB, *B, degRES, *RES1, *RES2, *RES3, e, n;
  MONTP_OPT2_AS_GENE  prime;

  printf("Testing univariate multiplication ... ");
  srand(getSeed());
  p=FP[(rand()%100)];
  EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
  degbound=2000;
  degA= rand()%degbound;
  degB= rand()%degbound;
  while (degA<2) degA = rand()%degbound;
  while (degB<2) degB = rand()%degbound;

  A = EX_RandomUniPolyCoeffsVec(degA, p);
  B = EX_RandomUniPolyCoeffsVec(degB, p);
  degRES =degA+degB;

  RES1 = (sfixn *)my_calloc(degRES+1, sizeof(sfixn));
  // plain univariate multiplicaton.
  EX_Mont_PlainMul_OPT2_AS_GENE(degRES, RES1, degA, A, degB, B, &prime);

  RES2 = (sfixn *)my_calloc(degRES+1, sizeof(sfixn));
  // TFT univariate multiplicaton.
  EX_Mont_TFTMul_OPT2_AS_GENE(RES2, degA, A, degB, B, &prime);

  RES3 = (sfixn *)my_calloc(degRES+1, sizeof(sfixn));
  // FFT univariate multiplicaton.
  e=logceiling(degRES+1);
  n=1<<e;
  EX_Mont_FFTMul_OPT2_AS_GENE(n, e, degRES, RES3, degA, A, degB, B, &prime);


  
  isOne1 = compareVec(degRES, RES1, RES2);
  isOne2 = compareVec(degRES, RES2, RES3); 





  if (isOne1 != 1){
    printf("TFT-based Univeraiate Multiplication breaks!\n");
    //fflush(stdout);
    my_free(A); my_free(B);
    my_free(RES1); my_free(RES2); my_free(RES3);
    #ifndef _mcompile_
    Throw 10000;
    #else
    MapleRaiseError(modpn_saved_kv, "TFT-based Univeraiate Multiplication breaks!");
    #endif
  }

  if (isOne2 != 1){
    printf("FFT-based Univeraiate Multiplication breaks!\n");
    //fflush(stdout);
    my_free(A); my_free(B);
    my_free(RES1); my_free(RES2); my_free(RES3);
    #ifndef _mcompile_
    Throw 10000;
    #else
    MapleRaiseError(modpn_saved_kv, "FFT-based Univeraiate Multiplication breaks!");
    #endif
  }

  my_free(A); my_free(B);
  my_free(RES1); my_free(RES2); my_free(RES3); 

  printf("done.\n");

}






void 
UniDiv(){  
  sfixn p;
  int isOne, i;
  sfixn degbound, degA, *A, degB, *B, degRES, *RES, degRem, *Rem, degQuo, *Quo;
  MONTP_OPT2_AS_GENE  prime;

  printf("Testing univariate division ... ");

  srand(getSeed());
  p=FP[(rand()%100)];
  EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
  degbound=10;
  degA= rand()%degbound;
  degB= rand()%degbound;
  while (degA<2) degA = rand()%degbound;
  while (degB>degA) degB = rand()%degbound;



  A = EX_RandomUniPolyCoeffsVec(degA, p);
  B = EX_RandomUniPolyCoeffsVec(degB, p);
  
  // the divisor needs to be monic when using fast division.
  B[degB] = 2000;

  Rem = EX_UniRem(&degRem, degA, A, degB, B, &prime);
  Quo = EX_UniQuo(&degQuo, degA, A, degB, B, &prime);

  degRES = degB+degQuo;

  RES = (sfixn *)my_calloc(degRES+1, sizeof(sfixn));

  EX_Mont_PlainMul_OPT2_AS_GENE(degRES, RES, degB, B, degQuo, Quo, &prime);


  for(i=0; i<=degRem; i++) RES[i]=AddMod(RES[i], Rem[i],p);

  isOne = compareVec(degA, A, RES);  
  

  if(isOne != 1){
    printf("Univariate Euclidean division breaks!\n");
    //fflush(stdout);
    my_free(A); my_free(B); my_free(RES);
    my_free(Rem); my_free(Quo);
    #ifndef _mcompile_
    Throw 10000;
    #else
    MapleRaiseError(modpn_saved_kv, "Univariate Euclidean division breaks!");	
    #endif
  }
  
  my_free(A); my_free(B); my_free(RES);
  my_free(Rem); my_free(Quo);

  printf("done.\n");

}






int ExtendedGcd(){
  
  sfixn p,dbound, *A, *B, *CX, *DX, *GX;
  sfixn  dA, dB, dCX, dDX, dGX;
  sfixn dAA, *AA, dBB, *BB, dGcd, *Gcd;
  int isOne;
  MONTP_OPT2_AS_GENE  prime;

  printf("Testing Extend Gcd ... ");

  srand(getSeed());
	
  p=FP[(rand()%100)];

//p= 9001;
  EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

  dbound=1000;

  dAA = rand()%dbound;
  while(dAA<1) dAA = rand()%dbound;
  dBB = rand()%dbound;
  while(dBB>dAA) dBB = rand()%dbound;
  dGcd = rand() %dbound;
  while(dGcd<1) dGcd = rand()%dbound;


  dA = dAA+dGcd;
  dB = dBB+dGcd;


  A = (sfixn *) my_calloc(dA+1, sizeof(sfixn));
  B = (sfixn *) my_calloc(dB+1, sizeof(sfixn));

  
  CX = (sfixn *) my_calloc(dA+1, sizeof(sfixn));
  DX = (sfixn *) my_calloc(dA+1, sizeof(sfixn));
  GX = (sfixn *) my_calloc(dA+1, sizeof(sfixn));


  AA = EX_RandomUniPolyCoeffsVecMonic(dAA, p);
  BB = EX_RandomUniPolyCoeffsVecMonic(dBB, p);  
  Gcd = EX_RandomUniPolyCoeffsVecMonic(dGcd, p);   
 
 

  EX_Mont_PlainMul_OPT2_AS_GENE(dA, A, dAA, AA, dGcd, Gcd, &prime);
  EX_Mont_PlainMul_OPT2_AS_GENE(dB, B, dBB, BB, dGcd, Gcd, &prime);

  ExGcd_Uni(CX, &dCX, DX, &dDX, GX, &dGX,  A, dA, B, dB, &prime);

	//
	printf("\n %d %d %d %d %d \n", dCX, dDX, dGX, dA, dB);
	//

  isOne = compareVec(dGcd, Gcd, GX); // THIS IS BAD, Just Use Maple
 
  

  if(isOne != 1){
    printf("Extend Euclidean  breaks!\n");
    //fflush(stdout);
  my_free(A);
  my_free(B);
  my_free(CX);
  my_free(DX);
  my_free(GX);
  my_free(AA);
  my_free(BB);
  my_free(Gcd);
  #ifndef _mcompile_
    Throw 10000;
    #else
    MapleRaiseError(modpn_saved_kv, "Extend Euclidean  breaks!");	
  #endif
  }



  my_free(A);
  my_free(B);
  my_free(CX);
  my_free(DX);
  my_free(GX);
  my_free(AA);
  my_free(BB);
  my_free(Gcd);

  printf("done.\n");


  return 0;
}










void
Reduce( ){
  preFFTRep *f1, *f2, *input, *output1, *output2;
  TriSet *ts;
  TriRevInvSet *tris;
  sfixn  N, dgbound, p;
  MONTP_OPT2_AS_GENE  prime;
  int isOne;

  printf("Testing Normalform ... ");

  srand(getSeed());
  p=FP[(rand()%100)];
  EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

  N=4;
  dgbound=4;
  ts = EX_initRandomTriSet(N, dgbound, &prime);


  tris = EX_getRevInvTriSet(N,  ts,  &prime);



  f1 = EX_randomPoly(N, BDS(ts), p);
  f2 = EX_randomPoly(N, BDS(ts), p);
  output1 = EX_InitOnePoly(N, BDS(ts));
  output2 = EX_InitOnePoly(N, BDS(ts));


  input = EX_EX_mulPoly(N, f1, f2, &prime);



  MultiMod_DF(N, output1, input, ts, tris, &prime);

  MultiMod_BULL(N, output2, input, ts, tris, &prime);

  isOne = EX_IsEqualPoly(output1, output2);

  if(isOne != 1){
    printf("MultiMod()  breaks!\n");
    //fflush(stdout);
    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);
    EX_FreeOnePoly(input);
    EX_FreeOnePoly(output1);
    EX_FreeOnePoly(output2);
    EX_freeTriSet(ts);
    EX_freeTriRevInvSet(tris);
   
  }


  EX_FreeOnePoly(f1);
  EX_FreeOnePoly(f2);
  EX_FreeOnePoly(input);
  EX_FreeOnePoly(output1);
  EX_FreeOnePoly(output2);
  EX_freeTriSet(ts);
  EX_freeTriRevInvSet(tris);
  
  printf("done.\n");

}






void Lifting(){ 

  int i, j;
  TriSet * ts;
  sfixn N=3, p;
  MONTP_OPT2_AS_GENE * pPtr=(MONTP_OPT2_AS_GENE *)my_malloc(sizeof(MONTP_OPT2_AS_GENE));
  sfixn *dgs=(sfixn *)my_calloc(3, sizeof(sfixn));
  sfixn Y;
  sfixn y0;
  sfixn *inDGS; sfixn *inSIZS; sfixn *inCOEFS;
  sfixn *GNS;
  int offset=0, totSiz=1;
  sfixn MDAS[39*3+43*3]={0,16,0,0,10,0,2,2,2,5,2,3,4,1,4,0,31,0,1,2,0,5,6,7,1,3,0,5,8,9,4,5,10,0,51,0,1,2,0,5,12,13,1,1,0,5,14,15,4,11,16,0,77,0,1,2,0,5,18,19,4,17,20,0,95,0,2,3,2,5,22,23,4,21,24,1,3,0,1,1,0,5,26,27,4,25,28,1,3,0,4,29,30,0,55,0,2,1,2,5,32,33,4,31,34,0,28,0,1,1,0,5,36,37,4,35,38,0,43,0,0,30,0,2,2,2,5,2,3,4,1,4,0,27,0,1,2,0,5,6,7,1,3,0,5,8,9,4,5,10,0,15,0,1,2,0,5,12,13,1,1,0,5,14,15,4,11,16,0,59,0,1,2,0,5,18,19,4,17,20,0,96,0,2,3,2,5,22,23,4,21,24,0,72,0,1,3,0,5,26,27,1,1,0,5,28,29,4,25,30,0,87,0,1,3,0,5,32,33,4,31,34,0,47,0,2,1,2,5,36,37,4,35,38,0,90,0,1,1,0,5,40,41,4,39,42};

  sfixn esti_dy=32;
  sfixn esti_siz=esti_dy*(dgs[1]+1)+esti_dy*(dgs[1]+1)*(dgs[2]+1);
  sfixn *outPDGVECS=(sfixn *) my_calloc(esti_siz, sizeof(sfixn));
  sfixn *outCOEFVECS=(sfixn *)my_calloc(esti_siz, sizeof(sfixn));
  int iter=0;

  printf("Testing lifting ... ");



  dgs[0]=1;
  dgs[1]=4;
  dgs[2]=1;
 
  srand(getSeed());
  p=FP[rand()%100];
  p = 469762049;
  
  EX_MontP_Init_OPT2_AS_GENE(pPtr,  p);


  // vec_slg = example_1_PolyVec_SLG();



  ts=(TriSet *)my_malloc(sizeof(TriSet));
  init_example_1_DenseTriSetForLifting_y0_10(N, dgs, ts, pPtr);

  


  Y=1;  y0=10;

  //final_ts = EX_UniNewtonLift(Y, y0, vec_slg, ts, N, pPtr);

  

  inDGS = (sfixn *) my_calloc(N*N,sizeof(sfixn));
  for(i=1; i<=N; i++){
    for(j=1; j<=N; j++){
      inDGS[(i-1)*N + j-1]=BUSZSI(ELEMI(ts, i), j);
    }
  }

  inSIZS = (sfixn *) my_calloc(N,sizeof(sfixn));
  for(i=1; i<=N; i++){
      inSIZS[i-1]=SIZ(ELEMI(ts, i));
  }



  for(i=0; i<N; i++){
    totSiz*=(inSIZS[i]);
  }

  inCOEFS= (sfixn *) my_calloc(totSiz,sizeof(sfixn));

  for(i=1; i<=N; i++){
    for(j=0; j<inSIZS[i-1]; j++){
      inCOEFS[offset+j]=(DAT(ELEMI(ts, i)))[j];
    }
    offset+=inSIZS[i-1];
  }



  GNS=(sfixn *)my_calloc(2, sizeof(sfixn));
  GNS[0]=39;
  GNS[1]=43;


  printf("before lifting\n");

  iter = NewtonLiftUniCN(outPDGVECS, outCOEFVECS, Y, y0, N, GNS, MDAS,  dgs, inDGS, inSIZS, inCOEFS, p);

  printf("after lifting ... ");

  my_free(outPDGVECS); my_free(outCOEFVECS);  
  my_free(GNS);  
  my_free(inDGS); my_free(inSIZS); my_free(inCOEFS);

  EX_freeTriSet(ts);

  my_free(dgs);
  my_free(pPtr);

  printf("done.\n");



}


int checkPlainUniMul(int argc, char *argv[])
{
	if(chdir("test_files/uni_mul")!= 0)
	{
		printf("Missing directory : test_files/uni_mul\n");	
		exit(0);
	}
	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	//n and m are the length of the polynomials
	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d",n, m, p);
	fclose(a);

	int i = system("maple -q testmul.mm");

	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *C = (sfixn *)my_calloc(n+m-1, sizeof(sfixn));


	a = fopen("PolyM.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}
		//
		//printf("%d ",M[i]);
		//

	}	
	fclose(a);
	//
	//printf("\n");
	//
	a = fopen("PolyN.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
		//
		//printf("%d ",N[i]);
		//
	}	
	fclose(a);
	//
	//printf("\n");
	//

	EX_Mont_PlainMul_OPT2_AS_GENE(n+m-2, C, n-1, N, m-1, M, &prime);
	
	a = fopen("PolyCmodpn.dat","w");
	for(i = 0; i < (n+m-1); ++i)
		fprintf(a,"%d ",C[i]);	
	fclose(a);
	
	my_free(M);
	my_free(N);
	my_free(C);
	
	return system("diff -b PolyC.dat PolyCmodpn.dat");

	
}


int checkFFTUniMul(int argc, char *argv[])
{
	if(chdir("test_files/uni_mul")!= 0)
	{
		printf("Missing directory : test_files/uni_mul\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	//n and m are the length of the polynomials
	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d",n, m, p);
	fclose(a);

	int i = system("maple -q testmul.mm");

	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *C = (sfixn *)my_calloc(n+m-1, sizeof(sfixn));


	a = fopen("PolyM.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}
		//
		//printf("%d ",M[i]);
		//

	}	
	fclose(a);
	//
	//printf("\n");
	//
	a = fopen("PolyN.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
		//
		//printf("%d ",N[i]);
		//
	}	
	fclose(a);
	//
	//printf("\n");
	//

	sfixn e=logceiling(n+m-1);
	sfixn nL = 1<<e;
	EX_Mont_FFTMul_OPT2_AS_GENE(nL, e, n+m-2, C, n-1, N, m-1, M, &prime);
	
	a = fopen("PolyCmodpn.dat","w");
	for(i = 0; i < (n+m-1); ++i)
		fprintf(a,"%d ",C[i]);	
	fclose(a);
	
	my_free(M);
	my_free(N);
	my_free(C);
	
	return system("diff -b PolyC.dat PolyCmodpn.dat");

	
}


int checkTFTUniMul(int argc, char *argv[])
{
	if(chdir("test_files/uni_mul")!= 0)
	{
		printf("Missing directory : test_files/uni_mul\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	//n and m are the length of the polynomials
	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d",n, m, p);
	fclose(a);

	int i = system("maple -q testmul.mm");

	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *C = (sfixn *)my_calloc(n+m-1, sizeof(sfixn));


	a = fopen("PolyM.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}
		//
		//printf("%d ",M[i]);
		//

	}	
	fclose(a);
	//
	//printf("\n");
	//
	a = fopen("PolyN.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
		//
		//printf("%d ",N[i]);
		//
	}	
	fclose(a);
	//
	//printf("\n");
	//

	EX_Mont_TFTMul_OPT2_AS_GENE(C, n-1, N, m-1, M, &prime);	

	a = fopen("PolyCmodpn.dat","w");
	for(i = 0; i < (n+m-1); ++i)
		fprintf(a,"%d ",C[i]);	
	fclose(a);
	
	my_free(M);
	my_free(N);
	my_free(C);
	
	return system("diff -b PolyC.dat PolyCmodpn.dat");

	
}

void benchmarkUniGCD(int argc, char *argv[]){
	sfixn p = 469762049, n = 100, m = 50, howMany = 10;
	if (argc > 2) p       = atoi (argv[2]);
	if (argc > 3) n       = atoi (argv[3]) +1;
	if (argc > 4) m       = atoi (argv[4]) + 1;
	if (argc > 5) howMany = atoi (argv[5]);
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *G;
	
	sfixn i, j;
	for(i = 0; i < m; ++i) M[i] = rand()%p;
	for(i = 0; i < n; ++i) N[i] = rand()%p;


	time_t t1,t2;
	(void) time(&t1);
	for(i = 0; i < howMany; ++i){ 
		G = EX_GCD_Uni_Wrapper(&j, N, n-1, M, m-1, &prime);	
		my_free(G);
	}
	(void) time(&t2);
	
	printf("CPU gcd: %d %d %d\n",n,m, (int) t2-t1);

	(void) time(&t1);
	for(i = 0; i < howMany; ++i){
		G = gcdCUDAStraight(&j, N, n-1, M, m-1, &prime);	
		my_free(G);
	}
	(void) time(&t2);
	
	printf("GPU gcd: %d %d %d\n",n,m, (int) t2-t1);


	my_free(M);
	my_free(N);
	

}

void benchmarkResultant(int argc, char *argv[]){
	sfixn p = 469762049, n1 = 100, n2 = 100,  m1 = 50, m2 = 50, howMany = 10;
	if (argc > 2) p       = atoi (argv[2]);
	if (argc > 3) n1       = atoi (argv[3]);
	if (argc > 4) n2       = atoi (argv[4]);
	if (argc > 5) m1       = atoi (argv[5]);
	if (argc > 6) m2       = atoi (argv[6]);
	if (argc > 7) howMany = atoi (argv[7]);
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn dgs1[3] = {0, n1, n2};
	sfixn dgs2[3] = {0, m1, m2};


        preFFTRep *f1 = EX_randomPoly(2, dgs1, p);
        preFFTRep *f2 = EX_randomPoly(2, dgs2, p);

	sfixn i;
	
	SCUBE *S;

	time_t t1,t2;
	(void) time(&t1);
	for(i = 0; i < howMany; ++i)  
	{
		S = EX_SubResultantChainSelect(0, f1, f2, 2, &prime);	
		EX_SCUBE_Free(S);	
	}
	(void) time(&t2);
	
	printf("CPU subresultant: %d %d %d %d %d\n",n1, n2, m1, m2, (int) t2-t1);

	(void) time(&t1);
	for(i = 0; i < howMany; ++i) 
	{
		S = EX_SubResultantChainSelect(1, f1, f2, 2, &prime);	
		EX_SCUBE_Free(S);	
	} 
	(void) time(&t2);
	
	printf("GPU subresultant: %d %d %d %d %d\n",n1,n2, m1, m2, (int) t2-t1);

	
        EX_FreeOnePoly(f1);
	EX_FreeOnePoly(f2);
}

void benchmarkUniDiv(int argc, char *argv[]){
	sfixn p = 469762049, n = 100, m = 50, howMany = 10;
	if (argc > 2) p       = atoi (argv[2]);
	if (argc > 3) n       = atoi (argv[3]) +1;
	if (argc > 4) m       = atoi (argv[4]) + 1;
	if (argc > 5) howMany = atoi (argv[5]);
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *Quo = (sfixn *)my_calloc(n - m + 1, sizeof(sfixn));
	
	sfixn i, j;
	for(i = 0; i < m; ++i) M[i] = rand()%p;
	for(i = 0; i < n; ++i) N[i] = rand()%p;

	time_t t1,t2;
	(void) time(&t1);
	for(i = 0; i < howMany; ++i)  Quo = EX_UniQuo(&j, n-1, N, m-1, M, &prime);	
	(void) time(&t2);
	
	printf("CPU div: %d %d %d\n",n,m, (int) t2-t1);

	(void) time(&t1);
	for(i = 0; i < howMany; ++i)  divCUDAStraight(&j, n-1, N, m-1, M, &prime);	
	(void) time(&t2);
	
	printf("GPU div: %d %d %d\n",n,m, (int) t2-t1);
		
		
	my_free(M);
	my_free(N);
	my_free(Quo);

}


void benchmarkUniIsDivided(int argc, char *argv[]){
	sfixn p = 469762049, n = 100, m = 50, howMany = 10;
	if (argc > 2) p       = atoi (argv[2]);
	if (argc > 3) n       = atoi (argv[3]) +1;
	if (argc > 4) m       = atoi (argv[4]) + 1;
	if (argc > 5) howMany = atoi (argv[5]);
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));

	
	sfixn i, j;
	for(i = 0; i < m; ++i) M[i] = rand()%p;
	for(i = 0; i < n; ++i) N[i] = rand()%p;
	if(M[m-1] == 0) M[m-1] = 1;

	time_t t1,t2;
	(void) time(&t1);
	for(i = 0; i < howMany; ++i)  j = is_divided_by(n-1, N, m-1, M, &prime);
	(void) time(&t2);
	
	printf("CPU isDividedBy: %d %d %d\n",n,m, (int) t2-t1);

	(void) time(&t1);
	for(i = 0; i < howMany; ++i)  j = is_divided_byCUDA(n-1, N, m-1, M, &prime);	
	(void) time(&t2);
	
	printf("GPU isDividedBy: %d %d %d\n",n,m, (int) t2-t1);
		
		
	my_free(M);
	my_free(N);
	

}




int checkUniDiv(int argc, char *argv[])
{
	if(chdir("test_files/uni_div")!= 0)
	{
		printf("Missing directory : test_files/uni_div\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	//n and m are the length of the polynomials
	if(m>n) 
	{
		printf("degree of first polynomial must be greater or equal to the second polynomial\n");
		exit(0);
	}
	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d",n, m, p);
	fclose(a);

	int i = system("maple -q testdiv.mm");
	
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *Quo = (sfixn *)my_calloc(n - m + 1, sizeof(sfixn));
	sfixn *Rem = (sfixn *)my_calloc(m-1, sizeof(sfixn));


	a = fopen("PolyM.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}
		//
		//printf("%d ",M[i]);
		//

	}	
	fclose(a);
	//
	//printf("\n");
	//
	a = fopen("PolyN.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
		//
		//printf("%d ",N[i]);
		//
	}	
	fclose(a);
	sfixn degRem, degQuo;
	Rem = EX_UniRem(&degRem, n-1, N, m-1, M, &prime);
	Quo = EX_UniQuo(&degQuo, n-1, N, m-1, M, &prime);

	
	a = fopen("PolyQuomodpn.dat","w");
	for(i = 0; i < (n-m+1); ++i)
	{
		if(i <= degQuo)	fprintf(a,"%d ",Quo[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);


	a = fopen("PolyRemmodpn.dat","w");
	for(i = 0; i < (m-1); ++i)
	{
		if(i <= degRem)	fprintf(a,"%d ",Rem[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);

	
	my_free(M);
	my_free(N);
	my_free(Rem);
	my_free(Quo);
		
	
	i = system("diff -b PolyQuomodpn.dat PolyQuo.dat");
	if(i != 0) printf("Quo has error\n");
	sfixn j = system("diff -b PolyRemmodpn.dat PolyRem.dat");
	if(j != 0) printf("Rem has error\n");
	return i+j;
}


int checkExtGCD(int argc, char *argv[])
{
	if(chdir("test_files/uni_exGCD")!= 0)
	{
		printf("Missing directory : test_files/uni_exGCD\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	sfixn h = m/2;
	if(n < m) h = n/2;
	if (argc > 5) h = atoi (argv[5]) + 1;
	if(h >= n && n <=m )
	{
		printf("the degree of common part must be less than the smaller polynomial\n");
		exit(0);
	}

	if(h >= m && m <=n )
	{
		printf("the degree of common part must be less than the smaller polynomial\n");
		exit(0);
	}

	//n and m are the length of the polynomials

	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d %d",p, n, m, h);
	fclose(a);

	int i = system("maple -q testextGCD.mm");
	
        sfixn big = n;
	if(m>n) big = m;

	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *G = (sfixn *)my_calloc(big, sizeof(sfixn));
	sfixn *s = (sfixn *)my_calloc(big, sizeof(sfixn));
	sfixn *t = (sfixn *)my_calloc(big, sizeof(sfixn));
	sfixn ds, dt, dG;

	//
	//printf("\n");
	//
	a = fopen("Poly1.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}
		//
		//printf("%d ",N[i]);
		//
	}	
	fclose(a);

	//
	//printf("\n");
	//
	a = fopen("Poly2.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
		//
		//printf("%d ",M[i]);
		//

	}	
	fclose(a);
	//
	//printf("\n");
	//
	

	/*
	printf("\n");
	for(i = 0; i < n; ++i)
		printf("%d ",N[i]);
	
	printf("\n");
	for(i = 0; i < m; ++i)
		printf("%d ",M[i]);
	printf("\n");

	*/
	//printf("\n%d %d %d %d %d\n",dG, ds, dt, n-1, m-1);
	//

	//


	ExGcd_Uni(s,  &ds,  t,  &dt,  G, &dG,    N, n-1,M, m-1,&prime);
	
	//
	//printf("\n%d %d %d %d %d\n",dG, ds, dt, n-1, m-1);
	//

	a = fopen("PolyGgpu.dat","w");
	for(i = 0; i < big; ++i)
	{
		if(i <= dG)	fprintf(a,"%d ",G[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);


	a = fopen("PolysGPU.dat","w");
	for(i = 0; i < big; ++i)
	{
		if(i <= ds)	fprintf(a,"%d ",s[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);

	a = fopen("PolySgpu.dat","w");
	for(i = 0; i < big; ++i)
	{
		if(i <= ds)	fprintf(a,"%d ",s[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);

	a = fopen("PolyTgpu.dat","w");
	for(i = 0; i < big; ++i)
	{
		if(i <= dt)	fprintf(a,"%d ",t[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);

	my_free(M);
	my_free(N);
	my_free(G);
	my_free(s);
	my_free(t);	

	i = system("diff -b PolyG.dat PolyGgpu.dat");
	sfixn j = system("diff -b Polys.dat PolySgpu.dat");
	sfixn k = system("diff -b Polyt.dat PolyTgpu.dat");

	if(i!=0){ printf("Error in GCD computation"); exit(0);}
	if(j!=0){ printf("Error in s  computation"); exit(0); }
	if(k!=0){ printf("Error in t computation"); exit(0); }

	
	return i+j+k;
}




int checkHGCD(int argc, char *argv[])
{
	if(chdir("test_files/uni_hgcd")!= 0)
	{
		printf("Missing directory : test_files/uni_hgcd\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	sfixn h = m/2;
	if(n < m) h = n/2;
	if (argc > 5) h = atoi (argv[5]) + 1;
	if(h >= n && n <=m )
	{
		printf("the degree of common part must be less than the smaller polynomial\n");
		exit(0);
	}

	if(h >= m && m <=n )
	{
		printf("the degree of common part must be less than the smaller polynomial\n");
		exit(0);
	}

	//n and m are the length of the polynomials

	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d %d",p, n, m, h);
	fclose(a);

	int i = system("maple -q testHGCD.mm");
	
        sfixn small = n;
	if(m<n) small = m;

	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *G;// = (sfixn *)my_calloc(small, sizeof(sfixn));

	a = fopen("Poly1.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}		
	}	
	fclose(a);

	a = fopen("Poly2.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
	}	
	fclose(a);

	sfixn dG;
	if(m >= n) G = HalfGCD(&dG, N, n-1, M, m-1, &prime);
	else       G = HalfGCD(&dG, M, m-1, N, n-1, &prime);
	sfixn inv = inverseMod(G[dG], p);
	for(i = 0; i <= dG; ++i) G[i] = MulMod(G[i], inv, p);

	a = fopen("PolyGgpu.dat","w");
	for(i = 0; i < small; ++i)
	{
		if(i <= dG)	fprintf(a,"%d ",G[i]);	
		else fprintf(a,"0 ");	
	}
	fclose(a);

	my_free(M);
	my_free(N);
	my_free(G);

	i = system("diff -b PolyG.dat PolyGgpu.dat");


	return i;
}


int checkUniquo(int argc, char *argv[])
{
	if(chdir("test_files/uni_uniquo")!= 0)
	{
		printf("Missing directory : test_files/uni_uniquo\n");	
		exit(0);
	}

	sfixn p = 469762049;
	if (argc > 2) p = atoi (argv[2]);
	sfixn n = 10;
	if (argc > 3) n = atoi (argv[3]) +1;
	sfixn m = 10;
	if (argc > 4) m = atoi (argv[4]) + 1;
	//n and m are the length of the polynomials
	if(m>n) 
	{
		printf("degree of first polynomial must be greater or equal to the second polynomial\n");
		exit(0);
	}
	FILE *a = fopen("MNP.dat","w");
	fprintf(a, "%d %d %d",n, m, p);
	fclose(a);

	int i = system("maple -q testdiv.mm");
	
	MONTP_OPT2_AS_GENE  prime;
	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);
	sfixn *M = (sfixn *)my_calloc(m, sizeof(sfixn));
	sfixn *N = (sfixn *)my_calloc(n, sizeof(sfixn));
	sfixn *Quo = (sfixn *)my_calloc(n - m + 1, sizeof(sfixn));



	a = fopen("PolyM.dat", "r");
	for(i = 0; i < m; ++i)
	{
		if(fscanf(a,"%d",&M[i]) == 0)
		{
			printf("Error in reading poly1\n");
			exit(0);
		}	
	}	
	fclose(a);
	a = fopen("PolyN.dat", "r");
	for(i = 0; i < n; ++i)
	{
		if(fscanf(a,"%d",&N[i]) == 0)
		{
			printf("Error in reading poly2\n");
			exit(0);
		}
	}	
	fclose(a);
	sfixn degRem, degQuo;
	fastQuo(n-m, Quo, n-1, N, m-1, M, &prime);
	
	
	a = fopen("PolyQuomodpn.dat","w");
	for(i = 0; i < (n-m+1); ++i)
	{
		fprintf(a,"%d ",Quo[i]);	
		//if(i <= degQuo)	fprintf(a,"%d ",Quo[i]);	
		//else fprintf(a,"0 ");	
	}
	fclose(a);


		
	my_free(M);
	my_free(N);	
	my_free(Quo);
		
	
	i = system("diff -b PolyQuomodpn.dat PolyQuo.dat");
	return i;
}






void BinSolveTest(int argc, char *argv[])
{
	if(chdir("test_files/bi_solve")!= 0)
	{
		printf("Missing directory : test_files/bi_solve\n");	
		exit(0);
	}


	sfixn p =  469762049;
	sfixn dgs1[3] = {0, 60, 60};
	sfixn dgs2[3] = {0, 60, 60};

	CUDA_TAG = 1;
	computingChoince myChoice = CPU;
	sfixn whatUwant = 0, i;
        computingChoince method = myChoice;
	FILE *a, *b, *c;

    if(argc > 2)
    {
	a = fopen("inputFile.txt","w");
	//fprintf(a,"/home/shaque/SysPool/Systems/%s/sys.txt ",argv[1]);
	//fprintf(a,"/home/shaque/SysPool/Systems/%s/sysinfo.txt",argv[1]);
	fprintf(a,"%s/sys.txt ",argv[2]);
	fprintf(a,"%s/sysinfo.txt",argv[2]);


	fclose(a);

	if(argc > 3) p = atoi(argv[3]);
	
	b = fopen("prime.txt","w");
	fprintf(b,"%d",p);
	fclose(b);

	i = system("maple -q binSolve.mpl");

	c = fopen("outputFile.txt","r");
	if(fscanf(c,"%d",&i) == 0) printf("invalid entries in outputFile.txt\n");
	fclose(c);
	//printf("%s\n",argv[1]);
	if(i == -1)
	{
		printf("check the polynomial system. It might not be a bivariate system..\n");
		return;
	}
	c = fopen("outputFile.txt","r");
	if(fscanf(c,"%d %d %d %d %d",&i, &dgs1[1], &dgs1[2], &dgs2[1], &dgs2[2]) == 0) printf("invalid entries in outputFile.txt\n");
	fclose(c);

    }		
    else
    {
	printf("no file is passed. returning...");
	return;
    }
	
    MONTP_OPT2_AS_GENE prime;
    MONTP_OPT2_AS_GENE *pPtr = &prime;
    EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

    preFFTRep *f1 = EX_randomPoly(2, dgs1, p);
    preFFTRep *f2 = EX_randomPoly(2, dgs2, p);
    preFFTRep *zerop = CreateZeroPoly();

    polyBothInHostDev *F1, *F2, *Zerop;
	a = fopen("Poly1.dat","r");
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f1, i)) == 0) printf("Polynomial is missing\n");
	}
	fclose(a);

	a = fopen("Poly2.dat","r");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f2, i))== 0) printf("Polynomial is missing\n");
	}
	fclose(a);

	//for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i) printf(" %d ",DATI(f1, i));
	//printf("\n");
	//for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i) printf(" %d ",DATI(f2, i));
	//printf("\n");
	//*/
   F1 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F1->status = -1;
   F1->polyHost = f1; 

   F2 = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   F2->status = -1;
   F2->polyHost = f2; 

   if(testP(f1, f2, p) == 0)
   {
	printf( "fail to solve  as the prime %d is not a Fourier prime\n", prime);        	            
	return;
   }
   
   


   Zerop = (polyBothInHostDev *) my_malloc(sizeof(polyBothInHostDev));
   Zerop->status = -1;
   Zerop->polyHost = zerop; 

   LinkedQueue *resQueue = EX_LinkedQueue_Init();

   clock_t t0, t1;
   double time0, time1;


   // printf("\nMYCHOICE: %d\n",myChoice);

    if(argc > 4)
    {
	if(atoi(argv[4]) == 1)     method = GPU;
	if(atoi(argv[4]) == 2)     method = GPUsmart;
    }		 
    printf("\nmethodm: %d\n", method);
    t0 = clock();
    ModularSolve2(method, F1, F2, Zerop, pPtr, resQueue);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
  
    FILE *outFile = fopen("rc.txt","w");
    fprintf(outFile,"p := %d;\n",p);
    FILE_LinkedQueue_Print(resQueue, outFile);
    fclose(outFile);
   EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);

    sfixn totalSize= totalSizeRC(resQueue);
    printf("\nTotal size: %d\n",totalSize);
    sfixn *compactRC = (sfixn*)malloc(sizeof(sfixn)*totalSize);
    makecompactRC(totalSize, compactRC, resQueue);
    for(i = 0; i < totalSize; ++i) printf(" %d ",compactRC[i]);
    printf("\n\n");

    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f  || \n ", time0);


    i = system("maple -q verifyWithMaple.mpl");

    
    EX_FreeOnePoly(zerop);
   
    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);

}


void checkResultant3(int argc, char *argv[])	
{
	if(chdir("test_files/resultant3")!= 0)
	{
		printf("Missing directory : test_files/resultant3\n");	
		exit(0);
	}

	sfixn p =  919601153;
	if(argc > 2) p = atoi(argv[2]);
    //sfixn p = 257;
	sfixn N = 3, d1 = 6, d2 = 6, d3 = 6;
	sfixn dgs1[4] = {0, 6, 6, 6};
	sfixn dgs2[4] = {0, 6, 6, 6};

	if(argc > 3) dgs1[1] = atoi(argv[3]);
	if(argc > 4) dgs1[2] = atoi(argv[4]);
	if(argc > 5) dgs1[3] = atoi(argv[5]);
	if(argc > 6) dgs2[1] = atoi(argv[6]);
	if(argc > 7) dgs2[2] = atoi(argv[7]);
	if(argc > 8) dgs2[3] = atoi(argv[8]);


	preFFTRep *F1;
    	preFFTRep *F2;
	preFFTRep *result = NULL;
	MONTP_OPT2_AS_GENE prime;
	MONTP_OPT2_AS_GENE *pPtr = &prime;
	double time0, time1;
	clock_t t0, t1;

	EX_MontP_Init_OPT2_AS_GENE(pPtr, p);
	/*
	if (argc == 2) {
        	d3 = d2 = d1 = atoi(argv[1]);
        	dgs1[1] = dgs2[1] = d1;
        	dgs1[2] = dgs2[2] = d2;
        	dgs1[3] = dgs2[3] = d3;
    	} 
	*/

    F1 = EX_randomPoly(N, dgs1, p);
    F2 = EX_randomPoly(N, dgs2, p);
    
    CUDA_TAG = 1;
    t0 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t0 = clock() - t0;
    
    FILE *outFile = fopen("prime.txt","w");
    fprintf(outFile,"p := %d:\n",p);
    fprintf(outFile,"R := PolynomialRing([c,b,a],p):\n");
    fclose(outFile);

    outFile = fopen("poly1.txt","w");
    fprintf(outFile,"F1 := ");	
    FILE_printPoly(F1, outFile);
    fprintf(outFile,":");	
    fclose(outFile);

    outFile = fopen("poly2.txt","w");
    fprintf(outFile,"F2 := ");	
    FILE_printPoly(F2, outFile);
    fprintf(outFile,":");	
    fclose(outFile);

    outFile = fopen("resultant1.txt","w");
    fprintf(outFile,"resultant1 := ");	
    FILE_printPoly(result, outFile);
    fprintf(outFile,":");	
    fclose(outFile);




    //EX_Poly_Print(F1);
    //EX_Poly_Print(F2);
    //EX_Poly_Print(result);
    EX_FreeOnePoly(result);
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;

    CUDA_TAG = 0;
    t1 = clock();
    result = EX_Resultant_Multi(F1, F2, N, pPtr);
    t1 = clock() - t1;

    outFile = fopen("resultant2.txt","w");
    fprintf(outFile,"resultant2 := ");	
    FILE_printPoly(result, outFile);
    fprintf(outFile,":");	
    fclose(outFile);


    EX_FreeOnePoly(result);
    time1 = (double)t1 / (double)CLOCKS_PER_SEC;
    
   // printf("d1 = %d, d2 = %d, d3 = %d, GPU = %.3f, CPU = %.3f\n", d1, d2, d3, time0, time1);

	sfixn i = system("maple -q resultant3.mpl");


}

	

void checkReduce(int argc, char *argv[])	
{
	if(chdir("test_files/uni_reduce")!= 0)
	{
		printf("Missing directory : test_files/uni_reduce\n");	
		exit(0);
	}

	sfixn p =  919601153;
	if(argc > 2) p = atoi(argv[2]);

	preFFTRep *f1, *f2, *input, *output1, *output2;
	TriSet *ts;
	TriRevInvSet *tris;
	sfixn  N, dgbound;
  	MONTP_OPT2_AS_GENE  prime;
  	
 
  	EX_MontP_Init_OPT2_AS_GENE(&prime,  p);

  	N=4;
	if(argc > 3) N = atoi(argv[3]);
  	dgbound=4;
	if(argc > 4) dgbound = atoi(argv[4]);

	ts = EX_initRandomTriSet(N, dgbound, &prime);
	tris = EX_getRevInvTriSet(N,  ts,  &prime);

	f1 = EX_randomPoly(N, BDS(ts), p);
  	f2 = EX_randomPoly(N, BDS(ts), p);
  	output1 = EX_InitOnePoly(N, BDS(ts));
  	output2 = EX_InitOnePoly(N, BDS(ts));

  	input = EX_EX_mulPoly(N, f1, f2, &prime);

  	MultiMod_DF(N, output1, input, ts, tris, &prime);
  	MultiMod_BULL(N, output2, input, ts, tris, &prime);

	// This code works with 7 variables: var[0] is unused
	char var[8] = {'x','a','b','c','d','e','f','g'};

	FILE *outFile = fopen("prime.txt","w");
    	fprintf(outFile,"p := %d:\n",p);
	fprintf(outFile,"R := PolynomialRing([");//c,b,a
	sfixn i;
	for(i = N; i > 1; --i)
		fprintf(outFile,"%c, ",var[i]);
	fprintf(outFile,"%c ",var[1]);

	fprintf(outFile,"],p):\n");
	fprintf(outFile,"N := %d: ",N);
    	fclose(outFile);


	outFile = fopen("poly.txt","w");
	fprintf(outFile,"f := ");	
    	FILE_printPoly(input, outFile);
    	fprintf(outFile,":");	
    	fclose(outFile);

	outFile = fopen("rc.txt","w");
	for(i = 1; i <= N; ++i)
	{
		fprintf(outFile,"rc%d := ",i);
    		FILE_printPoly(ts->elems[i], outFile);	

	    	fprintf(outFile,":\n");
	}
	fprintf(outFile,"rc := [");
	for(i = 1; i < N; ++i) fprintf(outFile,"rc%d, ",i);
	fprintf(outFile,"rc%d ",N);

	fprintf(outFile," ]:");

    	fclose(outFile);

	outFile = fopen("outputDF.txt","w");
	fprintf(outFile,"outputDF := ");
    	FILE_printPoly(output1, outFile);
    	fprintf(outFile,":");
    	fclose(outFile);
	
	outFile = fopen("outputBULL.txt","w");
	fprintf(outFile,"outputBULL := ");
    	FILE_printPoly(output2, outFile);
    	fprintf(outFile,":");
    	fclose(outFile);

	i = system("maple -q reduce.mm");



    	





  	

  	EX_FreeOnePoly(f1);
  	EX_FreeOnePoly(f2);
  	EX_FreeOnePoly(input);
	EX_FreeOnePoly(output1);
	EX_FreeOnePoly(output2);
	EX_freeTriSet(ts);
	EX_freeTriRevInvSet(tris);
  
}



extern int Interrupted;

int main(int argc, char *argv[]){
	
	sfixn problemID = 0;	
	
	if (argc > 1) problemID = atoi (argv[1]);
	
	if(problemID == 0)
	{

		if( checkPlainUniMul(argc, argv)!= 0) printf("Error in Univariate plain multiplication code\n");
		else printf("Univariate plain multiplication code PASS\n");
	}		
	if(problemID == 1)
	{

		if( checkFFTUniMul(argc, argv)!= 0) printf("Error in FFT based univariate multiplication code\n");
		else printf("FFT based univariate multiplication code PASS\n");
	}
	if(problemID == 2)
	{

		if( checkTFTUniMul(argc, argv)!= 0) printf("Error in TFT based univariate multiplication code\n");
		else printf("TFT based univariate multiplication code PASS\n");
	}
	if(problemID == 3)
	{

		if( checkUniDiv(argc, argv)!= 0) printf("Error in univariate division code\n");
		else printf("univariate division code PASS\n");
	}
	if(problemID == 4)
	{
		//Extended gcd does not work for small prime (may be non fourier prime and small degree)
		if( checkExtGCD(argc, argv)!= 0) printf("Error in extended GCD code\n");
		else printf("extended GCD code PASS\n");
	}
	if(problemID == 5)
	{

		//HalfGCD code has problem for these degrees 1000 500
		if( checkHGCD(argc, argv)!= 0) printf("Error in Half GCD code\n");
		else printf("Half GCD code PASS\n");
	}
	if(problemID == 6)
	{

		
		if( checkUniquo(argc, argv)!= 0) printf("Error in Uniquo code\n");
		else printf("Uniquo code PASS\n");
	}
	
	if(problemID == 7)
	{

		BinSolveTest(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}

	if(problemID == 8)
	{

		checkResultant3(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}

	if(problemID == 9)
	{

		checkReduce(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}
	if(problemID == 10)
	{

		benchmarkUniDiv(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}

	if(problemID == 11)
	{

		benchmarkUniGCD(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}
	if(problemID == 12)
	{

		benchmarkResultant(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}
	if(problemID == 13)
	{

		benchmarkUniIsDivided(argc, argv);	
		//if( checkReduce(argc, argv)!= 0) printf("Error in Reduce code\n");
		//else printf("Reduce code PASS\n");
	}
	
	


	
	
	

    //FILE *file;
    //file = fopen ("TEST-RESULT","a");
    //UniMul();
    //fprintf(file, "Univariate plain/FFT/TFT multiplication is correct.\n");
    //fclose(file);
    //if(Interrupted) return 0;
    //fflush(stdout);
    
    //file = fopen ("TEST-RESULT","a");
    //UniDiv();
    //fprintf(file, "Univariate plain/fast Euclidean division is correct.\n");
    //fclose(file);
    
    //if(Interrupted) return 0;
    //fflush(stdout);
    
   // file = fopen ("TEST-RESULT","a");
   // ExtendedGcd();
   // fprintf(file, "Univariate Extended Euclidean algorithm is correct.\n");
    //fclose(file);
    //if(Interrupted) return 0;
   // fflush(stdout);
    
    //file = fopen ("TEST-RESULT","a");
    //Reduce( );                           //  THIS IS a NORMAL FORM COMPUTATION
    //fprintf(file, "Univariate NormalForm is correct.\n");
   // fclose(file);
    //if(Interrupted) return 0;
   // fflush(stdout);
    
    //file = fopen ("TEST-RESULT","a");
   // Lifting();   // NOT FOR THE TEST SUITE
    //fprintf(file, "Univariate Lifting is correct.\n");
    //fclose(file);
    
    //if(Interrupted) return 0;
    //fflush(stdout); 
    //poly_tst();   // NOT FOR THE TEST SUITE
    //solve2_tst(); // NOT FOR THE TEST SUITE
    //random_dense_solve2(argc, argv); // NOT FOR THE TEST SUITE
    

	//solve for random polynomial created by Anis
    //newSolve2(argc, argv);  
	//solve for polynomial from any maple sorce file created by Anis
    //newBinarySolve2(argc, argv);	
	


    //hass_solve2(argc, argv); // NOT FOR THE TEST SUITE
    
    //uniquo_tst(argc, argv);  // WE MAKE IT INTO THE TEST SUITE
    //resultant3_tst(argc, argv); // WE MAKE IT INTO THE TEST SUITE

    return 0;
}
