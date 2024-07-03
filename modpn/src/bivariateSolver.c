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

extern sfixn FP[];
extern sfixn CUDA_TAG;

void newBinarySolve2(int argc, char** argv)
{
    sfixn p =  469762049;
    sfixn dgs1[3] = {0, 60, 60};
    sfixn dgs2[3] = {0, 60, 60};

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
	fscanf(c,"%d",&i);
	fclose(c);
	if(i == -1)
	{
		printf("check the polynomial system. It might not be a bivariate system..\n");
		return;
	}
	c = fopen("outputFile.txt","r");
	fscanf(c,"%d %d %d %d %d",&i, &dgs1[1], &dgs1[2], &dgs2[1], &dgs2[2]);
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
		fscanf(a,"%d",&DATI(f1, i));
	fclose(a);

	a = fopen("Poly2.dat","r");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
		fscanf(a,"%d",&DATI(f2, i));
	fclose(a);

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

    if(argc > 3) method = atoi(argv[3]);
	
    t0 = clock();
    ModularSolve2(method, F1, F2, Zerop, pPtr, resQueue);
    t0 = clock() - t0;
    time0 = (double)t0 / (double)CLOCKS_PER_SEC;
  
    FILE *outFile = fopen("rc.txt","w");
    fprintf(outFile,"p := %d;\n",p);
    FILE_LinkedQueue_Print(resQueue, outFile);
    fclose(outFile);
   EX_LinkedQueue_Print(resQueue, EX_RegularChain2_Print);
    EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
    printf(" total = %.3f  || \n ", time0);


    i = system("maple -q verifyWithMaple.mpl");

    EX_FreeOnePoly(f1);
    EX_FreeOnePoly(f2);
    EX_FreeOnePoly(zerop);





}



int main(int argc, char *argv[]){


    newBinarySolve2(argc, argv);	
	

    return 0;
}
