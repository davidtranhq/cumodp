/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/



/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */

#include "types.h"
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

void BinSolveTest(int argc, char *argv[])
{

	sfixn p =  469762049;
	sfixn dgs1[3] = {0, 60, 60};
	sfixn dgs2[3] = {0, 60, 60};

	CUDA_TAG = 1;
	computingChoince myChoice = CPU;
	sfixn whatUwant = 0, i;
  computingChoince method = myChoice;
	FILE *a, *b, *c;

  if(argc == 1)
  {
  	printf("Input directory name is missing ... returning...");
  	return;
  }

  char* poly1 = "Poly1.dat";
  char* poly1file = malloc(strlen(argv[1]) + strlen(poly1) + 2);
  // strcpy(poly1file, argv[1]);
  // strcat(poly1file, "/");
  // strcat(poly1file, poly1);

  sprintf(poly1file, "%s/%s", argv[1], poly1);

  char* poly2 = "Poly2.dat";
  char* poly2file = malloc(strlen(argv[1]) + strlen(poly2) + 2);
  // strcpy(poly2file, argv[1]);
  // strcat(poly2file, "/");
  // strcat(poly2file, poly2);
  sprintf(poly2file, "%s/%s", argv[1], poly2);

  char* degrees = "degrees.txt";
  char* degreesfile = malloc(strlen(argv[1]) + strlen(degrees) + 2);
  // strcpy(degreesfile, argv[1]);
  // strcat(degreesfile, "/");
  // strcat(degreesfile, degrees);
  sprintf(degreesfile, "%s/%s", argv[1], degrees);
  if(argc > 2)
  {
    p = atoi(argv[2]);
  }
  else
  {
    printf("Characteristic is missing ... . returning...");
    return;
  }
  c = fopen(degreesfile,"r");
  if(fscanf(c,"%d %d %d %d %d",&i, &dgs1[1], &dgs1[2], &dgs2[1], &dgs2[2]) == 0) 
    printf("invalid entries in the file degrees.txt\n");
  fclose(c);

  MONTP_OPT2_AS_GENE prime;
  MONTP_OPT2_AS_GENE *pPtr = &prime;
  EX_MontP_Init_OPT2_AS_GENE(pPtr, p);

  preFFTRep *f1 = EX_randomPoly(2, dgs1, p);
  preFFTRep *f2 = EX_randomPoly(2, dgs2, p);
  preFFTRep *zerop = CreateZeroPoly();

  polyBothInHostDev *F1, *F2, *Zerop;
  a = fopen(poly1file,"r");
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f1, i)) == 0) 
      printf("Polynomial is missing\n");
	}
	fclose(a);

	a = fopen(poly2file,"r");
	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i)
	{
		if(fscanf(a,"%d",&DATI(f2, i))== 0) 
      printf("Polynomial is missing\n");
	}
	fclose(a);

  /**************************************************/
	for(i = 0; i < ((dgs1[1] + 1)*(dgs1[2] + 1)); ++i) 
    printf(" %d ",DATI(f1, i));
	printf("\n");

	for(i = 0; i < ((dgs2[1] + 1)*(dgs2[2] + 1)); ++i) 
    printf(" %d ",DATI(f2, i));
	printf("\n");
  /**************************************************/
	
  F1 = (polyBothInHostDev *) my_malloc (sizeof(polyBothInHostDev));
  F1->status = -1;
  F1->polyHost = f1; 

  F2 = (polyBothInHostDev *) my_malloc (sizeof(polyBothInHostDev));
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

  if(argc > 3)
  {
    if(atoi(argv[3]) == 1)  
      method = GPU;
    if(atoi(argv[3]) == 2)     
      method = GPUsmart;
  }		 
  printf("\nmethod: %d \n", method);
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
  // printf("\nTotal size: %d\n",totalSize);
  sfixn *compactRC = (sfixn*)malloc(sizeof(sfixn)*totalSize);
  makecompactRC(totalSize, compactRC, resQueue);
  // for(i = 0; i < totalSize; ++i) printf(" %d ",compactRC[i]);
  // printf("\n\n");

  EX_LinkedQueue_Free(resQueue, EX_RegularChain2_Free);
  printf(" total = %.3f  sec. \n ", time0);

  EX_FreeOnePoly(zerop);
  EX_FreeOnePoly(f1);
  EX_FreeOnePoly(f2);

  free(poly1file);
  free(poly2file);
  free(degreesfile);
}

extern int Interrupted;

int main(int argc, char *argv[])
{
  BinSolveTest(argc, argv);	
	return 0;
}