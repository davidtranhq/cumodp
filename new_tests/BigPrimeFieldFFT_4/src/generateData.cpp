//============================================================================
// Name        : generateData.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include "dataTools.h"
#include "generateData.h"
using namespace std;

int main(int argc, char* argv[])
{
	/* arguments: Mode=0, seed=0, Size=1024, coefficientSize=8, range=9223372054034644992 */
	usfixn64 *vector;

	int inputVectorSize = 1024; //1k
	int coefficientSize = 16;
	usfixn64 range = 4611686087146864640;
	int mode = 0;
	int seed = 0;

	/* Mode:
	 * 0: all less than R-1 and random (default)
	 * 1: all zero except the last one (corner case)
	 * 2: all in vector equal to seed
	 */

	if (argc > 1)
		mode = atoi(argv[1]);

	if (argc > 2)
		seed = atoi(argv[2]);

	if (argc>3)
		inputVectorSize = atoi(argv[3]);
	
	if(argc>4)
		coefficientSize = atoi(argv[4]);

	if(argc>5)
		range = strtoul(argv[5], NULL, 0);

	printf("Mode=%d, Seed=%d, inputVectorSize=%d \n", mode, seed, inputVectorSize);
	int defaultVectorSize=1024;
	vector = (usfixn64*) malloc(sizeof(usfixn64) * defaultVectorSize * coefficientSize);
	// vector = (usfixn64*) malloc(sizeof(usfixn64) * inputVectorSize * coefficientSize);

	// usfixn64* vectorGenerator(usfixn64 * vector, usfixn64 &vectorSize, int &coeffientSize, usfixn64 range, int mode, int seed)
	// vector = vectorGenerator(vector, inputVectorSize, coefficientSize, range, mode, seed);
	// printVectorToFile(vector, inputVectorSize, coefficientSize, "vector", 0);

	// inputVectorSize should be multiple of defaultVectorSize (1024)

	if(mode!=4)
	for(int b=0; b<inputVectorSize/defaultVectorSize; b++)
		{
			// if ( (b+1)*defaultVectorSize > inputVectorSize )
			// 	{
			// 		vector = vectorGenerator(vector, inputVectorSize - b*defaultVectorSize, coefficientSize, range, mode, seed);
			// 		printVectorToFile(vector, inputVectorSize, coefficientSize, "vector", 1);
			// 		break;
			// 	}
			vector = vectorGenerator(vector, defaultVectorSize, coefficientSize, range, mode, seed);
			printVectorToFile(vector, defaultVectorSize, coefficientSize, "vector", 1);
		}
	// cout<<inputVectorSize<<endl;
	// void printVectorToFile(usfixn64 * vector, int &vectorSize, int &coefficientSize, char* fileName = "vector")
	
	if(mode==4)
	// for(int b=0; b<inputVectorSize/defaultVectorSize; b++)
		{
			// if ( (b+1)*defaultVectorSize > inputVectorSize )
			// 	{
			// 		vector = vectorGenerator(vector, inputVectorSize - b*defaultVectorSize, coefficientSize, range, mode, seed);
			// 		printVectorToFile(vector, inputVectorSize, coefficientSize, "vector", 1);
			// 		break;
			// 	}
			// mode=4;
			delete(vector);
			vector = (usfixn64*) malloc(sizeof(usfixn64) * inputVectorSize * coefficientSize);
			vector = vectorGenerator(vector, inputVectorSize, coefficientSize, range, mode, seed);
			printVectorToFile(vector, inputVectorSize, coefficientSize, "vector", 1);
		}	

	return 0;
}