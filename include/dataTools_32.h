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



#include <iostream>
#include <stdlib.h>
#include <stdio.h>
//#include "bigPrimeField.h"
using namespace std;

typedef unsigned int usfixn32Data;

/**********************************************/

//takes a row-major vector, with number of rows and columns as parameters
usfixn32Data* transposeMatrix(usfixn32Data *matrix, int nRows, int nColumns)
{
	usfixn32Data *tmp = (usfixn32Data *) malloc(nRows * nColumns * sizeof(usfixn32Data));
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nColumns; j++)
			tmp[j * nRows + i] = matrix[i * nColumns + j];
	free(matrix);
	return tmp;
}

/**********************************************/
usfixn32Data* blockPermutation(usfixn32Data *matrix, int nRows, int coefficientSize, int blockSize)
{
	usfixn32Data* permutatedMatrix = (usfixn32Data *) malloc(nRows * coefficientSize * sizeof(usfixn32Data));

//	blockSize = 32;
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < coefficientSize; j++)
		{
//			cout<<matrix[i*coefficientSize+j]<<endl;
			permutatedMatrix[(j + (i / blockSize) * coefficientSize) * blockSize + i % blockSize] = matrix[i * coefficientSize + j];
		}
	return permutatedMatrix;
}
/**********************************************/

usfixn32Data* inverseBlockPermutation(usfixn32Data *matrix, int nRows, int coefficientSize, int blockSize)
{
	usfixn32Data* permutatedMatrix = (usfixn32Data *) malloc(nRows * coefficientSize * sizeof(usfixn32Data));

//	blockSize = 32;
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < coefficientSize; j++)
		{
//			cout<<matrix[i*coefficientSize+j]<<endl;
			permutatedMatrix[i * coefficientSize + j] = matrix[(j + (i / blockSize) * coefficientSize) * blockSize + i % blockSize];
		}
	return permutatedMatrix;
}
/**********************************************/

void printVector(usfixn32Data * vector, usfixn32Data vectorSize, usfixn32Data &coefficientSize, char* vectorName = "v")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
			cout << vectorName << "[" << i << "," << j << "]=" << vector[i * coefficientSize + j] << endl;
		cout << endl;
	}
}

/**********************************************/

void printVectorToFile(usfixn32Data * vector, int &vectorSize, int &coefficientSize, char* fileName = "vector", int mode = 0)
{
	FILE * writeFile;
	if (mode == 0) //write only
		writeFile = fopen(fileName, "w");
	else if (mode == 1)
	{
		writeFile = fopen(fileName, "a");
//		cout << "opening file for appending" << endl;
	}
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
			fprintf(writeFile, "%lu\n", vector[i * coefficientSize + j]);
		fprintf(writeFile, "\n");
	}
	fclose(writeFile);
}

/**********************************************/
void printPermutatedVectorToFile(usfixn32Data * vector, int &vectorSize, int &coefficientSize, int blockSize, char* fileName = "vector")
{
	FILE * writeFile = fopen(fileName, "w");
	coefficientSize = 8;
//	for (int i = 0; i < vectorSize/32; i++)
//	{
//		for (int j = 0; j < coefficientSize; j++)
//		{//			fprintf(writeFile, "%lu\n", vector[i * coefficientSize + j]);
//			fprintf(writeFile, "%lu\n", vector[(i/32) * 8 + j]);
//		}
//		fprintf(writeFile, "\n");
//	}
	for (int b = 0; b < vectorSize / blockSize; b++)
	{
		for (int i = b * coefficientSize; i < (b + 1) * coefficientSize; i++)
//		for (int i = b * blockSize; i < (b + 1) * blockSize; i++)
		{
//			cout<<"i="<<i<<endl;
			for (int j = 0; j < blockSize; j++)
			{ //			fprintf(writeFile, "%lu\n", vector[i * coefficientSize + j]);
				fprintf(writeFile, "%lu\n", vector[i * blockSize + j]);
			}
			fprintf(writeFile, "\n");
		}
	}
	fclose(writeFile);
}

/**********************************************/

void printTwoVectors(usfixn32Data *vector1, usfixn32Data *vector2, int &vectorSize, int &coefficientSize, char* vectorName1 = "v1", char* vectorName2 = "v2")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
		{
			cout << vectorName1 << "[" << i << "," << j << "]=" << vector1[i * coefficientSize + j] << "\t";
			cout << vectorName2 << "[" << i << "," << j << "]=" << vector2[i * coefficientSize + j] << endl;
		}
		cout << "=================================" << endl;
	}
}

/**********************************************/

void printThreeVectors(usfixn32Data *vector1, usfixn32Data *vector2, usfixn32Data *vector3, int &vectorSize, int &coefficientSize, char* vectorName1 = "v1", char* vectorName2 =
		"v2", char* vectorName3 = "v3")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
		{
//			cout<<vectorName1<<"["<<i<<","<<j<<"]="<<vector1[i*coefficientSize+j]<<"\t";
//			cout<<vectorName2<<"["<<i<<","<<j<<"]="<<vector2[i*coefficientSize+j]<<"\t";
//			cout<<vectorName3<<"["<<i<<","<<j<<"]="<<vector3[i*coefficientSize+j]<<endl;
			cout << vector1[i * coefficientSize + j] << "\t";
			cout << vector2[i * coefficientSize + j] << "\t";
			cout << vector3[i * coefficientSize + j] << endl;
		}
		cout << "=================================" << endl;
	}
}
/**********************************************/

usfixn32Data* readVectorFromFile(usfixn32Data *vector, usfixn32Data &vectorSize, usfixn32Data &coefficientSize, char *fileName = NULL)
{
	FILE * readFile = fopen(fileName, "r");
	for (int i = 0; i < vectorSize; i++)
		for (int j = 0; j < coefficientSize; j++)
		{
			fscanf(readFile, "%lu", &vector[i * coefficientSize + j]);
		}
	return vector;
}
/**********************************************/
