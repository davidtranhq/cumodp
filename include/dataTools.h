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

typedef unsigned long long int usfixn64Data;
//typedef unsigned int usfixn64Data;

/**********************************************/

//takes a row-major vector, with number of rows and columns as parameters
usfixn64Data*
transposeMatrix (usfixn64Data *matrix, int nRows, int nColumns)
{
	usfixn64Data *tmp = (usfixn64Data *) malloc (
			nRows * nColumns * sizeof(usfixn64Data));
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nColumns; j++)
			tmp[j * nRows + i] = matrix[i * nColumns + j];
	free (matrix);
	return tmp;
}

/**********************************************/
usfixn64Data*
blockPermutation (usfixn64Data *matrix, int nRows, int coefficientSize,
									int blockSize)
{
	usfixn64Data* permutatedMatrix = (usfixn64Data *) malloc (
			nRows * coefficientSize * sizeof(usfixn64Data));

	int i, j;
//	blockSize = 32;
	for (i = 0; i < nRows; i++)
		for (j = 0; j < coefficientSize; j++)
		{
			// cout<<matrix[i*coefficientSize+j]<<endl;
			permutatedMatrix[(j + (i / blockSize) * coefficientSize) * blockSize
					+ i % blockSize] = matrix[i * coefficientSize + j];
		}
	return permutatedMatrix;
}
/**********************************************/
usfixn64Data*
simpleTranspose (usfixn64Data *matrix, int nRows, int nCols)
{
	usfixn64Data* permutatedMatrix = (usfixn64Data *) malloc (
			nRows * nCols * sizeof(usfixn64Data));

	int i, j, c;
//	blockSize = 32;
	usfixn64Data offset = 0;
	for (i = 0; i < nRows; i++)
		for (j = 0; j < nCols; j++)
		{
			// cout<<matrix[i*coefficientSize+j]<<endl;
			permutatedMatrix[i + j * nRows] = matrix[i * nCols + j];
		}
	return permutatedMatrix;
}
/**********************************************/

usfixn64Data*
inverseBlockPermutation (usfixn64Data *matrix, int nRows, int coefficientSize,
												 int blockSize)
{
	usfixn64Data* permutatedMatrix = (usfixn64Data *) malloc (
			nRows * coefficientSize * sizeof(usfixn64Data));

//	blockSize = 32;
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < coefficientSize; j++)
		{
//			cout<<matrix[i*coefficientSize+j]<<endl;
			permutatedMatrix[i * coefficientSize + j] = matrix[(j
					+ (i / blockSize) * coefficientSize) * blockSize + i % blockSize];
		}
	return permutatedMatrix;
}
/**********************************************/

void
printVector (usfixn64Data * vector, int vectorSize, int &coefficientSize,
						 char* vectorName = "v")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
			cout << vectorName << "[" << i << "," << j << "]="
					<< vector[i * coefficientSize + j] << endl;
		cout << endl;
	}
}

/**********************************************/

void
printVectorToFile (usfixn64Data * vector, int vectorSize, int coefficientSize,
									 char* fileName = "vector", int mode = 0)
{
	FILE * writeFile;

	if (mode == 0) //write only
		writeFile = fopen (fileName, "w");
	else if (mode == 1)
	{
		writeFile = fopen (fileName, "a");
//		cout << "opening file for appending" << endl;
	}

	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
			fprintf (writeFile, "%lu\n", vector[i * coefficientSize + j]);
		fprintf (writeFile, "\n");
	}
	fclose (writeFile);
}

/**********************************************/
void
printPermutatedVectorToFile (usfixn64Data * vector, int &vectorSize,
														 int &coefficientSize, int blockSize,
														 char* fileName = "vector")
{
	FILE * writeFile = fopen (fileName, "w");
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
				fprintf (writeFile, "%lu\n", vector[i * blockSize + j]);
			}
			fprintf (writeFile, "\n");
		}
	}
	fclose (writeFile);
}

/**********************************************/

void
printTwoVectors (usfixn64Data *vector1, usfixn64Data *vector2, int &vectorSize,
								 int &coefficientSize, char* vectorName1 = "v1",
								 char* vectorName2 = "v2")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
		{
			cout << vectorName1 << "[" << i << "," << j << "]="
					<< vector1[i * coefficientSize + j] << "\t";
			cout << vectorName2 << "[" << i << "," << j << "]="
					<< vector2[i * coefficientSize + j] << endl;
		}
		cout << "=================================" << endl;
	}
}

/**********************************************/

void
printThreeVectors (usfixn64Data *vector1, usfixn64Data *vector2,
									 usfixn64Data *vector3, int &vectorSize, int &coefficientSize,
									 char* vectorName1 = "v1", char* vectorName2 = "v2",
									 char* vectorName3 = "v3")
{
	for (int i = 0; i < vectorSize; i++)
	{
		for (int j = 0; j < coefficientSize; j++)
		{
			cout << vector1[i * coefficientSize + j] << "\t";
			cout << vector2[i * coefficientSize + j] << "\t";
			cout << vector3[i * coefficientSize + j] << endl;
		}
		cout << "=================================" << endl;
	}
}
/**********************************************/

usfixn64Data*
readVectorFromFile (usfixn64Data *vector, int vectorSize, int &coefficientSize,
										char *fileName = NULL)
{
	FILE * readFile = fopen (fileName, "r");
	for (int i = 0; i < vectorSize; i++)
		for (int j = 0; j < coefficientSize; j++)
		{
			fscanf (readFile, "%lu", &vector[i * coefficientSize + j]);
		}
	return vector;
}
/**********************************************/
