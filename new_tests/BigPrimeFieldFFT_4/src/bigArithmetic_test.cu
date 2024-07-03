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


//currently can have up to 4M elements in input vectors
#include <time.h>
#include <stdio.h>
#include <unistd.h>

/* data manipulation auxiliary functions*/
#include "dataTools.h"

#include "BigPrimeFieldFFT_4/bigArithmetic_host_P4.h"

int
host_big_prime_FFT_compute (usfixn64 * yVector, char * yDataPath,
														usfixn64 * xVector, char * xDataPath, usfixn32 e)
{
	data fd;

	float elapsedTime;

	usfixn64 *device_xVector, *device_yVector;
	usfixn64 *device_parameters;

	fd.operation = 14;
	fd.nIterations = 1;	//number of iterations to repeat
	fd.shuffle = 1;  //transpose the input matrix from Nx8 to 8xN
	fd.paddingMethod = 0; //no padding
	//K=16 for k=8
	//K=32 for k=16
	fd.K = 2 * COEFFICIENT_SIZE;
	fd.e = e;
	fd.inputVectorSize = pow (fd.K, fd.e); //1k elements

	/****************************************
	 * check if input vector is a power of K,
	 * if yes, e=log(N,K)
	 * if no, return -1
	 ***************************************/

	char projectPath[1024];

	//	int coefficientSize = 8;
	//	int strideSize = TN;
	//	int dynamicMemSize = coefficientSize * TN;
	//	int usingDynamicMem = 1;
	//	int nParameters = 8;
	//	usfixn64 BN = inputVectorSize / TN;

	fd.coefficientSize = COEFFICIENT_SIZE;
	fd.strideSize = TN;
	fd.nParameters = 16;
	fd.transposeBlockSize = 32; //achieving highest occupancy and least time for doing permutation (block by block transposition)

	char *cumodp_home;
	cumodp_home = getenv ("CUMODP_HOME");
	sprintf (projectPath, "%s/new_tests/BigPrimeFieldFFT_4", cumodp_home);

	fd.operation = 14;

	xVector = (usfixn64*) malloc (
			sizeof(usfixn64) * fd.inputVectorSize * fd.coefficientSize);
	yVector = (usfixn64*) malloc (
			sizeof(usfixn64) * fd.inputVectorSize * fd.coefficientSize);

	fd.parameters = (usfixn64*) malloc (sizeof(usfixn64) * fd.nParameters);
	memset (fd.parameters, 0x00, sizeof(usfixn64) * fd.nParameters);

	fd.permutationStride = fd.inputVectorSize;

	fd.parameters[0] = fd.operation;
	fd.parameters[1] = fd.nIterations;
	fd.parameters[2] = fd.paddingMethod;
	fd.parameters[3] = fd.dynamicMemSize;
	fd.parameters[4] = fd.strideSize;
	fd.parameters[5] = fd.permutationStride;
	fd.parameters[6] = fd.shuffle;
	fd.parameters[7] = fd.inputVectorSize;
	fd.parameters[8] = fd.coefficientSize;

	//instead, read x and y from files xData and yData
	xVector = readVectorFromFile (xVector, fd.inputVectorSize, fd.coefficientSize,
																xDataPath);
	xVector = simpleTranspose (xVector, fd.inputVectorSize, fd.coefficientSize);

	cudaEventCreate (&fd.startEvent);
	cudaEventCreate (&fd.stopEvent);

	cudaMalloc ((void **) &device_xVector,
							sizeof(usfixn64) * fd.coefficientSize * fd.inputVectorSize);
	cudaMalloc ((void **) &device_yVector,
							sizeof(usfixn64) * fd.coefficientSize * fd.inputVectorSize);

	cudaMalloc ((void **) &device_parameters, sizeof(usfixn64) * fd.nParameters);
	cudaMemcpy (device_xVector, xVector,
							sizeof(usfixn64) * fd.coefficientSize * fd.inputVectorSize,
							cudaMemcpyHostToDevice);
	cudaMemcpy (device_parameters, fd.parameters,
							sizeof(usfixn64) * fd.nParameters, cudaMemcpyHostToDevice);

	printf (
			"operation = %d, iteration = %d, shuffled = %d, paddingMethod = %d, input = %s, vectorSize = %u \n",
			fd.operation, fd.nIterations, fd.shuffle, fd.paddingMethod, 0,
			fd.inputVectorSize);

	fd.blockSize = BLOCK_SIZE;
	fd.device_xVector = device_xVector;
	fd.device_yVector = device_yVector;
	fd.device_parameters = device_parameters;

	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);

	host_init_device_pow_omega_16w (&fd, projectPath);

	host_init_fft32_shNo (&fd);

	host_init_multiplication_permutated (&fd);

	fd.paddingMethod = fd.e;
	host_fft_general (&fd);
//	host_fft32_permutated(&fd);
	/***********************************************************/
	cudaEventRecord (fd.stopEvent, 0);
	cudaDeviceSynchronize ();
	/***********************************************************/
	cudaError_t lastError = cudaGetLastError ();

	if (cudaSuccess != lastError)
	{
		printf ("Error! = %s \n", cudaGetErrorString (lastError));
		return 0;
	}
	/***********************************************************/

	cudaEventSynchronize (fd.stopEvent);
	cudaEventElapsedTime (&elapsedTime, fd.startEvent, fd.stopEvent);

	printf ("Elapsed-GPU-Time = %.3f\n", elapsedTime);
	cudaMemcpy (xVector, device_xVector,
							sizeof(usfixn64) * fd.coefficientSize * fd.inputVectorSize,
							cudaMemcpyDeviceToHost);

	cudaMemcpy (yVector, device_yVector,
							sizeof(usfixn64) * fd.coefficientSize * fd.inputVectorSize,
							cudaMemcpyDeviceToHost);

	xVector = simpleTranspose (xVector, fd.coefficientSize, fd.inputVectorSize);
//	yVector = simpleTranspose (yVector, fd.coefficientSize, fd.inputVectorSize);

	printVectorToFile (xVector, fd.inputVectorSize, fd.coefficientSize,
										 "xData_cuda");

	printVectorToFile (yVector, fd.inputVectorSize, fd.coefficientSize,
										 yDataPath);

	cudaFree (device_xVector);
	cudaFree (device_yVector);
	cudaFree (device_parameters);
	free (xVector);
	free (yVector);

}
/*************************************************************
 * Main Function
 *************************************************************/
int
main (int argc, char *argv[])
{
	//allocate pointers for input data
	usfixn64* xVector, *yVector;
	char x_path[1024], y_path[1024];

	sprintf (
			x_path,
			"/home/dave/Desktop/FFT/repository/modpn/cumodp/new_tests/BigPrimeFieldFFT_4/test/data/data_0_11_22_1m/xData");
	sprintf (y_path, "output.txt");
	usfixn32 e = 3; //FFT of size K^e
	host_big_prime_FFT_compute (yVector, y_path, xVector, x_path, e);

	return 0;
}

