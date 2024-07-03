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


/* basic constants of big prime field for P=r^8+1*/
//#include "bigPrimeFielld.h"
/* definition of all kernels (host and device invokers) + #include "bigArithmetic_.h" headers */
#include "BigPrimeFieldFFT_4/bigArithmeticKernels_P4.h"

int i;
void
host_twiddle_revised_permutated_16w (data*fd, short s);

/********************************************************/

//compute L_{k,m} for any "input size", k, and m are multiples of 32 (TILE_DIMENSION)
//computes input/(k*m) block permutations
void
host_transposition (data* fd, usfixn64 k, usfixn64 m)
{

	//checks
	if (k % TILE_DIMENSION != 0)
	{
		printf ("k in host_transposition is not a multiple of 32! return!;\n");
		return;
	}
	else if (m % TILE_DIMENSION != 0)
	{
		printf ("m in host_transposition is not a multiple of 32! return!;\n");
		return;
	}
	else if (fd->inputVectorSize % TILE_DIMENSION != 0)
	{
		printf (
			"inputVectorSize in host_transposition is not a multiple of 32! return!;\n");
		return;
	}

	dim3 blockDims (TILE_DIMENSION, TILE_DIMENSION, 1);
	dim3 gridDims (k / TILE_DIMENSION, m / TILE_DIMENSION, 1);

	/***************************************************/
	usfixn64 blockStride = k * m;
	usfixn64 nPermutationBlocks = fd->inputVectorSize / (blockStride);
//	printf ("npermutationBlocks = %lu, blockStride= %lu \n", nPermutationBlocks,
//					blockStride);
	/***************************************************/
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
//	printf ("k/32=%lu, m/32=%lu \n", gridDims.x, gridDims.y);
	kernel_p4_transposition <<<gridDims, blockDims>>> (fd->device_xVector,
		fd->device_yVector,
		nPermutationBlocks,
		blockStride,
		fd->device_parameters);
	host_cpy_d2d (fd->device_xVector, fd->device_yVector, fd->inputVectorSize,
		fd);
//	cudaMemcpy (fd->device_xVector, fd->device_yVector,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE,
//							cudaMemcpyDeviceToDevice);
}
/*************************************************************
 * nblocks = N / blockSize for add/sub/bigMult/shift
 * nblocks = N*8 / blockSize for mult_improved
 * nblcoks = (N/2)/blockSize for fft_base and fft_2
 *************************************************************/

void
computeGridDims (data* fd)
{
	if (fd->nBlocks <= MAX_DIM)
	{
		fd->gridSize.x = fd->nBlocks;
		fd->gridSize.y = 1;
	}
	else
	{
		fd->gridSize.x = MAX_DIM;
		fd->gridSize.y = (fd->nBlocks + MAX_DIM - 1) / MAX_DIM;
	}
}

///*************************************************************/
void
host_init_multiplication_permutated (data*fd)
{
	cudaMalloc ((void**) &fd->device_lVector,
		fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMalloc ((void**) &fd->device_hVector,
		fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMalloc ((void**) &fd->device_cVector,
		fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_lVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_hVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_cVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	//	cudaMemset (fd->device_signVector, 0x00, sizeof(usfixn64) * COEFFICIENT_SIZE);

	cudaMalloc ((void**) &fd->device_signVector,
		COEFFICIENT_SIZE * sizeof(usfixn32) * fd->inputVectorSize);

//	cudaMalloc ((void **) &fd->device_lVectorSub,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMalloc ((void **) &fd->device_hVectorSub,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMalloc ((void **) &fd->device_cVectorSub,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);

//	cudaMemset (fd->device_lVectorSub, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_hVectorSub, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_cVectorSub, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);

}
/*************************************************************/
void
host_init_device_pow_omega_16w (data*fd, char* projectPath)
{
	//needs fix
	//permutate host_pow_omega

	char host_pow_omega_file[1024];
	sprintf (host_pow_omega_file, "%s/src/powers_of_omega.dev/powers_of_omega_K%d", projectPath,
		fd->e);
//	printf ("%s\n", projectPath);
//	printf ("%s\n", host_pow_omega_file);

//	int J = fd->inputVectorSize / (fd->K); //J=N/2k = N/K
	int J = pow (fd->K, fd->e - 1);
//	J = 16;

	fd->host_powers_omega = (usfixn64*) malloc (
		J * fd->coefficientSize * sizeof(usfixn64));
	usfixn64 * tmpArray = (usfixn64*) malloc (
		J * COEFFICIENT_SIZE * sizeof(usfixn64));

	readVectorFromFile (fd->host_powers_omega, J * fd->coefficientSize,
		fd->coefficientSize, host_pow_omega_file);

//	printf (
//			"\nReading powers of omega from file for DFT_N = K^e, with K=%d, e=%d \n\n",
//			fd->K, fd->e);
//	for(short i=0;i<16;i++)
//		{
//			for(short j=0;j<8;j++)
////				host_pow_omega[i*8+j]=i*8+j;
//			printf("[%d][%d]=%lu\n", i,j,host_pow_omega[i*8+j]);
//		}

	fd->host_powers_omega = blockPermutation (fd->host_powers_omega, J,
		fd->coefficientSize, J);
	/* uncomment the following to verify if reading correct data */
//	host_pow_omega = blockPermutation (host_pow_omega, fd->coefficientSize, J, fd->coefficientSize);
//	printVectorToFile (host_pow_omega, J, fd->coefficientSize, "pw_omega", 0);
//	for(short i=0;i<8;i++)
//	{
//		for(short j=0;j<16;j++)
//		printf("[%d][%d]=%lu \t", i,j,host_pow_omega[i*16+j]);
//		printf("===\n");
//	}
//	printVector(host_pow_omega, J, fd->coefficientSize, "pow_omega");
//	for(i=0;i<128;i++)
//		host_pow_omega[i]=123123;
//	memcpy (tmpArray, fd->host_powers_omega,
//					J * fd->coefficientSize * sizeof(usfixn64));
	//this part needs more work, a lot more
	usfixn64 j = 0;
	usfixn64 d = 1;
	int l = J;
	usfixn64 offset = 0;
	short x;
	char fileName[1024];

	//filling device_powers_omega_K{i} for i in [2,e];
	for (short i = fd->e; i >= 2; i--)
	{
//		printf ("copy powers of omega for K%d \n", i);
//		printf ("l=%lu, d=%lu \n", l, d);
		cudaMalloc ((void**) &fd->device_powers_omega_K[i],
			l * fd->coefficientSize * sizeof(usfixn64));

		memset (tmpArray, 0, l * COEFFICIENT_SIZE * sizeof(usfixn64));

		if (i == fd->e)
			for (x = 0; x < fd->coefficientSize; x++)
				for (j = 0; j < l; j++)
				{
					//stride in permutated matrix = J (it is a row-major J*coefficientSize matrix)
					//stride in tmpArray = l (it is a row-major l*coefficientSize matrix
					tmpArray[j + x * l] = fd->host_powers_omega[j * d + x * J];
					if (tmpArray[j + x * l] >= R)
					{
						printf ("one omega is bigger than R!\n");
					}
				}
				else
					for (x = 0; x < fd->coefficientSize; x++)
						for (j = 0; j < l - 1; j++)
						{
					//stride in permutated matrix = J (it is a row-major J*coefficientSize matrix)
					//stride in tmpArray = l (it is a row-major l*coefficientSize matrix
							tmpArray[j + x * l] = fd->host_powers_omega[(j + 1) * d - 1 + x * J];
					//most important part!
							if (tmpArray[j + x * l] >= R)
							{
								printf ("a power of omega is bigger than R!\n");
							}
						}

		//for testing purposes
//		printVector (tmpArray, l, fd->coefficientSize, "tmpJ");
//		sprintf (fileName, "tmpJ_%d", i);
//		tmpArray = blockPermutation (tmpArray, fd->coefficientSize, l, fd->coefficientSize);
//		printVectorToFile (tmpArray, l, fd->coefficientSize, fileName, 0);

//		printf ("copied to dev_pow_omega %d \n", i);
						cudaMemcpy (fd->device_powers_omega_K[i], tmpArray,
							l * fd->coefficientSize * sizeof(usfixn64),
							cudaMemcpyHostToDevice);
						l /= fd->K;
						d *= fd->K;
					}

//	printf ("end of init pow omega\n");
//	cudaMalloc ((void**) &fd->device_powers_omega,
//							J * fd->coefficientSize * sizeof(usfixn64));
//	cudaMemcpy (fd->device_powers_omega, fd->host_powers_omega,
//							J * fd->coefficientSize * sizeof(usfixn64), cudaMemcpyHostToDevice);
				}
/*************************************************************/
				void
				host_init_fft32_shNo (data*fd)
				{
	//('twiddle,D4', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8])

					short host_shNo_FFT_32_r2[16] =
		{ 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8 }; //to support fft32 (5 rounds, 16 cases), round 2/5

	//('twiddle,D8', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 4, 4, 4, 4, 12, 12, 12, 12])
		short host_shNo_FFT_32_r3[16] =
		{ 0, 0, 0, 0, 8, 8, 8, 8, 4, 4, 4, 4, 12, 12, 12, 12 }; //to support fft32 (5 rounds, 16 cases), round 3/5

	//('twiddle,D16', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 4, 4, 12, 12, 2, 2, 10, 10, 6, 6, 14, 14])
		short host_shNo_FFT_32_r4[16] =
		{ 0, 0, 8, 8, 4, 4, 12, 12, 2, 2, 10, 10, 6, 6, 14, 14 }; //to support fft32 (5 rounds, 16 cases), round 4/5

	//('twiddle,D32', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15])
		short host_shNo_FFT_32_r5[16] =
		{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 }; //to support fft32 (5 rounds, 16 cases), round 5/5

		short host_shNoListGlobal[32] =
		{ 0, 15, 6, 21, 0, 15, 6, 21, -6, 9, 0, 15, -6, 9, 0, 15, -15, 0, -9, 6,
			-15, 0, -9, 6, -21, -6, -15, 0, -21, -6, -15, 0 };

			cudaMalloc ((void**) &fd->device_shNo_FFT_32_permutation, 32 * sizeof(short));
			cudaMalloc ((void**) &fd->device_shNo_FFT_32_r2, 16 * sizeof(short));
			cudaMalloc ((void**) &fd->device_shNo_FFT_32_r3, 16 * sizeof(short));
			cudaMalloc ((void**) &fd->device_shNo_FFT_32_r4, 16 * sizeof(short));
			cudaMalloc ((void**) &fd->device_shNo_FFT_32_r5, 16 * sizeof(short));

			cudaMemcpy (fd->device_shNo_FFT_32_r2, host_shNo_FFT_32_r2,
				16 * sizeof(short), cudaMemcpyHostToDevice);
			cudaMemcpy (fd->device_shNo_FFT_32_r3, host_shNo_FFT_32_r3,
				16 * sizeof(short), cudaMemcpyHostToDevice);
			cudaMemcpy (fd->device_shNo_FFT_32_r4, host_shNo_FFT_32_r4,
				16 * sizeof(short), cudaMemcpyHostToDevice);
			cudaMemcpy (fd->device_shNo_FFT_32_r5, host_shNo_FFT_32_r5,
				16 * sizeof(short), cudaMemcpyHostToDevice);
			cudaMemcpy (fd->device_shNo_FFT_32_permutation, host_shNoListGlobal,
				32 * sizeof(short), cudaMemcpyHostToDevice);
		}

		void
		host_twiddle (data* fd)
		{
//	if (fd->paddingMethod == 0)
			{
				host_twiddle_revised_permutated_16w (fd, fd->paddingMethod);
			}
//	if (fd->paddingMethod == 1)
//	{
//		host_twiddle_permutated_16w (fd, fd->e);
//	}
		}
/*************************************************************/
		void
		checkCudaError ()
		{
			cudaError_t cudaResult;
			cudaResult = cudaGetLastError ();
			if (cudaResult != cudaSuccess)
			{
				printf ("Error! = %s \n", cudaGetErrorString (cudaResult));
			}
		}
/*************************************************************/
		void
		host_addition_permutated (data* fd)
		{
//N threads for N coefficients
			fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
			computeGridDims (fd);

			printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
				fd->nBlocks, fd->blockSize, fd->permutationStride);
			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
			if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);

//#pragma unroll N_ITERATIONS
				for (i = 0; i < N_ITERATIONS; i++)
//			shuffled<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
					kernel_p4_addition_permutated <<<fd->gridSize, fd->blockSize>>> (
						fd->device_xVector, fd->device_yVector, fd->device_parameters);
			}
		}
/*************************************************************/

		void
		host_addition_plain (data* fd)
		{
			fd->dynamicMemSize = fd->coefficientSize * TN;
			fd->BN = fd->inputVectorSize / TN;
			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
			printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);
			if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
				for (i = 0; i < N_ITERATIONS; i++)
//			kernel_p4_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
					kernel_p4_addition_plain <<<fd->BN, TN>>> (fd->device_xVector,
						fd->device_yVector,
						fd->device_parameters);
			}
		}

/*************************************************************/
		void
		host_addition (data* fd)
		{
			if (fd->shuffle == 1)
			{
				host_addition_permutated (fd);
			}
			else if (fd->shuffle == 0)
			{
				host_addition_plain (fd);
			}
		}
/*************************************************************/

		void
		host_subtraction_permutated (data* fd)
		{
			fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
			computeGridDims (fd);
			printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
				fd->nBlocks, fd->blockSize, fd->permutationStride);

			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);

			if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);

//#pragma unroll N_ITERATIONS
				for (i = 0; i < N_ITERATIONS; i++)
					kernel_p4_subtraction_permutated <<<fd->gridSize, fd->blockSize>>> (
						fd->device_xVector, fd->device_yVector, fd->device_parameters);
			}
//
//	if (fd->paddingMethod == 1)
//	{
//		cudaEventRecord (fd->startEvent, 0);
////#pragma unroll N_ITERATIONS
//		for (i = 0; i < N_ITERATIONS; i++)
//			shuffled_sharedMem <<<fd->gridSize, fd->blockSize,
//														3 * fd->dynamicMemSize * sizeof(usfixn64)>>> (
//					fd->device_xVector, fd->device_yVector, fd->device_uVector,
//					fd->device_parameters);
//	}
		}

/*************************************************************/

		void
		host_subtraction_plain (data* fd)
		{
			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
			fd->dynamicMemSize = fd->coefficientSize * TN;
			fd->BN = fd->inputVectorSize / TN;
			printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);

			if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
				for (i = 0; i < N_ITERATIONS; i++)
					kernel_p4_subtraction_plain <<<fd->BN, TN>>> (fd->device_xVector,
						fd->device_yVector,
						fd->device_parameters);
//			kernel_p4_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
			}

		}

/*************************************************************/

		void
		host_subtraction (data* fd)
		{
			if (fd->shuffle == 1)
			{
				host_subtraction_permutated (fd);
			}
			else if (fd->shuffle == 0)
			{
				host_subtraction_plain (fd);
			}
		}
/*************************************************************/

		void
		host_cyclicShift_permutated (data* fd)
		{
			fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
			computeGridDims (fd);
			printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
				fd->nBlocks, fd->blockSize, fd->permutationStride);

			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);

//	if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);

//#pragma unroll N_ITERATIONS
//		fd->gridSize=128;
				for (i = 0; i < N_ITERATIONS; i++)
					kernel_p4_cyclicShift_permutated <<<fd->gridSize, fd->blockSize>>> (
						fd->device_xVector, fd->device_parameters);
			}

//  if (fd->paddingMethod == 1)
//    {
//      cudaEventRecord (fd->startEvent, 0);
////#pragma unroll N_ITERATIONS
//      for (i = 0; i < N_ITERATIONS; i++)
//	shuffled_sharedMem <<<fd->gridSize, fd->blockSize,
//			      3 * fd->dynamicMemSize * sizeof(usfixn64)>>> (
//	    fd->device_xVector, fd->device_yVector, fd->device_uVector,
//	    fd->device_parameters);
//    }
		}

/*************************************************************/

		void
		host_cyclicShift_plain (data* fd)
		{
			cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
			fd->dynamicMemSize = fd->coefficientSize * TN;
			fd->BN = fd->inputVectorSize / TN;
			printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);

			if (fd->paddingMethod == 0)
			{
				cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
				for (i = 0; i < N_ITERATIONS; i++)
					kernel_p4_cyclicShift_plain <<<fd->BN, TN>>> (fd->device_xVector,
						fd->device_parameters);
//			kernel_p4_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
			}

		}

/*************************************************************/

		void
		host_cyclicShift (data* fd)
		{
			if (fd->shuffle == 1)
			{
				host_cyclicShift_permutated (fd);
			}
			else if (fd->shuffle == 0)
			{
				host_cyclicShift_plain (fd);
			}
		}
/*************************************************************/
		void
		host_fft32_permutated (data* fd)
		{

			fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
			computeGridDims (fd);
//	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
//					fd->nBlocks, fd->blockSize, fd->permutationStride);

//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
			{
//		cudaDeviceSynchronize();
//		cudaEventRecord (fd->startEvent, 0); //uncomment for measuring time of fft32 alone
				kernel_p4_base_fft32_r1_DFT <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters);

				kernel_p4_base_fft32_r2_shift <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r2);
				kernel_p4_base_fft32_r2_DFT <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters);
				kernel_p4_base_fft32_r3_shift <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r3);
				kernel_p4_base_fft32_r3_DFT <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters);
				kernel_p4_base_fft32_r4_shift <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r4);
				kernel_p4_base_fft32_r4_DFT <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters);
				kernel_p4_base_fft32_r5_shift <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r5);
				kernel_p4_base_fft32_r5_DFT <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters);
				kernel_p4_base_fft32_r5_permutation <<<fd->gridSize, fd->blockSize>>> (
					fd->device_xVector, fd->device_parameters,
					fd->device_shNo_FFT_32_permutation);
			}
		}
/*************************************************************/
		void
		host_fft1k_permutated (data* fd)
		{
			fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
			computeGridDims (fd);
			printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
				fd->nBlocks, fd->blockSize, fd->permutationStride);

	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

	//	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);
			usfixn64 shmemsize = 0;
			if (fd->paddingMethod == 0)
			{
		//		cudaDeviceSynchronize();
//			cudaEventRecord (fd->startEvent, 0);
//		kernel_p4_base_fft32_r1_DFT <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters);
//
//		kernel_p4_base_fft32_r2_shift <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r2);
//		kernel_p4_base_fft32_r2_DFT <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters);
//		kernel_p4_base_fft32_r3_shift <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r3);
//		kernel_p4_base_fft32_r3_DFT <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters);
//		kernel_p4_base_fft32_r4_shift <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r4);
//		kernel_p4_base_fft32_r4_DFT <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters);
//		kernel_p4_base_fft32_r5_shift <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters, fd->device_shNo_FFT_32_r5);
//		kernel_p4_base_fft32_r5_DFT <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters);
//		kernel_p4_base_fft32_r5_permutation <<<fd->gridSize, fd->blockSize>>> (
//				fd->device_xVector, fd->device_parameters,
//				fd->device_shNo_FFT_32_permutation);
				host_transposition (fd, 32, 32);
				host_fft32_permutated (fd);
				host_twiddle_revised_permutated_16w (fd, 2);
				host_transposition (fd, 32, 32);
				host_fft32_permutated (fd);
				host_transposition (fd, 32, 32);
			}
		}
/*************************************************************/

		void
		host_fft1k_plain (data* fd)
		{

		}
/*************************************************************/

		void
		host_fft1k (data* fd)
		{
			if (fd->shuffle == 0)
			{
				host_fft1k_plain (fd);
			}
			else if (fd->shuffle == 1)
			{
				host_fft1k_permutated (fd);
			}
		}
/*************************************************************/
		void
		host_fft_general_plain (data*fd, short e)
		{
//	short i = 0;
//	short m = 1;
//
//	usfixn64 * tmpPtr;
//
//	/**********************************/
//	m = (1 << (8 * e / 2));
//
//	for (i = e; i > 2; i -= 2)
//	{
//		kernel_p4_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//		m /= 256;
//	}
//
//	/**********************************/
//	host_fft256_plain_for_64k(fd);
//	//twiddle
//	/**********************************/
//
//	m = 256;
//	for (i = 2; i < e; i += 2)
//	{
//		m *= 256;
//		//twiddle(K, K+2);
//		kernel_p4_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//
//		host_fft256_plain_for_64k(fd);
//
//		kernel_p4_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//	}

		}
/*************************************************************/
//twiddle D(K,K^{s-1}) for computing DFT_{K^s}
		void
		host_twiddle_revised_permutated_16w (data*fd, short s)
		{
	/*
	 * K=2*COEFFICIENT_SIZE
	 * for N=K^e, we read K^(e-1)-1 powers of omega from file:
	 * 		powers_of_omega_K*e*
	 * stride = N/K = K^(e-1)
	 * i = tid % stride
	 * j = tid /stride
	 * h = (i*j) / stride -> xVector[tid] = multPowR (xVector[tid],h);
	 * w = (i*j) % stride -> xVector[tid] = xVector[tid]*omega^w
	 */
//	usfixn32 height = (1 << e);
			usfixn64 l = pow (fd->K, s - 1);
//	printf ("s=%d \n", s);
//	kernel_p4_twiddle_general_permutated_allSteps<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l, fd->device_parameters);

//	host_init_multiplication_permutated (fd);
//	printf ("fd->K = %lu , l = %lu \n", fd->K, l);
			kernel_p4_twiddle_general_permutated_step1 <<<fd->inputVectorSize / BLOCK_SIZE,
			BLOCK_SIZE>>> (fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
				fd->device_parameters);
//#####################################
//	return;  problem is inside the first step! all other steps are correct.
//########################################
//		kernel_p4_twiddle_general_permutated_step2<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//		(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l, fd->device_parameters);

			fd->gridSize = 32;
//	mult_revised_8lhc_step1<<<fd->gridSize,fd->blockSize >>>
//	(fd->device_xVector, fd->device_yVector,fd->device_parameters,
//			fd->device_lVector, fd->device_hVector, fd->device_cVector,
//			fd->device_signVector);
//
//	mult_revised_8lhc_step2<<<fd->gridSize,fd->blockSize >>>
//	(fd->device_parameters,
//			fd->device_lVector, fd->device_hVector, fd->device_cVector,
//			fd->device_signVector);
//
//	mult_revised_8lhc_step3<<<fd->gridSize, fd->blockSize>>>(
//			fd->device_xVector, fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);

	/*******************/
//	kernel_p4_twiddle_general_permutated_step21 <<<fd->inputVectorSize / BLOCK_SIZE,
//	BLOCK_SIZE>>> (fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//								 fd->device_lVector, fd->device_hVector, fd->device_cVector,
//								 fd->device_signVector, fd->device_parameters);
			kernel_p4_twiddle_general_permutated_step21_lessReg <<<
			fd->inputVectorSize / BLOCK_SIZE,
			BLOCK_SIZE>>> (fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
				fd->device_lVector, fd->device_hVector, fd->device_cVector,
				fd->device_signVector, fd->device_parameters);

//	kernel_p4_twiddle_general_permutated_step21 <<<64,
//		BLOCK_SIZE>>> (fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//									 fd->device_lVector, fd->device_hVector, fd->device_cVector,
//									 fd->device_signVector, fd->device_parameters);

//	short nThreads=1;
//	kernel_p4_twiddle_general_permutated_step21_4threads <<<nThreads*fd->inputVectorSize / BLOCK_SIZE,
//		BLOCK_SIZE>>> (fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//									 fd->device_lVector, fd->device_hVector, fd->device_cVector,
//									 fd->device_signVector, fd->device_parameters,nThreads);

	//2nd and third step should not do anything for i0*j0=0;
	//also, check when i0*j0=0 then return

//	kernel_p4_twiddle_general_permutated_step22 <<<fd->gridSize, BLOCK_SIZE>>> (
//			fd->device_xVector, fd->K, l, fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);
//
//	kernel_p4_twiddle_general_permutated_step23 <<<fd->gridSize, BLOCK_SIZE>>> (
//			fd->device_xVector, fd->K, l, fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);

			kernel_p4_twiddle_general_permutated_step22 <<<fd->inputVectorSize / BLOCK_SIZE,
			BLOCK_SIZE>>> (fd->device_xVector, fd->K, l, fd->device_lVector,
				fd->device_hVector, fd->device_cVector, fd->device_signVector,
				fd->device_parameters);

			kernel_p4_twiddle_general_permutated_step23 <<<fd->inputVectorSize / BLOCK_SIZE,
			BLOCK_SIZE>>> (fd->device_xVector, fd->K, l, fd->device_lVector,
				fd->device_hVector, fd->device_cVector, fd->device_signVector,
				fd->device_parameters);

//
	/*******************/

//		mult_revised_8lhc_step2<<<fd->gridSize,fd->blockSize >>>
//		(fd->device_parameters,
//				fd->device_lVector, fd->device_hVector, fd->device_cVector,
//				fd->device_signVector);
//
//		mult_revised_8lhc_step3<<<fd->gridSize, fd->blockSize>>>(
//				fd->device_xVector, fd->device_lVector, fd->device_hVector,
//				fd->device_cVector, fd->device_signVector, fd->device_parameters);
//	kernel_p4_twiddle_general_permutated_step22<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//			fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);
//
//	kernel_p4_twiddle_general_permutated_step23<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//			fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);
//	kernel_p4_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, fd->device_parameters);
//	kernel_p4_twiddle_16_general_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, fd->device_pow_omega, height, fd->device_lVector,fd->device_hVector,
//			fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//	kernel_p4_twiddle_16_general_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, height, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//	kernel_p4_twiddle_16_general_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, height, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);
//	kernel_p4_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//	kernel_p4_twiddle_16_vectorized_permutated_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_pow_omega, fd->device_parameters);
		}

/*************************************************************/
		void
		host_fft_general_permutated (data*fd, short e = 0)
		{
//	base case fft is K=16, but it doesn't go below K^2=256, all will be through 256 functions
//	short K = 6;

//	short i = 0;
//	int m = 1;
	/************************************/
	/* computing L's: permutation(2,k) for( k=6; k-=2; k>2) */
//	m = 256 * 256;
//	m = (1<<(8*K/2));
//	m = 256;
//	for (i = K; i > 2; i -= 2)
//	for (i = e; i > 2; i -= 2)
//	{
////		host_permutation_256_permutated(fd, 256, m);
//		host_permutation_256_general_permutated(fd, i-1);//permutation L(2^(i-1),256);
////		m /= 256;
//	}
//	/************************************/
//	/* compute one DFT_256 on whole input vector*/
//	host_fft256_permutated_for_64k(fd);
//	/************************************/
//	m = 1;
//	for (i = 2; i < e; i += 2)
//	{
//		m *= 256;
//		//twiddle(K, K+2);
//		host_permutation_256_permutated(fd, 256, m);
//		host_fft256_permutated_for_64k(fd);
//		host_permutation_256_permutated(fd, 256, m);
//	}
//	e=4;
			cudaEventRecord (fd->startEvent, 0);
//	e = fd->e;
//	printf ("in general FFT for K^e for K=%d and e=%d\n", fd->K, e);

//	short m = e;
			short j = 0;
			usfixn64 stride;
			usfixn64 m;
			m = fd->inputVectorSize / 32;
/////***********************************************/
			for (i = 0; i < e - 1; i++)
			{
//		printf ("tranposing 32, m= %d \n", m);
				host_transposition (fd, 32, m);
				m /= 32;
//		printf ("===================================\n");
			}

			host_fft32_permutated (fd);
//	printf ("===================================\n");
			m = 32;
			for (i = 2; i <= e; i++)
			{
				host_twiddle_revised_permutated_16w (fd, i);
//		printf ("===================================\n");
				host_transposition (fd, m, 32);
//		printf ("===================================\n");
				host_fft32_permutated (fd);
//		printf ("===================================\n");
				host_transposition (fd, 32, m);
//		printf ("===================================\n");
				m *= 32;
			}

		}

/*************************************************************/
		void
		host_fft_general (data* fd)
		{

			if (fd->shuffle == 0)
			{
				host_fft_general_plain (fd, fd->paddingMethod);
			}
			else if (fd->shuffle == 1)
			{
				host_fft_general_permutated (fd, fd->paddingMethod);
			}
		}

/*************************************************************/

		void
		host_multiplication_permutated (data* fd)
		{
			fd->parameters[0] = fd->inputVectorSize;
			fd->parameters[1] = 0;
	fd->parameters[15] = 0; //Multiplication step from 0 to 7, increased inside kernel by tid=0;
	cudaMemcpy (fd->device_parameters, fd->parameters,
		sizeof(usfixn64) * fd->nParameters, cudaMemcpyHostToDevice);
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;

	computeGridDims (fd);

	cout << "nBlocks in shuffled Multiplication = " << fd->nBlocks << endl;

	host_init_multiplication_permutated (fd);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);

	cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
//	for (i = 0; i < fd->nIterations; i++)
	{
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_join2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_join2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
		//

//		newMult15_rAdd<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult15_rSub<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//		fd->nBlocks = ((8*fd->inputVectorSize)/fd->blockSize);//needs fix
//		computeGridDims(fd);
//		cudaEventRecord(fd->startEvent, 0);
//		newMult15_rAdd_8threads<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult15_rSub_8threads<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//
//		fd->nBlocks = ((fd->inputVectorSize)/fd->blockSize);//needs fix
//		computeGridDims(fd);

//		newMult15_rAdd_1thread<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult15_rSub_1thread<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//
//		newMult15_allSteps_1thread<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//		newMult15_allSteps_1thread_2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);

//		for (short k = 0; k < 8; k++)
//			newMult15_allSteps_1thread_2_k<<<fd->gridSize, fd->blockSize>>>(k,fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);

//		fd->nBlocks= ((fd->inputVectorSize)/fd->blockSize); //needs fix
//		computeGridDims(fd);

//		newMult15_join<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			newMult15_join_r1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//			newMult15_join_r2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//		newMult15_revised<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			mult_revised_singleThread<<<1, 1>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
//		if (fd->paddingMethod == 1)
//			mult_revised_singleThread_8lhc<<<1, 1>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
//
//			if(fd->paddingMethod==2)
//			{
////			mult_revised_v0<<<fd->gridSize,fd->blockSize >>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
//				mult_revised_step1<<<fd->gridSize,fd->blockSize >>>
//				(fd->device_xVector, fd->device_yVector,fd->device_parameters,
//						fd->device_lVector, fd->device_hVector, fd->device_cVector,
//						fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//				mult_revised_step2<<<fd->gridSize,fd->blockSize >>>
//				(fd->device_xVector, fd->device_yVector,fd->device_parameters,
//						fd->device_lVector, fd->device_hVector, fd->device_cVector,
//						fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub);
//
//				mult_revised_step3<<<fd->gridSize, fd->blockSize>>>(
//						fd->device_xVector, fd->device_lVector, fd->device_hVector,
//						fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub,
//						fd->device_cVectorSub, fd->device_parameters);
//			}

//		kernel_p4_lhc <<<1, 1>>> (fd->device_xVector, fd->device_yVector,
//													 fd->device_parameters);
//		return;
//		if (fd->paddingMethod == 0)
		{
			fd->gridSize = 32;
			printf ("blocksize = %d\n", fd->blockSize);
//			for(short i=0;i<COEFFICIENT_SIZE;i++)
			kernel_p4_mult_revised_8lhc_step1 <<<fd->gridSize, fd->blockSize>>> (
				fd->device_xVector, fd->device_yVector, fd->device_parameters,
				fd->device_lVector, fd->device_hVector, fd->device_cVector,
				fd->device_signVector);

//			mult_revised_8lhc_step1_parallel <<<8*fd->inputVectorSize/fd->blockSize, fd->blockSize>>> (
//								fd->device_xVector, fd->device_yVector, fd->device_parameters,
//								fd->device_lVector, fd->device_hVector, fd->device_cVector,
//								fd->device_signVector);

			kernel_p4_mult_revised_8lhc_step2 <<<fd->gridSize, fd->blockSize>>> (
				fd->device_parameters, fd->device_lVector, fd->device_hVector,
				fd->device_cVector, fd->device_signVector);

//			kernel_p4_one_shift_right <<<fd->inputVectorSize / fd->blockSize,
//																fd->blockSize>>> (fd->device_hVector,
//																									fd->device_parameters);
//			kernel_p4_one_shift_right <<<fd->inputVectorSize / fd->blockSize,
//																fd->blockSize>>> (fd->device_cVector,
//																									fd->device_parameters);
//			kernel_p4_one_shift_right <<<fd->inputVectorSize / fd->blockSize,
//																			fd->blockSize>>> (fd->device_cVector,
//																												fd->device_parameters);
//			mult_revised_8lhc_step3 <<<fd->gridSize, fd->blockSize>>> (
//					fd->device_xVector, fd->device_lVector, fd->device_hVector,
//					fd->device_cVector, fd->device_signVector, fd->device_parameters);
//
//			mult_revised_8lhc_step4 <<<fd->gridSize, fd->blockSize>>> (
//					fd->device_xVector, fd->device_lVector, fd->device_hVector,
//					fd->device_cVector, fd->device_signVector, fd->device_parameters);

			kernel_p4_mult_revised_8lhc_step3_v0 <<<fd->gridSize, fd->blockSize>>> (
				fd->device_xVector, fd->device_lVector, fd->device_hVector,
				fd->device_cVector, fd->device_signVector, fd->device_parameters);

//		newMult15_join<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			kernel_p4_lhc <<<1, 1>>> (fd->device_xVector, fd->device_yVector,
//														 fd->device_parameters);
		}
	}
}

/*************************************************************/

void
host_multiplication_plain (data* fd)
{
	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	fd->dynamicMemSize = fd->coefficientSize * TN;
	fd->BN = fd->inputVectorSize / TN;
	printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);

	if (fd->paddingMethod == 0)
	{
		cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
		for (i = 0; i < N_ITERATIONS; i++)
			kernel_p4_bigMult_plain <<<fd->BN, TN>>> (fd->device_xVector,
				fd->device_yVector,
				fd->device_parameters);
//			kernel_p4_bigMult_plain_2<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			kernel_p4_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
	}

//	if (fd->paddingMethod == 1)
//	{
//		cudaEventRecord(fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
//		for (i = 0; i < N_ITERATIONS; i++)
//			plain_sharedMem<<<fd->BN, TN, 3 * fd->dynamicMemSize * sizeof(usfixn64)>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
//		}

//#pragma unroll N_ITERATIONS
//	for (i = 0; i < fd->nIterations; i++)
//	{
//		newMult7<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
//		//					newMult2<<<nBlocks, blockDims, 2 * coefficientSize * sizeof(usfixn64) * blockSize / reductionStride>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					cudaMemcpy(device_parameters, parameters, sizeof(short) * nParameters, cudaMemcpyHostToDevice);
//		//					for (int j = 0; j < 3; j++)
//		//						newMult2<<<nBlocks, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult2<<<nBlocks / 2, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult4<<<nBlocks / 2, blockDims, 3 * blockSize * sizeof(usfixn64) / 2>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult4<<<nBlocks / 4, blockDims>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult4<<<nBlocks / 4, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult5<<<nBlocks / 4, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult6<<<nBlocks / 2, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult5<<<nBlocks / 8, blockDims>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult2<<<nBlocks / 8, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult2<<<nBlocks, blockDims, 2 * 3 * sizeof(usfixn64) * blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult3<<<nBlocks2, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult7<<<nBlocks, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult7<<<nBlocks, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					permutationStride = 16;
//		//					newMult11<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult7<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult9<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult8<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//					newMult10<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		//
//	}
}
/*************************************************************/

void
host_multiplication (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_multiplication_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_multiplication_permutated (fd);
	}
}
/*************************************************************/

void
host_reduction (data* fd)
{
	usfixn64 K = 2 * COEFFICIENT_SIZE;
	fd->blockSize = BLOCK_SIZE;
	dim3 blockDims (fd->blockSize, 1, 1);

	usfixn64 stride = usfixn64 (pow (1.0 * fd->K, fd->paddingMethod));
	usfixn64 nIterations = (fd->inputVectorSize + stride - 1) / stride;
	if (fd->inputVectorSize < stride)
	{
		printf ("inputVectorSize is less than K^e!  return! \n");
		return;
	}
	fd->nIterations = nIterations;
	dim3 gridDims (K * fd->inputVectorSize / fd->blockSize, 1, 1);
	if (gridDims.x > fd->maxGridDim)
		gridDims.x = fd->maxGridDim;
	//read precomputed powers of r from file
	//or precompute powers of radix here

	usfixn64* device_r_array[K];
	usfixn64* host_r_array;
	usfixn64 primeList[32] =
	{ 962592769, 957349889, 950009857, 943718401, 940572673, 938475521,
		935329793, 925892609, 924844033, 919601153, 918552577, 913309697,
		907018241, 899678209, 897581057, 883949569, 880803841, 862978049,
		850395137, 833617921, 824180737, 802160641, 800063489, 818937857,
		799014913, 786432001, 770703361, 754974721, 745537537, 740294657,
		718274561, 715128833 };
		host_r_array = (usfixn64*) malloc (K * COEFFICIENT_SIZE * sizeof(usfixn64));
		short i, j;

		for (i = 0; i < K * COEFFICIENT_SIZE; i++)
		{
			host_r_array[i] = R;
		}

		for (i = 0; i < K; i++)
		{
		//read host_r_array[i] from file for r=R and p=p[i]
			cudaMalloc ((void**) &device_r_array[i],
				COEFFICIENT_SIZE * sizeof(usfixn64));
			cudaMemcpy (device_r_array[i], host_r_array,
				COEFFICIENT_SIZE * sizeof(usfixn64),
				cudaMemcpyHostToDevice);
		}

		usfixn64 * device_r_array_all_primes;
		cudaMalloc ((void**) &device_r_array_all_primes,
			K * COEFFICIENT_SIZE * sizeof(usfixn64));
		cudaMemcpy (device_r_array_all_primes, host_r_array,
			K * COEFFICIENT_SIZE * sizeof(usfixn64), cudaMemcpyHostToDevice);

		usfixn64 * device_primeList;
		cudaMalloc ((void**) &device_primeList, K * sizeof(usfixn64));
		cudaMemcpy (device_primeList, primeList, K * sizeof(usfixn64),
			cudaMemcpyHostToDevice);
		cudaEventRecord (fd->startEvent, 0);

		cudaFree (fd->device_lVector);
		cudaFree (fd->device_hVector);

		cudaMalloc ((void**) &fd->device_lVector,
			K * fd->inputVectorSize * sizeof(usfixn64));
		kernel_p4_reduction_16 <<<gridDims, blockDims>>> (fd->device_xVector,
			fd->device_lVector,
			device_primeList,
			device_r_array_all_primes,
			fd->device_parameters);
//	kernel_p4_reduction_16_v2 <<<fd->inputVectorSize/fd->blockSize, blockDims>>> (fd->device_xVector,
//																									 fd->device_lVector,
//																									 fd->device_yVector,
//																									 device_primeList,
//																									 device_r_array_all_primes,
//																									 fd->device_parameters);

	}
/*************************************************************/
	void
	initializeData (int argc, char**argv, data* fd)
	{

	}
