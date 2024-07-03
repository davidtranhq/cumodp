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
#include "BigPrimeFieldFFT_3/bigArithmeticKernels_P3.h"

int i;

void
host_twiddle_revised_permutated (data*fd, short e);

/**********************************************/
int
computeDynamicMemSize (int paddingMethod, int coefficientSize)
{
	int dynamicMemSize = 0;
	//??
	return dynamicMemSize;
}
/**********************************************/

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

/**********************************************/
void
host_init_multiplication_permutated (data*fd)
{
	cudaMalloc ((void**) &fd->device_lVector,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMalloc ((void**) &fd->device_hVector,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMalloc ((void**) &fd->device_cVector,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMemset (fd->device_lVector, 0x00,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMemset (fd->device_hVector, 0x00,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
	cudaMemset (fd->device_cVector, 0x00,
							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);

	cudaMalloc ((void**) &fd->device_signVector,
							COEFFICIENT_SIZE * sizeof(usfixn32) * fd->inputVectorSize);
	cudaMemset (fd->device_signVector, 0x00, sizeof(usfixn64) * COEFFICIENT_SIZE);

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
/**********************************************/
void
host_init_device_pow_omega (data*fd, char* projectPath)
{
	//needs fix
	//permutate host_pow_omega

	char host_pow_omega_file[1024];
	sprintf (host_pow_omega_file, "%s/src/powers_of_omega_K%d", projectPath,
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
				}
		else
			for (x = 0; x < fd->coefficientSize; x++)
				for (j = 0; j < l - 1; j++)
				{
					//stride in permutated matrix = J (it is a row-major J*coefficientSize matrix)
					//stride in tmpArray = l (it is a row-major l*coefficientSize matrix
					tmpArray[j + x * l] = fd->host_powers_omega[(j + 1) * d - 1 + x * J];
					//most important part!
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

//	cudaMalloc ((void**) &fd->device_powers_omega,
//							J * fd->coefficientSize * sizeof(usfixn64));
//	cudaMemcpy (fd->device_powers_omega, fd->host_powers_omega,
//							J * fd->coefficientSize * sizeof(usfixn64), cudaMemcpyHostToDevice);
}
/**********************************************/
void
host_twiddle (data* fd)
{
//	if (fd->paddingMethod == 0)
	{
		host_twiddle_revised_permutated (fd, fd->paddingMethod);
	}
//	if (fd->paddingMethod == 1)
//	{
//		host_twiddle_permutated (fd, fd->e);
//	}
}
/**********************************************/
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

/**********************************************/
void
host_permutation_16_plain (data *fd)
{
	printf ("in plain permutation \n");

	usfixn64 stride = fd->paddingMethod;
	stride = 16;
	if (stride == 16)
	{
		cudaEventRecord (fd->startEvent);
		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
	}
	else if(stride==64)
	{
		cudaEventRecord(fd->startEvent);
		kernel_p3_permutation_64_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
	}
	else if(stride==256)
	{
//		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
	}
	else
	{
		printf("permutation of size = %d is not available! \n", fd->paddingMethod);
	}

}

/**********************************************/
void
host_permutation_16_permutated (data *fd)
{
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);

	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	kernel_p3_permutation_permutated<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//	kernel_p3_permutation_permutated<<<(fd->inputVectorSize/256), 16>>>(fd->device_xVector, fd->device_parameters);
	short l = 0;
	int baseBlockSize = 16;
//	for(l=0;l<3;l++)
	{
//		baseBlockSize <<= l;
//		baseB2 <<= l;
		cudaEventRecord (fd->startEvent);
//		printf ("baseB2=%d, baseBlockSize=%d \n", baseB2, baseBlockSize);
//		kernel_p3_permutation_16_permutated<<<(fd->inputVectorSize/baseB2),baseBlockSize,baseB2*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_parameters);
		if (fd->paddingMethod == 0)
			kernel_p3_permutation_16_permutated<<<(fd->inputVectorSize/256),16,256*sizeof(usfixn64) >>>
			(fd->device_xVector, fd->device_parameters);

		if (fd->paddingMethod == 1)
			kernel_p3_permutation_16_permutated_v1<<<64,256>>>
			(fd->device_xVector, fd->device_parameters);
		if (fd->paddingMethod == 2)
			kernel_p3_permutation_16_permutated_v2<<<64,256>>>
			(fd->device_xVector, fd->device_parameters);

//		cudaMemcpy(fd->device_xVector, fd->device_yVector, fd->inputVectorSize*COEFFICIENT_SIZE*sizeof(usfixn64), cudaMemcpyDeviceToDevice);

		}
	}
	/**********************************************/
void
host_permutation_16 (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_permutation_16_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_permutation_16_permutated (fd);
	}
}

/**********************************************/
void
host_permutation_256_plain (data *fd, int l, int m)
{
	printf ("in plain permutation \n");
//	usfixn64 stride = fd->paddingMethod;
//	stride = 256;
	cudaEventRecord (fd->startEvent);
//	printf("par[1]=%d, par[2]=%d\n", fd->parameters[1], fd->parameters[2]);
	kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters,l,m );
//		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
}
/**********************************************/

void
host_permutation_256_permutated (data *fd, int l, int m)
{
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

//	printf("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
//			fd->nBlocks, fd->blockSize, fd->permutationStride);

	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	kernel_p3_permutation_permutated<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//	kernel_p3_permutation_permutated<<<(fd->inputVectorSize/256), 16>>>(fd->device_xVector, fd->device_parameters);
//	short l = 0;
//	short m = 0;
	int baseBlockSize = 32;
//	int baseB2 = (32*32); //(baseBlockSize/256)*4096
	int shmemSize = (32 * 32);
//	m = fd->nIterations;
//	l = fd->paddingMethod;
	if (fd->inputVectorSize < l * m)
	{
		printf ("inputSize < l*m in permutation_256_permutated, l= %d, m=%d \n", l,
						m);
		return;
	}
	if (l % baseBlockSize != 0)
	{
		printf (
				"l=stride is not divisible by baseBlockSize=%d in permutation_256_permutated ! \n",
				baseBlockSize);
		return;
	}
	if (m % baseBlockSize != 0)
	{
		printf (
				"m=height is not divisible by baseBlockSize=%d in permutation_256_permutated! \n",
				baseBlockSize);
		return;
	}
//	for(l=0;l<3;l++)
//	{baseBlockSize<<=l;
//	baseB2 <<=l;
//	cudaEventRecord(fd->startEvent);
//		printf("baseB2=%d, baseBlockSize=%d \n", baseB2,baseBlockSize);
	kernel_p3_permutation_256_permutated<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, l, m);
//	kernel_p3_permutation_256_permutated<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);

//	kernel_p3_permutation_256_permutated<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//	int nBlocks=fd->inputVectorSize/(baseBlockSize*m);
//	printf("blocks = %d \n", nBlocks);
//	fd->shuffle=1;
//	usfixn64 *tmpPtr;
//	tmpPtr = fd->device_xVector;
//	fd->device_xVector = fd->device_yVector;
//	fd->device_yVector = tmpPtr;

	cudaMemcpy (fd->device_xVector, fd->device_yVector,
							COEFFICIENT_SIZE * fd->inputVectorSize * sizeof(usfixn64),
							cudaMemcpyDeviceToDevice);
	cudaMemset (fd->device_yVector, 0x00,
							COEFFICIENT_SIZE * fd->inputVectorSize * sizeof(usfixn64));

//	fd->device_xVector = fd->device_yVector;
}
/**********************************************/
void
host_permutation_256_permutated_2 (data *fd, int l, int m)
{
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);

	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	kernel_p3_permutation_permutated<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//	kernel_p3_permutation_permutated<<<(fd->inputVectorSize/256), 16>>>(fd->device_xVector, fd->device_parameters);
	int baseBlockSize = 256;
	int baseB2 = (32 * 32); //(baseBlockSize/256)*4096
	int shmemSize = (2 * 1024); //16K
//	m = fd->nIterations;
//	l = fd->paddingMethod;
	if (fd->inputVectorSize < l * m)
	{
		printf ("inputSize < l*m in permutation_256_permutated! \n");
		return;
	}
//	if(l%baseBlockSize!=0)
//	{
//			printf ("l=stride is not divisible by baseBlockSize=%d in permutation_256_permutated! \n",baseBlockSize);
//			return ;
//	}
//	if(m%baseBlockSize!=0)
//	{
//			printf ("m=height is not divisible by baseBlockSize=%d in permutation_256_permutated! \n",baseBlockSize);
//			return ;
//	}
//	for(l=0;l<3;l++)
//	{baseBlockSize<<=l;
//	baseB2 <<=l;
	cudaEventRecord (fd->startEvent);
//		printf("baseB2=%d, baseBlockSize=%d \n", baseB2,baseBlockSize);
	int nb = (fd->inputVectorSize / (l * m));
//	nb = ((nb+(baseBlockSize/l-1))/(baseBlockSize/l));
//	nb = 1;
//#########################
	/*problem in order of blocks and indexing of strides, the whole algo is correct 	*/
//##########################
//	kernel_p3_permutation_256_permutated_2<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
	kernel_p3_permutation_256_permutated_2<<<nb,baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
//	(fd->inputVectorSize/(baseBlockSize))
//	kernel_p3_permutation_256_permutated<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);

//	kernel_p3_permutation_256_permutated<<<(fd->inputVectorSize/(baseBlockSize*m)),baseBlockSize,shmemSize*sizeof(usfixn64) >>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
	int nBlocks = fd->inputVectorSize / (baseBlockSize);
	printf ("blocks = %d \n", nBlocks);
	fd->shuffle = 1;
//	fd->device_xVector = fd->device_yVector;
}
/**********************************************/
void
host_permutation_256 (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_permutation_256_plain (fd, fd->paddingMethod, fd->nIterations);
	}
	else if (fd->shuffle == 1)
	{
//		cudaEventRecord(fd->startEvent);
//		host_permutation_256_permutated(fd, fd->paddingMethod, fd->nIterations);
		host_permutation_256_permutated_2 (fd, fd->paddingMethod, fd->nIterations);
	}
}
/**********************************************/

void
host_permutation_256_general_permutated (data* fd, usfixn64 n = 8)
{

	usfixn16 prevBlockSize = fd->blockSize;
	fd->blockSize = 256;
//	n=16;
//	n = 8;
//	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1)
//			/ (fd->blockSize * fd->blockSize);
//	computeGridDims(fd);
//	usfixn16 n=1024;
//	n=256;
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1)
			/ (fd->blockSize * (1 << n));
//	printf("nBlocks = %d, n/16=%d\n",fd->nBlocks, n/16);
//	dim3 currentGrid (fd->nBlocks, n/16);
	dim3 currentGrid;
//	if((n-4)<=15)
	{
		currentGrid.x = (fd->nBlocks);
		currentGrid.y = (1 << (n - 4));
	}
//	else
//	{
//		currentGrid.x = (fd->nBlocks);
//		currentGrid.y = (1 << (15));
//	}

//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
//	kernel_p3_permutation_256_general_permutated_v0<<<currentGrid, fd->blockSize>>>(n, fd->device_xVector, fd->device_yVector, fd->device_parameters);
//	short nIt=256;
//	cudaEventRecord(fd->startEvent);
//	for(nIt;nIt>0;nIt--)
//	fd->device_xVector = fd->device_yVector;
	kernel_p3_permutation_256_general_permutated_v0<<<currentGrid, fd->blockSize>>>(n, fd->device_xVector, fd->device_yVector, fd->device_parameters);

//	usfixn64 * tmpPtr = fd->device_xVector;
//	fd->device_xVector = fd->device_yVector;
//	fd->device_yVector = tmpPtr;
	fd->blockSize = prevBlockSize;
	cudaMemcpy (fd->device_xVector, fd->device_yVector,
							COEFFICIENT_SIZE * fd->inputVectorSize * sizeof(usfixn64),
							cudaMemcpyDeviceToDevice);

//	kernel_p3_swap<<<fd->inputVectorSize/fd->blockSize, fd->blockSize>>> (fd->device_xVector, fd->device_yVector, fd->device_parameters);

//	kernel_p3_swap<<<(COEFFICIENT_SIZE*fd->inputVectorSize)/fd->blockSize, fd->blockSize>>> (fd->device_xVector, fd->device_yVector, fd->device_parameters);

}
/**********************************************/

void
host_permutation_256_general (data* fd)
{
	if (fd->shuffle == 0)
	{
//		host_permutation_256_general_plain(fd, fd->paddingMethod, fd->nIterations);
	}
	else if (fd->shuffle == 1)
	{
//		cudaEventRecord(fd->startEvent);
//		host_permutation_256_permutated(fd, fd->paddingMethod, fd->nIterations);
//			host_permutation_256_permutated(fd, fd->paddingMethod);
//		host_permutation_256_general_permutated(fd, 8); //for every 256=2^8 elements
		host_permutation_256_general_permutated (fd, fd->paddingMethod); //for every 256=2^8 elements
	}
}
/**********************************************/
void
host_permutationTwiddle_plain (data *fd)
{
	printf ("in plain permutation \n");
	cudaEventRecord (fd->startEvent);
	usfixn64 stride = fd->paddingMethod;
//	if (stride == 16)
//		kernel_p3_permutationTwiddle_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////	else if(stride==64)
////			kernel_p3_permutationTwiddle_64_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////	else if(stride==256)
////				{
////		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////				}
//		else
//		{
//			printf("permutation of size = %d is not available! \n", fd->paddingMethod);
//		}

	}

	/**********************************************/
void
host_permutationTwiddle_permutated (data *fd)
{
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);

	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	kernel_p3_permutation_permutated<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);

}
/**********************************************/
void
host_permutationTwiddle (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_permutationTwiddle_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_permutationTwiddle_permutated (fd);
	}
}
/**********************************************/

void
host_addition_permutated (data* fd)
{
//N threads for N coefficients
	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Addition:\tShuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);
	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
//	if (fd->paddingMethod == 0)
	{
		cudaEventRecord (fd->startEvent, 0);

//#pragma unroll N_ITERATIONS
		for (i = 0; i < N_ITERATIONS; i++)
//			shuffled<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
			kernel_p3_addition_permutated<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);

		}

		/*****************************/

	}
	/**********************************************/

void
host_addition_plain (data* fd)
{
	fd->dynamicMemSize = fd->coefficientSize * TN;
	fd->BN = fd->inputVectorSize / TN;
	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);
//	if (fd->paddingMethod == 0)
	{
		cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
		for (i = 0; i < N_ITERATIONS; i++)
//			kernel_p3_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
			kernel_p3_addition_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
		}

		/*****************************/
	}

	/**********************************************/
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
/**********************************************/

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
			kernel_p3_subtraction_permutated<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
		}

	}

	/**********************************************/

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
			kernel_p3_subtraction_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
//			kernel_p3_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
		}

	}

	/**********************************************/

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
/**********************************************/

void
host_cyclicShift_permutated (data* fd)
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
			kernel_p3_cyclicShift_permutated<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		}

	}

	/**********************************************/

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
			kernel_p3_cyclicShift_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
		}
	}

	/**********************************************/

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
/**********************************************/

void
host_fft2_permutated (data* fd)
{
//N/2 threads for N coefficients
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
			kernel_p3_fft2_permutated<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
		}
	}
	/**********************************************/
void
host_fft2_plain (data* fd)
{

}
/**********************************************/

void
host_fft2 (data* fd)
{
	if (fd->shuffle == 1)
	{
		host_fft2_permutated (fd);
	}
	else if (fd->shuffle == 0)
	{
		host_fft2_plain (fd);
	}
}

/**********************************************/
void
host_fft16_permutated (data* fd)
{
	usfixn64 tmpBlockSize;
	//needs to be fixed, when inputVectorSize=256, need only 128 threads
	//excessive threads will interfere

//	if (fd->inputVectorSize/2 <= fd->blockSize)
//	{
//		tmpBlockSize = fd->blockSize;
//		fd->blockSize /=2;
//	}
	fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

//	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
//					fd->nBlocks, fd->blockSize, fd->permutationStride);

//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	cudaStream_t stream0, stream1;
//	cudaStreamCreate (&stream0);
//	cudaStreamCreate(&stream1);

//	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);
//	usfixn64 shmemsize = 0;
//	if (fd->paddingMethod == 0)
	{
//		cudaDeviceSynchronize();
//		cudaEventRecord (fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
//		for (int i = 0; i < N_ITERATIONS; i++)
//			kernel_p3_base_fft16<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
///*******************************************/
//		kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//		/*****************************************/

//		{
//			kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_3_r1_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//
//			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
//			//			kernel_p3_base_fft16_3_r1_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
		kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
//		//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
		kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
//		//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//
		kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
		kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		cudaDeviceSynchronize ();
//
//			/*****************************************/
////			host_permutation_2_general_permutated(fd, 8);
////			usfixn16 n=1;
////			kernel_p3_permutation_2_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_xVector, fd->device_yVector, fd->device_parameters);
//////			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//////			host_permutation_2_reverse_general_permutated(fd, 8);
//////			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//////			cudaDeviceSynchronize();
//////			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_yVector, fd->device_parameters);
//////			host_permutation_2_reverse_general_permutated(fd, 8);
////			kernel_p3_permutation_2_reverse_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_yVector, fd->device_xVector, fd->device_parameters);
////				kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			/****************************************/
////			kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//						host_permutation_2_reverse_general_permutated(fd,0);
//		}
//			kernel_p3_fft16_all<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector);
//			kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector, const_device_parameters);
//			kernel_p3_fft16_plain_singleThread<<<1, 1>>>(fd->device_xVector, fd->device_yVector);
//		cudaDeviceSynchronize();
//			kernel_p3_strided_permutation4_multiple<<<fd->gridSize, fd->blockSize, fd->blockSize*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
	}

//	if (fd->paddingMethod == 1)
//	{
//		cudaEventRecord(fd->startEvent, 0);
////#pragma unroll N_ITERATIONS
//		for (int i = 0; i < N_ITERATIONS; i++)
//		{
////			kernel_p3_base_fft16<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector);
////			kernel_p3_base_fft16_1_shmem<<<fd->gridSize, fd->blockSize, 2*fd->blockSize *COEFFICIENT_SIZE*sizeof(usfixn64)>>>(fd->device_xVector, const_device_parameters);
////			kernel_p3_fft16_plain_singleThread<<<1, 1>>>(fd->device_xVector, fd->device_yVector);
////		cudaDeviceSynchronize();
////			kernel_p3_strided_permutation4_multiple<<<fd->gridSize, fd->blockSize, fd->blockSize*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//		}
//	}
}
/**********************************************/

void
host_fft16_plain (data* fd)
{
//every thread computes one coefficient
//N threads in total for vector of size N
	fd->BN = fd->inputVectorSize / TN;
	printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);
//
//			if (fd->paddingMethod == 0)
//			{
//				//warm-up kernel call
	cudaEventRecord (fd->startEvent, 0);
//	kernel_p3_fft16_plain<<<fd->BN, TN>>>(fd->device_xVector);
//	kernel_p3_fft16_plain_1<<<fd->BN, TN>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_fft16_plain_2<<<fd->BN, TN>>>(fd->device_xVector, fd->device_parameters);
//				FFT16_1<<<1,1>>>(fd->device_xVector, fd->device_yVector);
//				cudaEventRecord(fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
//				for (i = 0; i < N_ITERATIONS; i++)
////				FFT16_1<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
//				FFT16_1<<<1,1>>>(fd->device_xVector, fd->device_yVector);
//		}
}
/**********************************************/

void
host_fft16 (data* fd)
{
	if (fd->shuffle == 1)
	{
		host_fft16_permutated (fd);
	}
	else if (fd->shuffle == 0)
	{
		host_fft16_plain (fd);
	}
}
/**********************************************/
void
host_fft256_permutated (data* fd)
{
//	fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
//	computeGridDims(fd);
//
//	cudaEventRecord(fd->startEvent, 0);
//	kernel_p3_combined_fft256<<<fd->nBlocks, fd->blockSize,fd->blockSize*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_strided_permutation_multiple<<<fd->nBlocks, fd->blockSize, (fd->blockSize*sizeof(usfixn64))>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_strided_permutation_multiple<<<fd->nBlocks, fd->blockSize, (fd->blockSize*sizeof(usfixn64))>>>(fd->device_xVector, fd->device_parameters);

	fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);

//		cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

	cudaStream_t stream0, stream1;
	cudaStreamCreate (&stream0);
	cudaStreamCreate (&stream1);

	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);
	host_init_multiplication_permutated (fd);
//	if (fd->paddingMethod == 0)
	{
		cudaEventRecord (fd->startEvent, 0);
		//#pragma unroll N_ITERATIONS
		//		for (int i = 0; i < N_ITERATIONS; i++)
		//			kernel_p3_base_fft16<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		///*******************************************/
		//		kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		/*****************************************/

		//		fd->paddingMethod =16;
//				host_permutation_plain(fd);
		//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//				host_fft16_plain(fd);
		//		host_permutationTwiddle_plain(fd);
//				kernel_p3_permutationTwiddle_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//				host_fft16_plain(fd);
//				kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16>>>(fd->device_xVector, fd->device_parameters);

		//permutation
		//fft-16
		//twiddle
		//permutation
		//fft-16
		//permutation
//		for(short t=0;t<32;t++)
		{
//			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v2<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);

			host_p3_transposition (fd, 16, 16);
//			kernel_p3_permutation_16_permutated_v1<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			cudaMemcpy(fd->device_xVector, fd->device_yVector, fd->inputVectorSize*COEFFICIENT_SIZE*sizeof(usfixn64), cudaMemcpyDeviceToDevice);

			//FFT16-round 1
			kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

//		kernel_p3_twiddle_16_permutated<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);

//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_twiddle_16_vectorized_permutated<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_powers_omega_K[2], fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_powers_omega_K[2], fd->device_parameters);
//			usfixn64 l =16;
//			kernel_p3_twiddle_general_permutated_allSteps<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//					(fd->device_xVector, fd->device_powers_omega[fd->e], fd->K, l, fd->device_parameters);
//			printf("here before that\n");
//			host_twiddle_revised_permutated(fd, 2);
			//multiplication part of twiddle factors, called r2
			//after doing the multiplication-core, need to addition addition and subtraction, for positive and negative values
			//on the left side and right side of table of multiples; called r2_join1 and r2_join2
//			kernel_p3_twiddle_16_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVector,
//					fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);

			//			newMult15_join_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r3<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector,fd->device_lVector,
//					fd->device_hVector, fd->device_cVector,
//					fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//return;
//			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v2<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);
			host_p3_transposition (fd, 16, 16);
			host_twiddle_revised_permutated (fd, 2);

			//FFT16-round 2
			kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v1<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);
			host_p3_transposition (fd, 16, 16);
		}
	}

}
/**********************************************/

void
host_fft256_plain (data* fd)
{
//every thread computes one coefficient
//N threads in total for vector of size N
	fd->BN = fd->inputVectorSize / TN;
	printf ("Plain: BN=%d, TN=%d \n", fd->BN, TN);
//	computeGridDims(fd);
	if (fd->paddingMethod == 1)
	{
		cudaEventRecord (fd->startEvent, 0);
//		kernel_p3_fft256_plain<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_fft256_plain_2<<<fd->BN, TN>>>(fd->device_xVector, fd->device_parameters);
	}

	if (fd->paddingMethod == 0)
	{
		printf ("here in fft256 plain\n");
//		fd->paddingMethod =16;
		host_permutation_16_plain (fd);
//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		host_fft16_plain (fd);
//		host_permutationTwiddle_plain(fd);
//		kernel_p3_permutationTwiddle_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		host_fft16_plain (fd);
		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16>>>(fd->device_xVector, fd->device_parameters);
		printf ("host 256 plain v1 \n");

	}
}
/**********************************************/

void
host_fft256 (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_fft256_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft256_permutated (fd);
	}
}
/**********************************************/
void
host_fft32_permutated (data* fd)
{

}
/**********************************************/

void
host_fft32_plain (data* fd)
{

}
/**********************************************/

void
host_fft32 (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_fft32_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft32_permutated (fd);
	}

}
/**********************************************/
void
host_fft64_permutated (data* fd)
{

}
/**********************************************/

void
host_fft64_plain (data* fd)
{

}
/**********************************************/

void
host_fft64 (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_fft64_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft64_permutated (fd);
	}

}
/**********************************************/
void
host_fft1k_permutated (data* fd)
{

}
/**********************************************/

void
host_fft1k_plain (data* fd)
{

}
/**********************************************/

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
/**********************************************/

void
host_fft4k_permutated (data* fd)
{
}
/**********************************************/
void
host_fft4k_plain (data* fd)
{

}
/**********************************************/

void
host_fft4k (data* fd)
{
	if (fd->shuffle == 0)
	{
		host_fft4k_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft4k_permutated (fd);
	}

}

/**********************************************/

void
host_fft16k_permutated (data* fd)
{
}
/**********************************************/

void
host_fft16k_plain (data* fd)
{

}
/**********************************************/

void
host_fft16k (data* fd)
{

	if (fd->shuffle == 0)
	{
		host_fft16k_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft16k_permutated (fd);
	}
}
/**********************************************/
void
host_fft256_plain_for_64k (data*fd)
{
//every thread computes one coefficient
//N threads in total for vector of size N
	fd->BN = fd->inputVectorSize / TN;
//	printf("Plain: BN=%d, TN=%d \n", fd->BN, TN);
//	computeGridDims(fd);
	if (fd->paddingMethod == 0)
	{
//		cudaEventRecord(fd->startEvent, 0);
//		kernel_p3_fft256_plain<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		kernel_p3_fft256_plain_2<<<fd->BN, TN>>>(fd->device_xVector, fd->device_parameters);
	}

//	if (fd->paddingMethod == 1)
//	{
////		fd->paddingMethod =16;
//		host_permutation_16_plain(fd);
////		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//		host_fft16_plain(fd);
////		host_permutationTwiddle_plain(fd);
//		kernel_p3_permutationTwiddle_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
//		host_fft16_plain(fd);
//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
////		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16>>>(fd->device_xVector, fd->device_parameters);
//		printf("host 256 plain v1 \n");
//
//	}
}
/**********************************************/

void
host_fft256_permutated_for_64k (data* fd)
{
	fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

//	printf("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
//			fd->nBlocks, fd->blockSize, fd->permutationStride);

	//		cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

	cudaStream_t stream0, stream1;
	cudaStreamCreate (&stream0);
//	cudaStreamCreate(&stream1);

//	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);
	usfixn16 shmemsize = 0;
//	if (fd->paddingMethod == 0)
	{
//		cudaEventRecord(fd->startEvent, 0);
		//#pragma unroll N_ITERATIONS
		//		for (int i = 0; i < N_ITERATIONS; i++)
		//			kernel_p3_base_fft16<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//			kernel_p3_base_fft16_2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		///*******************************************/
		//		kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		////		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
		//		/*****************************************/

		//		fd->paddingMethod =16;
		//				host_permutation_plain(fd);
		//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//				host_fft16_plain(fd);
		//		host_permutationTwiddle_plain(fd);
		//				kernel_p3_permutationTwiddle_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//				host_fft16_plain(fd);
		//				kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
		//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16>>>(fd->device_xVector, fd->device_parameters);

		//permutation
		//fft-16
		//twiddle
		//permutation
		//fft-16
		//permutation
		//		for(short t=0;t<32;t++)
		{
			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v2<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);

//			shmemsize=0;
			//FFT16-round 1
			kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

			/*****************************************/
//			//			host_permutation_2_general_permutated(fd, 8);
//			usfixn16 n=1;
//			kernel_p3_permutation_2_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			//			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_yVector, fd->device_parameters);
//			//			host_permutation_2_reverse_general_permutated(fd, 8);
//			//			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//			//			cudaDeviceSynchronize();
//			//			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_yVector, fd->device_parameters);
//			//			host_permutation_2_reverse_general_permutated(fd, 8);
//			kernel_p3_permutation_2_reverse_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_yVector, fd->device_xVector, fd->device_parameters);
//						kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			/****************************************/

			kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

			//		kernel_p3_twiddle_16_permutated<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_twiddle_16_vectorized_permutated<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_pow_omega, fd->device_parameters);

			usfixn32 height=16;
//			kernel_p3_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_twiddle_16_general_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_pow_omega, height, fd->device_lVector,
//					fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			kernel_p3_twiddle_16_general_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, height, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//			kernel_p3_twiddle_16_general_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, height, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);

//			host_twiddle_permutated(fd, 4);
			host_twiddle_revised_permutated(fd,2);
//			kernel_p3_twiddle_16_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVector,
//					fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			newMult15_join_r1<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//			newMult15_join_r2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);

//			kernel_p3_twiddle_16_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_pow_omega, fd->device_lVector,
//					fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//			kernel_p3_twiddle_16_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);
////

//
//				kernel_p3_twiddle_16_general_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, height, fd->device_lVector,
//						fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//				kernel_p3_twiddle_16_general_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, height, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//				kernel_p3_twiddle_16_general_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, height, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);

			//		kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);

			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v2<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);

			shmemsize=0;
			//FFT16-round 2
			kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

			/*****************************************/
//			//			host_permutation_2_general_permutated(fd, 8);
////			usfixn16 n=1;
//			kernel_p3_permutation_2_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			//			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_base_fft16_3_r41_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_yVector, fd->device_parameters);
//			//			host_permutation_2_reverse_general_permutated(fd, 8);
//			//			kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
//			//			cudaDeviceSynchronize();
//			//			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
////			kernel_p3_base_fft16_3_r42_v1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_yVector, fd->device_parameters);
//			//			host_permutation_2_reverse_general_permutated(fd, 8);
//			kernel_p3_permutation_2_reverse_general_permutated_v0<<<(fd->inputVectorSize/BLOCK_SIZE),BLOCK_SIZE>>> (n,fd->device_yVector, fd->device_xVector, fd->device_parameters);
//				kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
			/****************************************/

			kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);

			//		kernel_p3_permutation_16_plain<<<1,1>>>(fd->device_xVector, fd->device_parameters);
			kernel_p3_permutation_16_permutated<<<fd->inputVectorSize/256,16,256*sizeof(usfixn64)>>>(fd->device_xVector, fd->device_parameters);
//			kernel_p3_permutation_16_permutated_v1<<<32,256>>>
//			(fd->device_xVector, fd->device_parameters);
		}
	}
}

/**********************************************/
void
host_fft64k_permutated (data* fd)
{

	fd->nBlocks = (fd->inputVectorSize / 2 + fd->blockSize - 1) / fd->blockSize;
	computeGridDims (fd);

	printf ("Shuffled: nBlocks = %d, blockSize= %d , permutationStride = %d \n",
					fd->nBlocks, fd->blockSize, fd->permutationStride);

	//		cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

//	cudaStream_t stream0, stream1;
//	cudaStreamCreate(&stream0);
//	cudaStreamCreate(&stream1);

//	usfixn64 *tmpPtr;
//	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);

//	host_init_multiplication_permutated (fd);
//	if (fd->paddingMethod == 0)
	{
		cudaEventRecord (fd->startEvent, 0);
//		host_permutation_256_permutated(fd, 256, 256);
		host_permutation_256_general_permutated (fd, 8);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//
		host_fft256_permutated_for_64k (fd);
		//twiddle
//		host_twiddle_permutated (fd, 4);
		host_twiddle_revised_permutated (fd, 3);
//		host_permutation_256_permutated(fd, 256, 256);
		host_permutation_256_general_permutated (fd, 8);

//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
		host_fft256_permutated_for_64k (fd);
//		host_permutation_256_permutated(fd, 256, 256);
		host_permutation_256_general_permutated (fd, 8);

//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
	}
}
/**********************************************/

void
host_fft64k_plain (data* fd)
{
	usfixn64 * tmpPtr;
	cudaEventRecord (fd->startEvent, 0);
	kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, 256);
	printf ("Plain fft 64k \n");

	cudaMemcpy (fd->device_xVector, fd->device_yVector,
							fd->inputVectorSize * COEFFICIENT_SIZE * sizeof(usfixn64),
							cudaMemcpyDeviceToDevice);
//	tmpPtr = fd->device_xVector;
//	fd->device_xVector = fd->device_yVector;
//	fd->device_yVector = tmpPtr;

	host_fft256_plain_for_64k (fd);
	//twiddle

	kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, 256);
	cudaMemcpy (fd->device_xVector, fd->device_yVector,
							fd->inputVectorSize * COEFFICIENT_SIZE * sizeof(usfixn64),
							cudaMemcpyDeviceToDevice);
//	tmpPtr = fd->device_xVector;
//	fd->device_xVector = fd->device_yVector;
//	fd->device_yVector = tmpPtr;

	host_fft256_plain_for_64k (fd);
	kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, 256);
	cudaMemcpy (fd->device_xVector, fd->device_yVector,
							fd->inputVectorSize * COEFFICIENT_SIZE * sizeof(usfixn64),
							cudaMemcpyDeviceToDevice);
//	tmpPtr = fd->device_xVector;
//	fd->device_xVector = fd->device_yVector;
//	fd->device_yVector = tmpPtr;

}
/**********************************************/

void
host_fft64k (data* fd)
{

	if (fd->shuffle == 0)
	{
		host_fft64k_plain (fd);
	}
	else if (fd->shuffle == 1)
	{
		host_fft64k_permutated (fd);
	}
}
/**********************************************/
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
//		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
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
//		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//
//		host_fft256_plain_for_64k(fd);
//
//		kernel_p3_permutation_256_plain<<<1,1>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters, 256, m);
//		tmpPtr = fd->device_xVector;
//		fd->device_xVector = fd->device_yVector;
//		fd->device_yVector = tmpPtr;
//	}

}

/**********************************************/

//twiddle D(K,K^{s-1}) for computing DFT_{K^s}
void
host_twiddle_revised_permutated (data*fd, short s)
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
//	kernel_p3_twiddle_general_permutated_allSteps<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l, fd->device_parameters);

//	host_init_multiplication_permutated (fd);
	kernel_p3_twiddle_general_permutated_step1<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l, fd->device_parameters);

//		kernel_p3_twiddle_general_permutated_step2<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
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

	kernel_p3_twiddle_general_permutated_step21<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
			fd->device_lVector, fd->device_hVector,
			fd->device_cVector, fd->device_signVector, fd->device_parameters);

	//2nd and third step shouldnot do anything for i0*j0=0;
	//also, check when i0*j0=0 then return

	kernel_p3_twiddle_general_permutated_step22<<<fd->gridSize , BLOCK_SIZE>>>
	(fd->device_xVector, fd->K, l,
			fd->device_lVector, fd->device_hVector,
			fd->device_cVector, fd->device_signVector, fd->device_parameters);

	kernel_p3_twiddle_general_permutated_step23<<<fd->gridSize , BLOCK_SIZE>>>
	(fd->device_xVector, fd->K, l,
			fd->device_lVector, fd->device_hVector,
			fd->device_cVector, fd->device_signVector, fd->device_parameters);

//
//		mult_revised_8lhc_step2<<<fd->gridSize,fd->blockSize >>>
//		(fd->device_parameters,
//				fd->device_lVector, fd->device_hVector, fd->device_cVector,
//				fd->device_signVector);
//
//		mult_revised_8lhc_step3<<<fd->gridSize, fd->blockSize>>>(
//				fd->device_xVector, fd->device_lVector, fd->device_hVector,
//				fd->device_cVector, fd->device_signVector, fd->device_parameters);

//	kernel_p3_twiddle_general_permutated_step22<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//			fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);
//
//	kernel_p3_twiddle_general_permutated_step23<<<fd->inputVectorSize/BLOCK_SIZE, BLOCK_SIZE>>>
//	(fd->device_xVector, fd->device_powers_omega_K[s], fd->K, l,
//			fd->device_lVector, fd->device_hVector,
//			fd->device_cVector, fd->device_signVector, fd->device_parameters);

//	kernel_p3_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, fd->device_parameters);
//	kernel_p3_twiddle_16_general_vectorized_permutated_r2_v1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, fd->device_pow_omega, height, fd->device_lVector,fd->device_hVector,
//			fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//	kernel_p3_twiddle_16_general_vectorized_permutated_r2_join1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, height, fd->device_lVector, fd->device_hVector, fd->device_cVector,fd->device_parameters);
//	kernel_p3_twiddle_16_general_vectorized_permutated_r2_join2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>
//			(fd->device_xVector, height, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub,fd->device_parameters);

//	kernel_p3_twiddle_16_vectorized_permutated_r1<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_parameters);
//	kernel_p3_twiddle_16_vectorized_permutated_r2<<<fd->inputVectorSize/BLOCK_SIZE,BLOCK_SIZE>>>(fd->device_xVector, fd->device_pow_omega, fd->device_parameters);
}

/**********************************************/

void
host_fft16_general (data*fd)
{

//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);

	cudaStream_t stream0, stream1;
	cudaStreamCreate (&stream0);
//	cudaStreamCreate(&stream1);

//	usfixn64 shmemsize = (COEFFICIENT_SIZE * sizeof(usfixn64) * fd->blockSize);
	usfixn16 shmemsize = 0;

	kernel_p3_base_fft16_3_r1<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	//		kernel_p3_base_fft16_3_r2<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r21<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r22<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	//		kernel_p3_base_fft16_3_r3<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r31<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r32<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	//		kernel_p3_base_fft16_3_r4<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r41<<<fd->gridSize, fd->blockSize,shmemsize, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r42<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
	kernel_p3_base_fft16_3_r5<<<fd->gridSize, fd->blockSize,0, stream0>>>(fd->device_xVector, fd->device_parameters);
}
/**********************************************/
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
		cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
	cudaEventRecord (fd->startEvent, 0);
	//	e = fd->e;
	//	printf ("in general FFT for K^e for K=%d and e=%d\n", fd->K, e);

	//	short m = e;
	short j = 0;
	usfixn64 stride;
	usfixn64 m;
	m = fd->inputVectorSize / (fd->K);
	/////***********************************************/
	for (i = 0; i < e - 1; i++)
	{
		//		printf ("tranposing 32, m= %d \n", m);
		host_p3_transposition (fd, fd->K, m);
		m /= fd->K;
		//		printf ("===================================\n");
	}

	host_fft16_permutated (fd);
	//	printf ("===================================\n");
	m = fd->K;
	for (i = 2; i <= e; i++)
	{
		host_twiddle_revised_permutated (fd, i);
		//		printf ("===================================\n");
		host_p3_transposition (fd, m, fd->K);
		//		printf ("===================================\n");
		host_fft16_permutated (fd);
		//		printf ("===================================\n");
		host_p3_transposition (fd, fd->K, m);
		//		printf ("===================================\n");
		m *= fd->K;
	}
}

/**********************************************/
void
host_fft_general (data* fd)
{

	if (fd->shuffle == 0)
	{
		host_fft_general_plain (fd, fd->paddingMethod);
	}
	else if (fd->shuffle == 1)
	{
		usfixn32 tmpPaddingMethod = fd->paddingMethod;
		host_fft_general_permutated (fd, fd->paddingMethod);
		fd->paddingMethod = tmpPaddingMethod;
	}
}

/**********************************************/

//void host_multiplication_permutated(data* fd)
//{
//	fd->parameters[0] = fd->inputVectorSize * fd->coefficientSize;
//	fd->parameters[1] = 0;
//	cudaMemcpy(fd->device_parameters, fd->parameters,
//			sizeof(usfixn64) * fd->nParameters, cudaMemcpyHostToDevice);
////	//	blockSize = 64;
////	//	blockSize = 128;
////	fd->blockSize = 256;
//	fd->nBlocks = (fd->coefficientSize * fd->inputVectorSize + fd->blockSize - 1)
//			/ fd->blockSize;
//
//	computeGridDims(fd);
//
//	cout << "nBlocks in shuffled Multiplication = " << fd->nBlocks << endl;
//
////			newMult<<< nBlocks, blockSize, 3*sizeof(usfixn64)*blockSize >>> (device_xVector, device_yVector, device_uVector, device_parameters);
////			newMult2<<< nBlocks, blockSize, 1*sizeof(usfixn64)*blockSize >>> (device_xVector, device_yVector, device_uVector, device_parameters);
////			for (int i=0; i<nIterations; i++)
////			cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
////			cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
////			newMult4<<< nBlocks/4, blockDims, 3*blockSize*sizeof(usfixn64)/2>>> (device_xVector, device_yVector, device_uVector, device_parameters);
////			newMult4<<< nBlocks/4, blockDims>>> (device_xVector, device_yVector, device_uVector, device_parameters);
////			newMult5<<< nBlocks/4, blockSize>>> (device_xVector, device_yVector, device_uVector, device_parameters);
//	cudaEventRecord(fd->startEvent, 0);
//#pragma unroll N_ITERATIONS
//	for (i = 0; i < fd->nIterations; i++)
//	{
////					newMult12<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//		newMult13<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
//	}
////	fd->permutationStride = fd->blockSize / fd->coefficientSize;
//}
//
///**********************************************/
void
host_multiplication_permutated (data* fd)
{
	fd->parameters[0] = fd->inputVectorSize;
	fd->parameters[1] = 0;
	fd->parameters[15] = 0;	//Multiplication step from 0 to 7, increased inside kernel by tid=0;
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

		if (fd->paddingMethod == 0)
		{
			fd->gridSize = 32;
			kernel_p3_mult_revised_8lhc_step1<<<fd->gridSize,fd->blockSize >>>
			(fd->device_xVector, fd->device_yVector,fd->device_parameters,
					fd->device_lVector, fd->device_hVector, fd->device_cVector,
					fd->device_signVector);

			kernel_p3_mult_revised_8lhc_step2<<<fd->gridSize,fd->blockSize >>>
			(fd->device_parameters,
					fd->device_lVector, fd->device_hVector, fd->device_cVector,
					fd->device_signVector);

			kernel_p3_mult_revised_8lhc_step3<<<fd->gridSize, fd->blockSize>>>(
					fd->device_xVector, fd->device_lVector, fd->device_hVector,
					fd->device_cVector, fd->device_signVector, fd->device_parameters);
		}
//		newMult15_join<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_lVector, fd->device_hVector, fd->device_cVector, fd->device_lVectorSub, fd->device_hVectorSub, fd->device_cVectorSub, fd->device_parameters);
//		kernel_p3_lhc<<<1, 1>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
	}
}

/**********************************************/

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
			kernel_p3_bigMult_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			kernel_p3_bigMult_plain_2<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);
//			kernel_p3_plain<<<fd->BN, TN>>>(fd->device_xVector, fd->device_yVector, fd->device_uVector, fd->device_parameters);
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
	/**********************************************/

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
/**********************************************/
//void
//host_multiplication_omega_plain (data* fd)
//{
//
//	kernel_p3_bigMult_plain<<<fd->inputVectorSize/fd->blockSize,fd->blockSize>>> (fd->device_xVector, fd->device_yVector, fd->device_parameters);
//	//to be filled with plain host_multiplication_pow_omega
////	fd->parameters[0] = fd->inputVectorSize;
////	fd->parameters[1] = 0;
////	fd->parameters[15] = 0;	//Multiplication step from 0 to 7, increased inside kernel by tid=0;
////	cudaMemcpy(fd->device_parameters, fd->parameters,
////			sizeof(usfixn64) * fd->nParameters, cudaMemcpyHostToDevice);
////	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
////
////	computeGridDims(fd);
////
////	usfixn64 * device_lVector, *device_hVector, *device_cVector;
////	cudaMalloc((void**) &device_lVector,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cudaMalloc((void**) &device_hVector,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cudaMalloc((void**) &device_cVector,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cudaMemset(device_lVector, 0x00,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cudaMemset(device_hVector, 0x00,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cudaMemset(device_cVector, 0x00,
////			fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
////	cout << "nBlocks in shuffled Multiplication = " << fd->nBlocks << endl;
//////			cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
//////			cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);
////	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
////	cudaEventRecord(fd->startEvent, 0);
////#pragma unroll N_ITERATIONS
////	for (i = 0; i < fd->nIterations; i++)
////	{
//////					newMult12<<<gridSize, blockSize>>>(device_xVector, device_yVector, device_uVector, device_parameters);
//////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////		newMult14_incStep<<<1,1>>>(fd->device_parameters);
////
////		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////
////		newMult14_join2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector,fd->device_parameters, device_lVector, device_hVector, device_cVector);
////	}
//}
/**********************************************/
//void
//host_multiplication_omega_permutated (data* fd)
//{
//	fd->parameters[0] = fd->inputVectorSize;
//	fd->parameters[1] = 0;
//	fd->parameters[15] = 0;	//Multiplication step from 0 to 7, increased inside kernel by tid=0;
//	cudaMemcpy (fd->device_parameters, fd->parameters,
//							sizeof(usfixn64) * fd->nParameters, cudaMemcpyHostToDevice);
//	fd->nBlocks = (fd->inputVectorSize + fd->blockSize - 1) / fd->blockSize;
//
//	computeGridDims (fd);
//
//	cudaMalloc ((void**) &fd->device_lVector,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMalloc ((void**) &fd->device_hVector,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMalloc ((void**) &fd->device_cVector,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_lVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_hVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cudaMemset (fd->device_cVector, 0x00,
//							fd->inputVectorSize * sizeof(usfixn64) * COEFFICIENT_SIZE);
//	cout << "nBlocks in shuffled Multiplication = " << fd->nBlocks << endl;
//
//	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
//	//	cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
//	cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
//
//	cudaEventRecord (fd->startEvent, 0);
////#pragma unroll N_ITERATIONS
//	for (i = 0; i < fd->nIterations; i++)
//	{
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_incStep<<<1,1>>>(fd->device_parameters);
//
//		newMult14<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector, fd->device_yVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//		newMult14_join2<<<fd->gridSize, fd->blockSize>>>(fd->device_xVector,fd->device_parameters, fd->device_lVector, fd->device_hVector, fd->device_cVector);
//	}
//}
//
/**********************************************/
//void
//host_multiplication_omega (data* fd)
//{
//	if (fd->shuffle == 0)
//	{
//		host_multiplication_omega_plain (fd);
//	}
//	else if (fd->shuffle == 1)
//	{
//		host_multiplication_omega_permutated (fd);
//	}
//}
/**********************************************/
void
host_mulLong_revised (data * fd)
{

	fd->nBlocks = fd->inputVectorSize / fd->blockSize;
	computeGridDims (fd);
	cudaEventRecord (fd->startEvent);
	kernel_p3_mulLong_revised<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);

}
/**********************************************/
void
host_mulLong_plain (data * fd)
{
	fd->nBlocks = fd->inputVectorSize / fd->blockSize;
	computeGridDims (fd);
	cudaEventRecord (fd->startEvent);
	kernel_p3_mulLong_plain<<<fd->nBlocks, fd->blockSize>>>(fd->device_xVector, fd->device_yVector, fd->device_parameters);

}
/**********************************************/
void
host_mulLong (data* fd)
{
	if (fd->shuffle == 0)
		host_mulLong_plain (fd);
	if (fd->shuffle == 1)
		host_mulLong_revised (fd);

}
/**********************************************/
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
	usfixn64 primeList[16] =
		{ 962592769, 957349889, 950009857, 943718401, 940572673, 938475521,
				935329793, 925892609, 924844033, 919601153, 918552577, 913309697,
				907018241, 899678209, 897581057, 883949569 };
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
								COEFFICIENT_SIZE * sizeof(usfixn64), cudaMemcpyHostToDevice);
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
	//	for (i = 0; i < K / 2; i++)
	//	{
	//		kernel_p3_reduction_16 <<<gridDims, blockDims>>> (fd->device_xVector,
	//																									 fd->device_lVector,
	//																									 primeList[i],
	//																									 device_r_array[i], i,
	//																									 fd->device_parameters);
	//	}
	//	for (i = K/2; i < K ; i++)
	//	{
	//		kernel_p3_reduction_16 <<<gridDims, blockDims>>> (fd->device_xVector,
	//																									 fd->device_hVector,
	//																									 primeList[i],
	//																									 device_r_array[i], i,
	//																									 fd->device_parameters);
	//	}
	cudaFree (fd->device_lVector);
	cudaFree (fd->device_hVector);

	cudaMalloc ((void**) &fd->device_lVector,
							K * fd->inputVectorSize * sizeof(usfixn64));
	kernel_p3_reduction_8<<<gridDims, blockDims>>> (fd->device_xVector,
			fd->device_lVector,
			device_primeList,
			device_r_array_all_primes,
			fd->device_parameters);

}
/**********************************************/
void
initializeData (int argc, char**argv, data* fd)
{

}
