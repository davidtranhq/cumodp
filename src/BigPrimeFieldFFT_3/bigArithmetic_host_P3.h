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



int i;

void
host_twiddle_revised_permutated (data*fd, short e);

/**********************************************/
int
computeDynamicMemSize (int paddingMethod, int coefficientSize)
;
/**********************************************/

void
computeGridDims (data* fd)
;

/**********************************************/
void
host_init_multiplication_permutated (data*fd)
;
/**********************************************/
void
host_init_device_pow_omega (data*fd, char* projectPath)
;
/**********************************************/
void
host_twiddle (data* fd)
;
/**********************************************/
void
checkCudaError ()
;

/**********************************************/
void
host_permutation_16_plain (data *fd)
;

/**********************************************/
void
host_permutation_16_permutated (data *fd)
;
	/**********************************************/
void
host_permutation_16 (data* fd)
;

/**********************************************/
void
host_permutation_256_plain (data *fd, int l, int m)
;
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
	;
	if (l % baseBlockSize != 0)
	;
	if (m % baseBlockSize != 0)
	;
//	for(l=0;l<3;l++)
//	;
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
	;
//	if(l%baseBlockSize!=0)
//	;
//	if(m%baseBlockSize!=0)
//	;
//	for(l=0;l<3;l++)
//	;
/**********************************************/
void
host_permutation_256 (data* fd)
;
/**********************************************/

void
host_permutation_256_general_permutated (data* fd, usfixn64 n = 8)
;
/**********************************************/

void
host_permutation_256_general (data* fd)
;
/**********************************************/
void
host_permutationTwiddle_plain (data *fd)
;

	/**********************************************/
void
host_permutationTwiddle_permutated (data *fd)
;
/**********************************************/
void
host_permutationTwiddle (data* fd)
;
/**********************************************/

void
host_addition_permutated (data* fd)
;
	/**********************************************/

void
host_addition_plain (data* fd)
;

	/**********************************************/
void
host_addition (data* fd)
;
/**********************************************/

void
host_subtraction_permutated (data* fd)
;

	/**********************************************/

void
host_subtraction_plain (data* fd)
;

	/**********************************************/

void
host_subtraction (data* fd)
;
/**********************************************/

void
host_cyclicShift_permutated (data* fd)
;

	/**********************************************/

void
host_cyclicShift_plain (data* fd)
;

	/**********************************************/

void
host_cyclicShift (data* fd)
;
/**********************************************/

void
host_fft2_permutated (data* fd)
;
	/**********************************************/
void
host_fft2_plain (data* fd)
;
/**********************************************/

void
host_fft2 (data* fd)
;

/**********************************************/
void
host_fft16_permutated (data* fd)
;
/**********************************************/

void
host_fft16_plain (data* fd)
;
/**********************************************/

void
host_fft16 (data* fd)
;
/**********************************************/
void
host_fft256_permutated (data* fd)
;
/**********************************************/

void
host_fft256_plain (data* fd)
;
/**********************************************/

void
host_fft256 (data* fd)
;
/**********************************************/
void
host_fft32_permutated (data* fd)
;
/**********************************************/

void
host_fft32_plain (data* fd)
;
/**********************************************/

void
host_fft32 (data* fd)
;
/**********************************************/
void
host_fft64_permutated (data* fd)
;
/**********************************************/

void
host_fft64_plain (data* fd)
;
/**********************************************/

void
host_fft64 (data* fd)
;
/**********************************************/
void
host_fft1k_permutated (data* fd)
;
/**********************************************/

void
host_fft1k_plain (data* fd)
;
/**********************************************/

void
host_fft1k (data* fd)
;
/**********************************************/

void
host_fft4k_permutated (data* fd)
;
/**********************************************/
void
host_fft4k_plain (data* fd)
;
/**********************************************/

void
host_fft4k (data* fd)
;

/**********************************************/

void
host_fft16k_permutated (data* fd)
;
/**********************************************/

void
host_fft16k_plain (data* fd)
;
/**********************************************/

void
host_fft16k (data* fd)
;
/**********************************************/
void
host_fft256_plain_for_64k (data*fd)
;
/**********************************************/

void
host_fft256_permutated_for_64k (data* fd)
;

/**********************************************/
void
host_fft64k_permutated (data* fd)
;
/**********************************************/

void
host_fft64k_plain (data* fd)
;
/**********************************************/

void
host_fft64k (data* fd)
;
/**********************************************/
void
host_fft_general_plain (data*fd, short e)
;

/**********************************************/

//twiddle D(K,K^;) for computing DFT_;
void
host_twiddle_revised_permutated (data*fd, short s)
;

/**********************************************/

void
host_fft16_general (data*fd)
;
/**********************************************/
void
host_fft_general_permutated (data*fd, short e = 0)
;

/**********************************************/
void
host_fft_general (data* fd)
;

/**********************************************/

//void host_multiplication_permutated(data* fd)
//;
//
///**********************************************/
void
host_multiplication_permutated (data* fd)
;

/**********************************************/

void
host_multiplication_plain (data* fd)
;
	/**********************************************/

void
host_multiplication (data* fd)
;
/**********************************************/
//void
//host_multiplication_omega_plain (data* fd)
//;
/**********************************************/
//void
//host_multiplication_omega_permutated (data* fd)
//;
//
/**********************************************/
//void
//host_multiplication_omega (data* fd)
//;
/**********************************************/
void
host_mulLong_revised (data * fd)
;
/**********************************************/
void
host_mulLong_plain (data * fd)
;
/**********************************************/
void
host_mulLong (data* fd)
;
/**********************************************/
void
host_reduction (data* fd)
;
/**********************************************/
void
initializeData (int argc, char**argv, data* fd)
;
