// for 16-word prime
/*
 * bigPrimeField.h
 */

#ifndef BIG_PRIME_FIELD_H_
#define BIG_PRIME_FIELD_H_

#include <stdio.h>
using namespace std;

typedef unsigned long long int usfixn64;
typedef unsigned int usfixn32;
typedef unsigned short usfixn16;
#define SQRTR 3037000502
#define R 4611686087146864640ULL

#define TILE_DIMENSION 32

//#define R ( (1<<63)+(1<<34) (x<<29+1)<<34)

//RC R complement  RC=2^64-R
//#define RC 9223372019674906624

//ULMAX=2^64-1 //equivalent to -1
#define ULMAX 18446744073709551615ULL

#define R32bit 2147483652
#define R_MINUS_ONE 4611686087146864639ULL
#define RC 13835057986562686976ULL
#define ZERO 0

#define GRID_SIZE 64

#define R0 0
#define R1 2147483652ULL
#define RC0 0
#define RC1 2147483644ULL
const usfixn64 STRIDE = 1048576ULL;

#define COEFFICIENT_SIZE 16
//#define STEPS 8
#define N_ITERATIONS 1
//#define MAX_DIM 32768 //2^15
#define MAX_DIM 1073741824
//#define UNROLL (#pragma unroll COEFFICIENT_SIZE)

#define MAX_INPUT_SIZE 1048576//1M
struct __align__(8) uConstShortArray8_align8
{
//	unsigned long long int i0, i1, i2, i3, i4, i5, i6, i7;
	short i[8];
};

struct __align__(8) uConstShortArray16_align8
{
//	unsigned long long int i0, i1, i2, i3, i4, i5, i6, i7;
	short i[16];
};

struct __align__(8) uConstArray8_align8
{
//	unsigned long long int i0, i1, i2, i3, i4, i5, i6, i7;
	usfixn64 i[8];
};

struct __align__(8) uConstArray16_align8
{
//	unsigned long long int i0, i1, i2, i3, i4, i5, i6, i7;
	usfixn64 i[16];
};

struct __align__(8) uvector8
{
	usfixn64 i0, i1, i2, i3, i4, i5, i6, i7;
};

struct __align__(8) uvector16
{
	usfixn64 i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
};

/*************************************************************/

typedef struct
{
	usfixn64* device_xVector;
	usfixn64* device_yVector;
	usfixn64* device_uVector;
	usfixn64* device_lVector;
	usfixn64* device_hVector;
	usfixn64* device_cVector;
	usfixn64* device_lVectorSub;
	usfixn64* device_hVectorSub;
	usfixn64* device_cVectorSub;
	usfixn64* device_parameters;
	usfixn64* device_pow_omega;
	usfixn64* device_powers_omega;
	usfixn64* device_powers_omega_K[10];
	//for N = K^e and omega^N=1
	//powers_KJ keeps powers of omega up to (K^J)-1
	//device_powers_omega_K[0] (null) and
	// device_powers_omega_K[1] (etta) are always equal to zero;

	usfixn64* parameters;
	usfixn32* device_signVector;
	usfixn64* host_powers_omega;
	usfixn32 K;
	usfixn64 e;

	usfixn64 maxGridDim = ((1 << 16) - 1);

	short * device_shNo_FFT_32_r2;
	short * device_shNo_FFT_32_r3;
	short * device_shNo_FFT_32_r4;
	short * device_shNo_FFT_32_r5;
	short * device_shNo_FFT_32_permutation;

	int operation = 0;
	int nIterations = 1;	//number of iterations to repeat
	int shuffle = 0;  //transpose the input matrix from Nx8 to 8xN
	int paddingMethod = 0;
	int dynamicMemSize;
	dim3 gridSize;
	cudaEvent_t startEvent;
	cudaEvent_t stopEvent;
	usfixn64 nBlocks;
	int blockSize;
	int BN;
	int inputVectorSize;
	int coefficientSize;
	int permutationStride;
	int nParameters;
	int strideSize;
	int transposeBlockSize;
} data;

/*************************************************************/

__device__              __inline__ usfixn64
getUvector8Element (const uvector8 & __restrict__ xm, const short en)
{
	if (en == 0)
		return xm.i0;
	if (en == 1)
		return xm.i1;
	if (en == 2)
		return xm.i2;
	if (en == 3)
		return xm.i3;
	if (en == 4)
		return xm.i4;
	if (en == 5)
		return xm.i5;
	if (en == 6)
		return xm.i6;
	if (en == 7)
		return xm.i7;
}

/*************************************************************/
__device__              __inline__ usfixn64
getUvector16Element (const uvector16 & __restrict__ xm, const short en)
{
	if (en == 0)
		return xm.i0;
	if (en == 1)
		return xm.i1;
	if (en == 2)
		return xm.i2;
	if (en == 3)
		return xm.i3;
	if (en == 4)
		return xm.i4;
	if (en == 5)
		return xm.i5;
	if (en == 6)
		return xm.i6;
	if (en == 7)
		return xm.i7;

	if (en == 8)
		return xm.i8;
	if (en == 9)
		return xm.i9;
	if (en == 10)
		return xm.i10;
	if (en == 11)
		return xm.i11;
	if (en == 12)
		return xm.i12;
	if (en == 13)
		return xm.i13;
	if (en == 14)
		return xm.i14;
	if (en == 15)
		return xm.i15;
}
/*************************************************************/
__device__ __inline__ void
setUvector8Element (uvector8 & __restrict__ xm, const short en,
										const usfixn64 value)
{
	if (en == 0)
	{
		xm.i0 = value;
		return;
	}

	if (en == 1)
	{
		xm.i1 = value;
		return;
	}

	if (en == 2)
	{
		xm.i2 = value;
		return;
	}

	if (en == 3)
	{
		xm.i3 = value;
		return;
	}

	if (en == 4)
	{
		xm.i4 = value;
		return;
	}

	if (en == 5)
	{
		xm.i5 = value;
		return;
	}

	if (en == 6)
	{
		xm.i6 = value;
		return;
	}

	if (en == 7)
	{
		xm.i7 = value;
		return;
	}
}
/*************************************************************/
__device__ __inline__ void
setUvector16Element (uvector16 & __restrict__ xm, const short en,
										 const usfixn64 value)
{
	if (en == 0)
	{
		xm.i0 = value;
		return;
	}

	if (en == 1)
	{
		xm.i1 = value;
		return;
	}

	if (en == 2)
	{
		xm.i2 = value;
		return;
	}

	if (en == 3)
	{
		xm.i3 = value;
		return;
	}

	if (en == 4)
	{
		xm.i4 = value;
		return;
	}

	if (en == 5)
	{
		xm.i5 = value;
		return;
	}

	if (en == 6)
	{
		xm.i6 = value;
		return;
	}

	if (en == 7)
	{
		xm.i7 = value;
		return;
	}
	if (en == 8)
	{
		xm.i8 = value;
		return;
	}

	if (en == 9)
	{
		xm.i9 = value;
		return;
	}

	if (en == 10)
	{
		xm.i10 = value;
		return;
	}

	if (en == 11)
	{
		xm.i11 = value;
		return;
	}

	if (en == 12)
	{
		xm.i12 = value;
		return;
	}

	if (en == 13)
	{
		xm.i13 = value;
		return;
	}

	if (en == 14)
	{
		xm.i14 = value;
		return;
	}

	if (en == 15)
	{
		xm.i15 = value;
		return;
	}
}
/*************************************************************/
__device__ __inline__ void
incUvector8Element (uvector8 & __restrict__ xm, const short en,
										const short incValue)
{
	if (en == 0)
	{
		xm.i0 += incValue;
		return;
	}

	if (en == 1)
	{
		xm.i1 += incValue;
		return;
	}

	if (en == 2)
	{
		xm.i2 += incValue;
		return;
	}

	if (en == 3)
	{
		xm.i3 += incValue;
		return;
	}

	if (en == 4)
	{
		xm.i4 += incValue;
		return;
	}

	if (en == 5)
	{
		xm.i5 += incValue;
		return;
	}

	if (en == 6)
	{
		xm.i6 += incValue;
		return;
	}

	if (en == 7)
	{
		xm.i7 += incValue;
		return;
	}
}
/*************************************************************/
__device__ __inline__ void
incUvector16Element (uvector16 & __restrict__ xm, const short en,
										 const short incValue)
{
	if (en == 0)
	{
		xm.i0 += incValue;
		return;
	}

	if (en == 1)
	{
		xm.i1 += incValue;
		return;
	}

	if (en == 2)
	{
		xm.i2 += incValue;
		return;
	}

	if (en == 3)
	{
		xm.i3 += incValue;
		return;
	}

	if (en == 4)
	{
		xm.i4 += incValue;
		return;
	}

	if (en == 5)
	{
		xm.i5 += incValue;
		return;
	}

	if (en == 6)
	{
		xm.i6 += incValue;
		return;
	}

	if (en == 7)
	{
		xm.i7 += incValue;
		return;
	}
	if (en == 8)
	{
		xm.i8 += incValue;
		return;
	}

	if (en == 9)
	{
		xm.i9 += incValue;
		return;
	}

	if (en == 10)
	{
		xm.i10 += incValue;
		return;
	}

	if (en == 11)
	{
		xm.i11 += incValue;
		return;
	}

	if (en == 12)
	{
		xm.i12 += incValue;
		return;
	}

	if (en == 13)
	{
		xm.i13 += incValue;
		return;
	}

	if (en == 14)
	{
		xm.i14 += incValue;
		return;
	}

	if (en == 15)
	{
		xm.i15 += incValue;
		return;
	}
}
/*************************************************************/
__device__ __inline__ void
decUvector8Element (uvector8 & __restrict__ xm, const short en,
										const short decValue)
{
	if (en == 0)
	{
		xm.i0 -= decValue;
		return;
	}

	if (en == 1)
	{
		xm.i1 -= decValue;
		return;
	}

	if (en == 2)
	{
		xm.i2 -= decValue;
		return;
	}

	if (en == 3)
	{
		xm.i3 -= decValue;
		return;
	}

	if (en == 4)
	{
		xm.i4 -= decValue;
		return;
	}

	if (en == 5)
	{
		xm.i5 -= decValue;
		return;
	}

	if (en == 6)
	{
		xm.i6 -= decValue;
		return;
	}

	if (en == 7)
	{
		xm.i7 -= decValue;
		return;
	}
}
/*******************************************/
//N the FFT size, where each element contains eight unsigned integers.
#define N 4096
//BN block number
#define BLOCK_DIM 512
//TN thread number in a block
#define THREAD_DIM 512
//input number count

#define TN 256

//#define BLOCK_SIZE 128
#define BLOCK_SIZE 256

/*************************************************************/
__device__           __inline__ usfixn64
get_tid (usfixn64& tid)
{
	//blockSize=128
	tid = blockIdx.x;
	if (blockDim.x == 32)
		tid <<= 5;
	if (blockDim.x == 64)
		tid <<= 6;
	if (blockDim.x == 128)
		tid <<= 7;
	if (blockDim.x == 256)
		tid <<= 8;
	if (blockDim.x == 512)
		tid <<= 9;
	if (blockDim.x == 1024)
		tid <<= 10;
	tid += threadIdx.x;
}

/*******************************************
 * FROM FFT
 *******************************************/

#define INC 256
//const int sharedMemSizePadding0 = 8 * TN;
//const int sharedMemSizePadding1 = 8 * TN + 28;
//const int sharedMemSizePadding2 = 8 * TN + 32;
//const int sharedMemSizePadding3 = 8 * TN + (TN * TN) / 16;
//const int sharedMemSizePadding4 = 8 * TN + (TN * TN - 4 * TN) / 32;

__constant__ short shNo_FFT_16_r5[16] =
	{ 0, 2, -2, 0, -7, -5, -9, -7, 7, 9, 5, 7, 0, 2, -2, 0 }; //to support fft16 (4 rounds, 8 cases), round 4/4

//('twiddle,D4', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8])

__constant__ short shNo_FFT_32_r2[16] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8 }; //to support fft32 (5 rounds, 16 cases), round 2/5

//('twiddle,D8', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 4, 4, 4, 4, 12, 12, 12, 12])
__constant__ short shNo_FFT_32_r3[16] =
	{ 0, 0, 0, 0, 8, 8, 8, 8, 4, 4, 4, 4, 12, 12, 12, 12 }; //to support fft32 (5 rounds, 16 cases), round 3/5

//('twiddle,D16', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 4, 4, 12, 12, 2, 2, 10, 10, 6, 6, 14, 14])
__constant__ short shNo_FFT_32_r4[16] =
	{ 0, 0, 8, 8, 4, 4, 12, 12, 2, 2, 10, 10, 6, 6, 14, 14 }; //to support fft32 (5 rounds, 16 cases), round 4/5

//('twiddle,D32', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15])
__constant__ short shNo_FFT_32_r5[16] =
	{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 }; //to support fft32 (5 rounds, 16 cases), round 5/5

__constant__ short shNo_FFT_32_r6[32] =
	{ 0, 15, 6, 21, 0, 15, 6, 21, -6, 9, 0, 15, -6, 9, 0, 15, -15, 0, -9, 6, -15,
			0, -9, 6, -21, -6, -15, 0, -21, -6, -15, 0 };

/***********************************************
 * product used in multiplication
 **********************************************/

typedef struct
{
	usfixn64 l, h;
	short c;
} product;

/***********************************************
 * uform used in multiplication
 **********************************************/

typedef struct
{
	usfixn64 u0;
	short u1;
} urform;

/**********************************************/
/**********************************************/
#endif /* BIG_PRIME_FIELD */
