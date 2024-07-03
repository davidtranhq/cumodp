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


/**
 * This file is to calculate determinant in finite field by
 * using the Gauss elemination method (no-pivot fraction-free).
 * This program just call one kernel function, while the determinant
 * calculation need call serveral times kernel function.
 * This program now only support to calculate the determinant of
 * matrix whose size is not greater than 1024, because we only use
 * one block. We can continuously to improve to support more bigger
 * matrix calculation.
 */

#include <time.h>
#include <stdio.h>

#define BASE_1 31

/************************************************/
__device__ __host__ __inline__
void
egcd (int x, int y, int *ao, int *bo, int *vo)
{
	int t, A, B, C, D, u, v, q;

	u = y;
	v = x;
	A = 1;
	B = 0;
	C = 0;
	D = 1;

	do
	{
		q = u / v;
		t = u;
		u = v;
		v = t - q * v;
		t = A;
		A = B;
		B = t - q * B;
		t = C;
		C = D;
		D = t - q * D;
	}
	while (v != 0);

	*ao = A;
	*bo = C;
	*vo = u;
}

/************************************************/
/*! \fn void egcdCPUFF(int x, int y, int *ao, int *bo, int *vo, int P)
 \brief This function computes extended gcd of integer x and y over finite field.

 \param x an integer.
 \param y an integer.
 \param ao integer pointer.
 \param bo integer pointer.
 \param vo integer pointer.
 \param P a prime.
 */

/************************************************/
__device__ __host__ __inline__
int
inv_mod (int n, int p)
{
	int a, b, v;
	egcd (n, p, &a, &b, &v);
	if (b < 0)
		b += p;
	return b % p;
}

/************************************************/
void
egcdCPUFF (int x, int y, int *ao, int *bo, int *vo, int P)
{
	int t, A, B, C, D, u, v, q;

	u = y;
	v = x;
	A = 1;
	B = 0;
	C = 0;
	D = 1;

	do
	{
		q = u / v;
		t = u;
		u = v;
		v = t - q * v;
		t = A;
		A = B;
		B = t - q * B;
		t = C;
		C = D;
		D = t - q * D;
	}
	while (v != 0);

	*ao = A;
	*bo = C;
	*vo = u;
}
/************************************************/
/*! \fn int inv_modCPUFF(int n, int P)
 \brief This function computes inverse of integer n over finite field.

 \param n an integer.
 \param P a prime.
 */
/************************************************/
int
inv_modCPUFF (int n, int P)
{
	int a, b, v;
	egcdCPUFF (n, P, &a, &b, &v, P);
	if (b < 0)
		b += P;
	return b % P;
}
/************************************************/
/*! \fn int mul_modCPUFF(int a, int b, int P)
 \brief This function computes product of integer a and b over finite field.

 \param a an integer.
 \param b an integer.
 \param P a prime.
 */
/************************************************/
int
mul_modCPUFF (int a, int b, int P)
{
	double ninv = 1 / (double) P;
	int q = (int) ((((double) a) * ((double) b)) * ninv);
	int res = a * b - q * P;
	res += (res >> BASE_1) & P;
	res -= P;
	res += (res >> BASE_1) & P;
	return res;
}

/************************************************/
/*! \fn int sub_modCPUFF(int a, int b, int P)
 \brief This function computes subtraction of integer a from b over finite field.

 \param a an integer.
 \param b an integer.
 \param P a prime.
 */
int
sub_modCPUFF (int a, int b, int P)
{
	int r = a - b;
	r += (r >> BASE_1) & P;
	return r;
}
/************************************************/
__device__ __host__ __inline__ int
add_mod (int a, int b, int p)
{
	int r = a + b;
	r -= p;
	r += (r >> BASE_1) & p;
	return r;
}
/************************************************/
__device__ __host__ __inline__ int
sub_mod (int a, int b, int p)
{
	int r = a - b;
	r += (r >> BASE_1) & p;
	return r;
}
/************************************************/
__device__ __host__ __inline__ int
neg_mod (int a, int p)
{
	int r = -a;
	r += (r >> BASE_1) & p;
	return r;
}
/************************************************/
__device__ __host__ __inline__ int
mul_mod (int a, int b, int n)
{
	double ninv = 1 / (double) n;
	int q = (int) ((((double) a) * ((double) b)) * ninv);
	int res = a * b - q * n;
	res += (res >> BASE_1) & n;
	res -= n;
	res += (res >> BASE_1) & n;
	return res;
}

/************************************************/
__device__ __host__ __inline__
int
pow_mod (int a, int ee, int n)
{
	int x, y;
	int e = ee;

	if (ee == 0)
		return 1;
	if (ee == 1)
		return a;

	x = 1;
	y = a;
	while (e)
	{
		if (e & 1)
			x = mul_mod (x, y, n);
		y = mul_mod (y, y, n);
		e = e >> 1;
	}

	if (ee < 0)
		x = inv_mod (x, n);
	return x;
}

/*! \fn __global__ void condensationFF(int *Aa, int *L, int n, int P)
 \brief This function computes one step of condensation.

 \param Aa an integer array that store the matrix before one condensation step.
 \param L the product of pivot elements
 \param n the dimension of matrix.
 \param P a prime.
 */
/************************************************/
__global__ void
condensationFF (int *Aa, int *lsd, int n, int P)
{

	int id = threadIdx.x; //0-127
	__shared__ int nonZeroNum;
	__shared__ int nonZeroIdx;
	int i = 0;
	int r = 0;
	int np = n - 1;
	int pos, pos1, currVal;
	nonZeroNum = 0;
	nonZeroIdx = 0;
	pos = 0;
	pos1 = id * n;
	for (r = 0; r < np; r++)
	{
		//find the nonzero number in the rth row
		if (id == 0)
		{
			nonZeroNum = 0;
			nonZeroIdx = -1;

			for (i = 0; i < n; i++)
			{
				if (Aa[pos + i] != 0)
				{
					nonZeroIdx = i;
					nonZeroNum = Aa[pos + i];
					lsd[r] = nonZeroNum;
					break;
				}
			}
		}
		__syncthreads ();
		if (nonZeroNum == 0)
			break; //there are no nonZeroNum in this row.
		if (id > r)
		{
			for (i = 0; i < nonZeroIdx; i++)
			{
				Aa[pos1 + i] = mul_mod ((-1) * Aa[pos1 + i], nonZeroNum, P);
			}
			currVal = Aa[pos1 + nonZeroIdx];
			for (i = nonZeroIdx + 1; i < n; i++)
			{
				Aa[pos1 + i] = sub_mod (mul_mod (nonZeroNum, Aa[pos1 + i], P),
																mul_mod (Aa[pos + i], currVal, P), P);
			}
			Aa[pos1 + nonZeroIdx] = 0;
		}
		__syncthreads ();
		pos += n;
	}
	if (id == 0)
	{
		pos = np * n;
		for (i = 0; i < n; i++)
		{
			if (Aa[pos + i] != 0)
			{
				lsd[np] = Aa[pos + i];
				break;
			}
		}
	}
}
/************************************************/
int
detFF (int *A, int n, int P)
{
	int *Aa, *ls, *lsd;
	int bm, bmv, i, j;
	float outerTime;

	cudaMalloc ((void **) &Aa, sizeof(int) * n * n);
	cudaMemcpy (Aa, A, sizeof(int) * n * n, cudaMemcpyHostToDevice);

	ls = (int *) malloc ((sizeof(int) * n));
	cudaMalloc ((void **) &lsd, sizeof(int) * n);
	cudaMemcpy (lsd, ls, sizeof(int) * n, cudaMemcpyHostToDevice);  //0

	cudaEvent_t start, stop;
	cudaEventCreate (&start);
	cudaEventCreate (&stop);
	cudaEventRecord (start, 0);

	//ls 0 
	condensationFF <<<1, n>>> (Aa, lsd, n, P);
	cudaThreadSynchronize ();
	cudaMemcpy (ls, lsd, sizeof(int) * n, cudaMemcpyDeviceToHost);
	cudaEventRecord (stop, 0);
	cudaEventSynchronize (stop);
	cudaEventElapsedTime (&outerTime, start, stop);

	j = 0;
	for (i = n - 3; i >= 0; i--)
	{
		j++;
		if (j > 1)
		{
			bm = ls[i];
			for (int k = 1; k < j; k++)
			{
				bm = mul_mod (bm, ls[i], P);
			}
		}
		else
		{
			bm = ls[i];
		}
		//bm=pow_mod(ls[i],j,P);
		bmv = inv_mod (bm, P);
		ls[n - 1] = mul_mod (ls[n - 1], bmv, P);
	}

	printf ("the determinant result is %d\n", ls[n - 1]);
	printf ("This program costs time %f ms\n", outerTime);
	cudaFree (Aa);
	cudaFree (lsd);
	free (ls);
	return 0;

}
/************************************************/
void
writefile (char *filename, int *A, int n, int P)
{
	int i, j, pos;
	FILE *fp2;
	char buf[1024];
	if ((fp2 = fopen (filename, "w")) == NULL)
	{
		printf ("fail to %s", filename);
		fclose (fp2);
		exit (-2);
	}
	memset (buf, 0, 1024);
	sprintf (buf, "P := %d:\n", P);
	fputs (buf, fp2);
	memset (buf, 0, 1024);
	sprintf (buf, "n := %d:\n", n);
	fputs (buf, fp2);
	memset (buf, 0, 1024);
	sprintf (buf, "[");
	fputs (buf, fp2);
	for (i = 0; i < n; i++)
	{
		memset (buf, 0, 1024);
		sprintf (buf, "[");
		fputs (buf, fp2);
		pos = i * n;
		for (j = 0; j < n; j++)
		{
			if (j < n - 1)
			{
				memset (buf, 0, 1024);
				sprintf (buf, "%d,", A[pos + j]);
				fputs (buf, fp2);
			}
			else
			{
				memset (buf, 0, 1024);
				sprintf (buf, "%d]", A[pos + j]);
				fputs (buf, fp2);
			}
		}
		if (i < n - 1)
		{
			memset (buf, 0, 1024);
			sprintf (buf, ",");
			fputs (buf, fp2);
		}
		else
		{
			memset (buf, 0, 1024);
			sprintf (buf, "];\n");
			fputs (buf, fp2);
		}
	}
	memset (buf, 0, 1024);
	sprintf (buf, "mat := Matrix(%) :\n");
	fputs (buf, fp2);
	memset (buf, 0, 1024);
	sprintf (buf, "with(LinearAlgebra) :\n");
	fputs (buf, fp2);
	fclose (fp2);
}
/************************************************/
int
main (int argc, char *argv[])
{
	int n = 128, i, P = 469762049;
	char fileName[1024];
	int *a;
	if (argc > 1)
		n = atoi (argv[1]);
	if (argc > 2)
		P = atoi (argv[2]);
	if (n > 1024)
	{
		printf (
				"The current program only supports to calculate the determinant of matrix with size not greater than 1024\n");
		exit (-1);
	}
	a = (int *) malloc ((sizeof(int) * n * n));

	for (i = 0; i < n * n; ++i)
	{
		a[i] = rand () % P;
	}
	memset (fileName, 0, sizeof(char) * 1024);
	sprintf (fileName, "determinant.input");
	writefile (fileName, a, n, P);
	printf ("The input matrix with size %d*%d is stored into %s\n", n, n,
					fileName);

	//	cout<<"Determinant (from Maple): ";
	//	i = system("maple -q mapleDet.mm");

	// the orinal array, 4
	int result = detFF (a, n, P);
	free (a);
	return 0;
}

