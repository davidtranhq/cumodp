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



#include "new_mixed_radix_rep.h"

/** Transpose polynominal from r*c to c*r
 * @transpose: Transposed polynominal
 * @polynominal: Original polynominal
 * @r: Original rows
 * @c: Original columns
 */
__global__ void transpose(sfixn* transpose, sfixn* polynominal, int r, int c) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < r * c) {
                transpose[idx] = polynominal[(idx%r)*c+idx/r];
        }
}


/** Compute elements needed in matrices
 ** m[i,j], where m^{-1}[i] = m[i,j] mod m[j] with i < j
 ** n[i,j] = m[j] - m[i,j]
 * @matrices: Array storing m[i,j] and n[i,j] continuously
 * @pnum: Number of primes
 * @m: Primes
 */
__global__ void initMRR(sfixn* matrices, int pnum, const int* m) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int row = pnum - 1;
	int size = row * row;
	if (idx < size) {
		int i = idx / row;
		int j = idx % row + 1;

		if (i < j) {
			matrices[idx] = inv_mod(m[i], m[j]);
			matrices[size+idx] = m[j] - matrices[idx];
		}
		else {
			matrices[idx] = 0;
		}
	}
}


/** Convert coefficients into Mixed Radix Representaion
 ** Compute ( .. ((X A_{0}) A_{1}) .. A_{s-1})
 * @mrr: Coefficients in MRR manner
 * @N: Size of original polynominal
 * @matrices: Array storing m[i,j] and n[i,j] continuously
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void MRR(sfixn* mrr, int N, sfixn* matrices, const int* m, int pnum) {
	int tid = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + tid;

	if (idx < N * pnum) {
		__shared__ sfixn mrrtmp[512], crttmp[512];
		mrrtmp[tid] = mrr[idx];

		int rows = pnum - 1;
		for (int k = 0; k < rows; k++) {
			crttmp[tid] = mrrtmp[tid];
			__syncthreads();

			// Compute X A_{k}
			int i = tid / pnum;
			int j = tid % pnum;
			if (k < j) {
				mrrtmp[tid] = crttmp[i*pnum+k] * matrices[rows*rows+k*rows+j-1] + crttmp[tid] * matrices[k*rows+j-1];
				mrrtmp[tid] = (int)mrrtmp[tid] % m[tid%pnum];
			}
		}

		__syncthreads();

		mrr[idx] = mrrtmp[tid];
	}
}


/** Compute the sum of signs changed within each thread block
 * @A: Array storing the sign changed
 * @N: Size of A
 */
__global__ void signsSum(int* A, int N, int TBS) {
	int tid = threadIdx.x;
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N/2) {
		__shared__ sfixn a[512];
		a[2*tid] = A[2*idx];
		a[2*tid+1] = A[2*idx+1];

		for (int i = TBS>>1; i > 0; i >>= 1) {
        		__syncthreads();
        	        if (tid < i) {
                		int ai = (N/(2*i))*(2*tid+1)-1;
                        	int bi = (N/(2*i))*(2*tid+2)-1;

                        	a[bi] += a[ai];
                	}
		}

		__syncthreads();

		A[2*idx] = a[2*tid];
		A[2*idx+1] = a[2*tid+1];
	}
}


/** Check the sign of a MRR coefficient
 ** If a = 0, return 0
 ** If a < 0, return -1
 ** Otherwise, return 1
 */
__device__ __host__ __inline__
int checkSign(sfixn* a, const int* m, int pnum) {
	int isZero = 0;
	for (int i = 0; i < pnum; i++) {
		if (!a[i]) { isZero++; }
	}
	if (isZero == pnum) { return 0; }

	for (int i = pnum-1; i > -1; i--) {
                int middle = (m[i]-1)/2;
                if (a[i] < middle) { return 1; }
                else if (a[i] > middle) { return -1; }
        }
	return 1;
}


/** Cout number of sign changed
 * @signs: If a coefficient's sign is changed, store 1, otherwise store 0;
 *	Then perform prefix sum on it, and the last element is the number of signs changed
 * @A: Polynominal
 * @N: Size of polynominal
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void signChange(int* signs, sfixn* A, int N, const int* m, int pnum) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N) {
		if (!idx) { signs[idx] = 0; }
		else {
			int ischange = checkSign(&A[(idx-1)*pnum], m, pnum) * checkSign(&A[idx*pnum], m, pnum);

			if (ischange < 0) { signs[idx] = 1; } 
			else { signs[idx] = 0; } 
		}
	}
}


/** Compare a and b
 ** If a > b, return 1
 ** If a < b, return -1
 ** Otherwise, return 0 
 */
__device__ int compareMRR(sfixn* a, sfixn* b, const int* m, int pnum) {
	int asign = checkSign(a, m, pnum);
	int bsign = checkSign(b, m, pnum);

	if ((asign >= 0 && bsign < 0) || (asign > 0 && bsign == 0)) { return 1; }
	else if ((asign < 0 && bsign >= 0) || (asign == 0 && bsign > 0)) {return -1; }
	else if (asign == 0 && bsign == 0) { return 0; }
	else {
		for (int i = pnum-1; i > -1; i--) {
			if (a[i] > b[i]) { return 1; }
			else if (a[i] < b[i]) { return -1; }
		}
		return 0;
	}
}


/** Compute the largest and smallest number among the coefficients within a thread block
 * @index: [0], store the index with smallest coefficient;
	 [N-1], store the index with largest coefficient
 * @A: Coefficients
 * @N: Size of A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Primes
 * @pnum: Number of primes
 */
__global__ void abslargestMRR(int* index, sfixn* A, int N, int TBS, const int* m, int pnum) {
	int tid = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx < N/2) {
		__shared__ sfixn itmp[512];
		itmp[2*tid] = 2*idx;
		itmp[2*tid+1] = 2*idx+1;

		for (int i = TBS>>1; i > 0; i >>= 1) {
        	        __syncthreads();
                	if (idx < i) {
				int offset = N/(2*i);
				// Compare to get largest coefficient
                        	int lai = itmp[offset*(2*idx+1)-1];
                        	int lbi = itmp[offset*(2*idx+2)-1];

				int lab = compareMRR(&A[lai*pnum], &A[lbi*pnum], m, pnum);

				// Compare to get smallest coefficient
				int sai = itmp[offset*2*idx];
				int sbi = itmp[offset*(2*idx+1)];

				int sab = compareMRR(&A[sai*pnum], &A[sbi*pnum], m, pnum);

				if (lab > 0) { itmp[offset*(2*tid+2)-1] = lai; }
				if (sab > 0) { itmp[offset*2*tid] = sbi; }
                	}
        	}

		__syncthreads();

		index[2*idx] = itmp[2*tid];
		index[2*idx+1] = itmp[2*tid+1];

	}
}
