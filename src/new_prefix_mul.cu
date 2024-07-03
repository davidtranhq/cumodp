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



#include "new_prefix_mul.h"

/** Initialize an array for factorial sequence
 * @A: Arrary storing 1 2 .. N .. 1 .. N
 * @N: Size of n!
 * @pnum: Number of primes
 */
__global__ void initArray(sfixn* A, int N, const int* m, int pnum) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N * pnum) {
		A[idx] = (idx % N + 1) % m[idx/N];
	}
}

/** Initialize an array for powers, with same value **/
__global__ void initArray(sfixn* A, int N, const int* m, int pnum, sfixn base) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        if (idx < N * pnum) {
		if (base < 0) 
			A[idx] = ((int)base % m[idx/N] + m[idx/N]) % m[idx/N];
		else
                	A[idx] = (int)base % m[idx/N];
        }
}


/** Compute prefix multiplication in each thread block
 * @A: Array {a_1, a_2, ..., a_B, a_{B+1}, ..., a_n}, and store the result {a_1, a_1*a_2, .., a_1*...a_B, a_{B+1}, ..., a_{x*B+1}*..*a_n}
 * @N: Size of array A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Prime
 */
__global__ void inclusivePreMul(sfixn* A, int N, int TBS, const int m) {
        int tid = threadIdx.x;
        int idx = blockIdx.x * blockDim.x + tid;

        __shared__ sfixn a[512], b[512];	// 2 * THREADBLOCK_SIZE

	if (idx < N/2) {
		int k = 2*tid;
		int l = 2*tid+1;
		int r = 2*idx;
		int s = 2*idx+1;

	       	// load input into shared memory
		a[k] = A[r];
	        a[l] = A[s];
	        b[k] = A[r];
	        b[l] = A[s];

		// build mul in place up the tree
	        for (int i = TBS>>1; i > 0; i >>= 1) {
	                __syncthreads();
	                if (tid < i) {
				int offset = TBS / (2 * i);
	                        int ai = offset*(2*tid+1)-1;
	                        int bi = offset*(2*tid+2)-1;

	               	        b[bi] = mul_mod(b[ai], b[bi], m);
		        }
		}

		// clear the last element
	       	if (tid == 0) { 
			b[TBS-1] = 1;
		}
		// traverse down tree & build scan
	        for (int i = 1; i < TBS; i *= 2) {
	                __syncthreads();
	                if (tid < i) {
				int offset = TBS / (2 * i);
	                        int ai = offset*(2*tid+1)-1;
	       	                int bi = offset*(2*tid+2)-1;

	       	                sfixn t = b[ai];
	                        b[ai] = b[bi];
                        	b[bi] = mul_mod(t, b[bi], m);
			}
       		}

	       	__syncthreads();

		// write results to device memory
		A[r] = mul_mod(a[k], b[k], m);
		A[s] = mul_mod(a[l], b[l], m);
	}
}


/** Update array
 * @A: Updated array
 * @B: Array {b_1, ..., b_{N/2*TBS}}, which is the last element of each thread block in A from the previous step
 * @N: Size of array A
 * @TBS: Number of coefficients dealt by a thread block
 * @m: Prime
 */
__global__ void exclusivePreMul(sfixn* A, sfixn* B, int N, int TBS, const int m) {
        int tid = threadIdx.x;
        int bid = blockIdx.x;
        int idx = bid * blockDim.x + tid;

	if (bid < N/TBS-1) {
	       	A[TBS+2*idx] = mul_mod(A[TBS+2*idx], B[bid], m);
       		A[TBS+2*idx+1] = mul_mod(A[TBS+2*idx+1], B[bid], m);
	}
}


/** Colloect the last element of each thread block of previous step 
 * @lasts: Last elements
 * @A: Original polynomial
 * @nb: # of thread blocks in previous step
 * @TBS: Number of coefficients dealt by a thread block during previous step
 */
__global__ void collectLast(sfixn* lasts, sfixn* A, int nb, int TBS) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nb) {
		lasts[idx] = A[TBS*(idx+1)-1];
	}
}


/** Compute the prefix multiplication
 * @A: Array {a_1, a_2, ..., a_n}, and store the result {a_1, a_1*a_2, .., a_1*..*a_n}
 * @N: Size of array A
 * @m: Primes
 * @pnum: Number of primes
 */
void prefixMul(sfixn* A, int N, const int* m, int pnum) {
	int threadsPerBlock = THREADBLOCK_SIZE;
	int TBS = (N > 2*threadsPerBlock) ? 2*threadsPerBlock : N;
	int blocksPerGrid = (N > TBS)? N / TBS : 1;
	size_t size_poly = N * pnum * sizeof(sfixn);
	int nb = N / TBS;
	size_t size_tmp = nb * sizeof(sfixn);

	sfixn* d_premul;
	cudaMalloc((void**)&d_premul, size_poly);
	sfixn* d_tmp;
	cudaMalloc((void**)&d_tmp, size_tmp);

	cudaMemcpy(d_premul, A, size_poly, cudaMemcpyHostToDevice);	

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (int i = 0; i < pnum; i++) {
	       	inclusivePreMul<<<blocksPerGrid, threadsPerBlock>>>(&d_premul[i*N], N, TBS, m[i]);

		if (nb > 1) {
			collectLast<<<1, threadsPerBlock>>>(d_tmp, &d_premul[i*N], nb, TBS);
			inclusivePreMul<<<1, threadsPerBlock>>>(d_tmp, nb, TBS, m[i]);
			exclusivePreMul<<<blocksPerGrid, threadsPerBlock>>>(&d_premul[i*N], d_tmp, N, TBS, m[i]);
		}
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float outerTime;
	cudaEventElapsedTime(&outerTime, start, stop);
	std::cout<<"Prefix Multiplication: "<<outerTime/1000.0<<std::endl;

	cudaMemcpy(A, d_premul, size_poly, cudaMemcpyDeviceToHost);
	
	cudaFree(d_premul);
	cudaFree(d_tmp);
}
