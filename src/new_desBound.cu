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



#include "new_desBound.h"

/** Compute P_{k, c} = 2^{k*n} P(\frac{x+c}{2^{k}})
 * @h_A: Input polynomial storing as ( f_0 .. f_d ) mod prime[0] .. ( f_0 .. f_d ) mod prime[np-1]
 * @h_B: Ouput polynominal after computing PKC
 * @N: Number of coefficients of input polynominal
 * @pnum: Number of primes
 * @k: Input k
 * @c: Input c
 */
void pkc_cpu(sfixn* h_A, sfixn* h_B, int N, const int* m, int pnum, int k, sfixn c) {
	int threadsPerBlock = THREADBLOCK_SIZE;
	int blocksPerGrid;
	int mn = N * pnum;
	size_t size_poly = mn * sizeof(sfixn);

	sfixn* d_A;
	cudaMalloc((void**)&d_A, size_poly);
	int* d_m;
	cudaMalloc((void**)&d_m, pnum*sizeof(int));

	cudaMemcpy(d_A, h_A, size_poly, cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, m, pnum*sizeof(int), cudaMemcpyHostToDevice);

	// Compute powers of 2^k
	sfixn basek = pow(2, k);
	sfixn* h_2powers = new sfixn[mn];

	sfixn* d_2powers;
	cudaMalloc((void**)&d_2powers, size_poly);

	blocksPerGrid = (mn > threadsPerBlock)? mn / threadsPerBlock : 1;
	initArray<<<blocksPerGrid, threadsPerBlock>>>(d_2powers, N, d_m, pnum, basek);
	cudaMemcpy(h_2powers, d_2powers, size_poly, cudaMemcpyDeviceToHost);
	prefixMul(h_2powers, N, m, pnum);
	cudaMemcpy(d_2powers, h_2powers, size_poly, cudaMemcpyHostToDevice);

	// Divide by a power of 2^k
	blocksPerGrid = (N > threadsPerBlock)? N / threadsPerBlock : 1;
	DivPow<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_2powers, N, d_m, pnum);
	cudaMemcpy(h_B, d_A, size_poly, cudaMemcpyDeviceToHost);


	// Compute Taylor Shift
	TayShi(h_B, N, c, m, pnum);

	// Muliply 2^{k*(n-1)}
	blocksPerGrid = (N > threadsPerBlock)? N / threadsPerBlock : 1;
	cudaMemcpy(d_A, h_B, size_poly, cudaMemcpyHostToDevice);
	MulPow<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_2powers, N, d_m, pnum);
	cudaMemcpy(h_B, d_A, size_poly, cudaMemcpyDeviceToHost);


	cudaFree(d_A);
	cudaFree(d_m);
	cudaFree(d_2powers);
}


/** Represent coefficients in MRR manner, compute the number of signs changed, and the largest positive number and the smallest negative number
 * @h_A: Input polynominal, storing as ( f_0 .. f_d ) mod prime[0] .. ( f_0 .. f_d ) mod prime[np-1]
 * @h_B: Ouput polynominal after computing PKC
 * @N: Number of coefficients of input polynominal
 * @pnum: Number of primes
 * @sc: Number of signs changed
 * @mrrPosNorm: Largest positive number represented in MRR manner
 * @mrrNegNorm: Smallest negative number represented in MRR manner
 */
void mrr_cpu(sfixn* h_A, sfixn* h_B, int N, const int* m, int pnum, int* sc, sfixn* mrrPosNorm, sfixn* mrrNegNorm) {
	int threadsPerBlock = THREADBLOCK_SIZE;
	int blocksPerGrid;


	size_t size_poly = N * pnum * sizeof(sfixn);
	int MN = 2 * (pnum - 1) * (pnum - 1);
	size_t size_matrices = MN * sizeof(sfixn);

	sfixn* d_A;
	cudaMalloc((void**)&d_A, size_poly);
	sfixn* d_transA;
	cudaMalloc((void**)&d_transA, size_poly);
	sfixn* d_B;
	cudaMalloc((void**)&d_B, size_poly);
	sfixn* d_transB;
	cudaMalloc((void**)&d_transB, size_poly);
	int* d_m;
	cudaMalloc((void**)&d_m, pnum*sizeof(int));
	sfixn* d_matricesA;
	cudaMalloc((void**)&d_matricesA, size_matrices);

	cudaMemcpy(d_A, h_A, size_poly, cudaMemcpyHostToDevice);
	cudaMemcpy(d_B, h_B, size_poly, cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, m, pnum*sizeof(int), cudaMemcpyHostToDevice);


	// Transpose polynominals
	blocksPerGrid = N * pnum / THREADBLOCK_SIZE;
	if (!blocksPerGrid) { blocksPerGrid = 1; }
	transpose<<<blocksPerGrid, threadsPerBlock>>>(d_transA, d_A, pnum, N);
	transpose<<<blocksPerGrid, threadsPerBlock>>>(d_transB, d_B, pnum, N);


	// Convert coefficients in MRR manner
	blocksPerGrid = MN / threadsPerBlock;
	if (!blocksPerGrid) { blocksPerGrid = 1; }
	initMRR<<<blocksPerGrid, threadsPerBlock>>>(d_matricesA, pnum, d_m);


	blocksPerGrid = N * pnum / threadsPerBlock;
	if (!blocksPerGrid) { blocksPerGrid = 1; }
//	MRR<<<blocksPerGrid, threadsPerBlock>>>(d_transA, N, d_matricesA, d_m, pnum);
	MRR<<<blocksPerGrid, threadsPerBlock>>>(d_transB, N, d_matricesA, d_m, pnum);

	cudaMemcpy(h_A, d_transA, size_poly, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_B, d_transB, size_poly, cudaMemcpyDeviceToHost);	


	// Compute number of signs changed
	int* d_signs;
	cudaMalloc((void**)&d_signs, N*sizeof(int));

	blocksPerGrid = (N > threadsPerBlock)? N / threadsPerBlock : 1;
	signChange<<<blocksPerGrid, threadsPerBlock>>>(d_signs, d_transB, N, d_m, pnum);

	int TBS = (N > 2*threadsPerBlock)? 2*threadsPerBlock : N;
	signsSum<<<blocksPerGrid, threadsPerBlock>>>(d_signs, N, TBS);


	int* h_signs = new int[N];
	cudaMemcpy(h_signs, d_signs, N*sizeof(int), cudaMemcpyDeviceToHost);
	sc[0] = h_signs[N-1];
	for (int i = 1; i < N/TBS; i++) {
		sc[0] += h_signs[N-1-i*TBS];
	}


	// Compute the largest and smallest coefficient
	abslargestMRR<<<blocksPerGrid, threadsPerBlock>>>(d_signs, d_transB, N, TBS, d_m, pnum);

	cudaMemcpy(h_signs, d_signs, N*sizeof(int), cudaMemcpyDeviceToHost);
/*	int smallest = h_signs[0];
	for (int i = 0; i < pnum; i++) {
                mrrNegNorm[i] = h_B[smallest*pnum+i];
        }

	int largest = h_signs[N-1];
	for (int i = 0; i < pnum; i++) {
                mrrPosNorm[i] = h_B[largest*pnum+i];
        }
*/

	cudaFree(d_A);
	cudaFree(d_transA);
	cudaFree(d_B);
	cudaFree(d_transB);
	cudaFree(d_m);
	cudaFree(d_matricesA);
	cudaFree(d_signs);
}
