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



#include "new_taylor_shift.h"

/** Compute factorial sequence (x+1)^n, with each coefficient n!/[k! * (n-k)!]
 * @A: Array storing (x+c)^{2^i}, 0 \leq i < \log[2](N) continuously
 * @facseq: Factorial sequence storing 1!, ..., (\frac{N}{2})!
 * @powersc: Powers of c
 * @N: Size of original  polynominal
 * @m: Primes
 * @pnum: Number of primes
 ** Regard to (x+c)^{i}, the leading coefficient is always 1, such that it is not stored
 */
__global__ void binomialCoef(sfixn* A, sfixn* facseq, sfixn* powersc, int N, const int* m, int pnum) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        if (idx < (N-1)*pnum) {
		int x = idx / (N-1);
		int y = idx % (N-1);

        	int l = (int)(log((float)(y+1))/log((float)2));
                int n = (int)pow((float)2, (float)l);
                int k = y - n + 1;
                sfixn elem;
                if (k) {
                	sfixn elem1 = mul_mod(facseq[x*N+k-1], facseq[x*N+n-k-1], m[x]);
                        elem = quo_mod(facseq[x*N+n-1], elem1, m[x]);
               	}
                else { elem = 1; }
                A[x*(N-1) + y] = mul_mod(powersc[x*N+n-k-1], elem, m[x]);
        }
}


/** Compute each multiplication of polynomials in each thread block
 * @Mgpu1: Polynomials
 * @Mgpu2: Polynomials after multiplication
 * @length_poly: Size of each polynomial
 * @poly_on_layer: # of polynomials doing multiplication
 * @threadsForAmul: = # of threads per thread block / mulInThreadBlock
 * @mulInThreadBlock: # of multiplications computed in a thread block
 * @p: Prime number
 */
__global__ void modifiedListPlainMulGpu(sfixn *Mgpu1, sfixn *Mgpu2, int length_poly, int poly_on_layer, int threadsForAmul, int mulInThreadBlock, int p) {
        __shared__ sfixn sM[2048];
        int mulID = threadIdx.x/threadsForAmul;

        if (mulID < poly_on_layer/2 && threadIdx.x < threadsForAmul*mulInThreadBlock) {
                int start_offset = blockIdx.x*length_poly*2;
                int j = start_offset + mulID*length_poly*2;
                int q = start_offset + mulID*(2*length_poly-1);

                int t = threadIdx.x/threadsForAmul;
                int u = threadIdx.x % threadsForAmul;

                int s = t*(4*length_poly-1);
                int k = s + length_poly;
                int l = k + length_poly;
                int c = l+u;
                int a, b;

                sM[s+u] = Mgpu1[j + u];
                __syncthreads();

                if (u != 2*length_poly-1) {
                        if (u < length_poly) {
                                a = s;
                                b = k + u;
                                sM[c] = mul_mod(sM[a], sM[b], p);
                                ++a; --b;
                                for(int i = 0; i < u; ++i, ++a, --b) {
                                        sM[c] = add_mod(mul_mod(sM[a], sM[b], p), sM[c], p);
                                }
                                Mgpu2[q+u] = sM[c];
                        }
                        else {
                                b = l - 1;
                                a = (u - length_poly) + 1 + s;
                                sM[c] = mul_mod(sM[a], sM[b], p);
                                ++a; --b;

                                int utmp = (2*length_poly-2) - u;
                                for(int i = 0; i < utmp; ++i, ++a, --b) {
                                        sM[c] = add_mod(mul_mod(sM[a], sM[b], p), sM[c], p);
                                }
                                Mgpu2[q+u] = sM[c];
                        }
                }
                else {
                        Mgpu2[q+u] = 0;
                }
        }
}


/** Compuse input polynomial for multiplication in the sub-tree of Taylor Shift
 * @mul: Output polynomial
 * @poly: Original polynomial
 * @N: Size of polynomial
 * @facseq: Factorial sequence
 * @base: Size of part of the polynomial that muplies (c+x)^base
 */
__global__ void composeMulPoly(sfixn* mul, sfixn* poly, int N, sfixn* facseq, int base) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N/2) {
		mul[2*idx-idx%base] = poly[base+2*idx-idx%base];
		mul[base+2*idx-idx%base] = facseq[base-1+idx%base];
	}
}


/** Addition of two polynomials
 * @A: Polynomial
 * @B: Polynomial
 * @N: Size of polynomial
 */
__global__ void addPoly(sfixn* A, sfixn* B, int N) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < N) {
		A[idx] += B[idx];
	}
}


/** Shift polynomial by c
 * @poly: Input polynomial
 * @N: Size of polynomial
 * @c: Shift by c
 * @m: Primes
 * @pnum: Number of primes
 */
void TayShi(sfixn* poly, int N, sfixn c, const int* m, int pnum) {
	int threadsPerBlock = THREADBLOCK_SIZE;
	int blocksPerGrid;
	int mn = N * pnum;
	size_t size_poly = mn * sizeof(sfixn);

        int* d_m;
        cudaMalloc((void**)&d_m, pnum*sizeof(int));
        cudaMemcpy(d_m, m, pnum*sizeof(int), cudaMemcpyHostToDevice);

	// Compute n!
        sfixn* h_premul = new sfixn[mn];

        sfixn* d_premul;
        cudaMalloc((void**)&d_premul, size_poly);

        blocksPerGrid = (mn > threadsPerBlock)? mn / threadsPerBlock : 1;
        initArray<<<blocksPerGrid, threadsPerBlock>>>(d_premul, N, d_m, pnum);
        cudaMemcpy(h_premul, d_premul, size_poly, cudaMemcpyDeviceToHost);
        prefixMul(h_premul, N, m, pnum);
	cudaMemcpy(d_premul, h_premul, size_poly, cudaMemcpyHostToDevice);

        // Compute powers of c
        sfixn* h_cpowers = new sfixn[mn];

        sfixn* d_cpowers;
        cudaMalloc((void**)&d_cpowers, size_poly);

        blocksPerGrid = (mn > threadsPerBlock)? mn / threadsPerBlock : 1;
        initArray<<<blocksPerGrid, threadsPerBlock>>>(d_cpowers, N, d_m, pnum, c);
        cudaMemcpy(h_cpowers, d_cpowers, size_poly, cudaMemcpyDeviceToHost);
        prefixMul(h_cpowers, N, m, pnum);
	cudaMemcpy(d_cpowers, h_cpowers, size_poly, cudaMemcpyHostToDevice);


	// Compute factorial sequence
	// (c+x)^{1} (c+x)^{2} (c+x)^{4} .. (c+x)^{2^i}
	sfixn* d_facseq;
	cudaMalloc((void**)&d_facseq, size_poly);

	blocksPerGrid = (N * pnum > threadsPerBlock)? N * pnum / threadsPerBlock : 1;
	binomialCoef<<<blocksPerGrid, threadsPerBlock>>>(d_facseq, d_premul, d_cpowers, N, d_m, pnum);


///***************************** TEST ONLY *******************************
// Factorial Sequence
sfixn* h_facseq = new sfixn[N*pnum];
cudaMemcpy(h_facseq, d_facseq, size_poly, cudaMemcpyDeviceToHost);
std::ofstream facseqfile("factorial_gpu.dat", std::ofstream::out);
for (int i = 0; i < pnum; i++) {
        for (int j = 0; j < N-1; j++) {
                facseqfile << h_facseq[i*(N-1)+j] << " ";
        }
        facseqfile << "\n";
}
facseqfile.close();
// ***********************************************************************/


	// Taylor Shift by c using plain multiplication
	sfixn* d_poly;
	cudaMalloc((void**)&d_poly, size_poly);
	sfixn* d_tmp;
	cudaMalloc((void**)&d_tmp, N*sizeof(sfixn));
	sfixn* d_mul;
	cudaMalloc((void**)&d_mul, N*sizeof(sfixn));

	cudaMemcpy(d_poly, poly, size_poly, cudaMemcpyHostToDevice);
	for (int i = 0; i < pnum; i++) {
		int plainsize = (N > 1024)? 1024 : N;
		for (int j = 1; j < plainsize; j *= 2) {
			threadsPerBlock = THREADBLOCK_SIZE;
			blocksPerGrid = (N > 2*threadsPerBlock)? N/(2*threadsPerBlock) : 1;
			composeMulPoly<<<blocksPerGrid, threadsPerBlock>>>(d_tmp, &d_poly[i*N], N, &d_facseq[i*(N-1)], j);

			threadsPerBlock = 2 * j;
			blocksPerGrid = (N > threadsPerBlock)? N / threadsPerBlock : 1;
			int npoly = threadsPerBlock / j;
			//threadsForAmul = threadsPerBlock / mulInThreadBlock

			modifiedListPlainMulGpu<<<blocksPerGrid, threadsPerBlock>>>(d_tmp, d_mul, j, npoly, threadsPerBlock, 1, m[i]);

			threadsPerBlock = THREADBLOCK_SIZE;
			blocksPerGrid = (N > threadsPerBlock)? N / threadsPerBlock : 1;
			addPoly<<<blocksPerGrid, threadsPerBlock>>>(&d_poly[i*N], d_mul, N);
		}
	}
	cudaMemcpy(poly, d_poly, size_poly, cudaMemcpyDeviceToHost);


	cudaFree(d_facseq);
	cudaFree(d_premul);
	cudaFree(d_cpowers);
	cudaFree(d_poly);
	cudaFree(d_tmp);
	cudaFree(d_mul);
}
