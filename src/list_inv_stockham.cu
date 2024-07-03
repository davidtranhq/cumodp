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



#include "list_inv_stockham.h"


/**
 * Scaling after doing the normal fft with invw
 * 
 * X: the array that needs to be scaled
 * invn: the scaling factor
 * m: the number of polynomials [ignore it]
 * p: the prime
 * pinv: inverse of p
 * length_layer: the bound of X
 * */
__global__ void list_inv_mul_ker(sfixn *X,sfixn invn,  sfixn m, sfixn p, double pinv, sfixn length_layer)
{
	sfixn tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < length_layer)
	{
		//sfixn v = X[tid];
		X[tid] = mul_mod(invn, X[tid], p, pinv);
	}
}

// [EXPIRED]
__global__ void inv_mul_ker(sfixn *X, sfixn invn, sfixn m, sfixn p, double pinv)
{
	sfixn tid = threadIdx.x + blockIdx.x * blockDim.x;

	sfixn start_offset = tid * 16;
	
	for (sfixn i = 0; i<16; i++)
		X[start_offset + tid] =  mul_mod(invn, X[tid], p, pinv);
}


/**
 * @X, input / output host array of length m * n 
 * @m,
 * @k, n = 2^k
 * @w, n-th primitive root of unity
 * @p, fourier prime number
 *
 * X will be filled by I_m @ DFT_n(X, w)
 */
void list_inv_stockham_host(sfixn *X, sfixn m, sfixn k, sfixn w, sfixn p) {
    sfixn *X_d;
	sfixn invw = inv_mod(w, p);
	sfixn n = (1 << k);
	sfixn invn = inv_mod(n, p);
    cudaMalloc((void**)&X_d, sizeof(sfixn) * (m << k));
    cudaMemcpy(X_d, X, sizeof(sfixn) * (m << k), cudaMemcpyHostToDevice);
    ///////////////////////////////////////
	list_stockham_dev(X_d, m, k, w, p);
    list_stockham_dev(X_d, m, k, invw, p);
	if( k >= 16 )
		inv_mul_ker<<<(n/512),512>>>(X_d, invn, m, p, (1/(double)p));
	else
		list_inv_mul_ker<<<n,m>>>(X_d, invn, m, p, (1/(double)p),m*n);
	//inv_mul_ker<<<1,1>>>(X_d, n, invn, m, p, (1/(double)p));
    ///////////////////////////////////////
    cudaMemcpy(X, X_d, sizeof(sfixn) * (m << k), cudaMemcpyDeviceToHost);
    cudaFree(X_d);

	//for (sfixn i = 0; i < (m << k); i++) X[i] = mul_mod(invn, X[i], p, 1/(double)p);
}

/////////////////////////////////////////////////////////////////////////////
//BEGIN:list_inv_stockham_tst
/////////////////////////////////////////////////////////////////////////////

void list_inv_stockham_tst(sfixn p, sfixn m, sfixn k)
{
	sfixn *X = (sfixn *)malloc((m<<k)*sizeof(sfixn));
	sfixn w = primitive_root(k,p);
	sfixn n = (1L << k);
	
	for(sfixn u = 0; u < m; u++)
		for(sfixn v = 0; v < n; v++)
			X[u*n +v] = v;
	#ifndef _mcompile_
	printf("Input: \n");
	printf("w = %d\n", w);
	for (sfixn i = 0; i < m; i++) print_vector(n, X + i * n);
	
	list_inv_stockham_host(X, m, k, w, p);

	printf("Output: \n");
	for (sfixn i = 0; i < m; i++) print_vector(n, X + i * n);

	for(sfixn u = 0; u < m; u++)
		for(sfixn v = 0; v < n; v++)
			if( X[u*n +v] != v){
				printf("FAIL at u= %d, v=%d\n",u,v);
				return;
			}
	printf("PASS");
	#endif

}
/////////////////////////////////////////////////////////////////////////////
//END:list_inv_stockham_tst
/////////////////////////////////////////////////////////////////////////////
