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



/*! \file condensationFiniteField.cu
\brief This file contains the code for computing determinant over finite field using CUDA.
*/
#include "condensationFiniteField.h"


using namespace std;


/*! \fn void egcdCPUFF(sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo, int P)
\brief This function computes extended gcd of integer x and y over finite field.

\param x an integer.
\param y an integer.
\param ao integer pointer.
\param bo integer pointer.
\param vo integer pointer.
\param P a prime.
*/		


void egcdCPUFF(sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo, int P) {
    sfixn t, A, B, C, D, u, v, q;

    u = y; v = x;
    A = 1; B = 0;
    C = 0; D = 1;

    do {
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
    } while (v != 0);

    *ao = A;
    *bo = C;
    *vo = u;
}


/*! \fn sfixn inv_modCPUFF(sfixn n, int P)
\brief This function computes inverse of integer n over finite field.

\param n an integer.
\param P a prime.
*/		

sfixn inv_modCPUFF(sfixn n, int P) {
    sfixn a, b, v;
    egcdCPUFF(n,  P, &a, &b, &v,P);
    if (b < 0) b += P;
    return b % P;
}



/*! \fn sfixn mul_modCPUFF(sfixn a, sfixn b, int P)
\brief This function computes product of integer a and b over finite field.

\param a an integer.
\param b an integer.
\param P a prime.
*/
sfixn mul_modCPUFF(sfixn a, sfixn b, int P) {
    double ninv = 1 / (double)P;
    sfixn q  = (sfixn) ((((double) a) * ((double) b)) * ninv);
    sfixn res = a * b - q * P;
    res += (res >> BASE_1) & P;
    res -=      P;
    res += (res >> BASE_1) & P;
    return res;
}


/*! \fn sfixn sub_modCPUFF(sfixn a, sfixn b, int P)
\brief This function computes subtraction of integer a from b over finite field.

\param a an integer.
\param b an integer.
\param P a prime.
*/
sfixn sub_modCPUFF(sfixn a, sfixn b, int P) {
    sfixn r = a - b;
    r += (r >> BASE_1) & P;
    return r;
}





/*! \fn __global__ void condensationFF(int *Aa, int *Bb, int *L, size_t n, int P,  int BK)
\brief This function computes one step of condensation.

\param Aa an integer array that store the matrix before one condensation step.
\param Bb an integer array that store the matrix after one condensation step.
\param L the product of pivot elements
\param n the dimension of intermediate matrix.
\param P a prime.
\param BK an integer that stores the size of tiles that a thread works.
*/


__global__ void condensationFF(int *Aa, int *Bb, int *L, size_t n, int P,  int BK)
{

        int k = (blockIdx.x*blockDim.x + threadIdx.x)*BK;

	if(L[0] >=0 && k<n*n)
	{
		int lValue = Aa[ (n+1)*L[0]  ];
		int firstRow;
		int p = k%n;
		int q = k/n;//floor((double) k/ (double)n);
		firstRow = Aa[(n+1)*(q+1) ];			
				
		int i, j, s, d, start, limitA, limitB;				
		limitB = n*n;
		limitA = limitB + 2*n + 1;
		for(s=1;s<BK;++s)
		{
			if((k+s)%n ==0)
				break;
		}
		if(q >= L[0])
		{	
			i = k+q+2+n;
                        start = L[0]*(n+1)+p+1;
			for(d=0; d < s && (d+i) < limitA && (k+d) < limitB  ;++d)
				Bb[k+d]= sub_mod(mul_mod(lValue,Aa[d+i], P), mul_mod( Aa[start+d], firstRow, P  ) , P);				
                }
		else
		{
			i = k+q+1;
                        for(d=0; d < s && (d+i) < limitA && (k+d) < limitB;  ++d)
				Bb[k+d]=  neg_mod(  mul_mod(lValue,Aa[i+d], P) , P);			
                }			 
                
		q = q+1;
		if(s < BK && q < n)
		{						
			if(q>=L[0])
                 	{
				firstRow = Aa[(n+1)*(q+1) ];			
				k = n*q;				
				i = (n+1)*L[0] + 1;
				d = (n+1)*(q+1) + 1;
				 
				for(j = s, start=0; j < BK; ++j, ++start)
					Bb[k + start] = sub_mod(mul_mod(lValue,Aa[d + start], P), mul_mod( Aa[i+start], firstRow , P  ), P);				
                	}
                	else
                	{
				k = n*q;
				d = (n+1)*(q) + 1;
				for(j = s, start=0; j < BK; ++j, ++start)
					Bb[k+start]= neg_mod(  mul_mod(lValue, Aa[ d + start ], P), P);				
                	}			
		}			
	}

}

/*! \fn __global__ void findingLFF(int *A,  int *B,  int *C, size_t n, int P)
\brief This function computes the pivot element and store the product of the pivot elements.

\param A an integer array that store the matrix before one condensation step.
\param B an integer array that store the matrix after one condensation step.
\param C the product of pivot elements
\param n the dimension of intermediate matrix.
\param P a prime.
*/

__global__ void findingLFF(int *A,  int *B,  int *C, size_t n, int P)
{
	
	if(B[0]>=0)
	{
		int i;
		for(i = 0; i < n; ++i)
		{
			if(A[i*n]!=0)
				break;
		}
		if(i<n)
		{
			B[0] = i;
			int x =  pow_mod(A[i*n],n-2, P);
			C[0] =  mul_mod(C[0], x, P);
		}
		else
		{
			B[0] = -1;
		}
	}
}


/*! \fn int detFF(int *A, int n, int P, int BK, int TH, int BS, int verbose)
<a name="detFFCPU"> </a>
\brief This is the CPU function to be called by external functions (i.e. C code).
It calls the kernels for conputing determinant
over a prime field.


\param A an integer array that store the matrix before any condensation step.
\param n the dimension of input matrix.
\param P a prime.
\param BK an integer that fix the tile size for each thread.
\param TH the thread number as suggested from calling function.
\param BS an integer that indicates the dimension of the matrix when condensation method should stop working.
Currently this option is not used to keep it independent. Previousely, this function depends on NTL library for when 
the code reached here.
\param verbose unused currently. It exists for future as a means to communicate the calling function.
*/

int detFF(int *A, int n, int P, int BK, int TH, int BS, int verbose)
{
        int *Aa, *Bb, *B, *lvalues, *ls,*lvaluesd, *lsd, *swap,i;	
	//n1 =n;
	
	int *BB;

	cudaMalloc((void **)&Aa, sizeof(int)*n*n);
        cudaMemcpy(Aa, A, sizeof(int)*n*n, cudaMemcpyHostToDevice);

	B = (int *)malloc(sizeof(int)*(n-1)*(n-1));
	for(i =0; i<(n-1)*(n-1); ++i) B[i] = 0;
	cudaMalloc((void **)&Bb, sizeof(int)*(n-1)*(n-1));
        cudaMemcpy(Bb, B, sizeof(int)*(n-1)*(n-1), cudaMemcpyHostToDevice);

	lvaluesd = (int *)malloc(sizeof(int));
	lsd = (int *)malloc(sizeof(int));
	lvaluesd[0]= 1;
	lsd[0] = 0;
	
        cudaMalloc((void **)&ls, sizeof(int));
 	cudaMemcpy(ls, lsd, sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void **)&lvalues, sizeof(int));
 	cudaMemcpy(lvalues, lvaluesd, sizeof(int), cudaMemcpyHostToDevice);
        
	dim3 threadsPerBlock(TH);
        dim3 numBlocksN;
//	initialization<<<1, 1>>>(ls,lvalues);
	
       cudaEvent_t start, stop;
	cudaEventCreate(&start);
         cudaEventCreate(&stop);

        cudaEventRecord(start, 0);


	while(n>=3)
	{
        	findingLFF<<<1, 1>>>(Aa, ls, lvalues,  n, P);

        
		--n;
		if(n<BK)
			BK = n;

        	numBlocksN=(ceil((double)(n*n)/(double)(TH*BK)));
        	condensationFF<<<numBlocksN, threadsPerBlock>>>(Aa,Bb,ls, n, P, BK);

		swap=Aa;
		Aa = Bb;
		Bb = swap;
	}


	 cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float outerTime;
        cudaEventElapsedTime(&outerTime, start, stop);
	 //printf("%d \t",n1);
	//if (verbose != 0)
	//{
	      //  printf("%f \t",outerTime/1000.0);
		//printf(" s\n");
	//}

	detFFtime = outerTime/1000.0;
	//long q=n;
	
	BB = (int *)malloc(sizeof(int)*n*n);
        cudaMemcpy(BB, Aa, sizeof(int)*n*n, cudaMemcpyDeviceToHost);
 	
	int temp1 = sub_modCPUFF(mul_modCPUFF(BB[0],BB[3],P), mul_modCPUFF(BB[1],BB[2],P),P);        

	//int temp2;
	free(BB);        
	
	BB = (int *)malloc(sizeof(int));
	cudaMemcpy(BB, lvalues, sizeof(int), cudaMemcpyDeviceToHost);
	
	int temp2 = mul_modCPUFF(inv_modCPUFF(BB[0],P), temp1, P);
	
	
	free(BB);

        cudaFree(Aa);
        cudaFree(Bb);
        cudaFree(lvalues);
        cudaFree(ls);
	
	free(B);
	//free(A);
	free(lsd);
	free(lvaluesd);
	return temp2;
	
}



