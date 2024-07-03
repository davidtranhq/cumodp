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



#include "opt_plain_mul.h"

using namespace std;

/*
This function creates a zero array.
@C: an array.
@n: length of C.
		
*/

__global__ void zeroC(sfixn *C, sfixn n)
{
	sfixn i = (blockIdx.x*blockDim.x + threadIdx.x)*BLOCK;
	if(i < n)
	{
		sfixn j = i + BLOCK;
		for(; i < j && i < n; ++i)
			C[i] = 0;	
	}
}

/*
This function shifts A to the left by BLOCK-1 coefficients.
It puts (BLOCK-1)  zeros at the front and end of A.
C is used as workspace and it is assumed that C is an exact copy of A
(That is the i-th slots of A and C contain the same value). 
@A, storage for first polynomial.
@C, temporary array, which is a copy of A.
@n, the degree of A.
*/

__global__ void copyA(sfixn *A, sfixn *C, sfixn n)
{
	sfixn i = (blockIdx.x*blockDim.x + threadIdx.x)*BLOCK;
	if(i < n)
	{
		sfixn j,k;
		if(i == 0)
		{
			j = BLOCK -1;
			for(; i < j; ++i)
				A[i] = 0; 
			i = n + BLOCK -1;
			j = i + BLOCK -1;
			for(; i < j; ++i)
				A[i] = 0;
			i = 0;
		}
		j = i + BLOCK;
		k = i + (BLOCK -1);
		for(; i < j && i < n; ++i, ++k)
			A[k] = C[i];		
	}					
}

/*
This function makes the final addition step of vectors into the result (namely CgpuFinal).
By vectors, we mean components of the array V in the HPCS 2012 paper.
@*Cgpu, the intermediate addition from all vectors. This is an 1-D arraye representing 2-D.
@CgpuFinal, the final additions
@c, the length of  one row of Cgpu if we visualize it as 2-D array.
@n,  the length of  CgpuFinal.
@w, which row of Cgpu stores all intermediate addition result if we visualize it as 2-D array.
@P, the prime number for finite field.
*/
				
__global__ void addFinal(sfixn *Cgpu, sfixn *CgpuFinal, sfixn c, sfixn n, sfixn w, sfixn P)
{
	sfixn i = (blockIdx.x*blockDim.x + threadIdx.x)*BLOCK;
	if(i < c)
	{
		
		sfixn k = i + w*BLOCK;
		i = i + c*w;
		sfixn j = i + BLOCK;
		
		for(; i < j && k < n; ++i, ++k )
			CgpuFinal[k] =  add_mod( Cgpu[i], CgpuFinal[k], P);
		
	}
}		


/*
This function makes all addition steps (except the last one) of vectors into the result (namely CgpuFinal)
By vectors, we mean components of the array V in the HPCS 2012 paper.
	*/
		
__global__ void add(sfixn *Cgpu, sfixn *CgpuFinal, sfixn c, sfixn n, sfixn w, sfixn i, sfixn j, sfixn k, sfixn P)
{
	sfixn virtualBlockID = blockIdx.x % w;
	sfixn rowNo = (blockIdx.x/ w);
	
	sfixn firstRow = i + rowNo*k;
	
	
	sfixn limitStart = (virtualBlockID*blockDim.x + threadIdx.x)*BLOCK;
	sfixn limitEnd, r1start, fStart, secondRow;	
	r1start = (firstRow*c) +  limitStart;
	limitEnd = limitStart + BLOCK;
	
	if(limitStart < j*BLOCK)
	{	
		/*
		rows to be added are not aligned. We are adding corresponding coefficients. So some 
		coefficients will be added directly to the CgpuFinal
		*/	
		if(limitEnd > j*BLOCK)	limitEnd = j*BLOCK;
		if(limitEnd > c) limitEnd = c;		
		fStart = limitStart + firstRow*BLOCK; 
		for(; limitStart < limitEnd; ++limitStart, ++r1start, ++fStart )
			CgpuFinal[fStart] = add_mod(Cgpu[r1start], CgpuFinal[fStart], P);				
	}		
	else
	{
		if(limitEnd > c) limitEnd = c;
		secondRow = firstRow + j;
		sfixn r2start = secondRow*c + limitStart - j*BLOCK;
		for(; limitStart < limitEnd; ++r1start, ++r2start,  ++limitStart)
			Cgpu[r2start] = add_mod(Cgpu[r1start], Cgpu[r2start], P);
	}
}

/*
This is the main kerenel for multiplication. More precisely, it does the multiplication
phasis described in the HPCS 2012 paper. Note that Cgpu in the code below corresponds
to V in the paper.
@A, the first polynomial.
@B, the  second polynomial.
@C, the 1-D array that we will treat as 2-D array. It stores all intermediate multiplication results.
@n, the degree of A.
@m,  the degree of B.
@c, the length of  one row of Cgpu if we visualize it as 2-D array.
@P, the prime for prime field.
@unitBlocks, number of thread block responsible for producing one row of C (if we visualize it as 1-D array).
*/
	
__global__ void mul(sfixn *A, sfixn *B, sfixn *C, sfixn n, sfixn m, sfixn c, sfixn P, sfixn unitBlocks)
{
	__shared__ sfixn sA[T];
	__shared__ sfixn sB[BLOCK];
	sfixn Res[BLOCK];
	/*
	The block are assigned from left to right. At end it starts assigning from left again. 
	virtualBlockID is the block id from left of the multiplication space.
	*/
	sfixn virtualBlockID = blockIdx.x % unitBlocks;

	sfixn i = (virtualBlockID*blockDim.x + threadIdx.x)*BLOCK;	
	sfixn j = i + BLOCK;
	sfixn k = threadIdx.x*BLOCK;
	sfixn l = k + BLOCK;
	
	sfixn whereInB = (blockIdx.x/ unitBlocks);
	/*
	Each thread copies a part from A and B
	*/
	for(; i < j && i < n; ++i, ++k) 
		sA[k] = A[i];
	for(; k < l ; ++k)
		sA[k] = 0;
	/*
	First thread of each block do some extra work. 
	Put some extra zeros in sB at both ends. If the length
	of B is not divisible by BLOCK then the extra coefficients
	are copied here. The same is true for A.
	
	*/
	/*if(threadIdx.x == 0)
	{
		i = whereInB*BLOCK;
		j = i + BLOCK;		
		k = 0;
		for(; i < j && i < m; ++i, ++k)	sB[k] = B[i];
		for(; k < BLOCK; ++k) sB[k] = 0;
		j = (virtualBlockID*blockDim.x + Tx)*BLOCK;
		for(i = Tx*BLOCK; i < T && j < n; ++i, ++j )	sA[i] = A[j];
		for(; i < T; ++i) sA[i] = 0;
	}
	__syncthreads();*/


        i = whereInB * BLOCK + threadIdx.x;
        if (i < m && threadIdx.x < BLOCK)
                sB[threadIdx.x] = B[i];
        else if (threadIdx.x < BLOCK)
                sB[threadIdx.x] = 0;

        i = Tx * BLOCK + threadIdx.x;
        j = (virtualBlockID * blockDim.x + Tx) * BLOCK + threadIdx.x;
        if (j < n && i < T)
                sA[i] = A[j];
        else if (i < T)
                sA[i] = 0;
        __syncthreads();


	l = threadIdx.x*BLOCK + BLOCK -1;
	whereInB = whereInB *c;
	sfixn s = (virtualBlockID*blockDim.x + threadIdx.x)*BLOCK;	

	/*
	Each thread is responsible for BLOCK coefficients and copied back to global C
	*/
	for(i = 0; i < BLOCK && s < c; ++i, ++s) 
	{		
		Res[i] = 0;	
		k = l + i;		
		for(j = 0; j <BLOCK; ++j, --k)
		{
			Res[i] = add_mod(Res[i] , mul_mod(sA[k], sB[j], P), P);			
		}
		C[s + whereInB] = Res[i];
	}		
}

/*! \fn void mulCUDA( sfixn *A,  sfixn *B, sfixn *C, sfixn n, sfixn m, sfixn Prime, sfixn verify)
<a name="multiCPU"> </a>
\brief This is the CPU function to be called by external functions (i.e. C code).
It calls the kernels for conputing univariate polynomial multiplication
over a prime field.
It computes C = A*B over a finite field.
@A, the first polynomial.
@B, the second polynomial.
@C, the polynomial that stores A*B.
@n, the degree of A.
@m, the degree of B. n>=m.
@Prime, prime number, the characterisitic of the finite field.
@verify, whether we want to check the result. it is 1 if we want to check the result.
*/
void mulCUDA( sfixn *A,  sfixn *B, sfixn *C, sfixn n, sfixn m, sfixn Prime, sfixn verify)
{
	sfixn i, j, k, l, w;
	sfixn *Agpu, *Bgpu, *Cgpu, *CgpuFinal;

	
	/* W.r.t. the HPCS 2012 Paper "Plain Polynomial Arithmetic on GPU"
           we have r = t = BLOCK and (m+1)/r = RowNo  
           log_2 of rowNoPerfect is the number of parallel steps in the addition phase */
      

	sfixn aRowSize = n + BLOCK -1;
	sfixn RowNo = (sfixn)ceil((double)m/(double)BLOCK);
	sfixn RowNo2Pow = (sfixn)ceil(log2((double)RowNo));
	sfixn rowNoPerfect = (sfixn)pow(2, (double)RowNo2Pow);
	

	cudaMalloc((void **)&Agpu, ( sizeof(sfixn)*(n + 2*(BLOCK-1)))  );
	cudaMalloc((void **)&Bgpu, sizeof(sfixn)*m);
	cudaMalloc((void **)&Cgpu, sizeof(sfixn)*(aRowSize*rowNoPerfect));

        cudaMemcpy(Bgpu, B, sizeof(sfixn)*m, cudaMemcpyHostToDevice);
        cudaMemcpy(Cgpu, A, sizeof(sfixn)*n, cudaMemcpyHostToDevice);

	/* total_blocksX is the total number of blocks */
	sfixn total_blocksX = (sfixn)ceil((double)(n) / (double)(BLOCK*Tx));
	dim3 threadsPerBlock(Tx); 
	dim3 numBlocks(total_blocksX);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
        cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	/*
	First we need to shift Agpu by left BLOCK-1 bits. And put BLOCK-1 zeros at the head
	and tail of A. So we do it with the help of Cgpu. We copy Agpu into Cgpu and then 
	shift Agpu to left. Later we make all entries of Cgpu zeros.
	*/
	copyA<<<numBlocks, threadsPerBlock>>>(Agpu, Cgpu,  n);	

        /* Next, we intialize all coefficients of C to zero */
	numBlocks.x = (sfixn)ceil( (double)(aRowSize*rowNoPerfect)/(double)(BLOCK*Tx) );	
	zeroC<<<numBlocks, threadsPerBlock>>>(Cgpu, aRowSize*rowNoPerfect);


        /* BLOCK*Tx is the number of coefficients per thread block 
           numBlocks.x is the number of thread blocks in the multiplication phasis */
	sfixn n1 = n + (BLOCK-1);
	total_blocksX = (sfixn)ceil((double)(n1) / (double)(BLOCK*Tx));
	numBlocks.x = total_blocksX* RowNo;

	n1 += (BLOCK-1);
	/*
	Each thread will compute BLOCK X BLOCK rectangle multiplication coefficients and stored
	into Cgpu.
	*/
	mul<<<numBlocks, threadsPerBlock>>>(Agpu, Bgpu, Cgpu,  n1, m, aRowSize, Prime, total_blocksX);	
	
	/* CgpuFinal will hold the result after the addition phasis completes
           where Cgpu holds the "partial sums" computed during the multiplication phasis */ 
	i =(n + m -1);
	cudaMalloc((void **)&CgpuFinal, sizeof(sfixn)*i );
	numBlocks.x = (sfixn)ceil((double)(i) / (double)(BLOCK*Tx));
	zeroC<<<numBlocks, threadsPerBlock>>>(CgpuFinal, i);
	/*
        For the general principle of the addition phasis, see the HPCS paper.
        Here's an example.
	Let we have 0, 1, 2, 3, 4, 5 ,6 ,7 rows of intermediate results. We have to add in log8 or 3 steps.
	First Step: 
		add row 0 and row 1 and store into row 1
		add row 2 and row 3 and store into row 3
		add row 4 and row 5 and store into row 5
		add row 6 and row 7 and store into row 7
	Second Step:
		add row 1 and row 3 and store into row 3
		add row 5 and row 7 and store into row 7
	Third step:
		add row 3 and row 7 and store into row 7
	
	this is done in add kernel

	Finally, addFinal copy row 7 to CgpuFinal


	*/
	i = 0;
	j = 1;
	k = 2;
	w = (sfixn)ceil((double)(aRowSize) / (double)(BLOCK*Tx));
	for(l = 0; l < RowNo2Pow; ++l)
	{		 
		numBlocks.x = w*((sfixn)pow(2, (double)RowNo2Pow-l-1));
		add<<<numBlocks, threadsPerBlock>>>(Cgpu, CgpuFinal, aRowSize, m+n-1, w, i, j, k, Prime);			

		i = i + j;
		j = 2*j;
		k = 2*k;
					
	}	

	numBlocks.x = w;
	addFinal<<<numBlocks, threadsPerBlock>>>(Cgpu, CgpuFinal, aRowSize, m+n-1, i, Prime);			


	cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float outerTime;
        cudaEventElapsedTime(&outerTime, start, stop);
	#ifndef _mcompile_
	if(verify != 1)
		cout<<n<<" "<<m<<" "<<outerTime/1000.0<<endl;

	#endif
	cudaMemcpy(C, CgpuFinal, sizeof(sfixn)*(m+n-1), cudaMemcpyDeviceToHost);
	
	/* 
		MULTIPLICATION RESULT IS IN CgpuFinal
	*/

	
       cudaFree(Agpu);
        cudaFree(Bgpu);
	
        cudaFree(Cgpu);
	cudaFree(CgpuFinal);
}

