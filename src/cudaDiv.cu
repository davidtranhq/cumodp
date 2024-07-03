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




/*! \file cudaDiv.cu
\brief This file contains the code for univariate polynomial division over finite field using CUDA.

The functions in this file use an array called status[0..4] to keep track 
some important information about the polynomials and the intermediate division
steps. The description of each entry of the array is given below. 


 status[0]: the first entry is the degree of the bigger polynomial, called A 
          (this property holds after executing the kernel called statusUpdate).

status[1]: the second entry is the degree of the second (or smaller) polynomial, called B.

status[2]: the third entry is the prime number (cardinality of the underlying field).

status[3]:  the fourth one is the number of required thread blocks in order to perform a "block of division steps".

status[4]: the fifth one indicates whether we are done or not -1 means iff we are not done.
     
*/

#include<iostream>
#include "cudaDiv.h"



/*! \fn __global__ void statusUpdate(sfixn *A,  sfixn *status)
\brief After s division steps, this function
checks whether we need to do any more step
to complete the division.

\param A the first (bigger) polynomial.
\param status the status array.
*/		
__global__ void statusUpdate(sfixn *A,  sfixn *status)
{
	if(status[4] == -1)
	{
		if( (status[0] - s) < status[1] )
		{
			//if( (status[0] - s) >= 0 ) status[4] = status[0] - s;
			//else 
			status[4] = -3;
		}
		else
		{	
			for(status[0] -= s; status[0] >= status[1]; --status[0])
			{
				if(A[ status[0] ] != 0)	
					break;
			}
			if( status[0] <  status[1])
			{
				status[4] = -4;
			}
		}
	}
}

/*! \fn __global__ void zero(sfixn *A, sfixn n) 
\brief This function makes all entries of A zero. This is used to initialize quotient and remainder array.

\param A an array
\param n the length of A.
*/
__global__ void zero(sfixn *A, sfixn n)
{
	sfixn k = (blockIdx.x*blockDim.x + threadIdx.x);
	if(k < n)
		A[k] = 0;
}

/*! \fn __global__ void	copyRem(sfixn *A, sfixn *R, sfixn *status)
\brief This function creates the remainder.
Once the division procedure is done,
the previously bigger polynomial becomes the remainder.
This function extracts the remainder polynomial from A and writes into R.
\param A the bigger polynomial.
\param R the remainder.
\param status the status array.
*/
__global__ void	copyRem(sfixn *A, sfixn *R, sfixn *status)
{
	sfixn k = (blockIdx.x*blockDim.x + threadIdx.x);
	if(k <= status[0] )
	{
		R[k] = A[k];
	}
}

/*! \fn __global__ void	reduceM(sfixn *B,  sfixn *status)
\brief Before  starting the computation,
the division function has to be sure that the 
leading coefficient of the small 
polynomial is nonzero.
This function finds out the first 
nonzero coefficient of the smaller polynomial, B.
\param B the smaller polynomial.
\param status the array that keeps the status.
*/
__global__ void	reduceM(sfixn *B,  sfixn *status)
{	
	for(; status[1] >= 0; --status[1])
	{
		if(B[ status[1] ] != 0)	
			break;
	}	
	if(status[1] < 0)
	{
		status[4] = -2;
	}
		
}

/*! \fn __global__ void divGPU(sfixn *A, sfixn *B, sfixn *status, sfixn *Q, sfixn *R)
\brief This function performs s division  steps in the Euclidean division of A by B.
\param A the first polynomial.
\param B the second polynomial.
\param status the status array.
\param Q quotient.
\param R remainder.
*/
__global__ void divGPU(sfixn *A, sfixn *B, sfixn *status, sfixn *Q, sfixn *R)
{

	if(status[4] == -1)
	{
		sfixn i, j;



		__shared__ sfixn headA;
		__shared__ sfixn headB;
		__shared__ sfixn startA,  startBind;

		__shared__ sfixn inv, p;
		__shared__ sfixn sA[t*s];
		__shared__ sfixn sB[s + t*s];
		__shared__ sfixn sAcommon[s];
		__shared__ sfixn sBcommon[s];

		//Deciding from which index a thread will copy to its shared memory.
		if(threadIdx.x == 0)
		{
					
			headA = status[0];
			headB = status[1];
			startA = headA -s - s*t*blockIdx.x;			
			startBind = headB  - s*t*blockIdx.x -1;
			p = status[2];
			inv = inv_mod(B[headB],p);		
				
		}
		__syncthreads();
	
		// Copynig the necessary coefficients into shared memory.
		if(headA >= headB && headB >= 0)
		{
			j = headA - threadIdx.x;	
			if(j >= 0)
				sAcommon[threadIdx.x] = A[j];
			else
				sAcommon[threadIdx.x] = 0;
			j = headB - threadIdx.x;
			if(j >= 0)
				sBcommon[threadIdx.x] = B[j];
			else
				sBcommon[threadIdx.x] = 0;	
		
			j = startA - t*threadIdx.x; 
			for(i = 0; i < t; ++i)
			{
				if((j - i) >= 0 )
					sA[threadIdx.x*t + i] = A[j-i];
				else
					sA[threadIdx.x*t + i] = 0;			
			}		
	
			j = startBind - (t+1)*threadIdx.x;
			for(i = 0; i <= t; ++i)
			{
				if((j - i) >= 0 )
					sB[threadIdx.x*(t+1) + i] = B[j-i];
				else
					sB[threadIdx.x*(t+1) + i] = 0;				
			}
			__syncthreads();
		// The division steps	
			sfixn factor;		
			for(i = 0; i < s && (headA-i) >= headB; ++i)
			{			
				factor = mul_mod(sAcommon[i], inv, p);

				if(blockIdx.x == 0 && threadIdx.x == 0)
					Q[ headA - i - headB ] = factor;
				if(i <= threadIdx.x )
					sAcommon[threadIdx.x ] = sub_mod(sAcommon[threadIdx.x ], mul_mod(sBcommon[threadIdx.x - i], factor, p), p);
	
				for(j = 0; j < t; ++j)
					sA[threadIdx.x*t + j ] = sub_mod(sA[threadIdx.x*t + j ], mul_mod(sB[s - 1 - i + threadIdx.x*t + j], factor, p), p);
						
				__syncthreads();	
			}		
		//check for the termination of division.
			if(i < s && blockIdx.x == 0 && threadIdx.x == 0)
			{
				status[4] = -5;
				status[0] -= s;		
				 
				for(status[3] = headB - 1; i < s &&  status[3] >= 0; ++i, --status[3])
				{
					R[status[3]] = sAcommon[i];							
				}

			}		
				
			//write back
			

			j = startA - t*threadIdx.x; 
			for(i = 0; i < t && (j - i) >= 0; ++i)
				A[j - i ] = sA[threadIdx.x*t + i];
			


			
		}
	}
	
}


/*!\fn __global__ void MulScalar( sfixn *Q, sfixn *A, sfixn n, sfixn factor, sfixn p)
\brief It simply multiplies each coefficient of A with factor and savies it into Q.
This function assumes that the degree of B is zero.
\param Q univariate polynomial.
\param A univariate polynomial.
\param n the length of A
\param factor an integer in field.
\param p the prime number of finite field.
*/
__global__ void MulScalar( sfixn *Q, sfixn *A, sfixn n, sfixn factor, sfixn p)
{
        sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n) Q[tid] = mul_mod(A[tid], factor, p);
}



/*! \fn float divPrimeField(sfixn *A, sfixn *B, sfixn *R, sfixn *Q, sfixn n, sfixn m, sfixn p)
<a name="divCPU"> </a>
\brief This is the CPU function to be called by external functions (i.e. C code).
It calls the kernels for conputing univariate polynomial division
over a prime field.
It computes Q and R such that A = BQ + R holds modulo the prime number p
and either R=0 holds or we have degree(R) < degree(B).
This function returns the time spent on the device and the allocation/de-allocation/transfer times.
\param A the first polynomial.
\param B the second polynomial.
\param Q quotient.
\param R remainder.
\param n the degree of polynomial in A is n-1.
\param m the degree of polynomial in B is m-1. Note that n >= m is assumed.
\param p the characterisitic of the prime field.
\return it will return the time (in second) that  it takes to compute the quatient and remainder of A and B.
*/
float divPrimeField(sfixn *A, sfixn *B, sfixn *R, sfixn *Q, sfixn n, sfixn m, sfixn p)
{

	sfixn numBlocksN, i;
	sfixn *Agpu, *Bgpu, *Qgpu;  

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
        cudaEventCreate(&stop);

	float outerTime;

	if(m == 1 )
	{
		R[0] = 0;
		if(B[0] == 0)
		{			
			//this is wrong as it is divide by zero.
			Q[0] = 0;		
		}
		else
		{
			i =  inv_mod(B[0],p);	
			if(n > 1)
			{
				cudaMalloc((void **)&Agpu, sizeof(sfixn)*n);
				cudaMalloc((void **)&Qgpu, sizeof(sfixn)*n);			

				cudaMemcpy(Agpu, A, sizeof(sfixn)*n, cudaMemcpyHostToDevice);
			
				cudaEventRecord(start, 0);
				numBlocksN = (sfixn)(ceil((double)( n)/(double)(T)));
				MulScalar<<<numBlocksN, T>>>(Qgpu, Agpu, n, i, p);
	
				cudaEventRecord(stop, 0);
        			cudaEventSynchronize(stop);
        		        cudaEventElapsedTime(&outerTime, start, stop);

				cudaMemcpy(Q, Qgpu,  sizeof(sfixn)*n, cudaMemcpyDeviceToHost);

				cudaFree(Agpu);
				cudaFree(Qgpu);
				return outerTime;
			}
			else Q[0] = mul_mod(A[0], i, p);

		}		
		return 0.0;
	}


	
	sfixn *Rgpu,*stateGpu;

	cudaMalloc((void **)&Agpu,     sizeof(sfixn)*n);
	cudaMalloc((void **)&Bgpu,     sizeof(sfixn)*m);
	cudaMalloc((void **)&Qgpu,     sizeof(sfixn)*(n - m + 1));
	cudaMalloc((void **)&Rgpu,     sizeof(sfixn)*(m - 1));
	cudaMalloc((void **)&stateGpu, sizeof(sfixn)*5);

        cudaMemcpy(Agpu, A, sizeof(sfixn)*n, cudaMemcpyHostToDevice);	
        cudaMemcpy(Bgpu, B, sizeof(sfixn)*(m), cudaMemcpyHostToDevice);        

	
	//Both A and B has at least one nonzero coefficient
       
	sfixn state[5]= { n-1 + s, m-1, p, 0, -1};
	cudaMemcpy(stateGpu, state, sizeof(sfixn)*5, cudaMemcpyHostToDevice);

	
	cudaEventRecord(start, 0);

	numBlocksN = (sfixn)(ceil((double)(n - m + 1)/(double)(T)));
	zero<<<numBlocksN, T>>>(Qgpu, (n - m + 1));

	numBlocksN = (sfixn)(ceil((double)(m-1.0)/(double)(T)));
	zero<<<numBlocksN, T>>>(Rgpu, (m - 1));

        // Ensures that the divisor (i.e. B) is non-zero
	reduceM<<<1, 1>>>(Bgpu, stateGpu);

	////
	//numBlocksN = (sfixn)(ceil((double)( m - 1.0)/(double)(T*t)));
	////

	numBlocksN = (sfixn)(ceil((double)( m )/(double)(T*t)));

	for(i = n; i>= m; i -= s )
	{
		statusUpdate<<<1, 1>>>(Agpu, stateGpu);
		divGPU<<<numBlocksN, T>>>(Agpu, Bgpu, stateGpu, Qgpu, Rgpu);
	}
	numBlocksN = (sfixn)(ceil((double)(m-1)/(double)(T)));
        // etxracts the remainder R from A
	copyRem<<<numBlocksN, T>>>(Agpu, Rgpu, stateGpu);


	cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        
        cudaEventElapsedTime(&outerTime, start, stop);

	

	cudaMemcpy(state, stateGpu,  sizeof(sfixn)*5, cudaMemcpyDeviceToHost);
	
		
	
	cudaMemcpy(Q, Qgpu,  sizeof(sfixn)*(n-m+1), cudaMemcpyDeviceToHost);
	
	
	cudaMemcpy(R, Rgpu,  sizeof(sfixn)*(m - 1), cudaMemcpyDeviceToHost);
	
	
	
	cudaFree(Agpu);
        cudaFree(Bgpu);
        cudaFree(stateGpu);
        cudaFree(Qgpu);
        cudaFree(Rgpu);
	return outerTime;

}





