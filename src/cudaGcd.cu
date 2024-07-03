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



/*! \file cudaGcd.cu
\brief This file contains the code for univariate polynomial GCD over finite field using CUDA.

The functions in this file use an array called status[0..5] to keep track 
some important information about the polynomials and the intermediate division
steps. The description of each entry of the array is given below. 

status[0] store the degree of polynomial A.

status[1] store the degree of polynomial B.

status[2] stores the prime.

status[3] stores the number of requires thread blocks

status[4] is very important.
status[4] == -1 means we can continue.
status[4] == -2 means both A and B  contains null polynomial.
status[4] == -3 means A becomes null polynomial.
status[4] == -4 means B becomes null polynomial.

status[5] tells us the which polynomial is treated as small polynomial in the next division steps.
status[5] = 0 means A is smaller.
status[5] = 1 means B is smaller.
*/



#include<iostream>
#include "cudaGcd.h"

using namespace std;


/*! \fn __global__ void	shrinkA(sfixn *A, sfixn *status)
\brief Shrinking the degree of A so that the leading coefficient is nonzero.
\param A the polynomial .
\param status The status array.
*/
__global__ void	shrinkA(sfixn *A, sfixn *status)
{
	for(; status[0] >= 0; --status[0])
	{
		if(A[ status[0] ] != 0)	
			break;
	}			
}


/*! \fn __global__ void	shrinkB(sfixn *B, sfixn *status)
\brief Shrinking the degree of B so that the leading coefficient is nonzero.
\param B the polynomial.
\param status The status array.
*/
__global__ void	shrinkB(sfixn *B, sfixn *status)
{
	for(; status[1] >= 0; --status[1])
	{
		if(B[ status[1] ] != 0)	
			break;
	}			
}

/*! \fn __global__ void MonicConversionA( sfixn *A, sfixn *status)
\brief It simply divide each coefficient of A by its leading coefficient.
\param A univariate polynomial.
\param status The status array.
*/
__global__ void MonicConversionA( sfixn *A, sfixn *status)
{
        sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
        if(tid < status[0]) A[tid] = mul_mod(inv_mod(A[ status[0] ],status[2]), A[tid], status[2]); 
}

/*! \fn __global__ void MonicConversionB( sfixn *B, sfixn *status)
\brief It simply divide each coefficient of B by its leading coefficient.
\param B univariate polynomial.
\param status The status array.
*/
__global__ void MonicConversionB( sfixn *B, sfixn *status)
{
        sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < status[1]) B[tid] = mul_mod(inv_mod(B[ status[1] ],status[2]),  B[tid], status[2]);
}


/*! \fn __global__ void	reduceMgcd(sfixn *B,  sfixn *status)
\brief Before starting any division steps,
We have to be sure that the leading cofficient 
of polynomial B is nonzero. This function ensure this.
status[1] store the degree of polynomial B.
status[2] stores the prime.
\param B the polynomial B.
\param status The status array.
*/
__global__ void	reduceMgcd(sfixn *B,  sfixn *status)
{
	if(status[4] == -1)
	{		
		for(; status[1] >= 0; --status[1])
		{
			if(B[ status[1] ] != 0)	
				break;
		}
	}		
}

/*! \fn __global__ void	reduceNgcd(sfixn *A,  sfixn *status)
\brief Before starting any division steps,
We have to be sure that the leading cofficient 
of polynomial A is nonzero. This function ensure this.
status[0] store the degree of polynomial A.
status[2] stores the prime.
\param A the polynomial A.
\param status The status array.
*/
__global__ void	reduceNgcd(sfixn *A,  sfixn *status)
{	
	if(status[4] == -1)
	{
		for(; status[0] >= 0; --status[0])
		{
			if(A[ status[0] ] != 0)	
				break;
		}
	}			
}

/*! \fn __global__ void	status3(sfixn *status)
\brief Compute status[3] that tells us how many thread block will be active during gcd computation.
\param status The status array.
*/
__global__ void	status3(sfixn *status)
{
	if(status[4] == -1)
	{
	     if(status[1] <  status[0]) status[3] = (status[1] + T)/(2*T) + 1;
	     else status[3] = (status[0] + T)/(2*T) + 1;
	}	
}

/*! \fn __global__ void	status4(sfixn *status)
\brief compute status[4].
\param status The status array.
*/
__global__ void	status4(sfixn *status)
{
	if(status[4] == -1)
	{
		if(status[0] <  0 && status[1] <  0)   status[4] = -2;
		if(status[0] <  0 && status[1] >= 0 )  status[4] = -3;
		if(status[0] >= 0 && status[1] <  0 )  status[4] = -4;
	}	
}

/*! \fn __global__ void	status5(sfixn *status)
\brief computes status[5] that  tells us the which polynomial is treated as small polynomial in the next division steps.
\param status The status array.
*/
__global__ void	status5(sfixn *status)
{
	if(status[4] == -1)
	{
	          if(status[0] <  status[1] && status[5] == 1) status[5] = 0;
             else if(status[0] >  status[1] && status[5] == 0) status[5] = 1;	     
	}	
}

/*! \fn __global__ void gcdGPU(sfixn *A, sfixn *B, sfixn *status)
\brief This is the main step of computing gcd.
\param A The first polynomial.
\param B The second polynomial.
\param status The status array.
*/
__global__ void gcdGPU(sfixn *A, sfixn *B, sfixn *status)
{
	
        // The second contion below checks whether a thread block is active or not
	if(status[4] == -1 &&  blockIdx.x < status[3])
	{
		sfixn i, j, k;
		

		__shared__ sfixn startAcom, startBcom, endAcom, endBcom;		
		__shared__ sfixn startAsh,  startBsh, endAsh, endBsh;
		__shared__ sfixn headA, headB, startA, startB;
		__shared__ sfixn selectedA, factor, p;
		__shared__ sfixn sA[3*T], sB[3*T], sAcommon[T], sBcommon[T];
		//from where each thread will copy.
		if(threadIdx.x == 0)
		{					
			headA = status[0]; headB     = status[1]; 
			p     = status[2]; selectedA = status[5];

			startA = status[0] - 2*T*blockIdx.x;
			startB = status[1] - 2*T*blockIdx.x;

			startAcom = 0; 	 startBcom = 0; 
			endAcom   = T;	 endBcom   = T;
			startAsh  = 0;   startBsh  = 0; 
			endAsh    = 3*T; endBsh    = 3*T;
		
		
			
					
		}
		__syncthreads();		
		// copying the coefficients from the global to shared memory.
		i = 3*threadIdx.x;
		sAcommon[threadIdx.x] = 0;
		sBcommon[threadIdx.x] = 0;
		sA[i] = 0;    sA[i+1] = 0;  sA[i+2] = 0;
		sB[i] = 0;    sB[i+1] = 0;  sB[i+2] = 0;
		k = threadIdx.x;			
		if( (headA - k) >= 0)	sAcommon[threadIdx.x] = A[headA - threadIdx.x];		
		if( (headB - k) >= 0) sBcommon[threadIdx.x] = B[headB - threadIdx.x];		
		
		j = startA - 3*threadIdx.x;		
		for(k = 0; k < 3 && (j-k) >= 0; ++k)
			sA[i + k] =  A[j - k];	
		
		j = startB - 3*threadIdx.x;		
		for(k = 0; k < 3 && (j-k) >= 0; ++k)
			sB[i + k] = B[j - k];
		
		__syncthreads();


		
		
	 	//begin gcd steps

		for(;selectedA >= 0;)
		{			
			for(;selectedA == 1;)
			{				
				if(threadIdx.x == 0)
					factor = mul_mod(sAcommon[startAcom], inv_mod(sBcommon[startBcom],p), p);
				__syncthreads();
				
				j = threadIdx.x + startBcom;
				i = threadIdx.x + startAcom;
				if(j < endBcom)
					sAcommon[i] =  sub_mod(sAcommon[i], mul_mod(sBcommon[j],factor,p), p);
				__syncthreads();

				j = (threadIdx.x*3) + startBsh;
				i = (threadIdx.x*3) + startAsh;
				for(k = 0; k < 3 && (j+k) < endBsh; ++k)
					sA[i+k] = sub_mod(sA[i+k], mul_mod(sB[j+k], factor, p), p);
				__syncthreads();
				
				if(threadIdx.x == 0)
				{
					++startAcom; --endBcom;	++startAsh; --endBsh; --startA; --headA;
					k = 0;
					if(startAcom < endAcom)
					{
						for(;sAcommon[startAcom] == 0; ++k, ++startAcom, --headA)
						{						
							if(startAcom >= endAcom || headA < 0) 
								break;						
						}
					}
					if(k > 0)
					{
						endBcom -= k;	startAsh += k; endBsh -= k; startA -= k; 
						
					}
					if( (startAcom >= endAcom)|| headA < 0 )
					{
						if(blockIdx.x == 0) status[5] = selectedA;
						 selectedA = -1;
					}
					else if(headA < headB) 
						selectedA = 0;										
				}						
				__syncthreads();				
			}							
			for(;selectedA == 0;)
			{
				if(threadIdx.x == 0)
					factor = mul_mod(sBcommon[startBcom], inv_mod(sAcommon[startAcom],p), p);
				__syncthreads();

				j = threadIdx.x + startAcom;
				i = threadIdx.x + startBcom;
				if(j < endAcom)
					sBcommon[i] = sub_mod(sBcommon[i], mul_mod(sAcommon[j],factor,p), p);
				__syncthreads();

				j = (threadIdx.x*3) + startAsh;
				i = (threadIdx.x*3) + startBsh;
				for(k = 0; k < 3 && (j+k) < endAsh; ++k)
					sB[i+k] = sub_mod(sB[i+k], mul_mod(sA[j+k], factor, p), p);
				__syncthreads();
			
				if(threadIdx.x == 0)
				{
					++startBcom; --endAcom;	++startBsh; --endAsh; --startB; --headB;
					k = 0;
					if(startBcom < endBcom)
					{
						for(;sBcommon[startBcom] == 0; ++k, ++startBcom, --headB)
						{						
							if(startBcom >= endBcom || headB < 0) 
								break;						
						}
					}
					if(k > 0)
					{
						endAcom -= k;	startBsh += k; endAsh -= k; startB -= k; 
						
					}
					if( (startBcom >= endBcom) || headB < 0 )
					{
						if(blockIdx.x == 0) status[5] = selectedA; 						 
						selectedA = -1;
					}
					else if(headB < headA) 
						selectedA = 1;
					
				}		
				__syncthreads();
			}
			
				
		}

		i = startA - threadIdx.x*2;
		if( i >= 0 )
		{
			j = startAsh + threadIdx.x*2;
			A[ i ] = sA[ j ];
			--i;
			if( i >= 0 )
				A[ i ] = sA[ j + 1 ];
			
		}
		
		i = startB - threadIdx.x*2;
		if( i >= 0 )
		{
			j = startBsh + threadIdx.x*2;
			B[ i ] = sB[ j ];
			--i;
			if( i >= 0 )
				B[ i ] = sB[ j + 1 ];				
		}
		
		if(threadIdx.x == 0 && blockIdx.x == 0)
		{
			status[0] = headA;
			status[1] = headB;						
		}

		
	}	
}


/*! \fn float gcdPrimeField(sfixn *Dgcd,sfixn *A, sfixn *B, sfixn *G, sfixn n, sfixn m, sfixn p, sfixn optionNormalize)
<a name="gcdCPU"> </a>
\brief This is the main cpu function that is called from outside 
for computing G = GCD(A,B) over a (small) finite field
\param Dgcd the degree of gcd will be passed through this pointer
\param A The first polynomial.
\param B The second polynomial.
\param G the gcd of A and B.
\param n the degree of A is n-1.
\param m the degree of B is m-1. n >= m
\param p the characteristic of the field
\return it will return the time (in second) that  it takes to compute the gcd of A and B.
*/
float gcdPrimeField(sfixn *Dgcd,sfixn *A, sfixn *B, sfixn *G, sfixn n, sfixn m, sfixn p, sfixn optionNormalize)
{
	sfixn numBlocksN;
	float outerTime;

	if(m <= 1 && n <= 1)
	{
		G[0] = 1;
		if(A[0] == 0 && B[0] == 0) G[0] = 0;		
		(*Dgcd) = 0;
		return 0.0;
	}
	
	
	sfixn *Agpu, *Bgpu, *stateGpu;

	cudaMalloc((void **)&Agpu,     sizeof(sfixn)*n);
	cudaMalloc((void **)&Bgpu,     sizeof(sfixn)*m);
	cudaMalloc((void **)&stateGpu, sizeof(sfixn)*6);

        cudaMemcpy(Agpu, A, sizeof(sfixn)*n, cudaMemcpyHostToDevice);	
        cudaMemcpy(Bgpu, B, sizeof(sfixn)*(m), cudaMemcpyHostToDevice);        

	
	//Both A and B have at least one nonzero coefficient
        // In the array state:
        // the first entry is the degree of the first poly 
        // the second entry is the degree of the second poly
        // the third entry is the prime
        // the fourth one is the number of requires thread blocks
        //  in order to perform a "block of division steps" 
        // the fifth one is defined before, see the kernel status4
        // the sixth one indicates which poly has the largest degree
	sfixn state[6]= { n-1, m-1, p, 0, -1, 1};
	if(m > n) state[5] = 0;
	cudaMemcpy(stateGpu, state, sizeof(sfixn)*6, cudaMemcpyHostToDevice);

	///
	//cout<<state[0]<<" "<<state[1]<<" "<<state[2]<<" "<<state[3]<<" "<<state[4]<<" "<<state[5]<<" "<<endl;
	///

	sfixn i, j;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
        cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	
	/*
	T = s, according to our HPCS 2012 paper.
	*/

	if(m == 1 && B[m-1] == 0)
	{
		if(optionNormalize == 1)
		{
			///
			//cout<<"monic is chosen"<<endl;
			///
			numBlocksN = (sfixn)(ceil((double)( n)/(double)(T)));
			///
			//cout<<"numBlocksN :"<<numBlocksN<<endl;
			///
			MonicConversionA<<<numBlocksN, T>>>(Agpu, stateGpu);
			///
			//cout<<"monic is done"<<endl;
			///			
			cudaMemcpy(G, Agpu, sizeof(sfixn)*n, cudaMemcpyDeviceToHost);
			G[n-1] = 1;		
		}
		else
			cudaMemcpy(G, Agpu, sizeof(sfixn)*n, cudaMemcpyDeviceToHost);

		cudaEventRecord(stop, 0);
	        cudaEventSynchronize(stop);        
	        cudaEventElapsedTime(&outerTime, start, stop);

		cudaFree(Agpu);
	        cudaFree(Bgpu);
	        cudaFree(stateGpu);
	
		(*Dgcd) = n-1;
		return outerTime;		
	}
	if(m == 1)
	{
		if(optionNormalize == 1) G[0] = 1;
		else G[0] = B[0];
		(*Dgcd) = 0;
		cudaFree(Agpu);
	        cudaFree(Bgpu);
	        cudaFree(stateGpu);
		return 0.0;
	}




	if(m <= n) 	j = m + T;
	else j = n + T;		
	numBlocksN = (sfixn)(ceil((double)( j)/(double)(T*2)));
	///
	//cout<<numBlocksN<<" "<<T<<endl;
	///
	
	for(i = n+m; i>= 0; i -= T )
	{
		///
		//cout<<i<<endl;
		///
		reduceNgcd<<<1, 1>>>(Agpu, stateGpu);
		reduceMgcd<<<1, 1>>>(Bgpu, stateGpu);

		status4<<<1, 1>>>(stateGpu);		
		status3<<<1, 1>>>(stateGpu);		
		status5<<<1, 1>>>(stateGpu);
		
		gcdGPU<<<numBlocksN, T>>>(Agpu, Bgpu, stateGpu);		
	}


	cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        
        cudaEventElapsedTime(&outerTime, start, stop);

	shrinkA<<<1, 1>>>(Agpu, stateGpu);
	shrinkB<<<1, 1>>>(Bgpu, stateGpu);
	
        // We need to copy the status vector on the host side
        // in order to know where the GCD will be stored on the
        // CPU side. 
	cudaMemcpy(state, stateGpu,  sizeof(sfixn)*6, cudaMemcpyDeviceToHost);
	
	////
	//cout<<state[0]<<" "<<state[1]<<" "<<state[2]<<" "<<state[3]<<" "<<state[4]<<" "<<state[5]<<" "<<endl;
	//cudaMemcpy(A, Agpu, sizeof(sfixn)*n, cudaMemcpyDeviceToHost);	
        //cudaMemcpy(B, Bgpu, sizeof(sfixn)*m, cudaMemcpyDeviceToHost);
	//for(i = 0; i < n; ++i) cout<<A[i]<<" ";
	//cout<<endl;  
	//for(i = 0; i < m; ++i) cout<<B[i]<<" ";
	//cout<<endl;  

	////


	j = -1;

        // Recall that cudaMemcpy is a CPU function. Thus we need to read
        // from state (not from stateGpu) where to take the GCD from
	if(state[0] >= 0)
	{
		j = state[0]+1;
			///
			//cout<<"in GPU: "<<j<<endl;
			///		
		if(optionNormalize == 1)
		{
			///
			//cout<<"monic is choosen"<<endl;
			///
			if(j>1)
			{				
				numBlocksN = (sfixn)(ceil((double)( j)/(double)(T)));
				MonicConversionA<<<numBlocksN, T>>>(Agpu, stateGpu);				
				cudaMemcpy(G, Agpu, sizeof(sfixn)*j, cudaMemcpyDeviceToHost);
			}
			G[j-1] = 1;
		}		
		else 
			cudaMemcpy(G, Agpu, sizeof(sfixn)*j, cudaMemcpyDeviceToHost);
	}
	else if(state[1] >= 0)
	{
		j = state[1] +1;
			///
			//cout<<"in GPU: "<<j<<endl;
			///

		if(optionNormalize == 1)
		{
			///
			//cout<<"monic is choosen"<<endl;
			///
			if(j>1)
			{
				numBlocksN = (sfixn)(ceil((double)( j)/(double)(T)));
				MonicConversionB<<<numBlocksN, T>>>(Bgpu, stateGpu);			
				cudaMemcpy(G, Bgpu, sizeof(sfixn)*j, cudaMemcpyDeviceToHost);
			}
			G[j-1] = 1;
		}		
		else 
			cudaMemcpy(G, Bgpu, sizeof(sfixn)*j, cudaMemcpyDeviceToHost);
	}
	
	cudaFree(Agpu);
        cudaFree(Bgpu);
        cudaFree(stateGpu);

	(*Dgcd) = j-1;
	return outerTime;
}


