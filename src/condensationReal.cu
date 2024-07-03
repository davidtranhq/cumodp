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



/*! \file condensationReal.cu
\brief This file contains the code for computing determinant over finite field using CUDA.
*/

#include "condensationReal.h"



using namespace std;

/*! \fn long double mulSeq(int n, long double *A)
\brief This function computes the product of a list of real numbers.

\param n an integer that stores the size of the list.
\param A an integer array that stores the list.
*/		

long double mulSeq(int n, long double *A)
{
	long double result = 1.0, swapping;
		

	int i, j, k,  flag;
	/*
	cout<<endl<<endl<<endl;
	for(i = 0; i < n; ++i)
		cout<<A[i]<<" ";
	cout<<endl;
	*/
	
	i = 0;	
	k = n-1;
	j = n;
	flag = 0;
	
	while(flag ==0)
	{
		flag = 1;

		while( (i < j)  && (A[i] >= 1.0 || A[i] <= -1.0))
			++i;
		
		//
		//cout<<i<<" ";
		//
		while( (k > -1) && (A[k] < 1.0 && A[k] > -1.0))
			--k;
		
		//
		//cout<<k<<" ";
		//
		j = k;
		if(i>k)
			break;
		else if(i < j && k > -1 &&(A[i] < 1.0 && A[i] > -1.0 ) && (A[k] >= 1.0 || A[k] <= -1.0) )		
		{
			swapping = A[i];
			A[i] = A[k];
			A[k] =swapping;		
			flag = 0;
		}
		//
		//cout<<j<<endl;
		//
	}
	/*
	cout<<j<<endl;
	for(i = 0; i < n; ++i)
		cout<<A[i]<<" ";
	cout<<endl;
	*/
	
	while( j > -1 && j < n-1)
	{
		//		
		//cout<<endl<<j<<" "<<n;
		//
		A[j] = A[j]*A[n-1];
		--n;
		if(A[j] < 1.0 && A[j] > -1.0)
			--j;
	}
	/*
	cout<<endl;
	for(i = 0; i < n; ++i)
		cout<<A[i]<<" ";
	cout<<endl;
	*/
	for(i = 0; i < n; ++i)
	{
		result*=A[i];
		//cout<<result<<" ";
	}
	
	return result;
}


/*! \fn void condensation(realV *Aa, realV *Bb, int n, int l, int BKReq, int BK,  int BKID, int THid)
\brief This function simulates the execution of condensation method of GPU on CPU. This is no longer in use.

\param Aa an real valued array that stores the matrix before one condensation method.
\param Bb an real valued array that stores the matrix after one condensation method.
\param n the dimension of the intermediate matrix.
\param l the index of pivot.
\param BKReq the required number of thread block.
\param BK the actual tile size for a thread.
\param BKID the block id.
\param THid the thread id.
*/		


void condensation(realV *Aa, realV *Bb, int n, int l, int BKReq, int BK,  int BKID, int THid)
{

	int k = (BKID *THno + THid)*BK;
	
	
	if(k<n*n)
	{
		//
		//cout<<k<<" ";
		//
		realV firstRow;
		int p = k%n;
		int q = k/n;
				
				
		int i, j, s, d, start, limitA, limitB;				
		limitB = n*n;
		limitA = limitB + 2*n + 1;
		for(s=1;s<BK;++s)
		{
			if((k+s)%n ==0)
				break;
		}
		//
		//cout<<l<<" "<<p<<" "<<q<<" "<<s<<" ";
		//
		///*
		if(q >= l)
		{
			firstRow = Aa[(n+1)*(q+1) ];
			//
			//cout<<firstRow<<" ";
			//		
			i = k+q+2+n;
                        start = l*(n+1)+p+1;
			
			for(d=0; d < s && (d+i) < limitA && (k+d) < limitB;++d)
			{
				Bb[k+d]= Aa[d+i] - firstRow*Aa[start+d];
				//
				//cout<<"("<<k+d<<" "<<d+i<<" "<<start+d<<")";
				//
			}
                }
		else
		{
			firstRow = Aa[(n+1)*(q) ];
			//
			//cout<<firstRow<<" ";
			//	
			i = k+q+1;
                        start = l*(n+1)+p+1;
                        for(d=0; d < s && (d+i) < limitA && (k+d) < limitB;  ++d)
			{
				Bb[k+d] = firstRow*Aa[start+d] - Aa[i+d];
				//
				//cout<<"("<<k+d<<" "<<d+i<<" "<<start+d<<")";
				//
			}
                }			 
                
		q = q+1;
		if(s < BK && q < n)
		{						
			if(q>=l)
                 	{
				
				firstRow = Aa[(n+1)*(q+1) ];
				//
				//cout<<firstRow<<" ";
				//			
				k = n*q;				
				i = (n+1)*l + 1;
				d = (n+1)*(q+1) + 1;
				 
				for(j = s, start=0; j < BK; ++j, ++start)
				{
					Bb[k + start] = Aa[d + start] - firstRow*Aa[i+start];		
					//
					//cout<<"("<<k+start<<" "<<d + start<<" "<<i+start<<")";
					//
				}
                	}
                	else
                	{
				firstRow = Aa[(n+1)*(q) ];
				//
				//cout<<firstRow<<" ";
				//
				k = n*q;
				d = (n+1)*(q) + 1;
				i = l*(n+1)+1;;
				for(j = s, start=0; j < BK; ++j, ++start)
				{
					Bb[k + start] =  firstRow*Aa[i+start] - Aa[d + start];	
					//
					//cout<<"("<<k+start<<" "<<d + start<<" "<<i+start<<")";
					//		
				}
                	}			
		}
		//
		//cout<<endl;
		//
		//*/			
	}
}



/*! \fn int det(realV *A, int n, int BS, int BK)
\brief This function is the top level function of condensation(). This is no longer in use.


\param A an integer array that store the matrix before any condensation step.
\param n the dimension of input matrix.
\param P a prime.
\param BS an integer that indicates the dimension of the matrix when condensation method should stop working.
Currently this option is not used to keep it independent. Previousely, this function depends on NTL library for when 
the code reached here.
\param BK an integer that fix the tile size for each thread.
*/

void det(realV *A, int n, int BS, int BK)
{
 	realV *temp, *B, *lvalues, closest=1.0;
	int n1=n;
	int l;
	int i, j, flag = -1;
	int BKReq;
	B = ( realV *)malloc(sizeof( realV)*(n-1)*(n-1));
	lvalues = ( realV *)malloc(sizeof( realV)*(n-1));		
	while(n>=3)
	{
		flag = -1;
		for(i = 0; i < n ; ++i)
		{
			if(A[i*n] != 0.0)
			{
				if(flag == -1)
				{
					closest = abs(A[i*n]-1.0);
					flag = i;
				}
				else if(abs(A[i*n]-1.0) < closest)
				{
					closest = abs(A[i*n]-1.0);
					flag = i;			
				}		
			}
		}
		lvalues[n-2] = A[n*flag];
		//
		//cout<<n<<" "<<flag<<" "<<A[n*flag]<<endl;
		//

		for(i = 1; i < n ; ++i)
			A[flag*n + i] = A[flag*n + i]/A[n*flag];
		l = flag;
		/*
		for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
                	{
                	        cout<<A[j*n+i]<<" ";
                	}
			cout<<endl;
        	}
		cout<<endl;
		 */
		
		--n;
		if(n<BK) BK = n;
		BKReq = ceil((double)(n*n)/(double) (THno*BK));
		//
		//cout<<BKReq<<endl;
		//
		for(i= 0; i < BKReq ; ++i)
		{	
			for(j=0; j<THno; ++j)
			{

				condensation(A, B, n, l,  BKReq, BK, i, j);
			}
		}
		
		//if(n>=BS){
		/*
		for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
                	{
                	        cout<<B[j*n+i]<<" ";
                	}
			cout<<endl;
        	}
		cout<<endl;
		cout<<endl;
		cout<<endl;
		//}
		*/
		
		temp = A;
		A = B;
		B = temp;	

	}
	lvalues[0] = A[0]*A[3] - A[1]*A[2];
	//cout<<closest<<" "<<endl;;
	
	/*
	for(i = 0; i < n1-2; ++i)
	{
		cout<<lvalues[i]<< " ";
		//closest = closest*lvalues[i];
	}
	*/
	//closest *= 
	//cout<<endl<<endl<<endl<<mulSeq(n1-1,lvalues)<<endl;

	long double *multi = (long double *)malloc(sizeof(long double)*(n1-1));
	for(i = 0; i<n1-1; ++i)
		multi[i] = lvalues[i];
	//long double aaa = mulSeq(n1-1,multi);
	cout<<endl<<"determinant from simulation GPU  on CPU method: "<<mulSeq(n1-1,multi)<<endl;
	free(multi);		


	//cout<<endl<<"determinant from simulation method: "<<mulSeq(n1-1,lvalues)<<endl;
	free(B);
	free(A)	;      
	
}



/*! \fn void direct(realV *A, int n)
\brief This is a CPU method of the implementation of condensation method.

\param A an integer array that store the matrix before any condensation step.
\param n the dimension of input matrix.
*/

void direct(realV *A, int n)
{
        realV *temp, *B, *lvalues, closest=1.0;
	int n1=n;
	int l;
	int i, j, flag = -1;
	B = ( realV *)malloc(sizeof( realV)*(n-1)*(n-1));
	lvalues = ( realV *)malloc(sizeof( realV)*(n-1));		
	while(n>=3)
	{
		flag = -1;
		for(i = 0; i < n ; ++i)
		{
			if(A[i*n] != 0.0)
			{
				if(flag == -1)
				{
					closest = abs(A[i*n]-1.0);
					flag = i;
				}
				else if(abs(A[i*n]-1.0) < closest)
				{
					closest = abs(A[i*n]-1.0);
					flag = i;			
				}		
			}
		}
		lvalues[n-2] = A[n*flag];
		//
		//cout<<n<<" "<<flag<<" "<<A[n*flag]<<endl;
		//

		for(i = 1; i < n ; ++i)
			A[flag*n + i] = A[flag*n + i]/A[n*flag];
		l = flag;
		/*
		for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
                	{
                	        cout<<A[j*n+i]<<" ";
                	}
			cout<<endl;
        	}
		cout<<endl;
		*/
		for(j= 0; j < l; ++j)
		{
			for(i = 0; i < (n-1); ++i)
			{
				//cout<<A[j*n]<<" "<<A[l*n + i +1]<<" "<<A[j*n+i+1]<<endl;
				B[j*(n-1)+i] = A[j*n]*A[l*n + i +1] - A[j*n+i+1];
			}
		}

		for(j= l; j < (n-1); ++j)
		{
			for(i = 0; i < (n-1); ++i)
			{
				B[j*(n-1)+i] =  A[n*(j+1) +i+1]- A[n*(j+1)]*A[l*n+i+1];
			}
		}
		--n;
		
		/*if(n>=BS){
		
		for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
                	{
                	        cout<<B[j*n+i]<<" ";
                	}
			cout<<endl;
        	}
		cout<<endl;
		cout<<endl;
		cout<<endl;
		//}*/
		
		
		temp = A;
		A = B;
		B = temp;	

	}
	lvalues[0] = A[0]*A[3] - A[1]*A[2];
	//cout<<closest<<" "<<endl;;
	
	/*
	for(i = 0; i < n1-2; ++i)
	{
		cout<<lvalues[i]<< " ";
		//closest = closest*lvalues[i];
	}
	*/
	//closest *= 
	//cout<<endl<<endl<<endl<<mulSeq(n1-1,lvalues)<<endl;

	long double *multi = (long double *)malloc(sizeof(long double)*(n1-1));
	for(i = 0; i<n1-1; ++i)
		multi[i] = lvalues[i];
	detDirect = mulSeq(n1-1,multi);
	//cout<<endl<<"determinant from direct method: "<<mulSeq(n1-1,multi)<<endl;
	free(multi);	



	//cout<<endl<<"determinant from direct method: "<<mulSeq(n1-1,lvalues)<<endl;
	free(B);
	free(A)	;
	
	//return aaa;
}

/*! \fn __global__ void find(realV *A, int *ls, realV *lvalues, int n, int temp)
\brief This function computes the pivot element and store the  pivot elements.

\param A an integer array that store the matrix before one condensation step.
\param ls an integer array that store the matrix after one condensation step.
\param lvalues the product of pivot elements
\param n the dimension of intermediate matrix.
\param temp an integer for intermediate computation.
*/


__global__ void find(realV *A, int *ls, realV *lvalues, int n, int temp)
{	
	if(ls[0]>=0)
	{

		realV closest;
        	__shared__ int min[THno];
		int flag;
        	int i = threadIdx.x;
        	min[i] = -1;
		flag = 0;

        	int j;
		int k = (i+1)*temp;        	

        	for(j = i*temp; j< n && j< k; ++j)
        	{
        	        if(A[j*n]!=0.0) 
        	        {
				if(flag == 0)
				{
					closest = abs(A[j*n]-1.0);  
					min[i] = j;
					flag = 1;
				}
				else if(abs(A[j*n]-1.0) < closest) 
				{
        	                	min[i] = j;
					closest = abs(A[j*n]-1.0);        		                
				}
        	        }
        	}
        	__syncthreads();
        	if(i==0)
        	{
        	        ls[0] = -1;
			for(j=0;j<THno;++j)
        	        {
				if(min[j] >=0 )
				{
					if(ls[0] == -1)
					{
						ls[0]=min[j];
						lvalues[n-2] =  abs(A[ls[0]*n]-1.0);
					}
					else if( abs(A[min[j]*n]-1.0) < lvalues[ n-2] )
					{
						ls[0]=min[j];
						lvalues[n-2] =  abs(A[ls[0]*n]-1.0);
					}
				}
        	        }
        	       if(ls[0] >= 0)
				lvalues[n-2] =  A[n*ls[0]];				
        	}
	}
}

/*! \fn __global__ void divideToOne(realV *A,  int *B,  realV *C, int n, int temp)
\brief This function makes the pivot element 1 and does the necessary adjustment.

\param A an integer array that store the matrix before one condensation step.
\param B an integer that works as an flag to indicate whether the computation needs to stop or not before n condensation steps.
\param C the list pivot elements
\param n the dimension of intermediate matrix.
\param temp an integer for intermediate computation.

*/
__global__ void divideToOne(realV *A,  int *B,  realV *C, int n, int temp)
{	
	int i = threadIdx.x *temp;
        int j;
	int k = (threadIdx.x+1)*temp;
	if(B[0]>=0 && i < n)
	{
		realV lValues = C[n-2];
		int ls= B[0]*n;
        	for(j = i; j< n && j< k; ++j)
        		A[ ls + j]  = A[ ls + j]/lValues;
	}
}

/*! \fn __global__ void condensGPU(realV *Aa, realV *Bb, int n, int *l, int BK)
\brief This function computes one step of condensation.

\param Aa an integer array that store the matrix before one condensation step.
\param Bb an integer array that store the matrix after one condensation step.
\param n the dimension of intermediate matrix.
\param l the list of pivot elements
\param BK an integer that stores the size of tiles that a thread works.
*/
__global__ void condensGPU(realV *Aa, realV *Bb, int n, int *l, int BK)
{

	int k = (blockIdx.x*blockDim.x + threadIdx.x)*BK;
	
	
	if(l[0]>=0&&k<n*n)
	{
		//
		//cout<<k<<" ";
		//
		realV firstRow;
		int p = k%n;
		int q = k/n;
				
				
		int i, j, s, d, start, limitA, limitB;				
		limitB = n*n;
		limitA = limitB + 2*n + 1;
		for(s=1;s<BK;++s)
		{
			if((k+s)%n ==0)
				break;
		}
		//
		//cout<<l<<" "<<p<<" "<<q<<" "<<s<<" ";
		//
		///*
		if(q >= l[0])
		{
			firstRow = Aa[(n+1)*(q+1) ];
			//
			//cout<<firstRow<<" ";
			//		
			i = k+q+2+n;
                        start = l[0]*(n+1)+p+1;
			
			for(d=0; d < s && (d+i) < limitA && (k+d) < limitB;++d)
			{
				Bb[k+d]= Aa[d+i] - firstRow*Aa[start+d];
				//
				//cout<<"("<<k+d<<" "<<d+i<<" "<<start+d<<")";
				//
			}
                }
		else
		{
			firstRow = Aa[(n+1)*(q) ];
			//
			//cout<<firstRow<<" ";
			//	
			i = k+q+1;
                        start = l[0]*(n+1)+p+1;
                        for(d=0; d < s && (d+i) < limitA && (k+d) < limitB;  ++d)
			{
				Bb[k+d] = firstRow*Aa[start+d] - Aa[i+d];
				//
				//cout<<"("<<k+d<<" "<<d+i<<" "<<start+d<<")";
				//
			}
                }			 
                
		q = q+1;
		if(s < BK && q < n)
		{						
			if(q>=l[0])
                 	{
				
				firstRow = Aa[(n+1)*(q+1) ];
				//
				//cout<<firstRow<<" ";
				//			
				k = n*q;				
				i = (n+1)*l[0] + 1;
				d = (n+1)*(q+1) + 1;
				 
				for(j = s, start=0; j < BK; ++j, ++start)
				{
					Bb[k + start] = Aa[d + start] - firstRow*Aa[i+start];		
					//
					//cout<<"("<<k+start<<" "<<d + start<<" "<<i+start<<")";
					//
				}
                	}
                	else
                	{
				firstRow = Aa[(n+1)*(q) ];
				//
				//cout<<firstRow<<" ";
				//
				k = n*q;
				d = (n+1)*(q) + 1;
				i = l[0]*(n+1)+1;;
				for(j = s, start=0; j < BK; ++j, ++start)
				{
					Bb[k + start] =  firstRow*Aa[i+start] - Aa[d + start];	
					//
					//cout<<"("<<k+start<<" "<<d + start<<" "<<i+start<<")";
					//		
				}
                	}			
		}
		//
		//cout<<endl;
		//
		//*/			
	}
}


/*! \fn void deter(realV *A, int n, int BS, int BK)
<a name="detRealCPU"> </a>
\brief This is the CPU function to be called by external functions (i.e. C code).
It calls the kernels for conputing determinant
for a real valued matrix.

\param A an integer array that store the input matrix.
\param n the dimension of input matrix.
\param BS an integer that indicates the dimension of the matrix when condensation method should stop working.
Currently this option is not used to keep it independent. Previousely, this function depends on NTL library for when 
the code reached here.
\param BK an integer that fix the tile size for each thread.
*/
void deter(realV *A, int n, int BS, int BK)
{
	int n1,  *ls, temp;
	realV *Aa, *Bb, *BB, *swap, *lvalues;
	n1 =n;
		//
		int i; //j;
		//
	cudaMalloc((void **)&Aa, sizeof(realV)*n*n);
        cudaMemcpy(Aa, A, sizeof(realV)*n*n, cudaMemcpyHostToDevice);

	cudaMalloc((void **)&Bb, sizeof(realV)*(n-1)*(n-1));
	cudaMalloc((void **)&lvalues, sizeof(realV)*(n-1));

	cudaMalloc((void **)&ls, sizeof(int));	
	int lsCPU[1];
	lsCPU[0] = 0;
	cudaMemcpy(ls, lsCPU, sizeof(int), cudaMemcpyHostToDevice);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

	while(n>=3)
	{
		

		/*
		BB = (realV *)malloc(sizeof(realV)*n*n);
        	cudaMemcpy(BB, Aa, sizeof(realV)*n*n, cudaMemcpyDeviceToHost);
 		int j;
        	for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
               		{
                		cout<<BB[(j)*n+i]<<" ";  
                	}
			cout<<endl;
               }
	       cout<< endl;
	       free(BB);
		*/
		temp = (int) ceil((double)n/(double)THno);

		find<<<1,THno>>>(Aa, ls, lvalues,  n, temp);
		/*
		BB = (realV *)malloc(sizeof(realV)*(n1-1));
        	cudaMemcpy(BB, lvalues, sizeof(realV)*(n1-1), cudaMemcpyDeviceToHost);
		for(i=0;i<n1-1;++i)
			cout<<BB[i]<<" ";
		cout<<endl;
		free(BB);
		cudaMemcpy(lsCPU, ls, sizeof(int), cudaMemcpyDeviceToHost);
		cout<<lsCPU[0]<<endl;
		*/
		divideToOne<<<1,THno>>>(Aa, ls, lvalues,  n, temp);
		
		--n;
		if(n<BK)
			BK = n;
		temp = (int)ceil((double)(n*n)/(double) (THno*BK));

		condensGPU<<< temp,THno>>>(Aa, Bb, n, ls, BK);
		swap = Aa;
		Aa= Bb;
		Bb = swap;

		/*
		BB = (realV *)malloc(sizeof(realV)*n*n);
        	cudaMemcpy(BB, Aa, sizeof(realV)*n*n, cudaMemcpyDeviceToHost);
 		
        	for(i=0;i<n;++i)
        	{
                	for(j=0;j<n;++j)
               		{
                		cout<<BB[(j)*n+i]<<" ";  
                	}
			cout<<endl;
               }
	       cout<< endl;
	       free(BB);
		*/
	}
	swap = (realV *)malloc(sizeof(realV)*(n*n));
	cudaMemcpy(swap, Aa, sizeof(realV)*(n*n), cudaMemcpyDeviceToHost);
	BB = (realV *)malloc(sizeof(realV)*(n1-1));
	cudaMemcpy(BB, lvalues, sizeof(realV)*(n1-1), cudaMemcpyDeviceToHost);
	BB[0] = swap[0]*swap[3] - swap[1]*swap[2];
	
	/*
	for(i = 0; i<n1-1; ++i)
		cout<<BB[i]<<" ";
	cout<<endl;
	*/
	long double *multi = (long double *)malloc(sizeof(long double)*(n1-1));
	for(i = 0; i<n1-1; ++i)
		multi[i] = BB[i];
	detCon = mulSeq(n1-1,multi);
	//cout<<endl<<"determinant from GPU method: "<<mulSeq(n1-1,multi)<<endl;

	cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float outerTime;
        cudaEventElapsedTime(&outerTime, start, stop);
	timeDet = outerTime/1000.0;
        //printf("%f \t",outerTime/1000.0);
	//printf(" s\n");
	
	free(multi);

	free(swap);
	
	free(BB);
	free(A);
	cudaFree(Aa);
	cudaFree(Bb);
	cudaFree(lvalues);
	cudaFree(ls);

	
	///return aaa;

}

