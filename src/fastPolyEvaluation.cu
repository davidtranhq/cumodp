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



/*! \file fastPolyEvaluation.cu
\brief This file contains the functions to evaluate an univariate polynomial for n points over finite field. 
*/


#include "fastPolyEvaluation.h"

struct PolyEvalSteps check; // this will hold all the intermediate results.

/*! \fn __global__ void copyMgpu(sfixn *dest, sfixn *source, sfixn length_poly)
\brief This copy list of polynomials in source to destination.
Each of the polynomial in source is of length_poly.
In dest all coefficients are copied except the most significant one.
\param dest a list of polynomial.
\param source a list of polynomial.
\param length_poly length of each poly  in source.
*/

__global__ void copyMgpu(sfixn *dest, sfixn *source, sfixn length_poly)
{
	dest[blockIdx.x *(length_poly-1) + threadIdx.x] = source[blockIdx.x *length_poly + threadIdx.x];	
}


/*! \fn __global__ void allZero( sfixn *X, sfixn n)
\brief It simply make X a zero vector
\param X an array.
\param n length of X.
*/
__global__ void allZero( sfixn *X, sfixn n)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n) X[tid] = 0;
}

/*! \fn __global__ void allNeg( sfixn *X, sfixn n, sfixn p)
\brief It simply negate X.
\param X an array.
\param n length of X.
\param p the prime for finite field.
*/
__global__ void allNeg( sfixn *X, sfixn n, sfixn p)
{
        sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
        if(tid < n) X[tid] = neg_mod(X[tid],p);
}

/*! \fn __global__ void pointAdd2(sfixn *dest, sfixn *source, sfixn l,  sfixn n, sfixn p)
\brief This function adds two coefficients of source and stores in dest.
\param dest list of polynomials.
\param source  list of polynomials.
\param l length of polynomial in source.  
\param n the length of dest or source is 2*n
\param p the prime for finite field.

*/
__global__ void pointAdd2(sfixn *dest, sfixn *source, sfixn l,  sfixn n, sfixn p)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
	{
		sfixn a = tid/l;
		sfixn b = tid%l;
		dest[l*(2*a+1) +b ] = add_mod(source[l*(2*a+1) +b], source[l*2*a + b], p);	
	}		
}

/*! \fn __global__ void  zeroInbetween(sfixn *X, sfixn *Y, sfixn n, sfixn l)
\brief This function puts zero in between the list of polynomials in X and Y.
\param X a list of polynomials.
\param Y a list of polynomials.
\param n total lenght of X or Y.
\param l length of each polynomial in X and Y.
*/

__global__ void  zeroInbetween(sfixn *X, sfixn *Y, sfixn n, sfixn l)
{
	sfixn tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid < n)
	{
		sfixn a = tid/l;
		sfixn b = tid%l;
		X[l*(2*a+1) + b ] = 0;
		Y[l*(2*a+1) + b ] = 0;
	}
}

/*! \fn __global__ void pointMul(sfixn *dest, sfixn *source, sfixn n, sfixn p)
\brief This function does pointwise multiplication between source and dest array and dest will store the point wise multiplications.
\param dest a list of polynomials.
\param source a list of polynomials.
\param n total lenght of source or destination.
\param p prime for prime field.
*/
__global__ void pointMul(sfixn *dest, sfixn *source, sfixn n, sfixn p)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
		dest[tid] = mul_mod(source[tid],dest[tid], p);	
}

/*! \fn __global__ void scalarMul(sfixn *A, sfixn ninv, sfixn L, sfixn p)
\brief  This function multiplies each entry of A with ninv over finite field.
\param A a list of polynomials.
\param ninv a value over prime field.
\param L total length of A.
\param p prime for prime field.
*/

__global__ void scalarMul(sfixn *A, sfixn ninv, sfixn L, sfixn p)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < L)
		A[tid] = mul_mod(A[tid], ninv, p);	
}


/*! \fn __global__ void pointAdd(sfixn *dest, sfixn *source, sfixn n, sfixn p)
\brief This function does pointwise addition between source and dest array and dest will store the point wise multiplications.
\param dest a list of polynomials.
\param source a list of polynomials.
\param n total length of dest or source.
\param p prime for prime field.
*/


__global__ void pointAdd(sfixn *dest, sfixn *source, sfixn n, sfixn p)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
		dest[tid] = add_mod(source[tid],dest[tid], p);	
}

/*! \fn __global__ void listPlainMulGpu( sfixn *Mgpu1, sfixn *Mgpu2 , sfixn length_poly, sfixn poly_on_layer, sfixn threadsForAmul, sfixn mulInThreadBlock, sfixn p)
\brief This function does a list of plain multiplications.

'sfixn *Mgpu1' is a list of 'poly_on_layer' polynomials. 
'poly_on_layer' is an even number.
Each polynomials has length of 'length_poly'.

'sfixn *Mgpu2' stores the list multiplications of consecutive pairs
i.e multiplication of poly i and poly i+1 in Mgpu1 is stored in poly i/2_{th}
slot of Mgpu2.

The length of a poly in Mgpu2 is '2*length_poly - 1'.
We dedicate 2*length_poly threads for a poly multiplication.
so the number of threads responsible for one poly multiplication
'threadsForAmul = 2*length_poly'.

The multiplication is done using shared memory.
The same shared memory can not be accesses by
two different thread blocks. So a thread block
is responsible for doing a number of poly multiplication.
But the reverse is not true. That means two or more 
thread block can not do a single poly multiplication together.
Thats why the number of multiplications done by one thread block
'mulInThreadBlock = (sfixn)floor((double)Tmul/(double)threadsForAmul)'
Where 'Tmul' is the number of threads in a thread block. We keep it as 512.

These constraint us to the limitation of our poly multiplications.
First 2*length_poly <= Tmul.
poly_on_layer/2 <= maximum number of thread block.

\param Mgpu1 list of polynomials
\param Mgpu2 list of polynomials after multiplications.
\param length_poly length of each poly in  Mgpu1.
\param poly_on_layer no of poly on the Mgpu1.
\param threadsForAmul no of threads required for one multiplication.
\param mulInThreadBlock the no of multiplications in a thread block. It must be at least 1.
\param p prime for prime field.
*/
__global__ void listPlainMulGpu( sfixn *Mgpu1, sfixn *Mgpu2 , sfixn length_poly, sfixn poly_on_layer, sfixn threadsForAmul, sfixn mulInThreadBlock, sfixn p)
{

	__shared__ sfixn sM[2*Tmul];
	/*
	sM is the shared memory where the all the coefficients and intermediate multiplications results
	are stored. For each multiplication it reserve 4*length_poly -1 spaces.
	mulID is the multiplication ID. It refers to the poly in Mgpu2 on which it will work.
	mulID must be less than (poly_on_layer/2).
	*/	
	sfixn mulID= ((threadIdx.x/threadsForAmul) + blockIdx.x*mulInThreadBlock);
	
	if( mulID < (poly_on_layer/2) && threadIdx.x < threadsForAmul*mulInThreadBlock)
	{
		/*
		The next 10 lines of code copy the polynomials in Mgpu1 from global memory to shared memory.
		Each thread is responsible of copying one coefficient.
		A thread will copy a coefficient from Mgpu1[( mulID* length_poly*2)...( mulID* length_poly*2) + length_poly*2 -1] 
		j+u gives the right index of the coefficient in Mgpu1.

		In sM, the coefficients are stored at the lower part.
		t will find the right (4*length_poly-1) spaced slot for it.
		s gives the start index of its right slot.
		s+u gives right position for the index.

		
		*/

		sfixn j = ( mulID* length_poly*2);		
		sfixn q = ( mulID*(2*length_poly-1));

		sfixn t = (threadIdx.x/threadsForAmul);
		sfixn u = threadIdx.x % threadsForAmul;

		sfixn s = t*(4*length_poly-1);
		sfixn k = s + length_poly;
		sfixn l = k + length_poly;
		sfixn c = l+u;
		sfixn a, b, i;

		sM[s+u] = Mgpu1[j + u];
		__syncthreads();
		
		if(u != (2*length_poly-1) )
		{
			/*
			In the multiplication space, the half of the leading coefficients 
			are computed differently than the last half. Here the computation of 
			first half are shown. the last half is shown in else statement.
			In both cases sM[c] is the cofficient on which this thread will work on.
			sM[a] is the coefficient of one poly.
			sM[b] is the coefficient of the other poly.
			*/
			if(u < length_poly)
			{
				a = s;
				b = k + u;			
				sM[c] =  mul_mod(sM[a],sM[b],p);
				++a; --b;
				for(i = 0; i < u; ++i, ++a, --b)
					sM[c] =  add_mod(mul_mod(sM[a],sM[b],p),sM[c] ,p);
				

				Mgpu2[q+u] = sM[c];				
			}
			else
			{
				b = l - 1;
				a = (u - length_poly) + 1 + s;
				sM[c] =  mul_mod(sM[a],sM[b],p);
				++a; --b;
				
				sfixn tempU = u;
				u = (2*length_poly-2) - u;
				for(i = 0; i < u; ++i, ++a, --b)
					sM[c] =  add_mod(mul_mod(sM[a],sM[b],p),sM[c] ,p);
				

				Mgpu2[q+tempU] = sM[c];			
			}	
		}
	}		
}

/*! \fn __global__ void listPolyinv( sfixn *Mgpu, sfixn *invMgpu,  sfixn poly_on_layer, sfixn prime)
\brief This kernel  computes the
inverse of a list of poly.
Each thread is responsible for 
computing the inverse of one polynomial.
The length of each polynomial is 257.
It first copy the polynomial into the 
shared memory. Then it computes the 
inverse of it (mod x^256) using Newton 
iteration using lain arithmetic.
The number of threads Tinv
in a thread block is small, 
as we a need enough spaces for 
storing all polynomials and intermediate 
results into shared memory.
We fix Tinv = 16.
Each thread is working 
with a polynomial.
For details please read Chapter 10 of the 
PhD thesis of Sardar Haque (Computer Science, UWO, 2013).
\param Mgpu list of polynomials.
\param invMgpu  list of computed inverse polynomials.
\param poly_on_layer no of poly in Mgpu.
\param prime prime for prime field.

*/
__global__ void listPolyinv( sfixn *Mgpu, sfixn *invMgpu,  sfixn poly_on_layer, sfixn prime)
{
	sfixn PolyId = blockIdx.x*blockDim.x + threadIdx.x;	
	
	if(PolyId < poly_on_layer)
	{
		const sfixn  SIZE = 256;
		const sfixn r = 8;
		sfixn i, j;
		sfixn f[SIZE+1], g[SIZE], gtemp[SIZE], gtemp2[SIZE];


		j = (PolyId+1)*(SIZE+1) -1;		
		for(i = 0; i <= SIZE; ++i, --j)
			f[i] = Mgpu[j];

		g[0] = 1;
		gtemp[0] = 1;
		gtemp2[0] = 1;
		for(i = 1; i < SIZE; ++i)		
		{
			g[i] = 0;
			gtemp[i] = 0;
			gtemp2[i] = 0;
		}
		sfixn  p, q, start = 1, end = 2;	
		for(i = 1; i <= r; ++i)
		{
			for(j = start; j < end; ++j)
			{
				for(p = 0, q = j; p < SIZE && q >= 0; ++p, --q)
					gtemp[j] = add_mod(mul_mod(g[p], f[q], prime), gtemp[j] ,prime); 
				gtemp[j] = neg_mod(gtemp[j], prime);	
			} 

			for(j = start; j < end; ++j)
				for(p = 0, q = j; p < SIZE && q >= start; ++p, --q)
					g[j] = add_mod(mul_mod(gtemp[q],gtemp2[p], prime),g[j], prime); 			


			for(j = start; j < end; ++j)
				gtemp2[j] = g[j];
			
			start = start*2;
			end = end*2;
		}		
		for(j = 0, i = PolyId*(SIZE); j < SIZE; ++i, ++j)
			invMgpu[i] = g[j];
	}
}



/*! \fn __global__ void listReversePoly(sfixn *revMgpu, sfixn *Mgpu, sfixn length_poly, sfixn poly_on_layer)
\brief This function reverse the list of polynomials in Mgpu to revMgpu.
\param revMgpu list of computed reverse polynomials. 
\param Mgpu list of polynomials.
\param length_poly length of each polynomial in Mgpu.
\param poly_on_layer no of polynomials in Mgpu.
*/
__global__ void listReversePoly(sfixn *revMgpu, sfixn *Mgpu, sfixn length_poly, sfixn poly_on_layer)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
	{
		sfixn polyID = tid/length_poly;
		sfixn offset= tid%length_poly;
		revMgpu[(polyID+1)*length_poly - offset -1  ] = Mgpu[polyID*length_poly + offset];
	}
}

/*! \fn __global__ void listCpLdZeroPoly(sfixn *B, sfixn *A, sfixn length_poly, sfixn poly_on_layer)
\brief This function copy A to B, which is a list of  polynomials except the leading coefficient.
\param B list of polynomials.
\param A list of polynomials.
\param length_poly  length of each polynomial in A.
\param poly_on_layer no of poly in A.

*/

__global__ void listCpLdZeroPoly(sfixn *B, sfixn *A, sfixn length_poly, sfixn poly_on_layer)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
        {
                sfixn polyID = tid/length_poly;
                sfixn offset= tid%length_poly;
		if(offset != length_poly -1)
                	B[polyID*length_poly + offset]=  A[polyID*length_poly + offset];
		else
			 B[polyID*length_poly + offset]= 0;
        }

}


/*! \fn __global__ void listPolyDegInc( sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly, sfixn poly_on_layer, sfixn newLength)
\brief Each polynomial is same as corresponding poly in Mgpu
Except each of its degree is equal to newDegree.
The extra coefficients are padded by zero.

(poly_on_layer*newLength]) < 2^24

\param Mgpu list of polynomials.
\param extDegMgpu list of polynomials with extended length.
\param length_poly length of each polynomials in Mgpu.
\param poly_on_layer no of polynomials in Mgpu.
\param newLength new length of extDegMgpu.
*/
__global__ void listPolyDegInc( sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly, sfixn poly_on_layer, sfixn newLength)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*newLength)
	{
		sfixn polyID = tid/newLength;
		sfixn offset= tid%newLength;		
		if(offset < length_poly) 
			extDegMgpu[polyID*newLength + offset] = Mgpu[ polyID*length_poly + offset];
		else 
			extDegMgpu[polyID*newLength + offset] = 0;		
	}
}


/*! \fn __global__ void listCpUpperCuda(sfixn *dest, sfixn *source, sfixn n, sfixn l)
\brief The lower part (lower l coefficients) of each poly in source 
are copied to the zero part of the corresponding poly in destination.
\paramdest list of polynomials.
\param source list of polynomials.
\param n the number of integers copied from source.
\param l each polynomial of dest is of length l.
*/
__global__ void listCpUpperCuda(sfixn *dest, sfixn *source, sfixn n, sfixn l)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
	{
		sfixn x = tid/l;
		x = (x*l) + tid;
		dest[x + l] = source[x];	
	}
}


/*! \fn __global__ void listCpLowerCuda(sfixn *dest, sfixn *source, sfixn n, sfixn l)
\brief The lower part (lower l coefficients) of each poly in source 
are copied to the corresponding poly in destination.
\param dest list of polynomials.
\param source list of polynomials.
\param n the number of integers copied from source.
\param l each polynomial of dest is of length l.

*/
__global__ void listCpLowerCuda(sfixn *dest, sfixn *source, sfixn n, sfixn l)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
	{
		sfixn x = tid/l;
		sfixn y = tid%l;
		dest[tid] = source[x*2*l + y];	
		//dest[tid] = x;
	}
}





/*! \fn __global__ void list2wayCp(sfixn *dest, sfixn *source, sfixn l, sfixn totalLength, sfixn p)
\brief The lower l coeffs of each poly in source and
the upper l coeffs of each poly in dest will be added
and stored in the upper l coeffs of dest.

\param dest list of polynomials.
\param source list of polynomials.
\param l the no of coefficients from each polynomial need to be added.
\param totalLength total no of coefficients need to be added.
\param p the prime number of finite field.
*/
__global__ void list2wayCp(sfixn *dest, sfixn *source, sfixn l, sfixn totalLength, sfixn p)
{
        sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
        if(tid < totalLength)
        {
                sfixn L = 2*l;

                sfixn x = tid/l;
                sfixn y = tid%l;
                dest[x*L+l+y] = add_mod(dest[x*L+l+y] , source[x*L +y], p);
        }
}


/*! \fn __global__ void	leavesSubproductTree(sfixn *M1gpu, sfixn *Mgpu, sfixn numPoints, sfixn rightSubtree)
\brief The code is tested for numpoints = 2^24

\param M1gpu holds the numPoints*2 poits x_0, x_1...
\param Mgpu it will create numPoints polynomials (X - x_i). 
\param numPoints total number of points.
\param rightSubtree if rightSubtree == 0, it works with first half of Mgpu. else it works with the last half of Mgpu. 

*/

__global__ void	leavesSubproductTree(sfixn *M1gpu, sfixn *Mgpu, sfixn numPoints, sfixn rightSubtree)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < numPoints)
	{
		M1gpu[tid*2] = -Mgpu[rightSubtree*numPoints + tid ];
		M1gpu[tid*2+1] = 1;		
	}
}

/*! \fn void subProductTree(sfixn k, sfixn p)
\brief  The function for creating subproduct tree.
\param k the number of levels.
\param p the prime number for finite field.
*/

void subProductTree(sfixn k, sfixn p)
{
	sfixn polyLengthCurrent = 2; // the length of poly at leaves.
	sfixn polyOnLayerCurrent = 1L << (k-1); //the number of poly in the leaves of each subtree. 
	sfixn polyLengthNext, polyOnLayerNext;
	sfixn threadsForAmul, mulInThreadBlock, blockNo;
	sfixn L = 1L << (k-1); // spaces required for a level in any subtree from plainMulLimit
	sfixn l;
	sfixn w, winv, ninv;
	sfixn *Al, *Bl, *Cl; 
	cudaMalloc((void **)&Al,    sizeof(sfixn)*L); // temporary storage required for doing FFT based multiplication.
	cudaMalloc((void **)&Bl,    sizeof(sfixn)*L); // temporary storage required for doing FFT based multiplication.			
	sfixn *Ar, *Br, *Cr; 
	cudaMalloc((void **)&Ar,    sizeof(sfixn)*L); // temporary storage required for doing FFT based multiplication.
	cudaMalloc((void **)&Br,    sizeof(sfixn)*L); // temporary storage required for doing FFT based multiplication.			


	for(sfixn i = 1; i < k; ++i)
	{
		if(i <= plainMulLimit)
		{
			polyOnLayerNext = polyOnLayerCurrent/2; // at the immediate upper level poly number will be reduced by half.
			polyLengthNext = 2*polyLengthCurrent -1; // at the immediate upper level poly length will be double minus 1.

			threadsForAmul = 2*polyLengthCurrent; // how many threads are necessary for each plain multiplication.
                        mulInThreadBlock = (sfixn)floor((double)Tmul/(double)threadsForAmul); // how many multiplications that a thread block can do.
                        blockNo = (sfixn)ceil( ((double) polyOnLayerCurrent/(double) mulInThreadBlock)*0.5  );			

			cudaMalloc((void **)&check.Ml[i], sizeof(sfixn)*polyOnLayerNext*polyLengthNext); // allocating space for the immediate upper level of left  subtree
			cudaThreadSynchronize();
			cudaMalloc((void **)&check.Mr[i], sizeof(sfixn)*polyOnLayerNext*polyLengthNext); // allocating space for the immediate upper level of right subtree
			cudaThreadSynchronize();
			listPlainMulGpu<<<blockNo, Tmul>>>(check.Ml[i-1], check.Ml[i], polyLengthCurrent, polyOnLayerCurrent, threadsForAmul, mulInThreadBlock, p);
			cudaThreadSynchronize(); // creating the immediate upper level of left  subtree
			listPlainMulGpu<<<blockNo, Tmul>>>(check.Mr[i-1], check.Mr[i], polyLengthCurrent, polyOnLayerCurrent, threadsForAmul, mulInThreadBlock, p);
			cudaThreadSynchronize(); // creating the immediate upper level of right subtree

			if(i == plainMulLimit)
			{				
				// We need the subinverse tree at level plainMulLimit directly.
				// From this level we will not store the leading coefficient of subproduct tree.
				// So poly in supproduct tree created at this level should also need to be modified. 
				//We will use the pointer to the next level as temporary storage.
				cudaMalloc((void **)&check.Ml[i+1], sizeof(sfixn)*(polyOnLayerNext)*(polyLengthNext-1));
				cudaThreadSynchronize(); 
				cudaMalloc((void **)&check.Mr[i+1], sizeof(sfixn)*(polyOnLayerNext)*(polyLengthNext-1));
				cudaThreadSynchronize(); 

				cudaMalloc((void **)&check.InvMl[i], sizeof(sfixn)*(polyOnLayerNext)*(polyLengthNext-1));
                                cudaThreadSynchronize();// allocating space for subinverse tree of left  subtree for this level
				cudaMalloc((void **)&check.InvMr[i], sizeof(sfixn)*(polyOnLayerNext)*(polyLengthNext-1));
                                cudaThreadSynchronize();// allocating space for subinverse tree of right subtree for this level

				blockNo = (sfixn)ceil((double)(polyOnLayerNext)/(double)(Tinv));
				listPolyinv<<<blockNo, Tinv>>>(check.Ml[i], check.InvMl[i], polyOnLayerNext, p);
				cudaThreadSynchronize();// creating subinverse tree of left  subtree for this level
				listPolyinv<<<blockNo, Tinv>>>(check.Mr[i], check.InvMr[i], polyOnLayerNext, p);
				cudaThreadSynchronize();// creating subinverse tree of right subtree for this level
				
				// Now we will store the poly in subproduct tree at this level according to our new data structure. i.e no leading coefficient

				copyMgpu<<<polyOnLayerNext ,(polyLengthNext -1)>>>(check.Ml[i+1], check.Ml[i], polyLengthNext);
				copyMgpu<<<polyOnLayerNext ,(polyLengthNext -1)>>>(check.Mr[i+1], check.Mr[i], polyLengthNext);
				cudaThreadSynchronize();

				cudaFree(check.Ml[i]); cudaFree(check.Mr[i]); // free the temporary spaces				
				check.Ml[i] = check.Ml[i+1]; check.Mr[i] = check.Mr[i+1]; // adjust the pointer	
							
			}			
			polyLengthCurrent = polyLengthNext;
			polyOnLayerCurrent = polyOnLayerNext;

		}
		else
		{

			l = 1L << (i-1); // length of each poly at (i-1)th level.

			Cl = check.Ml[i-1]; // Cl points the poly of level (i-1)th of left  subtree
			Cr = check.Mr[i-1]; // Cr points the poly of level (i-1)th of right subtree

			cudaMemcpy(Al,   Cl,     sizeof(sfixn)*L,     cudaMemcpyDeviceToDevice); // Al contains the poly of level (i-1)th of left  subtree
			cudaMemcpy(Ar,   Cr,     sizeof(sfixn)*L,     cudaMemcpyDeviceToDevice); // Ar contains the poly of level (i-1)th of right subtree
			cudaMemcpy(Bl, &(Cl[l]), sizeof(sfixn)*(L-l), cudaMemcpyDeviceToDevice); // Bl contains the poly of level (i-1)th of left  subtree excluding the first poly	
			cudaMemcpy(Br, &(Cr[l]), sizeof(sfixn)*(L-l), cudaMemcpyDeviceToDevice); // Bl contains the poly of level (i-1)th of right subtree excluding the first poly		
			// Every alternate poly is making 0 poly. By this the poly need to be multiplied are placed in the same indices. This is done for FFT. Now length of the poly 
			// becomes double. 
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(Tmax*2)), Tmax>>>(Al, Bl, L/2, l ); 
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(Tmax*2)), Tmax>>>(Ar, Br, L/2, l );
			// Allocating space for the poly in the immediate upper level
			cudaMalloc((void **)&check.Ml[i], sizeof(sfixn)*L);
			cudaMalloc((void **)&check.Mr[i], sizeof(sfixn)*L);
			// initializing check.Ml[i] and check.Mr[i] as zero array
			allZero<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>(check.Ml[i], L);	
			allZero<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>(check.Mr[i], L);
			// As we are not storing the leading coefficients, we need to add the coefficients with it.
			pointAdd2<<<(sfixn)ceil((double)(L/2)/(double)Tmax) , Tmax>>>(check.Ml[i] , check.Ml[i-1], l, L/2, p);	
			pointAdd2<<<(sfixn)ceil((double)(L/2)/(double)Tmax) , Tmax>>>(check.Mr[i] , check.Mr[i-1], l, L/2, p);	

			w = primitive_root(i, p); // primitive root of unity of 2^i
			l = 1L << (k-i -1); // no of poly in Al, Bl, Ar, Br
			// computing FFT
			list_stockham_dev(Al, l, i, w, p);
			list_stockham_dev(Bl, l, i, w, p);
			list_stockham_dev(Ar, l, i, w, p);
			list_stockham_dev(Br, l, i, w, p);
			//pointwise multiplication
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( Al, Bl, L, p);
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( Ar, Br, L, p);
			// inverse FFT
			winv = inv_mod(w, p);
			list_stockham_dev(Al, l, i, winv, p);
			list_stockham_dev(Ar, l, i, winv, p);

			w = (1L << i);
			ninv = inv_mod(w, p);			
			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( Al, ninv, L, p);
			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( Ar, ninv, L, p);
		
			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( check.Ml[i], Al, L, p);	
			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax) , Tmax>>>( check.Mr[i], Ar, L, p);						
		
		}
	}
	cudaFree(Al);
	cudaFree(Bl);	
	cudaFree(Ar);
	cudaFree(Br);	
}


/*! \fn void subInvTree(sfixn k, sfixn p)
\brief  The function for creating subinverse tree.
Please read Chapter 10 of the PhD thesis of Sardar Haque (Computer Science, UWO, 2013)
for the definition and implementation details of subinverse tree
\param k the number of levels.
\param p the prime number for finite field.
*/
void subInvTree(sfixn k, sfixn p)
{
	sfixn *Al, *Bl, *Cl, *Dl, *El, *Fl, *Gl;
	sfixn *Ar, *Br, *Cr, *Dr, *Er, *Fr, *Gr; 
	sfixn L = 1L << (k-1);     // the number of coefficients representing the polys of any subtree of subinverse tree of one level. 
	sfixn l = sizeof(sfixn)*L; // the length of the array storing the polys of any subtree of subinverse tree of one level. 
	sfixn j;
	cudaMalloc((void **)&Al, l);	cudaMalloc((void **)&Bl, l);	cudaMalloc((void **)&Cl, l);
	cudaMalloc((void **)&Ar, l);	cudaMalloc((void **)&Br, l);	cudaMalloc((void **)&Cr, l);
	l = l*2;
	cudaMalloc((void **)&Dl, l);	cudaMalloc((void **)&El, l);	cudaMalloc((void **)&Fl, l);	cudaMalloc((void **)&Gl, l);
	cudaMalloc((void **)&Dr, l);	cudaMalloc((void **)&Er, l);	cudaMalloc((void **)&Fr, l);	cudaMalloc((void **)&Gr, l);

	
	sfixn blockNo1 = (sfixn)ceil((double)(L)/(double)(Tmax));
	sfixn blockNo2 = (sfixn)ceil((double)(L*2)/(double)(Tmax));
	sfixn w, winv, ninv;
	for(sfixn i =  plainMulLimit; i < (k-1); ++i)
	{
		l =  1L << (i); // the length of each poly of subinverse tree at level i
		j = 1L << (k - i - 1); // the number of polys of any subtree of subinverse tree at level i
		/*
		rev(M)*inv_i(M) = 1 mod x^{2^i}
		We need inv_i(M) such that rev(M)*inv_i(M) = 1 mod x^{2^{i+1}}
		We need to compute rev(M)*inv_i(M) mod x^{2^{i+1}}
		We have, M' = M without leading term.
		I = rev(inv_i(M)) without leading term.
		M' and inv_i(M) has same degree.
		M*rev(inv_i(M)) = rev(conv(M', rev(inv_i(M))) + I) trailing by (2^i - 1) zeros then 1
		*/		
		//rev(inv_i(M))
		listReversePoly<<<blockNo1,  Tmax>>>(Al, check.InvMl[i], l, j); // reversing polys in level i of left  subtree of subinverse tree
		listReversePoly<<<blockNo1,  Tmax>>>(Ar, check.InvMr[i], l, j); // reversing polys in level i of right subtree of subinverse tree
		//I
		listCpLdZeroPoly<<<blockNo1, Tmax>>>(Bl, Al ,l , j); // copying Al to Bl except the leading coefficient which is set 0 in Bl
		listCpLdZeroPoly<<<blockNo1, Tmax>>>(Br, Ar ,l , j); // copying Ar to Br except the leading coefficient which is set 0 in Br
		//M'
		cudaMemcpy(Cl, check.Ml[i], sizeof(sfixn)*L, cudaMemcpyDeviceToDevice); // copying the left  subtree of subproduct tree into Cl
		cudaMemcpy(Cr, check.Mr[i], sizeof(sfixn)*L, cudaMemcpyDeviceToDevice); // copying the right subtree of subproduct tree into Cl
		// conv(M', rev(inv_i(M)))
		w =  primitive_root(i, p);
		list_stockham_dev(Al, j, i, w, p);
		list_stockham_dev(Ar, j, i, w, p);
		list_stockham_dev(Cl, j, i, w, p);
		list_stockham_dev(Cr, j, i, w, p);

		pointMul<<<blockNo1 , Tmax>>>( Al, Cl, L, p);
		pointMul<<<blockNo1 , Tmax>>>( Ar, Cr, L, p);
		winv = inv_mod(w, p);	
		list_stockham_dev(Al, j, i, winv, p);
		list_stockham_dev(Ar, j, i, winv, p);

		ninv = inv_mod(l, p);
		scalarMul<<<blockNo1, Tmax>>>(Al, ninv, L, p);
		scalarMul<<<blockNo1, Tmax>>>(Ar, ninv, L, p);
		// adding I
		pointAdd<<<blockNo1,  Tmax>>> (Al, Bl, L, p);
		pointAdd<<<blockNo1,  Tmax>>> (Ar, Br, L, p);
		//reversing
		listReversePoly<<<blockNo1, Tmax>>>(Bl, Al, l , j);
		listReversePoly<<<blockNo1, Tmax>>>(Br, Ar, l , j);
		//setting negative to all coefficients 		
		allNeg<<<blockNo1, Tmax >>>(Bl, L, p);
		allNeg<<<blockNo1, Tmax >>>(Br, L, p);
		// making inv_i(M) bigger poly by padding zeros
		listPolyDegInc<<<blockNo2, Tmax>>>(check.InvMl[i], Dl, l,   j, 2*l);
		listPolyDegInc<<<blockNo2, Tmax>>>(check.InvMr[i], Dr, l,   j, 2*l);
		// making -(M*rev(inv_i(M))) bigger by padding zeros
		listPolyDegInc<<<blockNo2, Tmax>>>(Bl, El, l,   j, 2*l);
		listPolyDegInc<<<blockNo2, Tmax>>>(Br, Er, l,   j, 2*l);
		cudaMemcpy(Fl,  Dl,  sizeof(sfixn)*L*2,  cudaMemcpyDeviceToDevice);
		cudaMemcpy(Fr,  Dr,  sizeof(sfixn)*L*2,  cudaMemcpyDeviceToDevice);
		// -(M*rev(inv_i(M))) * inv_i(M)
		w =  primitive_root(i+1, p);
		list_stockham_dev(El, j, i+1, w, p);
		list_stockham_dev(Er, j, i+1, w, p);
		list_stockham_dev(Dl, j, i+1, w, p);		
		list_stockham_dev(Dr, j, i+1, w, p);
		pointMul<<<blockNo2, Tmax>>>( El, Dl, 2*L, p);
		pointMul<<<blockNo2, Tmax>>>( Er, Dr, 2*L, p);
		winv = inv_mod(w, p);
		list_stockham_dev(El, j, i+1, winv, p);
		list_stockham_dev(Er, j, i+1, winv, p);
		ninv = inv_mod(l*2, p);
		scalarMul<<<blockNo2, Tmax>>>( El, ninv, L*2, p);
		scalarMul<<<blockNo2, Tmax>>>( Er, ninv, L*2, p);
		// copy the upper part as the lower part is known from inv_i(M) or F
		listCpUpperCuda<<<blockNo1 , Tmax>>>( Fl, El, L, l);
		listCpUpperCuda<<<blockNo1 , Tmax>>>( Fr, Er, L, l);
		//Now F contains inv_i(M) with more precision
		//Multiply pair of inv_i(M) to get inv_{i+1}(M)
		cudaMemcpy(Gl, &(Fl[2*l]), sizeof(sfixn)*(2*L-2*l), cudaMemcpyDeviceToDevice);
		cudaMemcpy(Gr, &(Fr[2*l]), sizeof(sfixn)*(2*L-2*l), cudaMemcpyDeviceToDevice);
		zeroInbetween<<<blockNo1, Tmax>>>(Fl, Gl, L, 2*l );
		zeroInbetween<<<blockNo1, Tmax>>>(Fr, Gr, L, 2*l );		
		j = j/2;
		w =  primitive_root(i+2, p);
		list_stockham_dev(Fl, j, i+2, w,   p);
		list_stockham_dev(Fr, j, i+2, w,   p);
		list_stockham_dev(Gl, j, i+2, w,   p);
		list_stockham_dev(Gr, j, i+2, w,   p);
		pointMul<<<blockNo2 , Tmax>>>( Fl, Gl, L*2, p);
		pointMul<<<blockNo2 , Tmax>>>( Fr, Gr, L*2, p);
		winv = inv_mod(w, p);
		list_stockham_dev(Fl, j, i+2, winv, p);
		list_stockham_dev(Fr, j, i+2, winv, p);
		ninv = inv_mod(l*4, p);
		scalarMul<<<blockNo2 , Tmax>>>( Fl, ninv, 2*L, p);
		scalarMul<<<blockNo2 , Tmax>>>( Fr, ninv, 2*L, p);
		cudaMalloc((void **)&check.InvMl[i+1], sizeof(sfixn)*L);
		cudaMalloc((void **)&check.InvMr[i+1], sizeof(sfixn)*L);
		listCpLowerCuda<<<blockNo1 , Tmax>>>( check.InvMl[i+1], Fl, L, 2*l);
		listCpLowerCuda<<<blockNo1 , Tmax>>>( check.InvMr[i+1], Fr, L, 2*l);
		cudaThreadSynchronize();		
	}
	cudaFree(Al);	cudaFree(Bl);	cudaFree(Cl);	cudaFree(Dl);	cudaFree(El);	cudaFree(Fl);	cudaFree(Gl);
	cudaFree(Ar);	cudaFree(Br);	cudaFree(Cr);	cudaFree(Dr);	cudaFree(Er);	cudaFree(Fr);	cudaFree(Gr);
}


/*! \fn __global__ void EvalistReversePolyDec( sfixn *Mgpu, sfixn *revMgpu, sfixn length_poly, sfixn poly_on_layer)
\brief This function reverse each poly in Mgpu and decrease its degree by putting zeros in its significant coefficients.	

\param Mgpu is a list of poly_on_layer/2 polynomials of length lengthFpoly of type F.
\param revMgpu is the reverse of fGPU except the leading lengthFpoly/2 coefficients are zero.
\param length_poly length of each poly in Mgpu.
\param poly_on_layer no of poly in Mgpu.

*/
__global__ void EvalistReversePolyDec( sfixn *Mgpu, sfixn *revMgpu, sfixn length_poly, sfixn poly_on_layer)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
	{
		sfixn polyID = tid/length_poly;
		sfixn offset= tid%length_poly;
		if(offset >= (length_poly/2) )
			revMgpu[(polyID+1)*length_poly - offset -1  ] = Mgpu[tid];
		else
			revMgpu[(polyID+1)*length_poly - offset -1  ] = 0;
	}
}


/*! \fn __global__ void EvalistPolyDegInc( sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly, sfixn poly_on_layer, sfixn newLength)
\brief MinvGPU is the array of poly_on_layer polynomoals of Minverse type. Each of their length is lengthMinvPoly.
"extMinvGPU" is the list of "poly_on_layer" polynomials of lengthMinvPoly*2 length.
These two polynomials are same except padded by zeros.

\param Mgpu is the array of poly_on_layer polynomoals of Minverse type.
\param extDegMgpu is the list of "poly_on_layer" polynomials of lengthMinvPoly*2 length.
\param length_poly length of each poly in Mgpu.
\param poly_on_layer no of poly.
\param newLength new length of the polynomials.

*/
__global__ void EvalistPolyDegInc( sfixn *Mgpu, sfixn *extDegMgpu, sfixn length_poly, sfixn poly_on_layer, sfixn newLength)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*newLength)
	{
		sfixn polyID = tid/newLength;
		sfixn offset= tid%newLength;		
		if(offset < length_poly) 
			extDegMgpu[tid] = Mgpu[ polyID*length_poly + offset];
		else 
			extDegMgpu[tid] = 0;		
	}
}

/*! \fn __global__ void EvalistZeroBetweenRev(sfixn *revMzero, sfixn *M, sfixn length_poly, sfixn poly_on_layer)
\brief This function reverse the poly M into revMzero and makes the lower half of revMzero as list of zero entries.
\param revMzero list of poly.
\param M list of poly.
\param length_poly  length of each poly in M.
\param poly_on_layer no of poly.

*/

__global__ void EvalistZeroBetweenRev(sfixn *revMzero, sfixn *M, sfixn length_poly, sfixn poly_on_layer)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
	{
		revMzero[tid] = 0;

		sfixn polyID = tid/length_poly;
		sfixn offset= tid%length_poly;
		
		if(offset < (length_poly/2))
			revMzero[tid] = M[ polyID*length_poly + ( ((length_poly/2) -1) - offset)  ];		
	}
}

/*! \fn __global__ void subtractLowerSingle(sfixn *U, sfixn *F, sfixn *G, sfixn n, sfixn p)
\brief Both G and F are of 2n length
This kerenel subtract the lower n coefficient and stored in U
\param U list of poly.
\param F list of poly.
\param  G  list of poly.

*/
__global__ void subtractLowerSingle(sfixn *U, sfixn *F, sfixn *G, sfixn n, sfixn p)
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < n)
		U[tid] = sub_mod(F[tid], G[tid], p);
}

/*! \fn __global__ void pointwiseMulHalf(sfixn *A, sfixn *B, sfixn length_poly, sfixn poly_on_layer, sfixn p )
\brief This function does  pointwise multiplcation between
the lower half of A and B.
\param A  list of poly.
\param B list of poly.
\param length_poly length of each poly.
\param poly_on_layer no of poly in A.
\param p prime number for finite field.
*/
__global__ void pointwiseMulHalf(sfixn *A, sfixn *B, sfixn length_poly, sfixn poly_on_layer, sfixn p )
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
	{
		sfixn polyID = tid/length_poly;
		sfixn offset= tid%length_poly;
		A[tid] = mul_mod(A[tid], B[ (polyID/2)*length_poly + offset], p);
	}
}

/*! \fn __global__ void subtractLower(sfixn *R, sfixn *A, sfixn *BQ, sfixn length_poly, sfixn poly_on_layer, sfixn p )
\brief This function does  pointwise subtraction between
the lower half of A and BQ.
\param R  list of poly.
\param A  list of poly.
\param BQ list of poly.
\param length_poly length of each poly.
\param poly_on_layer no of poly in A.
\param p prime number for finite field.
*/
__global__ void subtractLower(sfixn *R, sfixn *A, sfixn *BQ, sfixn length_poly, sfixn poly_on_layer, sfixn p )
{
	sfixn tid= blockIdx.x*blockDim.x + threadIdx.x;
	if(tid < poly_on_layer*length_poly)
	{
		sfixn polyIDinv = tid/(length_poly/2);
		sfixn polyIDf   = polyIDinv/2;
		sfixn offset   = tid%(length_poly/2);

		R[tid] = sub_mod(A[ (polyIDf*length_poly )+ offset],  BQ[ (polyIDinv*length_poly) + offset],p);
	}
}

/*! \fn __global__ void listPlainCudaDiv(sfixn *M, sfixn *F, sfixn start, sfixn length, sfixn threadsPerDiv, sfixn DivPerBlock, sfixn polyNum, sfixn P)
\brief This function does a list of plain division.
F = F mod M.
\param M list of poly.
\param F list of poly.
\param start it is 0 (to support if starting point is not 0).
\param length lenght of each poly in M.
\param threadsPerDiv no of threads required for a division.
\param DivPerBlock no of division per thread block.
\param polyNum number of poly in M.
\param P prime number for prime field.
*/
__global__ void listPlainCudaDiv(sfixn *M, sfixn *F, sfixn start, sfixn length, sfixn threadsPerDiv, sfixn DivPerBlock, sfixn polyNum, sfixn P)
{
	__shared__ sfixn sM[Tdiv*W];
	__shared__ sfixn sF[Tdiv*2*W];
	__shared__ sfixn invM[Tdiv];
	__shared__ sfixn oneF[Tdiv];
	__shared__ sfixn mID[Tdiv];
	__shared__ sfixn fID[Tdiv];
	__shared__ sfixn startF[Tdiv];
	
	sfixn i, j, k, l, s, t, polyID;
	polyID = ((threadIdx.x/threadsPerDiv) + blockIdx.x*DivPerBlock);
	if( polyID < polyNum && threadIdx.x < threadsPerDiv*DivPerBlock)
	{
		if( (threadIdx.x %threadsPerDiv) == (threadsPerDiv-1))
		{
			i = threadIdx.x;				
			mID[i] = (i/ threadsPerDiv) + (blockIdx.x*DivPerBlock);
			fID[i] = mID[i]/ 2;
	
			mID[i] = start + (mID[i]+1)*length -1;
			fID[i] = (fID[i]+1)*(2*length-2) - 1;
			invM[i] = inv_mod(M[mID[i]], P);			
		}
		else
		{
			i = threadIdx.x/threadsPerDiv;
			i = (i*threadsPerDiv) + threadsPerDiv -1;
		}
		__syncthreads();
		
		j = threadIdx.x;
		k = i - j;
		t = W-1;
		l = mID[i] - k*W;
		s = l - W;
		for(; l > s && l >= 0; --l, --t ) 
			sM[j*W + t] = M[l];	
	
		l = fID[i] - (k*W*2);
		s = l - 2*W;
		t = 2*W-1;	
		for(; l > s && l >= 0; --l, --t ) 
			sF[j*2*W +t ] = F[l];

		if(i == j)
		{
			fID[i] = 2*i*W + 2*W -1;
			mID[i] = fID[i] - length + 1;
			while( sF[ fID[ i ] ] == 0 &&  fID[ i ] > mID[i] ) --fID[ i ];

			oneF[i] = mul_mod(invM[i], sF[fID[i]], P);			
		}
		
		__syncthreads();


		//sfixn temp1 =0 , temp2 = 4;
		while(fID[ i ] > mID[i]) //while(temp1 < temp2)
		{
			//--temp2;

			l = fID[ i ] - k*W;
			s = l - W;
			if( (fID[ i ] - length) > s )	s = fID[ i ] - length;
			t = i*W + W -1 -k*W;
			for(; l > s; --l, --t)
				sF[l] = sub_mod(sF[l], mul_mod(sM[t], oneF[i], P), P);


			if(i == j)
			{
				--fID[ i ];
				while( sF[ fID[ i ] ] == 0 &&  fID[ i ] > mID[i] ) --fID[ i ];
				oneF[i] = mul_mod(invM[i], sF[fID[i]], P);
			
			}
			__syncthreads();
		}

		if(i == j)
		{
			fID[i] = mID[i] - length +1;
			startF[i] =  ( (polyID/2) +1)*(2*length-2) - 1;
			
			if(polyID%2 == 0)
				startF[i] = startF[i] - length +1;
			
		}
		__syncthreads();

	


		

		l = mID[ i ] - k*W;
		t = l - W;
		s =  startF[i] - k*W;
		if(t < fID[i]) t = fID[i];
		for(; l >  t; --l, --s)
			F[s] = sF[l];
		
			
	}		
}


const sfixn singDivLength = (1 << plainMulLimit);

/*! \fn __global__ void PlainCudaDiv(sfixn *M, sfixn *F, sfixn m, sfixn f, sfixn p)
\brief This is another implementation of plain division.
when the value of k is small s.t we
need to do plain arithmetic from the beninning.
we have a situation where we have one F polynomial 
and we need to divide it with one poly of M type.
this code does it.
\param M one of the top poly of subproduct tree.
\param F the polynomial need to be evaluated.
\param m length of M.
\param f length of F.
\param p prime number for finite field.
*/
__global__ void PlainCudaDiv(sfixn *M, sfixn *F, sfixn m, sfixn f, sfixn p)
{
	sfixn sM[singDivLength];
	sfixn sF[singDivLength];
	sfixn i;

	for(i = 0; i < m; ++i)    sM[i] = M[i];
	for(i = 0; i < f; ++i)    sF[i] = F[i];
	
	sfixn inv = inv_mod(sM[m-1], p);
	sfixn factor, j;
	
	for(i = f; i >= m; --i)
	{
		factor = mul_mod(inv, sF[i-1], p);
		for(j = 0; j < m; ++j)
			sF[i - 1 - j] = sub_mod(sF[i - 1 - j], mul_mod(sM[m-1-j], factor, p),  p);
	}
		
	for(i = 0; i < (m-1); ++i)
		F[i] = sF[i];
}

/*! \fn void fastEvaluationLow(sfixn *F, sfixn k, sfixn p, sfixn flag)
\brief This function does the polynomial evaluation
if the value of k is small s.t
all computations are done by plain polynomial arithmetics.
\param F the polynomial need to be evaluated.
\param points the list of points for which F need to be evaluated.
\param k the height of subproduct tree.
\param p prime number for finite field.
\param flag what we want to  do for checking the result with maplesoft.

*/

void fastEvaluationLow(sfixn *F, sfixn k, sfixn p, sfixn flag)
{

	sfixn L = (1 << k);
	sfixn length, no, threadsForAdiv, divInThreadBlock, blockNo;
	
	cudaMalloc((void **)&check.fL[k+1], sizeof(sfixn)*L);
        cudaMalloc((void **)&check.fR[k+1], sizeof(sfixn)*L);

	cudaMemcpy(check.fL[k+1], F, sizeof(sfixn)*L, cudaMemcpyHostToDevice); 
	cudaMemcpy(check.fR[k+1], F, sizeof(sfixn)*L, cudaMemcpyHostToDevice);	

	PlainCudaDiv<<<1,1>>>(check.Ml[k-1], check.fL[k+1], (L/2)+1, L, p);
	PlainCudaDiv<<<1,1>>>(check.Mr[k-1], check.fR[k+1], (L/2)+1, L, p);

	cudaMalloc((void **)&check.fL[k], sizeof(sfixn)*(L/2));
	cudaMalloc((void **)&check.fR[k], sizeof(sfixn)*(L/2));

	cudaMemcpy( check.fL[k], check.fL[k+1], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);
	cudaMemcpy( check.fR[k], check.fR[k+1], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);

	for(sfixn i =  k-2; i >= 0; --i )
	{	
		 //We took M[i],  F[i+2] and produce F[i+1]
                length = (1L<<i) + 1;
                no     = 1L<<(k-i-1);
		threadsForAdiv = (sfixn)ceil((double)length/(double)W);
		divInThreadBlock = (sfixn)floor((double)Tdiv/(double)threadsForAdiv);
		blockNo = (sfixn)ceil((double)no/(double)divInThreadBlock);

		cudaMalloc((void **)&check.fL[i+1], sizeof(sfixn)*(L/2));
                cudaMalloc((void **)&check.fR[i+1], sizeof(sfixn)*(L/2));

		cudaMemcpy( check.fL[i+1], check.fL[i+2], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);
	        cudaMemcpy( check.fR[i+1], check.fR[i+2], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);
						
		listPlainCudaDiv<<<blockNo,Tdiv>>>(check.Ml[i], check.fL[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);
		listPlainCudaDiv<<<blockNo,Tdiv>>>(check.Mr[i], check.fR[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);

		if(flag != 1)
		{
			cudaFree(check.Ml[i]);
			cudaFree(check.Mr[i]);

			cudaFree(check.fL[i+2]);
			cudaFree(check.fR[i+2]);
		}

	}

}

/*! \fn void fastEvaluation(sfixn *F, sfixn k, sfixn p,  sfixn flag)
\brief This function does the polynomial evaluation
if the value of k is not small. 
we are relying on FFT based arithmetic from the 
top level (k) to plainMulLimit level.
\param F the polynomial need to be evaluated.
\param points the list of points for which F need to be evaluated.
\param k the height of subproduct tree.
\param p prime number for finite field.
\param flag what we want to  do for checking the result with maplesoft.

*/

void fastEvaluation(sfixn *F, sfixn k, sfixn p,  sfixn flag)
{
	sfixn i, ninv, L, blockNo, w, winv;
	L = (1 << k);
	sfixn *AL, *AR, *BL, *BR, *CL, *CR, *DL, *DR, *Fgpu, *revFgpu;
	cudaMalloc((void **)&AL,      sizeof(sfixn)*L);
	cudaMalloc((void **)&AR,      sizeof(sfixn)*L);
	cudaMalloc((void **)&BL,      sizeof(sfixn)*L);
	cudaMalloc((void **)&BR,      sizeof(sfixn)*L);
	cudaMalloc((void **)&CL,      sizeof(sfixn)*L);
        cudaMalloc((void **)&CR,      sizeof(sfixn)*L);
	cudaMalloc((void **)&DL,      sizeof(sfixn)*(L/2));
        cudaMalloc((void **)&DR,      sizeof(sfixn)*(L/2));

	cudaMalloc((void **)&Fgpu,    sizeof(sfixn)*L);
	cudaMalloc((void **)&revFgpu, sizeof(sfixn)*L);

	cudaMemcpy(Fgpu, F, sizeof(sfixn)*L, cudaMemcpyHostToDevice); 

	blockNo = (sfixn)ceil((double)(L)/(double)(Tmax));

	EvalistPolyDegInc<<<blockNo, Tmax>>>(check.InvMl[k-1], AL, L/2, 1, L);
	EvalistPolyDegInc<<<blockNo, Tmax>>>(check.InvMr[k-1], AR, L/2, 1, L);

	EvalistReversePolyDec<<<blockNo, Tmax>>>(Fgpu, revFgpu, L, 1);
	
	w =  primitive_root(k, p);
	list_stockham_dev(revFgpu, 1, k, w, p);
	list_stockham_dev(AL,      1, k, w, p);
	list_stockham_dev(AR,      1, k, w, p);

	pointMul<<<blockNo , Tmax>>>( AL, revFgpu, L, p);
	pointMul<<<blockNo , Tmax>>>( AR, revFgpu, L, p);
	winv = inv_mod(w, p);

	list_stockham_dev(AL,   1, k, winv, p);
	list_stockham_dev(AR,   1, k, winv, p);
	ninv = inv_mod(L, p);
	scalarMul<<<blockNo, Tmax>>>( AL, ninv, L, p);
	scalarMul<<<blockNo, Tmax>>>( AR, ninv, L, p);

	EvalistZeroBetweenRev<<<blockNo, Tmax>>>(BL, AL, L, 1);
	EvalistZeroBetweenRev<<<blockNo, Tmax>>>(BR, AR, L, 1);


 	EvalistPolyDegInc<<<blockNo, Tmax>>>(check.Ml[k-1], CL, L/2, 1, L);
        EvalistPolyDegInc<<<blockNo, Tmax>>>(check.Mr[k-1], CR, L/2, 1, L);



	w =  primitive_root(k, p);
        list_stockham_dev(CL, 1, k, w, p);
        list_stockham_dev(CR, 1, k, w, p);
        list_stockham_dev(BL, 1, k, w, p);
        list_stockham_dev(BR, 1, k, w, p);


        pointMul<<<blockNo , Tmax>>>( BL, CL, L, p);
        pointMul<<<blockNo , Tmax>>>( BR, CR, L, p);
        winv = inv_mod(w, p);

        list_stockham_dev(BL, 1, k, winv, p);
        list_stockham_dev(BR, 1, k, winv, p);
        ninv = inv_mod(L, p);
        scalarMul<<<blockNo, Tmax>>>( BL, ninv, L, p);
        scalarMul<<<blockNo, Tmax>>>( BR, ninv, L, p);



	blockNo = (sfixn)ceil((double)(L)/(double)(Tmax * 2.0));

	cudaMalloc((void **)&check.fL[k], sizeof(sfixn)*(L/2));
	cudaMalloc((void **)&check.fR[k], sizeof(sfixn)*(L/2));

	subtractLowerSingle<<<blockNo, Tmax>>>(check.fL[k] , Fgpu, BL, L/2, p);
	subtractLowerSingle<<<blockNo, Tmax>>>(check.fR[k], Fgpu, BR, L/2, p);

	cudaFree(Fgpu);
        cudaFree(revFgpu);
	if(flag != 1)
	{
		cudaFree(check.Ml[k-1]);
		cudaFree(check.Mr[k-1]);
		cudaFree(check.InvMl[k-1]);
		cudaFree(check.InvMr[k-1]);
	}
			

	sfixn length, no;
	for(i = k-2; i >=  plainMulLimit; --i)
	{
		//We took M[i], invM[i], F[i+2] and produce F[i+1]
		length = 1L<<i;
		no     = 1L<<(k-i-1); 

		blockNo = (sfixn)ceil((double)(L)/(double)(Tmax));
		EvalistPolyDegInc<<<blockNo, Tmax>>>(check.InvMl[i], AL, length, no, length*2);
	        EvalistPolyDegInc<<<blockNo, Tmax>>>(check.InvMr[i], AR, length, no, length*2);

 		blockNo = (sfixn)ceil((double)(L)/(double)(Tmax*2.0));
		EvalistReversePolyDec<<<blockNo, Tmax>>>(check.fL[i+2], DL, length*2, no/2);
		EvalistReversePolyDec<<<blockNo, Tmax>>>(check.fR[i+2], DR, length*2, no/2);		

		w =  primitive_root(i+1, p);
	        list_stockham_dev(DL, no/2, i+1, w, p);
	        list_stockham_dev(DR, no/2, i+1, w, p);	
	        list_stockham_dev(AL, no,   i+1, w, p);
        	list_stockham_dev(AR, no,   i+1, w, p);

		blockNo = (sfixn)ceil((double)(L)/(double)(Tmax));
		pointwiseMulHalf<<<blockNo , Tmax>>>( AL, DL, length*2, no, p);
                pointwiseMulHalf<<<blockNo , Tmax>>>( AR, DR, length*2, no, p);
	        winv = inv_mod(w, p);

        	list_stockham_dev(AL, no, i+1, winv, p);
	        list_stockham_dev(AR, no, i+1, winv, p);
        	ninv = inv_mod(length*2, p);
	        scalarMul<<<blockNo, Tmax>>>( AL, ninv, L, p);
        	scalarMul<<<blockNo, Tmax>>>( AR, ninv, L, p);

		EvalistZeroBetweenRev<<<blockNo, Tmax>>>(BL, AL, length*2, no);
	        EvalistZeroBetweenRev<<<blockNo, Tmax>>>(BR, AR, length*2, no);

	        EvalistPolyDegInc<<<blockNo, Tmax>>>(check.Ml[i], CL, length, no, length*2);
	        EvalistPolyDegInc<<<blockNo, Tmax>>>(check.Mr[i], CR, length, no, length*2);
	
	        list_stockham_dev(CL, no, i+1, w, p);
        	list_stockham_dev(CR, no, i+1, w, p);
	        list_stockham_dev(BL, no, i+1, w, p);
	        list_stockham_dev(BR, no, i+1, w, p);

	        pointMul<<<blockNo , Tmax>>>( BL, CL, L, p);
	        pointMul<<<blockNo , Tmax>>>( BR, CR, L, p);

        	list_stockham_dev(BL, no, i+1, winv, p);
	        list_stockham_dev(BR, no, i+1, winv, p);
	       
	        scalarMul<<<blockNo, Tmax>>>( BL, ninv, L, p);
	        scalarMul<<<blockNo, Tmax>>>( BR, ninv, L, p);

		blockNo = (sfixn)ceil((double)(L)/(double)(Tmax * 2.0));

	        cudaMalloc((void **)&check.fL[i+1], sizeof(sfixn)*(L/2));
	        cudaMalloc((void **)&check.fR[i+1], sizeof(sfixn)*(L/2));

		subtractLower<<<blockNo, Tmax>>>(check.fL[i+1], check.fL[i+2], BL, length*2, no/2, p );
		subtractLower<<<blockNo, Tmax>>>(check.fR[i+1], check.fR[i+2], BR, length*2, no/2, p );

		//We took M[i], invM[i], F[i+2] and produce F[i+1]
		if(flag != 1)
		{
			cudaFree(check.Ml[i]);
			cudaFree(check.Mr[i]);

			cudaFree(check.InvMl[i]);
			cudaFree(check.InvMr[i]);

			cudaFree(check.fL[i+2]);
			cudaFree(check.fR[i+2]);
		}
	}
	
	sfixn  threadsForAdiv, divInThreadBlock;
	for( i =  plainMulLimit-1; i >= 0; --i )
	{	
		 //We took M[i],  F[i+2] and produce F[i+1]
                length = (1L<<i) + 1;
                no     = 1L<<(k-i-1);
		threadsForAdiv = (sfixn)ceil((double)length/(double)W);
		divInThreadBlock = (sfixn)floor((double)Tdiv/(double)threadsForAdiv);
		blockNo = (sfixn)ceil((double)no/(double)divInThreadBlock);

		cudaMalloc((void **)&check.fL[i+1], sizeof(sfixn)*(L/2));
                cudaMalloc((void **)&check.fR[i+1], sizeof(sfixn)*(L/2));

		cudaMemcpy( check.fL[i+1], check.fL[i+2], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);
	        cudaMemcpy( check.fR[i+1], check.fR[i+2], sizeof(sfixn)*(L/2), cudaMemcpyDeviceToDevice);
						
		listPlainCudaDiv<<<blockNo,Tdiv>>>(check.Ml[i], check.fL[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);
		listPlainCudaDiv<<<blockNo,Tdiv>>>(check.Mr[i], check.fR[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);

		if(flag != 1)
		{
			cudaFree(check.Ml[i]);
			cudaFree(check.Mr[i]);

			cudaFree(check.InvMl[i]);
			cudaFree(check.InvMr[i]);

			cudaFree(check.fL[i+2]);
			cudaFree(check.fR[i+2]);
		}

	}


	cudaFree(AL);
	cudaFree(AR);
	cudaFree(BL);
	cudaFree(BR);
	cudaFree(CL);
        cudaFree(CR);
	cudaFree(DL);
        cudaFree(DR);	
}


/*! \fn struct PolyEvalSteps fastEvaluation(sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag)
<a name="evaCPU"> </a>
\brief This function creates subproduct tree and subinverse tree based on n points and then evaluate the polynomial on those points.
\param F: the polynomial need to be evaluated.
\param points: the list of points for which F need to be evaluated.
\param k: the height of subproduct tree.
\param p: prime number for finite field.
\param flag: what we want to  do for checking the result with MAPLE.
\return evaluation of the polynomial on the points. 
*/

struct PolyEvalSteps fastEvaluation(sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag)
{
	//declaring the struct PolyEvalSteps check which holds the intermediate results.
	

	sfixn numPoints = (1 << k);
	sfixn *pointsGPU;
	cudaMalloc((void **)&pointsGPU,   sizeof(sfixn)*numPoints); // allocating device memory for points
	cudaMalloc((void **)&check.Ml[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of left  subproduct tree
	cudaMalloc((void **)&check.Mr[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of right subproduct tree
	cudaMemcpy(pointsGPU, points, sizeof(sfixn)*numPoints, cudaMemcpyHostToDevice); //coping the points in device memory
	cudaThreadSynchronize();

	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Ml[0], pointsGPU, numPoints/2, 0); // creating leaves of left  subproduct tree
	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Mr[0], pointsGPU, numPoints/2, 1); // creating leaves of right subproduct tree
	cudaThreadSynchronize();
	
	subProductTree( k, p); // creating  subproduct tree
	subInvTree(k, p);      // creating  subinverse tree
	if(k > plainMulLimit)	
		fastEvaluation( F, k, p, flag); // evaluating the polynomial F on points
	else
                fastEvaluationLow( F, k, p, flag); // evaluating the polynomial F on points

	cudaThreadSynchronize();

	cudaFree(pointsGPU);
	return check;
}


/*! \fn struct PolyEvalSteps onlySubproductTree(sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag)
\brief This function returns the subproduct tree based on n points.
\param F the polynomial need to be evaluated.
\param points the list of points for which F need to be evaluated.
\param k the height of subproduct tree.
\param p prime number for finite field.
\param flag what we want to  do for checking the result with MAPLE.
\return subproduct tree
*/
struct PolyEvalSteps onlySubproductTree(sfixn *F, sfixn *points, sfixn k, sfixn p, sfixn flag)
{
	//declaring the struct PolyEvalSteps check which holds the intermediate results.
	

	sfixn numPoints = (1 << k);
	sfixn *pointsGPU;
	cudaMalloc((void **)&pointsGPU,   sizeof(sfixn)*numPoints); // allocating device memory for points
	cudaMalloc((void **)&check.Ml[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of left  subproduct tree
	cudaMalloc((void **)&check.Mr[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of right subproduct tree
	cudaMemcpy(pointsGPU, points, sizeof(sfixn)*numPoints, cudaMemcpyHostToDevice); //coping the points in device memory
	cudaThreadSynchronize();

	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Ml[0], pointsGPU, numPoints/2, 0); // creating leaves of left  subproduct tree
	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Mr[0], pointsGPU, numPoints/2, 1); // creating leaves of right subproduct tree
	cudaThreadSynchronize();
	
	subProductTree( k, p); // creating  subproduct tree
	
	cudaThreadSynchronize();

	cudaFree(pointsGPU);
	return check;
}

/*! \fn struct PolyEvalSteps onlyTrees(sfixn *points, sfixn k, sfixn p, sfixn flag)
\brief This function returns the subproduct tree and subinverse tree based on n points.
\param F the polynomial need to be evaluated.
\param points the list of points for which F need to be evaluated.
\param k the height of subproduct tree.
\param p prime number for finite field.
\param flag what we want to  do for checking the result with MAPLE.
\return subproduct tree and subinverse tree
*/
struct PolyEvalSteps onlyTrees(sfixn *points, sfixn k, sfixn p, sfixn flag)
{
	//declaring the struct PolyEvalSteps check which holds the intermediate results.
	

	sfixn numPoints = (1 << k);
	sfixn *pointsGPU;
	cudaMalloc((void **)&pointsGPU,   sizeof(sfixn)*numPoints); // allocating device memory for points
	cudaMalloc((void **)&check.Ml[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of left  subproduct tree
	cudaMalloc((void **)&check.Mr[0], sizeof(sfixn)*numPoints); // allocating device memory for leaves of right subproduct tree
	cudaMemcpy(pointsGPU, points, sizeof(sfixn)*numPoints, cudaMemcpyHostToDevice); //coping the points in device memory
	cudaThreadSynchronize();

	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Ml[0], pointsGPU, numPoints/2, 0); // creating leaves of left  subproduct tree
	leavesSubproductTree<<<(sfixn)ceil((double)( numPoints/2.0)/(double)Tmax) , Tmax>>>(check.Mr[0], pointsGPU, numPoints/2, 1); // creating leaves of right subproduct tree
	cudaThreadSynchronize();
	
	subProductTree( k, p); // creating  subproduct tree

	subInvTree(k, p);      // creating  subinverse tree

	cudaThreadSynchronize();

	cudaFree(pointsGPU);
	return check;
}





