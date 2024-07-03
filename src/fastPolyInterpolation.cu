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




#include "fastPolyInterpolation.h"
#include<iostream>
#include<fstream>

using namespace std;

static int SFIXN_LENGTH = sizeof(sfixn);

struct PolyEvalSteps subproductSteps;

sfixn fullArraySize;// will be assigned to = n * SFIXN_LENGTH
sfixn halfArraySize;// will be assigned to = n * SFIXN_LENGTH / 2
sfixn defaultBlocksDimension;// will be assigned to = (sfixn) ceil((double) n / (double) Tmax)
sfixn n;// will be assigned to = 1 << k;

struct PolyInterpolateSteps interpolateSteps; // this will hold all the intermediate results.

/**
 * Construct the interpolation tree from existing leaves up to level (k-1) using existing subproduct-tree.
 * We use plain multiplication up to level 8, and after that we use FFT for multiplications. 
 * 
 * Here every pair of operations in each level of the trees are depicted:
 * _________________________________I_____________________________________
 * |...................|    I_{i}   |   I_{i+1}    |.....................|
 *
 * _________________________________M_____________________________________
 * |...................|    M_{i}   |   M_{i+1}    |.....................|
 * 
 * _________________________________R_____________________________________
 * |...................|I_{i}*M_{i+1}+I_{i+1}*M_{i}|.....................|
 * 
 *
 @param k depth of sub product tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void linearCombination(sfixn k, sfixn p){

	sfixn *subproduct, *subproductShifted, *interpolated, *interpolatedShifted;
	cudaMalloc((void **) &subproduct, halfArraySize);
	cudaMalloc((void **) &subproductShifted, halfArraySize);
	cudaMalloc((void **) &interpolatedShifted, halfArraySize);
	cudaMalloc((void **) &interpolated, halfArraySize);

	sfixn polyCount; // Number of the polynomials in the trees
	sfixn subPolyLength; // length of each poly at (i-1)th level in the subproduct-tree

	sfixn w, winv, ninv; // w = primitive root of unity of 2^i
	
	sfixn L = 1L << (k - 1);// n / 2;

	sfixn i;
	for (i = 1; i < k; i++) {
		
		// Allocate memories for results in CUDA
		cudaMalloc((void **)&interpolateSteps.lefts[i], halfArraySize);
		cudaThreadSynchronize();
		cudaMalloc((void **)&interpolateSteps.rights[i], halfArraySize);
		cudaThreadSynchronize();

		if (i <= plainMulLimit) {// We do the plain multiplication

			subPolyLength = 1 << (i - 1); 
			polyCount = 1 << (k - i); 

			sfixn threadsForAnOperation = 2 * (subPolyLength + 1);// No. of threads required for an operation.
			sfixn operationsInABlock = (sfixn)floor((double)Tmax/(double)threadsForAnOperation);// No. operations done by 1 block
			
			// Compute the result in the left interpolation-tree.			
			plainCrossMultiplications<<<(sfixn)ceil((double)(polyCount/2)/(double)operationsInABlock), Tmax>>>(
							interpolateSteps.lefts[i], subproductSteps.Ml[i - 1], interpolateSteps.lefts[i - 1],
							polyCount, threadsForAnOperation, operationsInABlock, p);
			cudaThreadSynchronize();
			
			// Compute the result in the right interpolation-tree.			
			plainCrossMultiplications<<<(sfixn)ceil((double)(polyCount/2)/(double)operationsInABlock), Tmax>>>(
							interpolateSteps.rights[i], subproductSteps.Mr[i - 1], interpolateSteps.rights[i - 1],
							polyCount, threadsForAnOperation, operationsInABlock, p);
			cudaThreadSynchronize();
		
		} else {// We do the FFT-based multiplication
	
			subPolyLength = 1L << (i - 1);
			polyCount = 1L << (k - i - 1);

			w = primitive_root(i, p); 
			winv = inv_mod(w, p);
			ninv = inv_mod((1L << i), p);
	
			//************************************************************************
			//****** COMPUTING POLYNOMIAL IN THE LEFT INTERPOLATION TREE *************
			//************************************************************************
 
			sfixn shiftedSize = halfArraySize - SFIXN_LENGTH * subPolyLength;
			cudaMemcpy(subproduct, subproductSteps.Ml[i - 1], halfArraySize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(interpolated, interpolateSteps.lefts[i - 1], halfArraySize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(subproductShifted, &subproductSteps.Ml[i - 1][subPolyLength], shiftedSize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(interpolatedShifted, &interpolateSteps.lefts[i - 1][subPolyLength],  shiftedSize, cudaMemcpyDeviceToDevice);
		
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(2*Tmax)), Tmax>>>(subproduct, interpolatedShifted, L / 2, subPolyLength ); 
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(2*Tmax)), Tmax>>>(subproductShifted, interpolated, L / 2, subPolyLength ); 

			allZero<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.lefts[i], L);

			// Since we are not storing the leading coefficient in subproduct-tree's polynomials, we have to add interpolation polynomials...
			pointAdd2<<<(sfixn)ceil((double)(L/2)/(double)Tmax), Tmax>>>(interpolateSteps.lefts[i] , interpolateSteps.lefts[i-1], subPolyLength, L/2, p);

			//computing FFT
			list_stockham_dev(subproduct, polyCount, i, w, p);
			list_stockham_dev(interpolatedShifted, polyCount, i, w, p);
			list_stockham_dev(interpolated, polyCount, i, w, p);
			list_stockham_dev(subproductShifted, polyCount, i, w, p);

			//pointwise multiplication
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproduct, interpolatedShifted, L, p);
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproductShifted, interpolated, L, p);

			// inverse FFT
			list_stockham_dev(subproduct, polyCount, i, winv, p);
			list_stockham_dev(subproductShifted, polyCount, i, winv, p);

			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproduct, ninv, L, p);
			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproductShifted, ninv, L, p);

			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.lefts[i], subproduct, L, p);
			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.lefts[i], subproductShifted, L, p);	
		
			//************************************************************************
			//****** COMPUTING POLYNOMIAL IN THE RIGHT INTERPOLATION TREE ************
			//************************************************************************

			cudaMemcpy(subproduct, subproductSteps.Mr[i - 1], halfArraySize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(interpolated, interpolateSteps.rights[i - 1], halfArraySize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(subproductShifted, &subproductSteps.Mr[i - 1][subPolyLength], shiftedSize, cudaMemcpyDeviceToDevice);
			cudaMemcpy(interpolatedShifted, &interpolateSteps.rights[i - 1][subPolyLength], shiftedSize, cudaMemcpyDeviceToDevice);
			
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(2*Tmax)), Tmax>>>(subproduct, interpolatedShifted, L / 2, subPolyLength ); 
			zeroInbetween<<<(sfixn)ceil((double)L/(double)(2*Tmax)), Tmax>>>(subproductShifted, interpolated, L / 2, subPolyLength ); 

			allZero<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.rights[i], L);

			// Since we are not storing the leading coefficient in subproduct-tree's polynomials, we have to add interpolation polynomials...
			pointAdd2<<<(sfixn)ceil((double)(L/2)/(double)Tmax), Tmax>>>(interpolateSteps.rights[i] , interpolateSteps.rights[i-1], subPolyLength, L/2, p);
			
			// computing FFT
			list_stockham_dev(subproduct, polyCount, i, w, p);
			list_stockham_dev(interpolatedShifted, polyCount, i, w, p);
			list_stockham_dev(interpolated, polyCount, i, w, p);
			list_stockham_dev(subproductShifted, polyCount, i, w, p);

			//pointwise multiplication
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproduct, interpolatedShifted, L, p);
			pointMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproductShifted, interpolated, L, p);

			// inverse FFT
			list_stockham_dev(subproduct, polyCount, i, winv, p);
			list_stockham_dev(subproductShifted, polyCount, i, winv, p);

			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproduct, ninv, L, p);
			scalarMul<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(subproductShifted, ninv, L, p);
		
			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.rights[i], subproduct, L, p);
			pointAdd<<<(sfixn)ceil((double)L/(double)Tmax), Tmax>>>(interpolateSteps.rights[i], subproductShifted, L, p);

		}
	}

	//Free CUDA-allocated memories...
	cudaFree(subproduct);
	cudaFree(subproductShifted);
	cudaFree(interpolated);
	cudaFree(interpolatedShifted);
}

/**
 * Cross multiplication of two lists of polynomials (M & I) in plain mode.
 * Every operation are as this picture:
 * 
 * _________________________________I_____________________________________
 * |...................|    I_{i}   |   I_{i+1}    |.....................|
 *
 * _________________________________M_____________________________________
 * |...................|    M_{i}   |   M_{i+1}    |.....................|
 * 
 * _________________________________R_____________________________________
 * |...................|I_{i}*M_{i+1}+I_{i+1}*M_{i}|.....................|
 * 
 * First, we calculate the indexes and copy the required polynomials from global to the shared-memory.
 * Every thread is computing one coefficient of the result which involves at most polyLength multiplications
 * and additions. So the number of the threads required for an operation is the degree of the each 
 * result polynomial. 'threadsForAnOperation' specifies that which is equal to 2*(polyLength+1) in
 * which the 'polyLength' is the length of each polynomial in the subproduct-tree.
 * After computing the coefficient we wright it back to the correct index in the global memory.
 * Since we are storing polynomials to the shared memory, every result polynomial should be computed in one
 * block. But every block may be responsible for multiple operations. 'operationsInABlock' specifies the
 * number of operations which will be done by a single block.
 * 
 * These are some cheat-equation between the mentioned variables:
 *	i = level No.	k = log2(n)	#ThreadsPerBlock = 512
 * 	polyLength = 2 ^ (i - 1)
 *	polyCount = 2 ^ (k - i)
 *	threadsForAnOperation = 2 * (polyLength + 1)
 *	operationsInABlock = #ThreadsPerBlock/threadsForAnOperation
 *
 * Constraints with using this approach:
 *	Tmul >= 2* polyLength
 *	#Blocks >= polyCount/2
 *
 @param R result, a pointer to list of polynomials at (i+1)-th level in the interpolation tree.
 @param M pointer to list of polynomials at i-th level in the subproduct tree.
 @param I pointer to list of polynomials at i-th level in the interpolation tree.
 @param polyCount number of polynomial at the level.
 @param threadsForAnOperation no of threads required for one multiplication. threadsForAnOperation = 2 * polyLength
 @param operationsInABlock the no of multiplications in a thread block. It must be at least 1.
 @param p the prime number for finite field.
 */
__global__ void plainCrossMultiplications(sfixn *R, sfixn *M, sfixn *I, sfixn polyCount, sfixn threadsForAnOperation, sfixn operationsInABlock, sfixn p){
	
	//______Tmul/2__________Tmul/2___________Tmul__________ 
	//|mmmmmmmmmmmmmmm|iiiiiiiiiiiiiii|rrrrrrrrrrrrrrrrrrr|
	__shared__ sfixn shared[2 * Tmul]; // Tmul = Tmax = #ThreadsPerBlock = 512

	sfixn operationId = threadIdx.x / threadsForAnOperation + blockIdx.x * operationsInABlock;
	if(operationId < polyCount / 2 && threadIdx.x < threadsForAnOperation * operationsInABlock) { //Check if this is a valid operation!
		sfixn index = threadIdx.x % threadsForAnOperation; // The index that is going to be computed by a thread.
	
		sfixn polyLength = threadsForAnOperation / 2;

		// Computing start indexes in the shared memory...
		sfixn mStart = (threadIdx.x / threadsForAnOperation) * polyLength;
		sfixn iStart = mStart + Tmul / 2; // (threadIdx.x / threadsForAnOperation) * polyLength + Tmul / 2
		sfixn rStart = 2 * mStart + Tmul; // 2 * (threadIdx.x / threadsForAnOperation) * polyLength + Tmul

		// IDs of which polynomials we are multiplying...
		sfixn mPolyId = operationId * 2;
		sfixn iPolyId = mPolyId + 1; //operationId * 2 + 1;
	
		sfixn offset, a, b;
		
		// Copy the required polynomials from global memory to the shared memory
		if(index >= polyLength) {
			// This half of the threads are responsible for copying interpolation-tree's polynomials
			offset = (polyLength-1) * iPolyId;// Length of polynomial in interpolation-tree is 1 less. 		
		        if((index-polyLength)<polyLength-1 )
				shared[iStart + index - polyLength] = I[offset + index - polyLength];
			else
				shared[iStart + index - polyLength] = 0;
			b = index - polyLength + 1; a = polyLength - 1;
		} else {
			// This half of the threads are responsible for copying subproduct-tree's polynomials
			offset = polyLength * mPolyId;		
			shared[mStart + index] = M[offset + index];
			b = 0; a = index;
		}

		__syncthreads();

		// Here we do the multiplication... Compute only one coefficient of the result
		shared[rStart + index ] = 0;
		while(a >= 0 && b < (polyLength - 1)) {
			shared[rStart + index] = add_mod(shared[rStart + index], mul_mod(shared[mStart + a], shared[iStart + b], p), p);
			a--; b++;
    		}

		__syncthreads();
	
		// Adjust polynomials ID for the next multiplication...
		++mPolyId; --iPolyId;

		// Copy the required polynomials from global memory to the shared memory
		if(index >= polyLength) {
			// This half of the threads are responsible for copying interpolation-tree's polynomials
			offset = (polyLength-1) * iPolyId;		
			if((index-polyLength)<polyLength-1 )
				shared[iStart + index - polyLength] = I[offset + index - polyLength];
			else
				shared[iStart + index - polyLength] = 0;
			b = index - polyLength + 1; a = polyLength - 1;
		} else {
			// This half of the threads are responsible for copying subproduct-tree's polynomials
			offset = polyLength * mPolyId;		
			shared[mStart + index] = M[offset + index];
			b = 0; a = index;
		}

		__syncthreads();

		// Here we do the multiplication... Compute only one coefficient of the result
		while(a >= 0 && b < (polyLength - 1)) {
			shared[rStart + index] = add_mod(shared[rStart + index], mul_mod(shared[mStart + a], shared[iStart + b], p), p);
			a--; b++;
		}
		
		__syncthreads();
		
		// Copy the result from shared-memory back to the global memory...
		offset = operationId * (2 * polyLength - 2);
		if(index < (2 * polyLength - 2))
			R[offset + index] = shared[rStart + index];
	        	
		__syncthreads();
	}
}

/**
 * Compute intermediate coefficients for computing interpolation of the polynomial:
 * ci=vi*si in which vi=evaulate(ui) & si=1/(u1-ui)...(u_{i-1}-ui)(u_{i+1}-ui)...(un-ui).
 *
 @param c intermediate coefficients for computing interpolation: c[i]=v[i]*s[i]
 @param v result of evaluation steps of f.
 @param sinv result of evaluation for derivation of m=(x-u0)...(x-un-1)
 @param n dimension of the polynomial.
 @param p the prime number for finite field.
 */
__global__ void computeC(sfixn *c, sfixn* v, sfixn *sinv, sfixn n, sfixn p) {

	sfixn tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < n)
		c[tid] = mul_mod(v[tid], inv_mod(sinv[tid], p), p); //ci=vi*si in which si=1/sinv[i]=1/m(ui).

}

/**
 * Computes the derivation of a polynomial.
 *
 @param f polynomial's coefficients.
 @param n the polynomial's dimension.
 @param p the prime number for finite field.
 */
__global__ void derivate(sfixn *f, sfixn n, sfixn p) { 
	sfixn tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < n) {
		sfixn tmp = f[tid + 1];
		__syncthreads();
		f[tid] = mul_mod(tmp, tid + 1, p);
	}
}


/**
 * Set an element of an array with a value.
 *
 @param arr The array.
 @param idx ID of the element.
 @param val Value of the element to be set.
 */
__global__ void change_elem(sfixn *arr, sfixn idx, sfixn val) {
    arr[idx] = val;
}

/**
 * Add portion1 and portion2 to the leading coefficients of the R.
 *
 @param R the polynomial.
 @param portion1 first polynomial that is going to be added to the R.
 @param portion2 second polynomial that is going to be added to the R.
 @param n the portion polynomial's dimension.
 @param p the prime number for finite field.
 */
__global__ void addToLeading(sfixn *R, sfixn *portion1, sfixn *portion2, sfixn n, sfixn p){

        sfixn tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid < n)
			R[tid + n] = add_mod(portion1[tid], portion2[tid], p);
}


/**
 * Evaluate F on the points. 
 * It uses the existing subproduct & subinverse trees that are stored in subproductSteps data-structure
 * It should be called when k is less than plainMulLimit. (otherwise use fastEvaluationWithExistingTrees)
 *
 * This code is the same as the one implemented in fastPolyEvaluation (By Sardar Haque). We rewrite it because:
 * 1. We already computed subproduct & subinverse tree (stored in subproductSteps). 
 * 2. There was some problems and conflicts for using global variable called 'subproductSteps' in fastPolyEvaluation.
 * 3. Memory allocations are removed for some intermediate polynomials. ( such as fRs and fLs) 
 *
 @param F the polynomial which is going to be evaluated.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void fastEvaluationWithExistingTreesLow(sfixn *F, sfixn k, sfixn p) {
	sfixn length, no, threadsForAdiv, divInThreadBlock, blockNo;
	
	cudaMalloc((void **)&subproductSteps.fL[k+1], fullArraySize);
        cudaMalloc((void **)&subproductSteps.fR[k+1], fullArraySize);

	cudaMemcpy(subproductSteps.fL[k+1], F, fullArraySize, cudaMemcpyDeviceToDevice); 
	cudaMemcpy(subproductSteps.fR[k+1], F, fullArraySize, cudaMemcpyDeviceToDevice);	

	PlainCudaDiv<<<1,1>>>(subproductSteps.Ml[k-1], subproductSteps.fL[k+1], (n/2)+1, n, p);
	PlainCudaDiv<<<1,1>>>(subproductSteps.Mr[k-1], subproductSteps.fR[k+1], (n/2)+1, n, p);

	cudaMalloc((void **)&subproductSteps.fL[k], halfArraySize);
	cudaMalloc((void **)&subproductSteps.fR[k], halfArraySize);
	
	cudaMemcpy( subproductSteps.fL[k], subproductSteps.fL[k+1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy( subproductSteps.fR[k], subproductSteps.fR[k+1], halfArraySize, cudaMemcpyDeviceToDevice);

	for(sfixn i =  k-2; i >= 0; --i ) {	

		//We took M[i],  F[i+2] and produce F[i+1]
        	length = (1L<<i) + 1;
        	no     = 1L<<(k-i-1);
		threadsForAdiv = (sfixn)ceil((double)length/(double)W);
		divInThreadBlock = (sfixn)floor((double)Tdiv/(double)threadsForAdiv);


		blockNo = (sfixn)ceil((double)no/(double)divInThreadBlock);
		
		cudaMalloc((void **)&subproductSteps.fL[i+1], halfArraySize);
                cudaMalloc((void **)&subproductSteps.fR[i+1], halfArraySize);

		cudaMemcpy( subproductSteps.fL[i+1], subproductSteps.fL[i+2], halfArraySize, cudaMemcpyDeviceToDevice);
	    	cudaMemcpy( subproductSteps.fR[i+1], subproductSteps.fR[i+2], halfArraySize, cudaMemcpyDeviceToDevice);
						
		listPlainCudaDiv<<<blockNo,Tdiv>>>(subproductSteps.Ml[i], subproductSteps.fL[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);
		listPlainCudaDiv<<<blockNo,Tdiv>>>(subproductSteps.Mr[i], subproductSteps.fR[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);

		cudaFree(subproductSteps.fL[i+2]);
		cudaFree(subproductSteps.fR[i+2]);

	}

}

/**
 * Evaluate F on the points. 
 * It uses the existing subproduct & subinverse trees that are stored in subproductSteps data-structure
 * It should be called when k is greater than plainMulLimit. (otherwise use fastEvaluationWithExistingTreesLow)
 *
 * This code is the same as the one implemented in fastPolyEvaluation (By Sardar Haque). We rewrite it because:
 * 1. We already computed subproduct & subinverse tree (stored in subproductSteps). 
 * 2. There was some problems and conflicts for using global variable called 'subproductSteps' in fastPolyEvaluation.
 * 3. Memory allocations are removed for some intermediate polynomials. ( such as fRs and fLs) 
 *
 @param Fgpu the polynomial which is going to be evaluated.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void fastEvaluationWithExistingTrees(sfixn *Fgpu, sfixn k, sfixn p){
	sfixn i, ninv, blockNo, w, winv;
	sfixn *AL, *AR, *BL, *BR, *CL, *CR, *DL, *DR, *revFgpu;
	cudaMalloc((void **)&AL, fullArraySize);
	cudaMalloc((void **)&AR, fullArraySize);
	cudaMalloc((void **)&BL, fullArraySize);
	cudaMalloc((void **)&BR, fullArraySize);
	cudaMalloc((void **)&CL, fullArraySize);
	cudaMalloc((void **)&CR, fullArraySize);
	cudaMalloc((void **)&DL, halfArraySize);
	cudaMalloc((void **)&DR, halfArraySize);

	cudaMalloc((void **)&revFgpu, fullArraySize);

	EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.InvMl[k-1], AL, n/2, 1, n);
	EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.InvMr[k-1], AR, n/2, 1, n);

	EvalistReversePolyDec<<<defaultBlocksDimension, Tmax>>>(Fgpu, revFgpu, n, 1);
	
	w =  primitive_root(k, p);
	list_stockham_dev(revFgpu, 1, k, w, p);
	list_stockham_dev(AL,      1, k, w, p);
	list_stockham_dev(AR,      1, k, w, p);

	pointMul<<<defaultBlocksDimension , Tmax>>>( AL, revFgpu, n, p);
	pointMul<<<defaultBlocksDimension , Tmax>>>( AR, revFgpu, n, p);
	winv = inv_mod(w, p);

	list_stockham_dev(AL,   1, k, winv, p);
	list_stockham_dev(AR,   1, k, winv, p);
	ninv = inv_mod(n, p);
	scalarMul<<<defaultBlocksDimension, Tmax>>>( AL, ninv, n, p);
	scalarMul<<<defaultBlocksDimension, Tmax>>>( AR, ninv, n, p);

	EvalistZeroBetweenRev<<<defaultBlocksDimension, Tmax>>>(BL, AL, n, 1);
	EvalistZeroBetweenRev<<<defaultBlocksDimension, Tmax>>>(BR, AR, n, 1);


 	EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.Ml[k-1], CL, n/2, 1, n);
        EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.Mr[k-1], CR, n/2, 1, n);

        list_stockham_dev(CL, 1, k, w, p);
        list_stockham_dev(CR, 1, k, w, p);
        list_stockham_dev(BL, 1, k, w, p);
        list_stockham_dev(BR, 1, k, w, p);

        pointMul<<<defaultBlocksDimension , Tmax>>>( BL, CL, n, p);
        pointMul<<<defaultBlocksDimension , Tmax>>>( BR, CR, n, p);

        list_stockham_dev(BL, 1, k, winv, p);
        list_stockham_dev(BR, 1, k, winv, p);
        scalarMul<<<defaultBlocksDimension, Tmax>>>( BL, ninv, n, p);
        scalarMul<<<defaultBlocksDimension, Tmax>>>( BR, ninv, n, p);

	blockNo = (sfixn)ceil((double)(n)/(double)(Tmax * 2.0));

	cudaMalloc((void **)&subproductSteps.fL[k], halfArraySize);
	cudaMalloc((void **)&subproductSteps.fR[k], halfArraySize);

	subtractLowerSingle<<<blockNo, Tmax>>>(subproductSteps.fL[k] , Fgpu, BL, n/2, p);
	subtractLowerSingle<<<blockNo, Tmax>>>(subproductSteps.fR[k], Fgpu, BR, n/2, p);

        cudaFree(revFgpu);
	cudaFree(subproductSteps.InvMl[k-1]);
	cudaFree(subproductSteps.InvMr[k-1]);
			
	sfixn length, no;
	for(i = k-2; i >=  plainMulLimit; --i) {

		//We took M[i], invM[i], F[i+2] and produce F[i+1]
		length = 1L<<i;
		no     = 1L<<(k-i-1); 

		EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.InvMl[i], AL, length, no, length*2);
	        EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.InvMr[i], AR, length, no, length*2);

		EvalistReversePolyDec<<<blockNo, Tmax>>>(subproductSteps.fL[i+2], DL, length*2, no/2);
		EvalistReversePolyDec<<<blockNo, Tmax>>>(subproductSteps.fR[i+2], DR, length*2, no/2);		

		w =  primitive_root(i+1, p);
	        list_stockham_dev(DL, no/2, i+1, w, p);
	        list_stockham_dev(DR, no/2, i+1, w, p);	
	        list_stockham_dev(AL, no,   i+1, w, p);
        	list_stockham_dev(AR, no,   i+1, w, p);

		pointwiseMulHalf<<<defaultBlocksDimension , Tmax>>>( AL, DL, length*2, no, p);
                pointwiseMulHalf<<<defaultBlocksDimension , Tmax>>>( AR, DR, length*2, no, p);
	        winv = inv_mod(w, p);

        	list_stockham_dev(AL, no, i+1, winv, p);
	        list_stockham_dev(AR, no, i+1, winv, p);
        	ninv = inv_mod(length*2, p);
	        scalarMul<<<defaultBlocksDimension, Tmax>>>( AL, ninv, n, p);
        	scalarMul<<<defaultBlocksDimension, Tmax>>>( AR, ninv, n, p);

		EvalistZeroBetweenRev<<<defaultBlocksDimension, Tmax>>>(BL, AL, length*2, no);
	        EvalistZeroBetweenRev<<<defaultBlocksDimension, Tmax>>>(BR, AR, length*2, no);

	        EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.Ml[i], CL, length, no, length*2);
	        EvalistPolyDegInc<<<defaultBlocksDimension, Tmax>>>(subproductSteps.Mr[i], CR, length, no, length*2);
	
	        list_stockham_dev(CL, no, i+1, w, p);
        	list_stockham_dev(CR, no, i+1, w, p);
	        list_stockham_dev(BL, no, i+1, w, p);
	        list_stockham_dev(BR, no, i+1, w, p);

	        pointMul<<<defaultBlocksDimension , Tmax>>>( BL, CL, n, p);
	        pointMul<<<defaultBlocksDimension , Tmax>>>( BR, CR, n, p);

        	list_stockham_dev(BL, no, i+1, winv, p);
	        list_stockham_dev(BR, no, i+1, winv, p);
	       
	        scalarMul<<<defaultBlocksDimension, Tmax>>>( BL, ninv, n, p);
	        scalarMul<<<defaultBlocksDimension, Tmax>>>( BR, ninv, n, p);

	        cudaMalloc((void **)&subproductSteps.fL[i+1], halfArraySize);
	        cudaMalloc((void **)&subproductSteps.fR[i+1], halfArraySize);

		subtractLower<<<blockNo, Tmax>>>(subproductSteps.fL[i+1], subproductSteps.fL[i+2], BL, length*2, no/2, p );
		subtractLower<<<blockNo, Tmax>>>(subproductSteps.fR[i+1], subproductSteps.fR[i+2], BR, length*2, no/2, p );

		//We took M[i], invM[i], F[i+2] and produce F[i+1]
		cudaFree(subproductSteps.InvMl[i]);
		cudaFree(subproductSteps.InvMr[i]);

		cudaFree(subproductSteps.fL[i+2]);
		cudaFree(subproductSteps.fR[i+2]);
	}
	
	sfixn  threadsForAdiv, divInThreadBlock;
	for( i =  plainMulLimit-1; i >= 0; --i ) {
	
		 //We took M[i],  F[i+2] and produce F[i+1]
                length = (1L<<i) + 1;
                no     = 1L<<(k-i-1);
		threadsForAdiv = (sfixn)ceil((double)length/(double)W);
		divInThreadBlock = (sfixn)floor((double)Tdiv/(double)threadsForAdiv);
		blockNo = (sfixn)ceil((double)no/(double)divInThreadBlock);

		cudaMalloc((void **)&subproductSteps.fL[i+1], halfArraySize);
                cudaMalloc((void **)&subproductSteps.fR[i+1], halfArraySize);

		cudaMemcpy( subproductSteps.fL[i+1], subproductSteps.fL[i+2], halfArraySize, cudaMemcpyDeviceToDevice);
	        cudaMemcpy( subproductSteps.fR[i+1], subproductSteps.fR[i+2], halfArraySize, cudaMemcpyDeviceToDevice);
						
		listPlainCudaDiv<<<blockNo,Tdiv>>>(subproductSteps.Ml[i], subproductSteps.fL[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);
		listPlainCudaDiv<<<blockNo,Tdiv>>>(subproductSteps.Mr[i], subproductSteps.fR[i+1], 0, length, threadsForAdiv, divInThreadBlock, no, p);

		cudaFree(subproductSteps.InvMl[i]);
		cudaFree(subproductSteps.InvMr[i]);

		cudaFree(subproductSteps.fL[i+2]);
		cudaFree(subproductSteps.fR[i+2]);

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

/**
 * Evaluate F on the points using the existing subproduct & subinverse trees.
 *
 @param F the polynomial which is going to be evaluated.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void fastEvaluate(sfixn *F, sfixn k, sfixn p){
	if(k > plainMulLimit)	
		fastEvaluationWithExistingTrees(F, k, p);
	else
        	fastEvaluationWithExistingTreesLow(F, k, p);
}

/**
 * Compute m=(x-u1)(x-u2)...(x-un) by multiplying the roots of left & right subproduct trees. 
 * This function is used only for small K (Less than 8). (It uses plain multiplications)
 *
 @param m the pointer of which the m is going to be stored in.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void generateMForSmallK(sfixn *m, sfixn k, sfixn p){
	
	// m contains subproductSteps.Ml[k-1] & subproductSteps.Mr[k-1]
	sfixn *mTmp;
	cudaMalloc((void **) &mTmp, fullArraySize + 2 * SFIXN_LENGTH);
	cudaMemcpy(mTmp, subproductSteps.Ml[k-1], halfArraySize + SFIXN_LENGTH, cudaMemcpyDeviceToDevice);
	cudaMemcpy(&mTmp[n/2+1], subproductSteps.Mr[k-1], halfArraySize + SFIXN_LENGTH, cudaMemcpyDeviceToDevice);	

	sfixn subPolyLength = (1 << (k - 1)) + 1;
	sfixn polyCount = 2;

	sfixn threadsForAmul = 2 * subPolyLength;
	sfixn mulInThreadBlock = (sfixn)floor((double)Tmax/(double)threadsForAmul);

	listPlainMulGpu<<<(sfixn)ceil((double)(polyCount/2)/(double)mulInThreadBlock), Tmul>>>(
						mTmp, m, subPolyLength, 2, threadsForAmul, mulInThreadBlock, p);

	cudaThreadSynchronize();

	cudaFree(mTmp);
}

/**
 * Compute m=(x-u1)(x-u2)...(x-un) by multiplying the roots of left & right subproduct trees. 
 * This function is used only for Big K (more than 8). (It uses FFT-based multiplications)
 *
 @param m the pointer that the m is going to be stored in.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void generateMForBigK(sfixn *m, sfixn k, sfixn p){
	sfixn *ml, *mr;

	cudaMalloc((void **) &ml, fullArraySize);
	cudaMalloc((void **) &mr, fullArraySize);
	
	allZero<<<defaultBlocksDimension, Tmax>>>(m, n);
	allZero<<<defaultBlocksDimension, Tmax>>>(ml, n);
	allZero<<<defaultBlocksDimension, Tmax>>>(mr, n);

	cudaMemcpy(ml, subproductSteps.Ml[k-1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(mr, subproductSteps.Mr[k-1], halfArraySize, cudaMemcpyDeviceToDevice);

	addToLeading<<<(sfixn) ceil((double) (n/2) / (double) Tmax), Tmax>>>(
			m, subproductSteps.Ml[k-1], subproductSteps.Mr[k-1], n/2, p);

	sfixn w = primitive_root(k, p);

	list_stockham_dev(ml, 1, k, w, p);
	list_stockham_dev(mr, 1, k, w, p);

	pointMul<<<defaultBlocksDimension, Tmax>>>(ml, mr, n, p);

	sfixn winv = inv_mod(w, p);
	list_stockham_dev(ml, 1, k, winv, p);

	w = 1L << k;

	sfixn ninv = inv_mod(w, p);
	scalarMul<<<defaultBlocksDimension, Tmax>>>(ml, ninv, n, p);

	pointAdd<<<defaultBlocksDimension, Tmax>>>(m, ml, n, p);

	cudaFree(ml); 
	cudaFree(mr);
}

/**
 * Compute m=(x-u1)(x-u2)...(x-un) by multiplying the roots of left & right subproduct trees. 
 *
 @param m the pointer that the m is going to be stored in.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void generateM(sfixn *m, sfixn k, sfixn p) {
	if(k > plainMulLimit)	
		generateMForBigK(m, k, p);
	else
        	generateMForSmallK(m, k, p);
}

/**
 * Computes c[i]s (ci=vi*si). C[i]s are the as the coefficients in lagrange interpolation.
 * Later we initialise leaves of interpolation-tree by c[i]s...
 *
 @param c results (the lagrange coefficients).
 @param v evaluated results for corresponding points.
 @param k the number of levels, which is equal to log(n).
 @param p the prime number for finite field.
 */
void computeCs(sfixn *c, sfixn *v, sfixn k, sfixn p) {
	
	// Compute m = subproductSteps.Ml[0] * subproductSteps.Mr[0]
	// We store m in c (c=m)
	generateM(c, k, p);
	
	//c = mPrime (derivation of m=(x-u0)...(x-un))
	derivate<<<defaultBlocksDimension, Tmax>>>(c, n - 1, p);

	// Evaluate mPrime (the derivation of m) at the evaluation points, 
	change_elem<<<1,1>>>(c, n - 1, n);// Since the last coefficient is 1, we assign the degree to the last derivation's coefficient

	// We use this function since we already built the sub-product & sub-inverse trees.
	fastEvaluate(c, k, p);

	sfixn *sinv;// sinv contains the evaluation results of mPrime
	cudaMalloc((void**) &sinv, fullArraySize);
	cudaMemcpy(sinv, subproductSteps.fL[1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(&sinv[n / 2], subproductSteps.fR[1], halfArraySize, cudaMemcpyDeviceToDevice);

	// Call the kernel in which we compute every ci.
	computeC<<<defaultBlocksDimension, Tmax>>>(c, v, sinv, n, p);
	
	cudaFree(sinv);// Free sinv, since we have computed c we don't need them anymore
}

/**
 * Compute interpolation result: I=(I(k-1)left) * (M(k-1)right) + (I(k-1)right) * (M(k-1)left)
 * This function is used only for Small K (less than 8). (It uses plain multiplications)
 *
 @param I Interpolation Result.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void computeIForSmallK(sfixn *I, sfixn k, sfixn p){
	sfixn *m, *i;

	// m contains subproductSteps.Ml[k-1] & subproductSteps.Mr[k-1]
	cudaMalloc((void **) &m, fullArraySize + 2 * SFIXN_LENGTH);
	cudaMemcpy(m, subproductSteps.Ml[k-1], halfArraySize + SFIXN_LENGTH, cudaMemcpyDeviceToDevice);
	cudaMemcpy(&m[n/2+1], subproductSteps.Mr[k-1], halfArraySize + SFIXN_LENGTH, cudaMemcpyDeviceToDevice);	

	// i contains interpolateSteps.lefts[k-1] & interpolateSteps.rights[k-1]
	cudaMalloc((void **) &i, fullArraySize);
	cudaMemcpy(i, interpolateSteps.lefts[k-1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(&i[n/2], interpolateSteps.rights[k-1], halfArraySize, cudaMemcpyDeviceToDevice);	

	sfixn subPolyLength = 1 << (k - 1);
	sfixn polyCount = 2;

	sfixn threadsForAnOperation = 2 * (subPolyLength + 1);
	sfixn operationsInABlock = (sfixn)floor((double)Tmax/(double)threadsForAnOperation);

	plainCrossMultiplications<<<(sfixn)ceil((double)(polyCount/2)/(double)operationsInABlock), Tmax>>>(
								I, m, i, polyCount, threadsForAnOperation, operationsInABlock, p);
	cudaThreadSynchronize();

	cudaFree(m);
	cudaFree(i);
}

/**
 * Compute interpolation result: I=(I(k-1)left) * (M(k-1)right) + (I(k-1)right) * (M(k-1)left)
 * This function is used only for Big K (more than 8). (It uses FFT-based multiplications)
 *
 @param I Interpolation Result.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 */
void computeIForBigK(sfixn *I, sfixn k, sfixn p){
	sfixn *ml, *mr;

	cudaMalloc((void **) &ml, fullArraySize);
	cudaMalloc((void **) &mr, fullArraySize);

	allZero<<<defaultBlocksDimension, Tmax>>>(I, n);
	allZero<<<defaultBlocksDimension, Tmax>>>(ml, n);
	allZero<<<defaultBlocksDimension, Tmax>>>(mr, n);

	cudaMemcpy(ml, subproductSteps.Ml[k-1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(mr, interpolateSteps.rights[k-1], halfArraySize, cudaMemcpyDeviceToDevice);

	addToLeading<<<(sfixn) ceil((double) (n/2) / (double) Tmax), Tmax>>>(
				I, interpolateSteps.rights[k-1], interpolateSteps.lefts[k-1], n/2, p);

	sfixn w = primitive_root(k, p);

	list_stockham_dev(ml, 1, k, w, p);
	list_stockham_dev(mr, 1, k, w, p);

	pointMul<<<defaultBlocksDimension, Tmax>>>(ml, mr, n, p);

	sfixn winv = inv_mod(w, p);
	list_stockham_dev(ml, 1, k, winv, p);

	w = 1L << k;

	sfixn ninv = inv_mod(w, p);
	scalarMul<<<defaultBlocksDimension, Tmax>>>(ml, ninv, n, p);

	pointAdd<<<defaultBlocksDimension, Tmax>>>(I, ml, n, p);

	allZero<<<defaultBlocksDimension, Tmax>>>(ml, n);
	allZero<<<defaultBlocksDimension, Tmax>>>(mr, n);

	cudaMemcpy(ml, subproductSteps.Mr[k-1], halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(mr, interpolateSteps.lefts[k-1], halfArraySize, cudaMemcpyDeviceToDevice);

	w = primitive_root(k, p);

	list_stockham_dev(ml, 1, k, w, p);
	list_stockham_dev(mr, 1, k, w, p);

	pointMul<<<defaultBlocksDimension, Tmax>>>(ml, mr, n, p);

	winv = inv_mod(w, p);
	list_stockham_dev(ml, 1, k, winv, p);

	w = 1L << k;

	ninv = inv_mod(w, p);
	scalarMul<<<defaultBlocksDimension, Tmax>>>(ml, ninv, n, p);

	pointAdd<<<defaultBlocksDimension, Tmax>>>(I, ml, n, p);

	cudaFree(ml);
	cudaFree(mr);
}

/**
 * Compute interpolation result: I=(I(k-1)left) * (M(k-1)right) + (I(k-1)right) * (M(k-1)left)
 *
 @param I Interpolation Result.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 @param flag set to 1 to not delete the intermedate results.
 */
void computeI(sfixn *I, sfixn k, sfixn p, sfixn flag){
	
	if(k > plainMulLimit)	
		computeIForBigK(I, k, p);
	else
        	computeIForSmallK(I, k, p);

	cudaFree(subproductSteps.Ml[k-1]);
	cudaFree(subproductSteps.Mr[k-1]);

	// If we are not comparing the results we free allocated memories...	
	if(flag != 1){
		cudaFree(interpolateSteps.lefts[k-1]);
		cudaFree(interpolateSteps.rights[k-1]);
	}
}

/**
 *
 * Generate the unique polynomial with existing subproduct tree and leaves of
 * the interpolation tree called c_{i} = v_{i}*s_{i}, 
 *     si=1/(u1-ui)...(u_{i-1}-ui)(u_{i+1}-ui)...(un-ui)
 * It will construct the interpolation tree, then generate the final polynomial
 *
 @param c ci, leaves of the interpolation tree. (lagrange coefficients)
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 @param flag set to 1 to not delete the intermedate results.
 */
struct PolyInterpolateSteps interpolate(sfixn *c, sfixn k, sfixn p, sfixn flag) {
	
	//Copy c_{i}s to the leaves of the interpolateTree
	cudaMalloc((void **) &interpolateSteps.lefts[0], halfArraySize);
	cudaMalloc((void **) &interpolateSteps.rights[0], halfArraySize);
	cudaMemcpy(interpolateSteps.lefts[0], c, halfArraySize, cudaMemcpyDeviceToDevice);
	cudaMemcpy(interpolateSteps.rights[0], &c[n / 2], halfArraySize, cudaMemcpyDeviceToDevice);

	// Compute the polynomial interpolation by constructing the interpolation-tree
	linearCombination(k, p);

	// Free unrequired CUDA memories...
	int i;
	for(i=0; i<(k-1); i++){
		cudaFree(subproductSteps.Ml[i]);
		cudaFree(subproductSteps.Mr[i]);
	}	

	// If we are not comparing the results we free allocated memories...
	if(flag != 1){
		for(i=0; i<(k-1); i++){
			cudaFree(interpolateSteps.lefts[i]);
			cudaFree(interpolateSteps.rights[i]);
		}
		cudaFree(interpolateSteps.c);
	}

	cudaMalloc((void **) &interpolateSteps.I, fullArraySize);
	computeI(interpolateSteps.I, k, p, flag);
	return interpolateSteps;
}

/**
 *
 * Generate the unique polynomial in which the evaluated values of the polynomial on the
 * points (u) are unique values 'v'. (generate f such that f(u_{i}) = v_{i})
 *
 @param v evaluated values at the points.
 @param u evaluation points.
 @param k depth of the tree which is equal to log(n).
 @param p the prime number for finite field.
 @param flag set to 1 to not delete the intermedate results.
 */
struct PolyInterpolateSteps fastInterpolate(sfixn *v, sfixn *u, sfixn k, sfixn p, sfixn flag) {

	// Setting global variables...
	n = 1 << k;
	defaultBlocksDimension = (sfixn) ceil((double) n / (double) Tmax);
	fullArraySize = n * SFIXN_LENGTH;
	halfArraySize = fullArraySize / 2;

	subproductSteps = onlyTrees(u, k, p, flag);

	sfixn *vGPU; // v contains the evaluation results. v_{i}=f(u_{i})
	cudaMalloc((void**) &vGPU, fullArraySize);
	cudaMemcpy(vGPU, v, fullArraySize, cudaMemcpyHostToDevice);

	//compute ci=vi*si in which vi=evaluate(ui) & si=1/(u1-ui)...(u_{i-1}-ui)(u_{i+1}-ui)...(un-ui).
	cudaMalloc((void**) &interpolateSteps.c, fullArraySize);
	computeCs(interpolateSteps.c, vGPU, k, p);
	
	cudaFree(vGPU);// Free v, since we have computed c we don't need them anymore

	return interpolate(interpolateSteps.c, k, p, flag);
}
