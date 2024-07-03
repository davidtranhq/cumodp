/*******************************************************************
* This code is intended for for computing remainder i.e f rem m    *
* The degree of m is T*W-1 that is 511 and the degree of f is 1020
********************************************************************/
#include "list_plain_division.h"



using namespace std;


void divCPU(sfixn *M, sfixn *F, sfixn n, sfixn length, sfixn start, sfixn P,sfixn sn);


__global__ void list_divCUDA(sfixn *M, sfixn *F, sfixn start, sfixn length, sfixn threadsPerDiv, sfixn DivPerBlock, sfixn polyNum, sfixn P)
{
	__shared__ sfixn sM[T*W];
	__shared__ sfixn sF[T*2*W];
	__shared__ sfixn invM[T];
	__shared__ sfixn oneF[T];
	__shared__ sfixn mID[T];
	__shared__ sfixn fID[T];
	__shared__ sfixn startF[T];
	
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
void list_div(sfixn *M, sfixn *F, sfixn n, sfixn m, sfixn  start_offset, sfixn length_poly, sfixn poly_on_layer, sfixn p )
{

	/**************************************************
	* Number of threads responsible for one division  *
	**************************************************/	
	sfixn threadsForAdiv = (sfixn)ceil((double)length_poly / (double)W);

	/******************************************
	* Number of divisions in one thread block *
	******************************************/
	sfixn divInThreadBlock = (sfixn)floor((double)T/(double)threadsForAdiv);

	/****************************
	* Number of blocks required *
	****************************/
	sfixn blockNo = (sfixn)ceil((double)poly_on_layer/(double) divInThreadBlock );


	//cout<<"threadsForAdiv: "<<threadsForAdiv<<" divInThreadBlock: "<<divInThreadBlock<<" blockNo: " <<blockNo<<endl;

	sfixn *Mgpu, *Fgpu;

	cudaMalloc((void **)&Mgpu, sizeof(sfixn)*n);
	cudaMalloc((void **)&Fgpu, sizeof(sfixn)*m);
    cudaMemcpy(Mgpu, M, sizeof(sfixn)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(Fgpu, F, sizeof(sfixn)*m, cudaMemcpyHostToDevice);
	
	sfixn Ftemp[3000];
	list_divCUDA<<<blockNo, T>>>(Mgpu, Fgpu, start_offset, length_poly, threadsForAdiv, divInThreadBlock, poly_on_layer, p);	

	cudaMemcpy(Ftemp, Fgpu, sizeof(sfixn)*3000, cudaMemcpyDeviceToHost);
	//cout<<endl<<"F: "<<endl;

	sfixn i, j = 0; 

	sfixn *FF=(sfixn *)malloc(length_poly*sizeof(sfixn));
    for(j=0; j < poly_on_layer; ++j) // for(j=0; j < threadsForAdiv; ++j)
	{
		//cout<<endl;
		for(i = 0; i < (length_poly - 1); ++i)//for(i = 0; i < 18; ++i)
			FF[i]=Ftemp[j*(length_poly-1)+i];
			//cout<<Ftemp[k++]<<" ";	
		//printf("Rem_CUDA[%d] := ",j);print_poly(length_poly-2, FF, 'x');printf(";\n");
		
	}
	free(FF);
	//cout<<endl;
	
	for(j=0; j < poly_on_layer; ++j)
	{
		divCPU(M, F, j, length_poly, start_offset,p,j);
	}
		
	

	cudaFree(Mgpu);
    cudaFree(Fgpu);
}
///////////////////////////////////////////////////////////////////////////
//BEGIN:list_plain_division_tst
///////////////////////////////////////////////////////////////////////////

void divCPU(sfixn *M, sfixn *F, sfixn n, sfixn length, sfixn start, sfixn P,sfixn sn)
{
	//cout<<"r["<<sn<<"] := Rem(";

	//sfixn t;
	sfixn k = ((n/2) + 1)*(2*length-2) -1;
	sfixn l = k - (2*length -2);
	k = k+1;
	l = l+1;
	
	//cout<<F[l];
	l=l+1;
	//t =1;
	//for(; l<k; ++l,++t)
	//	cout<<"+"<<F[l]<<"*x^"<<t;
	//cout<<",";


	sfixn i = start + (n+1)*length -1;
	sfixn j = i - length;
	j = j +1;
	i = i+1;

	
	//cout<<M[j];
	j = j +1;

	//t=1;
	
}

void list_plain_division_tst(sfixn start_offset, sfixn length_poly, sfixn poly_on_layer, sfixn p)
{
#ifndef _mcompile_
	/**************************************
	* poly_on_layer MUST BE DIVISIBLE BY 2 *
	***************************************/
	sfixn n = 3000, m = 3000, i;
	
	sfixn *M = new sfixn[n];
	sfixn *F = new sfixn[m];
	
	//cout<<endl;
	for(i = 0; i < n; ++i)
	{	
		M[i] = rand()% p;
		if(i % length_poly == 0) cout<<endl;
		if( (i % length_poly) == (length_poly-1) ) 
		{
			if(M[i] == 0)
			{
				while(M[i] ==0)
					M[i] = rand()% p;
			
			}
		}
		cout<<M[i]<<" ";
	}
	cout<<endl;		
	for(i = 0; i < m; ++i)	
	{
		F[i] = rand()% p;
		if(i % (2*length_poly -2) == 0) cout<<endl;
		cout<<F[i]<<" ";
	}
	cout<<"done"<<endl;
		
	list_div(M ,F, n, m, start_offset, length_poly, poly_on_layer, p );
	
//double multiplications =
//(double)(poly_on_layer)*(double)(length_poly)*((double)(length_poly)
//-1.0);
//double additions =
//(double)(poly_on_layer)*(double)(length_poly)*((double)(length_poly)
//-1.0);
//cout<<"Multiplications: "<<multiplications<<endl;
//cout<<"additions: "<<additions<<endl;


	delete [] M;
	delete [] F;
#endif	
}
///////////////////////////////////////////////////////////////////////////
//END:list_plain_division_tst
///////////////////////////////////////////////////////////////////////////
/*
int main(int argc, char *argv[])
{
	int   p= 7;
	int start_offset = 0, length_poly = 4, poly_on_layer = 4;
	

	if (argc > 1) start_offset = atoi(argv[1]);
	if (argc > 2) length_poly = atoi(argv[2]);
    if (argc > 3) poly_on_layer = atoi(argv[3]);
	if (argc > 4) p = atoi(argv[4]);
	
	list_plain_division_tst(start_offset, length_poly, poly_on_layer,p);
	
	return 0;
}
*/	
	
