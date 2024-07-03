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


#include<iostream>
#include <ctime>
#include<cmath>
#include<fstream>
#include "fastPolyInterpolation.h"
#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "printing.h"
#include "cudautils.h"
#include "types.h"

using namespace std;

struct PolyEvalSteps subproductSteps;

void generateRandomNumbers(sfixn *a, sfixn n, sfixn p){
	for(sfixn i = 0; i < n; ++i)
		a[i] = rand()%p;
}

void generateDifferentNumbers(sfixn *a, sfixn n, sfixn p){
	for(sfixn i = 0; i < n; ++i)
		a[i] = i % p;
}

void generateDifferentRandomNumbers(sfixn *a, sfixn n, sfixn p){
	for(sfixn i = 0; i < n; ++i){
		sfixn tmp = rand() % p;
		for(sfixn j=0; j<i; j++){
			if(tmp==a[j]){
				tmp = rand() % p;
				j=0;				
			}
		}
		a[i] = tmp;
	}
}

void readFromFile(sfixn *a, sfixn n, char *name) {
	ifstream ifs( name, ifstream::in  ); 
	for(int i = 0; i < n; ++i)
		ifs>>a[i];
	ifs.close();
}

void writeToFileGPU(sfixn *f, char *name, int k, sfixn n) {
	ofstream ofs( name, ofstream::out  );

	sfixn *fGPU = new sfixn[n];
	cudaMemcpy(fGPU, f, sizeof(sfixn)*n, cudaMemcpyDeviceToHost);

	int j;
	for(j = 0; j < n; ++j)
		ofs<<fGPU[j]<<" ";

	delete [] fGPU;
	ofs.close();
}


void writeToFileGPU(sfixn *f, char *name, int k) {
	ofstream ofs( name, ofstream::out  );

	sfixn n = 1 << k;		
	sfixn *fGPU = new sfixn[n];
	cudaMemcpy(fGPU, f, sizeof(sfixn)*n, cudaMemcpyDeviceToHost);

	int j;
	for(j = 0; j < n; ++j) 
		ofs<<fGPU[j]<<" ";

	delete [] fGPU;
	ofs.close();	
}

void writeToFileCPU(sfixn *f, char *name, int k) {
	ofstream ofs( name, ofstream::out  );

	sfixn n = 1 << k;		
	int j;
	for(j = 0; j < n; ++j) 
		ofs<<f[j]<<" ";

	ofs.close();	
}

void writeInterpolationTree(struct PolyInterpolateSteps interpolationSteps, char *name, int k){
	ofstream ofs( name, ofstream::out);
	
	int i, j, l, length, no;
	length = 2;
	no = 1 << (k-2);		
	sfixn *T = new sfixn[length*no ];
	for(i = 1; i < k; ++i )
	{
		cudaMemcpy(T, interpolationSteps.lefts[i], sizeof(sfixn)*(no) * length, cudaMemcpyDeviceToHost);
		for(j = 0; j < no; ++j)
		{
			for(l = 0; l < length; ++l)
			{
				ofs<<T[j*length + l]<<" ";
			}
			ofs<<endl;
		}

		cudaMemcpy(T, interpolationSteps.rights[i], sizeof(sfixn)*(no)*length, cudaMemcpyDeviceToHost);
		for(j = 0; j < no; ++j)
		{
			for(l = 0; l < length; ++l)
			{
				ofs<<T[j*length + l]<<" ";
			}
			ofs<<endl;
		}	
		
		no = no/2;
		length = 2*length; 
		
	}
	delete [] T;
	ofs.close();	
}

int main(int argc, char *argv[]) {
	int k = 10;
	sfixn p = 469762049;
	int verify = 0;
	if (argc > 1) k = atoi(argv[1]);
	if (argc > 2) p = atoi(argv[2]);
	if (argc > 3) verify = atoi(argv[3]);

	// Writing the value of k and p into a file	
	ofstream ofs1( "KP.dat", ofstream::out  );	
	ofs1<<k<<" "<<p;
	ofs1.close();
	
	sfixn n = 1 << k;

	sfixn *points = new sfixn[n];
	sfixn *evaled = new sfixn[n];

	struct PolyInterpolateSteps interpolateSteps;	

	if(verify == 1){

		system("maple -q fastInterpolation.mm");

		readFromFile(evaled, n, "value.dat"); // getting the polynomial F to be evaluated.
		readFromFile(points, n, "Points.dat"); // getting the values of x for which we need to evaluate F.

		cudaEvent_t start, stop;
                cudaEventCreate(&start);
        cudaEventCreate(&stop);
                cudaEventRecord(start, 0);

		interpolateSteps = fastInterpolate(evaled, points, k, p, 1);
	
		cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float outerTime;
        cudaEventElapsedTime(&outerTime, start, stop);

                cout<<k<<"  "<<outerTime/1000.0<<endl;

		writeToFileGPU(interpolateSteps.I, "final.dat", k);
		writeToFileGPU(interpolateSteps.c, "cGPU.dat", k);
		writeInterpolationTree(interpolateSteps, "IGPU.dat", k);

		int testI = system("diff -b PolyI.dat IGPU.dat");
	 	int testF = system("diff -b finalMaple.dat final.dat");
		int testC = system("diff -b c.dat cGPU.dat");
		if(testC!= 0) cout<<"Coefficients are not correct..."<<endl;
		else if(testI!= 0) cout<<"Interpolation-trees are not equivalant..."<<endl;
		else if(testF!=0) cout<<"Final result is not Correct !!!"<<endl;

		int i;
		for(i=0; i<k; i++){
			cudaFree(interpolateSteps.lefts[i]);
			cudaFree(interpolateSteps.rights[i]);
		}
		cudaFree(interpolateSteps.c);

	} else{
		generateRandomNumbers(evaled, n, p);
		generateDifferentNumbers(points, n, p);
		//generateDifferentRandomNumbers(points, n, p);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
        cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		interpolateSteps = fastInterpolate(evaled, points, k, p, 0);

		cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float outerTime;
        cudaEventElapsedTime(&outerTime, start, stop);

		cout<<k<<"  "<<outerTime/1000.0<<endl;

		writeToFileGPU(interpolateSteps.I, "final.dat", k);
	}

	// Free the allocated CUDA memory for the result (whether we are verifying the outputs or not)
	cudaFree(interpolateSteps.I);
	
	delete [] points;
	delete [] evaled;
	return 0;
}
