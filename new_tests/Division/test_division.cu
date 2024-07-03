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
#include "cudaDiv.h"
#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "printing.h"
#include "cudautils.h"
#include "types.h"


using namespace std;

/*
It should be called in the following sequence

the degree of bigger polynomial, n-1.
the degree of smaller polynomial, m-1. n >= m.
prime number for prime field, p.
whatUwant = 1 for checking with maple.

*/		


int main(int argc, char *argv[])
{
	int n = 8192, m = 1024, i, whatUwant = 0;
	sfixn p = 2147483647;

	 if(argc > 1) n = atoi(argv[1]);
	 if(argc > 2) m = atoi(argv[2]);
	 if(argc > 3) p = atoi(argv[3]);
	 if(argc > 4) whatUwant = atoi(argv[4]);


	sfixn *A = new sfixn[n];
	sfixn *B = new sfixn[m];

	if(whatUwant == 1)
	{
		ofstream ofs1( "PNM.dat", ofstream::out  );	
		ofs1<<p<<" "<<n<<" "<<m;
	        ofs1.close();
		i = system("maple -q maple_division.mm");

		
		ifstream ifs1( "Poly1.dat", ifstream::in  ); 
		for(i = 0; i < n; ++i)	
			ifs1>>A[i];
		ifs1.close();
		
		ifstream ifs2( "Poly2.dat", ifstream::in  ); 
		for(i = 0; i < m; ++i)	
			ifs2>>B[i];
		ifs2.close();
		

	}
	else
	{		
		for(i = 0; i < n; ++i)	
			A[i] =  rand()% p;
		for(i = 0; i < m; ++i)	
			B[i] =  rand()% p;
	}
	

	
	
	
	sfixn *Q = new sfixn [ n-m +1 ];
	sfixn *R = new sfixn [ m - 1 ];
	float outerTime = divPrimeField(A, B, R, Q, n, m, p);
	
	if(whatUwant == 1)
	{
		int rError, qError;
		ofstream ofsQ( "QuoGPU.dat", ofstream::out  );
		for(i = 0;  i < (n - m + 1); ++i)
			ofsQ<<Q[i]<<" ";
        	ofsQ.close();

	
		ofstream ofs( "RemGPU.dat", ofstream::out  );
		for(i = 0;  i < (m - 1); ++i)
			ofs<<R[i]<<" ";
        	ofs.close();

		qError = system("diff -b QuoGPU.dat Quo.dat");
		if(qError != 0) cout<<"ERROR"<<endl;

		rError = system("diff -b RemGPU.dat Rem.dat");
		if(rError != 0) cout<<"ERROR"<<endl;

	        if(qError == 0 && rError == 0) cout<<"PASS"<<endl;

	}
	else
		cout<<n<<" "<<m<<" "<<outerTime/1000.0<<endl;

	
	
	
	delete [] A;
	delete [] B;
	delete [] Q;
	delete [] R;
	return 0;
}





