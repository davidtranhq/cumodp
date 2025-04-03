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


#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include "cumodp_simple.h"
#include "defines.h"
#include "inlines.h"
#include "cudautils.h"
#include "types.h"
#include "opt_plain_mul.h"
#include "stockham.h"
// #include "two_convolution_poly_mul.h"

using namespace std;

/************************************************/
void
getMpoly (sfixn *F, int numPoints, int p, int verify)
{
	if (verify == 1)
	{
		ifstream ifs ("PolyM.dat", ifstream::in);
		for (int i = 0; i < numPoints; ++i)
			ifs >> F[i];
		ifs.close ();
	}
	else
	{
		for (int i = 0; i < numPoints; ++i)
			F[i] = rand () % p;
	}
}

/************************************************/
void
getNpoly (sfixn *F, int numPoints, int p, int verify)
{
	if (verify == 1)
	{
		ifstream ifs ("PolyN.dat", ifstream::in);
		for (int i = 0; i < numPoints; ++i)
			ifs >> F[i];
		ifs.close ();
	}
	else
	{
		for (int i = 0; i < numPoints; ++i)
			F[i] = rand () % p;
	}
}

/************************************************/
void
printCpoly (sfixn *C, int c, const std::string& output_file)
{
	ofstream ofs (output_file, ofstream::out);
	for (int i = 0; i < c; ++i)
		ofs << C[i] << " ";
	ofs.close ();
}

/*
void
printMpzPoly(const UnivariateMPZPolynomial& poly, const std::string& output_file)
{
	ofstream ofs (output_file, ofstream::out);
	for (auto x : poly)
		ofs << x << " ";
	ofs.close ();
}
	*/

/*
 IMPORTANT
 n >= m
 verify = 1, if you want to test with Maple.
 n, the degree of first poly.
 m, the degree of second poly.
 p, prime for finite field.

 */

/************************************************/
int
main (int argc, char *argv[])
{
	int n = 4096, m = 4096, p = 469762049, verify = 1;
	if (argc > 1)
		n = atoi (argv[1]);
	if (argc > 2)
		m = atoi (argv[2]);
	if (argc > 3)
		p = atoi (argv[3]);
	if (argc > 4)
		verify = atoi (argv[4]);

	if (verify == 1)
	{
		// Writing the value of m, n and p into a file	
		ofstream ofs1 ("MNP.dat", ofstream::out);
		ofs1 << n << " " << m << " " << p;
		ofs1.close ();

		/*
		 First call maple to create
		 i)  two random polynomials M and N of length (m+1) and (n+1) over finite field of prime p.
		 ii) the multiplication of M and N and write the result in polyC.dat

		 */
		int i = system ("maple -q maple_PolyMult.mm");
	}

	sfixn *M = new sfixn[m];
	sfixn *N = new sfixn[n];
	sfixn *C = new sfixn[m + n - 1];
	sfixn *D = new sfixn[m + n - 1];
	sfixn *E = new sfixn[m + n - 1];

	getMpoly (M, m, p, verify); // getting 1st poly
	getNpoly (N, n, p, verify); // getting 2nd poly
	// mulCUDA (M, N, C, m, n, p, verify);

    // printCpoly(C, m + n, "PolyCgpu.dat");

	// It has a segamentation fault at the end
	stockham_poly_mul_host(m+n, D, m-1, M, n-1, N, p);
    printCpoly(D, m + n, "PolyCgpu.dat");

	// UnivariateMPZPolynomial a(M, M + m);
	// UnivariateMPZPolynomial b(N, N + n);
	// UnivariateMPZPolynomial c = two_convolution_poly_mul(a, b);
	// printMpzPoly(c, "PolyCgpu.dat");

	/*
	 int j;
	 cout<<endl;
	 for(j = 0; j < m; ++j)
	 cout<<M[j]<<" ";
	 cout<<endl;
	 for(j = 0; j < n; ++j)
	 cout<<N[j]<<" ";
	 cout<<endl;
	 for(j = 0; j < m+n-1; ++j)
	 cout<<C[j]<<" ";
	 cout<<endl;

	 */

	if (verify == 1)
	{
		cout << "Testing polynomial multiplication with " << n << " " << m << " "
				<< p << endl;
		printCpoly (C, m + n - 1, "output.dat");
		int testI = system ("diff -b PolyC.dat PolyCgpu.dat");
		if (testI != 0)
			cout << "ERROR" << endl;
		else
			cout << "PASS" << endl;

	}

	delete[] M;
	delete[] N;
	delete[] C;
	delete[] D;
	return 0;
}

