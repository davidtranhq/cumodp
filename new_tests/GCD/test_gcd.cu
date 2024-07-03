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
#include "cudaGcd.h"
#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "printing.h"
#include "cudautils.h"
#include "types.h"

using namespace std;

/************************************************/
sfixn
add_modCPU (sfixn a, sfixn b, sfixn P)
{
	sfixn r = a + b;
	r -= P;
	r += (r >> BASE_1) & P;
	return r;
}

/************************************************/
sfixn
sub_modCPU (sfixn a, sfixn b, int P)
{
	sfixn r = a - b;
	r += (r >> BASE_1) & P;
	return r;
}

/************************************************/
sfixn
neg_modCPU (sfixn a, int P)
{
	sfixn r = -a;
	r += (r >> BASE_1) & P;
	return r;
}

/************************************************/
void
egcdCPU (sfixn x, sfixn y, sfixn *ao, sfixn *bo, sfixn *vo, int P)
{
	sfixn t, A, B, C, D, u, v, q;

	u = y;
	v = x;
	A = 1;
	B = 0;
	C = 0;
	D = 1;

	do
	{
		q = u / v;
		t = u;
		u = v;
		v = t - q * v;
		t = A;
		A = B;
		B = t - q * B;
		t = C;
		C = D;
		D = t - q * D;
	}
	while (v != 0);

	*ao = A;
	*bo = C;
	*vo = u;
}

/************************************************/
sfixn
inv_modCPU (sfixn n, int P)
{
	sfixn a, b, v;
	egcdCPU (n, P, &a, &b, &v, P);
	if (b < 0)
		b += P;
	return b % P;
}

/************************************************/
sfixn
mul_modCPU (sfixn a, sfixn b, int P)
{
	double ninv = 1 / (double) P;
	sfixn q = (sfixn) ((((double) a) * ((double) b)) * ninv);
	sfixn res = a * b - q * P;
	res += (res >> BASE_1) & P;
	res -= P;
	res += (res >> BASE_1) & P;
	return res;
}


/************************************************/
sfixn
pow_modCPU (sfixn a, sfixn ee, int P)
{
	sfixn x, y;
	usfixn e;

	if (ee == 0)
		return 1;
	if (ee == 1)
		return a;
	if (ee < 0)
		e = -((usfixn) ee);
	else
		e = ee;

	x = 1;
	y = a;
	while (e)
	{
		if (e & 1)
			x = mul_modCPU (x, y, P);
		y = mul_modCPU (y, y, P);
		e = e >> 1;
	}

	if (ee < 0)
		x = inv_modCPU (x, P);
	return x;
}

/*
 It should be called in the following sequence
 prime number for prime field, p.
 the degree of bigger polynomial, n-1.
 the degree of smaller polynomial, m-1. n >= m.
 the desired length of the gcd, h.
 whatUwant = 1 for checking with maple.

 */

/************************************************/
int
main (int argc, char *argv[])
{
	int n = 21, m = 14, whatUwant = 0, h = 1;
	int i, j, k;
	sfixn p = 469762049;
	sfixn monic = 1;

	if (argc > 1)
		p = atoi (argv[1]);
//if (argc > 1) whatUwant = atoi(argv[1]);
//if (argc > 1) whatUwant = atoi(argv[1]);
	if (argc > 2)
		n = atoi (argv[2]);
	if (argc > 3)
		m = atoi (argv[3]);
	if (argc > 4)
		h = atoi (argv[4]);
	if (argc > 5)
		whatUwant = atoi (argv[5]);
//if (argc > 5) p = atoi(argv[5]);
//if (argc > 5) p = atoi(argv[5]);
	if (argc > 6)
		monic = atoi (argv[6]);

	/*
	 if(m > n)
	 {
	 i = m;
	 m = n;
	 n = i;
	 }
	 */
	if (h > m)
		h = 1;

	sfixn *A = new int[n];
	sfixn *B = new int[m];
	sfixn *G = new int[m];

//A[0] = 0; A[1] = 0; A[2] = 0; A[3] =  200940312; A[4] =  382537351; A[5] =  292010474; A[6] = 459302886; A[7] =  451721134; A[8] =  433342116;
//A[9] =  330082711; A[10] =  175374780; A[11] =  467965863; A[12] =  320548451; A[13] =  285646418; A[14] =  388799201; A[15] =  458548079; 
//A[16] =  469651169; A[17] =  167988; A[18] = 249360; A[19] =  101400; A[20] =  12000;

//B[0] = 0; B[1] = 0; B[2] = 0; B[3] = 469516289; B[4] = 469669121; B[5] = 468323585; B[6] = 467486912 ; B[7] = 468359449; B[8] = 469317158; B[9] = 469757193 ; B[10] = 4062; 
//B[11] =  6444 ; B[12] = 5016; B[13] = 1440;
	if (whatUwant == 1)
	{
//ofstream ofs1( "PNM.dat", ofstream::out  );	
//ofs1<<p<<" "<<n<<" "<<m<<" "<<h;
//ofs1.close();	
//i = system("maple -q gcd.mm");

		ifstream ifs1 ("Poly1.dat", ifstream::in);
		for (i = 0; i < n; ++i)
			ifs1 >> A[i];
		ifs1.close ();

		ifstream ifs2 ("Poly2.dat", ifstream::in);
		for (i = 0; i < m; ++i)
			ifs2 >> B[i];
		ifs2.close ();

	}
	else
	{
		for (i = 0; i < n; ++i)
			A[i] = rand () % p;
		for (i = 0; i < m; ++i)
			B[i] = rand () % p;
	}

//for(i = 0; i < m; ++i)
//	G[i] = 0;
//////
	sfixn Dgcd;
//cout<<endl;
//for(i = 0; i < n; ++i) cout<<A[i]<<" ";
//cout<<endl;
//for(i = 0; i < m; ++i) cout<<B[i]<<" ";
//cout<<endl;
/////
	float outerTime = gcdPrimeField (&Dgcd, A, B, G, n, m, p, monic);
	cout << n << " " << m << " " << Dgcd << " " << outerTime / 1000.0 << endl;

	/*if(whatUwant == 1)
	 {

	 if(monic == 0)
	 {
	 ///
	 //cout<<"making monic"<<" "<<G[Dgcd]<<" "<<Dgcd<<endl;
	 ///
	 k = inv_modCPU( G[Dgcd], p);
	 for(j = 0; j <= Dgcd; ++j)
	 G[j] = mul_modCPU(G[j], k, p);
	 }
	 ofstream ofs2( "GCDcuda.dat", ofstream::out  );
	 for(i = 0; i <= Dgcd; ++i)
	 ofs2<<G[i]<<" ";
	 ofs2.close();

	 i = system("diff -b GCD.dat GCDcuda.dat");
	 if(i != 0) cout<<"ERROR"<<endl;
	 else cout<<"PASS"<<endl;

	 }
	 else
	 cout<<n<<" "<<m<<" "<<Dgcd<<" "<<outerTime/1000.0<<endl;
	 */

	delete[] A;
	delete[] B;
	delete[] G;

	return 0;
}

