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


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <time.h>
#include "condensationReal.h"

using namespace std;

/************************************************/
void
writefile (char *filename, realV *A, int n)
{
	int i, j;
	ofstream ofs (filename, ofstream::out);
	ofs << "n := " << n << " :" << endl << "[";
	for (i = 0; i < n; ++i)
	{
		ofs << "[";
		for (j = 0; j < n; ++j)
		{
			ofs << A[j * n + i] << " ";
			if (j == n - 1)
				ofs << "]";
			else
				ofs << ",";
		}
		ofs << endl;
		if (i == n - 1)
			ofs << "];" << endl;
		else
			ofs << ",";

	}
	ofs << endl << "mat := Matrix(%) :";
	ofs.close ();
}

/************************************************/
int
main (int argc, char *argv[])
{
	int n = 10, i, BK = 4, BS = 3, test = 0;

	if (argc > 1)
		n = atoi (argv[1]);
	if (argc > 2)
		BK = atoi (argv[2]);
	if (argc > 3)
		BS = atoi (argv[3]);
	if (argc > 4)
		test = atoi (argv[4]);

	realV *a, *b, *c; // closest;
	if (argc > 1)
		n = atoi (argv[1]);

	a = (realV *) malloc (sizeof(realV) * n * n);
	b = (realV *) malloc (sizeof(realV) * n * n);
	c = (realV *) malloc (sizeof(realV) * n * n);
	//cout<<n<<"\t";

	for (i = 0; i < n * n; ++i)
	{

		a[i] = (realV) (rand () % 5); //(realV)rand()/rand();
		if (rand () % 2 == 0)
			a[i] = -a[i];
		b[i] = a[i];
		c[i] = a[i];
	}
	a[0] = a[0] * 1.1;
	b[0] = a[0];
	c[0] = a[0];
	//writefile("determinant.input", a, n);
	/*

	 clock_t start, end;
	 double elapsed;
	 start = clock();

	 int j;
	 for(j = 0; j < 10; ++j)
	 { */
	if (test == 1)
		direct (b, n);

	//det(a, n,  BS,  BK);
	deter (c, n, BS, BK);
	if (test == 1)
	{
		cout << detDirect << " " << detCon << endl;

	}
	else
		cout << "Took " << timeDet << " s";
	cout << endl;

	return 0;
}

