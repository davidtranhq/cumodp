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
#include "condensationFiniteField.h"

using namespace std;

/************************************************/
void
writefile (char *filename, int *A, int n, int P)
{
	int i, j;
	ofstream ofs (filename, ofstream::out);
	ofs << "P := " << P << " :" << endl;
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
	int n = 10, i, P = 469762049;
	int BK = 4, TH = 128;
	int BS = 8;
	int *a;
	int verbose = 0;
	//int inputMatrix = 0;
	if (argc > 1)
		n = atoi (argv[1]);
	if (argc > 2)
		P = atoi (argv[2]);
	if (argc > 3)
		BK = atoi (argv[3]);
	if (argc > 4)
		TH = atoi (argv[4]);
	if (argc > 5)
		BS = atoi (argv[5]);
	if (argc > 6)
		verbose = atoi (argv[6]);
	//if (argc > 7) inputMatrix = atoi(argv[7]);

	//for(n = 1100; n<=4000;n=n+100){
	a = (int *) malloc ((sizeof(int) * n * n));

	//int j = 0;
	for (i = 0; i < n * n; ++i)
		a[i] = rand () % P;
	if (verbose == 1)
	{
		writefile ("determinant.input", a, n, P);
		cout << "Determinant (from Maple): ";

		i = system ("maple -q maple_determinantFF.mm");
	}
	sfixn det = detFF (a, n, P, BK, TH, BS, verbose);
	cout << "Determinant (from GPU): " << endl << det << endl;
	cout << detFFtime << " s" << endl;
	free (a);
	return 0;
}

