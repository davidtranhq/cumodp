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
#include "fastPolyEvaluation.h"
#include "cumodp_simple.h"
#include "defines.h"
#include "rdr_poly.h"
#include "inlines.h"
#include "printing.h"
#include "cudautils.h"
#include "types.h"

using namespace std;


/************************************************/
void
getFpoly (sfixn *F, int numPoints, int p, int verify)
{
	if (verify == 1)
	{
		ifstream ifs ("PolyF.dat", ifstream::in);
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
getPoints (sfixn *M, int numPoints, int p, int verify)
{
	if (verify == 1)
	{
		ifstream ifs ("Points.dat", ifstream::in);
		for (int i = 0; i < numPoints; ++i)
			ifs >> M[i];
		ifs.close ();
	}
	else
	{
		for (int i = 0; i < numPoints; ++i)
			M[i] = rand () % p;

	}
}

/************************************************/
void
printFtree (sfixn *F, struct PolyEvalSteps check, int k)
{
	ofstream ofs ("PolyFgpu.dat", ofstream::out);
	int i, j, l, length, no;
	length = 1L << (k);

	for (i = 0; i < length; ++i)
		ofs << F[i] << " ";
	ofs << endl;

	length = 1L << (k - 1);
	no = 1;
	sfixn *T = new sfixn[length * no];
	for (j = k; j >= 1; --j)
	{
		cudaMemcpy (T, check.fL[j], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (i = 0; i < no; ++i)
		{
			for (l = 0; l < length; ++l)
				ofs << T[i * length + l] << " ";
			ofs << endl;
		}
		cudaMemcpy (T, check.fR[j], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (i = 0; i < no; ++i)
		{
			for (l = 0; l < length; ++l)
				ofs << T[i * length + l] << " ";
			ofs << endl;
		}

		length = length / 2;
		no = no * 2;
	}
	delete[] T;
	ofs.close ();
}

/************************************************/
void
printSubInvTree (struct PolyEvalSteps check, int k)
{
	ofstream ofs ("PolyMinvMgpu.dat", ofstream::out);
	int j, length, no;
	int i = plainMulLimit;
	length = 1 << (i);
	no = 1 << (k - i - 1);
	sfixn *T = new sfixn[length * no];
	for (; i < k; ++i)
	{
		length = 1 << (i);
		no = 1 << (k - i - 1);

		cudaMemcpy (T, check.InvMl[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < length * no; ++j)
		{
			ofs << T[j] << " ";
			if (j % length == (length - 1))
				ofs << endl;
		}
		cudaMemcpy (T, check.InvMr[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < length * no; ++j)
		{
			ofs << T[j] << " ";
			if (j % length == (length - 1))
				ofs << endl;
		}
	}
	delete[] T;
	ofs.close ();
}

void
printSubProductTree (struct PolyEvalSteps check, int k)
{
	ofstream ofs ("PolyMgpu.dat", ofstream::out);

	int i, j, l, length, no;
	length = 2;
	no = 1 << (k - 1);
	sfixn *T = new sfixn[length * no];
	for (i = 0; i < plainMulLimit; ++i)
	{
		cudaMemcpy (T, check.Ml[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < no; ++j)
		{
			for (l = 0; l < length; ++l)
			{
				ofs << T[j * length + l] << " ";
			}
			ofs << endl;
		}

		cudaMemcpy (T, check.Mr[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < no; ++j)
		{
			for (l = 0; l < length; ++l)
			{
				ofs << T[j * length + l] << " ";
			}
			ofs << endl;
		}

		no = no / 2;
		length = 2 * length - 1;

	}
	length = 1 << (i);
	no = 1 << (k - i - 1);
	for (; i < k; ++i)
	{
		length = 1 << (i);
		no = 1 << (k - i - 1);
		cudaMemcpy (T, check.Ml[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < no; ++j)
		{
			for (l = 0; l < length; ++l)
			{
				ofs << T[j * length + l] << " ";
			}
			ofs << "1" << endl;
		}
		cudaMemcpy (T, check.Mr[i], sizeof(sfixn) * no * length,
								cudaMemcpyDeviceToHost);
		for (j = 0; j < no; ++j)
		{
			for (l = 0; l < length; ++l)
			{
				ofs << T[j * length + l] << " ";
			}
			ofs << "1" << endl;
		}
	}
	delete[] T;
	ofs.close ();
}


/************************************************/
int
main (int argc, char *argv[])
{
	/*
	 We are working on prime filed (mod p).
	 We have 2^k number of points to evaluate.
	 */
	int k = 10;
	sfixn p = 469762049;
	int verify = 0;
	int i;
	if (argc > 1)
		k = atoi (argv[1]);
	if (argc > 2)
		p = atoi (argv[2]);
	if (argc > 3)
		verify = atoi (argv[3]);
	int numPoints = 1 << k;

	if (verify == 1)
	{
		// Writing the value of k and p into a file	
		ofstream ofs1 ("KP.dat", ofstream::out);
		ofs1 << k << " " << p;
		ofs1.close ();

		/*
		 First call maple to create
		 i)   random 2^k points x_1, x_2,..x_{2^k}. Each x_i in field of prime p
		 ii)  a random polynomial F of length 2^k over finite field of prime p.
		 iii) subproduct tree with  (X-x_1)..(X-x_{2^k})
		 iv)  subinverse tree of the subproduct tree
		 v)   Fast evaluation
		 */
		i = system ("maple -q maple_fastEvaluation.mm");
	}
	sfixn *points = new sfixn[numPoints];
	sfixn *F = new sfixn[numPoints];
	getPoints (points, numPoints, p, verify); // getting the values of x for which we need to evaluate F.
	getFpoly (F, numPoints, p, verify); // getting the polynomial F to be evaluated.

	cudaEvent_t start, stop;
	cudaEventCreate (&start);
	cudaEventCreate (&stop);
	cudaEventRecord (start, 0);

	struct PolyEvalSteps check = fastEvaluation (F, points, k, p, verify);
	// the last parameter is a flag. 
	// 1 means we need the output at each level of subproduct tree, subinverse tree and evaluation.
	// 0 means we will need only the evaluation.
	// points will have the evaluation value for F.

	cudaEventRecord (stop, 0);
	cudaEventSynchronize (stop);
	float outerTime;
	cudaEventElapsedTime (&outerTime, start, stop);
	if (verify != 1)
		cout << k << "  " << outerTime / 1000.0 << endl;

	if (verify == 1)
	{
		int testI, testJ, testK;
		cout
				<< "Testing the outcome of Fast evaluation with Maple for every steps. with k = "
				<< k << endl;
		printSubProductTree (check, k);
		testI = system ("diff -b PolyM.dat PolyMgpu.dat");
		if (testI != 0)
			cout << " Error in subproduct tree" << endl;
		printSubInvTree (check, k);
		if (k > plainMulLimit)
		{
			testJ = system ("diff -b PolyMinv.dat PolyMinvMgpu.dat");
			if (testJ != 0)
				cout << " Error in subinverse tree" << endl;
		}

		printFtree (F, check, k);
		testK = system ("diff -b PolyF.dat PolyFgpu.dat");
		if (testK != 0)
			cout << " Error in evaluation step" << endl;

		for (i = 0; i < k; ++i)
		{
			cudaFree (check.Ml[i]);
			cudaFree (check.Mr[i]);
		}

		for (i = plainMulLimit; i < k; ++i)
		{
			cudaFree (check.InvMl[i]);
			cudaFree (check.InvMr[i]);
		}

		for (i = k; i > 1; --i)
		{
			cudaFree (check.fL[i]);
			cudaFree (check.fR[i]);
		}
		if (testI == 0 && testK == 0 && testJ == 0)
			cout << "PASS" << endl;
		else
			cout << "FAIL" << endl;

	}

	/* 
	 #########################################################
	 check.fL[1] and check.fR[1] contain the evaluation result
	 #########################################################
	 */

	cudaFree (check.fL[1]);
	cudaFree (check.fR[1]);

	delete[] points;
	delete[] F;
	return 0;
}

