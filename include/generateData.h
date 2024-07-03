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



#include <random>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
typedef unsigned long long int usfixn64;
// typedef unsigned int usfixn64;
using namespace std;

usfixn64* vectorGenerator(usfixn64 * vector, int vectorSize, int &coeffientSize, usfixn64 &range, int &mode, int &seed)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> rDistribution(0.0, 1.0);
	generator.seed(seed);
	double generated;

	switch (mode)
	{
	case 0:
		cout << "all less than R-1 and random (default)" << endl;
		break;
	case 1:
		cout << "all zero except the last one equal to R (corner case)" << endl;
		break;
	case 2:
		cout << "all equal to seed = " << seed << endl;
		break;

	}

//	for(int b=0; b<vectorSize/1024; b++)
	for (int i = 0; i < vectorSize * coeffientSize; i++)
//		for (int i = b*1024*coeffientSize; i < (b+1)*1024* coeffientSize; i++)
	{

		generated = rDistribution(generator);
		switch (mode)
		{
		case 0: /* 0: all less than R-1 and random (default)*/
//			if (i % 8 != 7) //all digits except the biggest one
			vector[i] = generated * (range - 1);
//				cout<<vector[i]<<endl;
//			else
//				vector[i] = generated*(R-1);
			break;

		case 1: /* 1: all zero except the last one equal to R (corner case)*/
			//all digits except the biggest one
			if (i % coeffientSize != coeffientSize-1) 
				vector[i] = 0;
			else
				vector[i] = range;
			break;

		case 2: /* 2: all equal to seed*/
			vector[i] = seed;
			break;

		case 3:
			vector[i] = i;
			break;
		case 4:
			vector[i] = i/coeffientSize;
			break;
		}
	}
	return vector;
}

/**********************************************/
