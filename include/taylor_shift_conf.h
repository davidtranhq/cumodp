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



#ifndef _TAYLOR_SHIFT_CONF_H_
#define _TAYLOR_SHIFT_CONF_H_

// Libraries :
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
using namespace std;


// Number of threads per block (size of a block) :
#define NB_THREADS 512

// Number of threads per block (Chinese Remainder) :
#define T_CHN 32

//#define MAX_LEVEL 25

typedef int sfixn;


const sfixn Tmul = 512;
//const int BASE_1 = 31;

// Debuging flags
//#define DEBUG 0

#endif // _TAYLOR_SHIFT_CONF_H_
