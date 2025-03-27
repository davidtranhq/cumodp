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
//#include "cumodp_simple.h"
//#include "defines.h"
//#include "inlines.h"
//#include "cudautils.h"
//#include "types.h"
#include "opt_plain_mul.h"
#include "stockham.h"

void test_poly_mul();

int main (int argc, char *argv[])
{
    int n = 0, m = 0;
    std::cout << "Enter the number of terms in the first polynomial: ";
    std::cin >> n;
    std::cout << "Enter the number of terms in the second polynomial: ";
    std::cin >> m;

    std::cout << "Enter the " << n << " terms of the first polynomial: ";
    sfixn *N = new sfixn[n];
    for (int i = 0; i < n; ++i)
        std::cin >> N[i];

    std::cout << "Enter the " << m << " terms of the second polynomial: ";
    sfixn *M = new sfixn[m];
    for (int i = 0; i < m; ++i)
        std::cin >> M[i];

    std::cout << "Enter the characteristic of the field to multiply over: ";
    int p = 0;
    std::cin >> p;

    sfixn *C = new sfixn[m + n - 1];
    // Note: if the polynomial is too small then the stockham_poly_mul_host will not perform
    // any work (it requires a minimum polynomial degree of DFTBASESIZE (or whatever its called ...))
    stockham_poly_mul_host(m + n, C, n - 1, N, m - 1, M, p);
    std::cout << "Product:\n";
    for (int i = 0; i < m + n - 1; ++i)
        std::cout << C[i] << ' ';
    std::cout << '\n';
	delete[] M;
	delete[] N;
	delete[] C;
	return 0;
}

