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



#ifndef _TAYLOR_SHIFT_CPU_H_
#define _TAYLOR_SHIFT_CPU_H_
#include "taylor_shift_conf.h"


// error messages if there is a lack of arguments to make the program
  void error_message(sfixn m);
  void error_message2(sfixn m);


// computes the nomber of blocks
  sfixn number_of_blocks(sfixn n);


// stocks a file in an array
  void stock_file_in_array(char* filename, sfixn n, sfixn* & a);


// stockes the array of Newton's coefficients in a file
  void stock_array_in_file(const char *name_file, sfixn *T, sfixn size);


// computes the number of lines of a file
  sfixn size_file(char* filename);


// display of an array
  void display_array(sfixn *T, sfixn size);


// addition of two arrays
  void add_arrays(sfixn *res, sfixn *T1, sfixn *T2, sfixn size, sfixn p);


// Horner's method to compute g(x) = f(x+1) (equivalent to Shaw & Traub's method for a=1)
  void horner_shift_CPU(sfixn *Polynomial, sfixn *Polynomial_shift, sfixn n, sfixn p);


#endif // _TAYLOR_SHIFT_CPU_H_
