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




#ifndef _RDR_POLY_H_
#define _RDR_POLY_H_
	
#include <cstring>
#include <iostream>
#include <cstdio>
#include "defines.h"

/**
 * Recursive dense multivariate polynomials over finite fields
 * @n      : the number of variables,
 * @ns     : array of size n, ns[i] is the size in the i-th variable,
 * @sz     : sz = ns[0] * ... * ns[n-1] is the size of coeffs,
 * @coeffs : the coefficient array of size sz.
 *
 * Example:  Let x < y < z and F = 1 + x + y^2 + z + x * y + y * z + z^2.
 * First, we set and compute
 *
 * ns[0] = 2, ns[1] = 3, ns[2] = 3 and sz = 2 * 3 * 3 = 18.
 *
 * Polynomial F can be written recursively as
 *
 * F = (1 + x + x * y + y^2) + (1 + y) * z + z^2
 *   = ((1 + 1 * x) + (0 + 1 * x) * y + (1 + 0 * x) * y^2) + 
 *     ((1 + 0 * x) + (1 + 0 * x) * y + (0 + 0 * x) * y^2) * z + 
 *     ((1 + 0 * x) + (0 + 0 * x) * y + (0 + 0 * x) * y^2) * z^2.  
 *
 * The coefficient vector is: [1 1 0 1 1 0, 1 0 1 0 0 0, 1 0 0 0 0 0].     
 ***/
struct rdr_poly {
    // ctors and dtor
    rdr_poly(sfixn _n, sfixn *_ns);
    rdr_poly(sfixn _n, sfixn *_ns, sfixn *_coeffs, bool host = true); 
    // specialized for bivariate polynomials
    rdr_poly(sfixn nx, sfixn ny);
    rdr_poly(sfixn nx, sfixn ny, const sfixn *_coeffs, bool host = true); 

    rdr_poly(const rdr_poly &other);            // deep copy used
    rdr_poly& operator=(const rdr_poly &other); // empty, not implemented

    ~rdr_poly();

    // methods
    sfixn pdeg(sfixn) const;
    sfixn tdeg() const;
    void info() const;
    bool is_zero() const;
    void print_to_maple(const char * msg = "") const;
    int num_of_vars() const { return n; }
    int partial_size(int i) const { return ns[i]; }   
    // fields
    sfixn n;
    sfixn *ns; 
    sfixn sz;
    sfixn *coeffs;   

private: 
    sfixn pdeg_inner(sfixn, const sfixn *) const;
    sfixn tdeg_inner(const sfixn *) const;
    sfixn shrink_deg(sfixn deg, const sfixn *F, sfixn coefSz) const;
    void print_inner(sfixn k, const sfixn *data, const sfixn *accum, 
                     const char* vars) const;
};

#endif // END_OF_FILE
