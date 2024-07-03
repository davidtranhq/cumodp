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


#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>



NTL_CLIENT

int main(int argc, char *argv[])
{
  int deg, prime, deg_2;

  if(argc > 1) deg = atoi(argv[1]);
  if(argc > 2) deg_2 = atoi(argv[2]);
  if(argc > 3) prime = atoi(argv[3]);

  zz_p::init(prime);
  zz_pX a, b, c, d, e;
  zz_p deg_p;

  random(a, deg);
  random(b, deg_2);

  double u = GetTime();
  
   div(c, a, b);
 
 u = GetTime() - u;
  cout<<deg<<" "<<deg_2<<" "<<u<<endl;
 
  return 0;
}
