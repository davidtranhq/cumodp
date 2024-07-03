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


/************************************************/
void
host_cudaMalloc (usfixn* ptr, usfixn numWords, int nBytes_per_word)
{
	cudaMalloc ((void**) &ptr, numWords * nBytes_per_word);
}
/************************************************/
void
host_cudaMemset_to_zero (usfixn* ptr, usfixn numWords, int nBytes_per_word)
{
	cudaMemset ((void**) &ptr, 0x00, numWords * nBytes_per_word);
}

/************************************************/
//read poly from file
//read poly from array
//read matrix from file
//read matrix from array
//print poly to file
//print poly
//print matrix to file
//print matrix

//methods for measuring time
//


/************************************************/
