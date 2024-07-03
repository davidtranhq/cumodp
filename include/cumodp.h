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



#ifndef _CU_MODP_H_
#define _CU_MODP_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

//////////////////////
// CUDA Device Query
//////////////////////
int is_cuda_enabled();

int num_of_cuda_devices(); 

int is_double_float_enabled();

unsigned int global_memory_in_bytes();

float global_memory_in_megabytes();

int can_call_fftmul_uni(sfixn df, sfixn dg);

int can_call_subres2(sfixn dx, sfixn d1, sfixn d2);

int can_call_subres3(sfixn dx, sfixn dy, sfixn d1, sfixn d2);

///////////////////////////////
//  FFT & FFT Multiplications
///////////////////////////////
cumodp_err 
cumodp_fft_uni(sfixn *X, sfixn n, sfixn k, sfixn w, sfixn p);

cumodp_err 
cumodp_invfft_uni(sfixn *X, sfixn n, sfixn k, sfixn w, sfixn p);

cumodp_err
cumodp_fft_bivariate(sfixn *X, sfixn em, sfixn wm, sfixn en, 
    sfixn wn, sfixn p);

cumodp_err
cumodp_invfft_bivariate(sfixn *X, sfixn em, sfixn wm, sfixn en, 
    sfixn wn, sfixn p);

cumodp_err
cumodp_fftmul_uni(sfixn dh, sfixn *H, sfixn df, const sfixn *F, 
    sfixn dg, const sfixn *G, sfixn p);

///////////////////////////////////////////
// subresultant chain construction for CPU
///////////////////////////////////////////
sfixn
cumodp_subres_chain2_fine(sfixn *S, sfixn B, sfixn w, 
    sfixn npx, sfixn npy, const sfixn *P, 
    sfixn nqx, sfixn nqy, const sfixn *Q, sfixn p);

sfixn
cumodp_subres_chain2_coarse(sfixn *S, sfixn B, sfixn w, 
    sfixn npx, sfixn npy, const sfixn *P, 
    sfixn nqx, sfixn nqy, const sfixn *Q, sfixn p);

sfixn
cumodp_subres_chain3_fine(sfixn *S, sfixn Bx, sfixn By, sfixn wx, sfixn wy,
    sfixn npx, sfixn npy, sfixn npz, const sfixn *P, 
    sfixn nqx, sfixn nqy, sfixn nqz, const sfixn *Q, sfixn p); 

sfixn
cumodp_subres_chain3_coarse(sfixn *S, sfixn Bx, sfixn By, sfixn wx, sfixn wy,
    sfixn npx, sfixn npy, sfixn npz, const sfixn *P, 
    sfixn nqx, sfixn nqy, sfixn nqz, const sfixn *Q, sfixn p);



///////////////////////////////////////////////
/////////////GCD function//////////////////////
///////////////////////////////////////////////
//float gcdPrimeField(sfixn *A, sfixn *B, sfixn *G, sfixn n, sfixn m, sfixn p);
float gcdPrimeField(sfixn *Dgcd, sfixn *A, sfixn *B, sfixn *G, sfixn n, sfixn m, sfixn p, sfixn optionNormalize);
float divPrimeField(sfixn *A, sfixn *B, sfixn *R, sfixn *Q, sfixn n, sfixn m, sfixn p);




///////////////////////////////////////////
// subresultant chain construction for GPU
///////////////////////////////////////////

void *init_cuda_scube(sfixn N, const sfixn *sz_p, const sfixn *sz_q, sfixn fp);

void free_cuda_scube(void *scb);

void print_cuda_scube(const void *scb);

cumodp_err 
build_cuda_scube(void *scb, const sfixn *sz_p, const sfixn *P,
    const sfixn *sz_q, const sfixn *Q);

const sfixn *interp_subres_coeff2(sfixn *nx, void *scb, 
    sfixn i, sfixn j);

const sfixn *interp_subres_coeff3(sfixn *nx, sfixn *ny, void *scb, 
    sfixn i, sfixn j);

#ifdef __cplusplus
}
#endif

#endif
