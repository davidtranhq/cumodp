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



#include "cumodp.h"
#include "stockham.h"
#include "list_stockham.h"
#include "subres.h"
#include "cudautils.h"
#include "inlines.h"
#include "scube.h"
#include <cstdio>

/**
 * Interface for modpn C library
 * 
 * (1) Device Query
 * (2) FFT & Multiplication
 * (3) Subresultant Chain 
 *  
 * Created at Aug 24, 2010, WP
 *
 **/

///////////////////
// Device queries 
///////////////////
#ifdef __cplusplus
extern "C" 
#endif
int is_cuda_enabled() {
    int count;
    cudaGetDeviceCount(&count);
    return (count > 0);
}

#ifdef __cplusplus
extern "C" 
#endif
int num_of_cuda_devices() {
    int count;
    cudaGetDeviceCount(&count);
    return count;
}

#ifdef __cplusplus
extern "C" 
#endif
int is_double_float_enabled() {
    cudaDeviceProp deviceProp;
    // only support a single cuda device, device 0
    cudaGetDeviceProperties(&deviceProp, 0);
    // Returns 9999 for both major & minor fields, 
    // if no CUDA capable devices are present.
    // If major is 1, then devices with minor >= 3 
    // support double float computations.
    // If major is 2, then all these devices do.
    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
        return 0;
    } else if (deviceProp.major == 1){
        return (deviceProp.minor >= 3);
    } 
    return (deviceProp.major > 1);
}

#ifdef __cplusplus
extern "C" 
#endif
unsigned int global_memory_in_bytes() {
    cudaDeviceProp deviceProp;
    // only support a single cuda device, device 0
    cudaGetDeviceProperties(&deviceProp, 0);
    return deviceProp.totalGlobalMem;
}

#ifdef __cplusplus
extern "C" 
#endif
float global_memory_in_megabytes() {
    cudaDeviceProp deviceProp;
    // only support a single cuda device, device 0
    cudaGetDeviceProperties(&deviceProp, 0);
    return (float)(deviceProp.totalGlobalMem >>20);
}

#ifdef __cplusplus
extern "C"
#endif
int can_call_fftmul_uni(sfixn df, sfixn dg) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    int c0 = (deviceProp.major == 9999 && deviceProp.minor == 9999);
    int c1 = (deviceProp.major == 1 && deviceProp.minor < 3);
    if (c0 || c1) return 0;
    // mem is in bytes
    size_t mem = deviceProp.totalGlobalMem;
    sfixn B = df + dg + 1;
    sfixn b = ceiling_log2(B);
    return (floor_log2(mem/4) > b + 2);
}

/* mem-usage in words */
size_t maxmem_subres(sfixn B, sfixn d1, sfixn d2) {
    if (d1 < d2) return maxmem_subres(B, d2, d1);
    // evaluation memusage + scube size + workspace
    size_t N1 = B * (1 + d1);
    size_t N2 = B * (1 + d2);
    size_t N3 = B * d2 * (d2 + 1) / 2;
    size_t N4 = 2 * B * d1 + B + 2;
    return N1 + N2 + N3 + N4;
}

#ifdef __cplusplus
extern "C"
#endif
int can_call_subres2(sfixn dx, sfixn d1, sfixn d2) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    int c0 = (deviceProp.major == 9999 && deviceProp.minor == 9999);
    int c1 = (deviceProp.major == 1 && deviceProp.minor < 3);
    if (c0 || c1) return 0;
    // mem is in bytes
    size_t mem = deviceProp.totalGlobalMem;
    sfixn B = dx * (d1 + d2) + 1;
    sfixn b = ceiling_log2(B);
    B = (sfixn)1 << b;
    return (mem / 4 > 1.3 * maxmem_subres(B, d1, d2));
}

#ifdef __cplusplus
extern "C"
#endif
int can_call_subres3(sfixn dx, sfixn dy, sfixn d1, sfixn d2) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    int c0 = (deviceProp.major == 9999 && deviceProp.minor == 9999);
    int c1 = (deviceProp.major == 1 && deviceProp.minor < 3);
    if (c0 || c1) return 0;
    // mem is in bytes
    size_t mem = deviceProp.totalGlobalMem;
    sfixn Bx = dx * (d1 + d2) + 1;
    sfixn By = dy * (d1 + d2) + 1;
    sfixn bx = ceiling_log2(Bx);
    sfixn by = ceiling_log2(By);
    Bx = ((sfixn)1 << bx);
    By = ((sfixn)1 << by);
    return (mem / 4 > 1.3 * maxmem_subres(Bx * By, d1, d2));
}

//////////////////////////////////////////////////
// FFT and FFT based polynomial multiplications
//////////////////////////////////////////////////
#ifdef __cplusplus
extern "C" 
#endif
cumodp_err 
cumodp_fft_uni(sfixn *X, sfixn n, sfixn k, sfixn w, sfixn p)
{
    // in-place 1d fft
    stockham_host(X, n, k, w, p);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        if (DEBUG) 
	fprintf(stderr, "cumodp_fft_uni: %s.\n", cudaGetErrorString(err));
        return CUMODP_FFT_ERROR;
    }
    return CUMODP_SUCCESS;
}

#ifdef __cplusplus
extern "C" 
#endif
cumodp_err 
cumodp_invfft_uni(sfixn *X, sfixn n, sfixn k, sfixn w, sfixn p)
{
    // in-place 1d inverse fft
    inverse_stockham_host(X, n, k, w, p);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        if (DEBUG) fprintf(stderr, 
            "cumodp_invfft_uni: %s.\n", cudaGetErrorString(err));
        return CUMODP_FFT_ERROR;
    }
    return CUMODP_SUCCESS;
}

#ifdef __cplusplus
extern "C" 
#endif
cumodp_err
cumodp_fft_bivariate(sfixn *X, sfixn em, sfixn wm, sfixn en,
    sfixn wn, sfixn p) 
{
    bivariate_stockham_host(X, em, wm, en, wn, p);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        if (DEBUG) fprintf(stderr, 
            "cumodp_fft_bivariate: %s.\n", cudaGetErrorString(err));
        return CUMODP_FFT_ERROR;
    }
    return CUMODP_SUCCESS;
}

#ifdef __cplusplus
extern "C" 
#endif
cumodp_err
cumodp_invfft_bivariate(sfixn *X, sfixn em, sfixn wm, sfixn en, 
    sfixn wn, sfixn p)
{
    inverse_bivariate_stockham_host(X, em, wm, en, wn, p);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        if (DEBUG) fprintf(stderr, 
            "cumodp_invfft_bivariate: %s.\n", cudaGetErrorString(err));
        return CUMODP_FFT_ERROR;
    }
    return CUMODP_SUCCESS;
}

#ifdef __cplusplus
extern "C" 
#endif
cumodp_err 
cumodp_fftmul_uni(sfixn dh, sfixn *H, sfixn df, const sfixn *F, 
    sfixn dg, const sfixn *G, sfixn p) 
{
    stockham_poly_mul_host(dh, H, df, F, dg, G, p);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        if (DEBUG) fprintf(stderr, 
            "cumodp_fftmul_uni: %s.\n", cudaGetErrorString(err));
        return CUMODP_FFT_ERROR;
    }
    return CUMODP_SUCCESS;
}

////////////////////////////////////
// Subresultant chain constructions
////////////////////////////////////
#ifdef __cplusplus
extern "C" 
#endif
sfixn cumodp_subres_chain2_coarse(sfixn *S, sfixn B, sfixn w, 
    sfixn npx, sfixn npy, const sfixn *P, 
    sfixn nqx, sfixn nqy, const sfixn *Q, sfixn p) 
{
    assert(npy < nqy);
    return subres_chain2_coarse_host(S, B, w, npx, npy, P, nqx, nqy, Q, p);   
}

#ifdef __cplusplus
extern "C" 
#endif
sfixn cumodp_subres_chain3_coarse(sfixn *S, sfixn Bx, sfixn By, 
    sfixn wx, sfixn wy, sfixn npx, sfixn npy, sfixn npz, const sfixn *P, 
    sfixn nqx, sfixn nqy, sfixn nqz, const sfixn *Q, sfixn p) 
{
    assert(npz >= nqz);
    return subres_chain3_coarse_host(S, Bx, By, wx, wy, npx, npy, npz, P, 
        nqx, nqy, nqz, Q, p);
}

#ifdef __cplusplus
extern "C" 
#endif
sfixn cumodp_subres_chain2_fine(sfixn *S, sfixn B, sfixn w, 
    sfixn npx, sfixn npy, const sfixn *P, 
    sfixn nqx, sfixn nqy, const sfixn *Q, sfixn p) 
{
    assert(npy >= nqy);
    return subres_chain2_fine_host(S, B, w, npx, npy, P, nqx, nqy, Q, p);
}

#ifdef __cplusplus
extern "C" 
#endif
sfixn 
cumodp_subres_chain3_fine(sfixn *S, sfixn Bx, sfixn By, sfixn wx, sfixn wy,
    sfixn npx, sfixn npy, sfixn npz, const sfixn *P, 
    sfixn nqx, sfixn nqy, sfixn nqz, const sfixn *Q, sfixn p) 
{
    assert(npz >= nqz);
    return subres_chain3_fine_host(S, Bx, By, wx, wy, npx, npy, npz, P, 
        nqx, nqy, nqz, Q, p);
}

///////////////////////////////////////////////////////////////////////////////
// Callback functions for keeping scube inside GPU
///////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C" 
#endif
void *init_cuda_scube(sfixn N, const sfixn *sz_p, const sfixn *sz_q, sfixn fp) 
{   
    if (!is_double_float_enabled()) return NULL;

    scube_t *scb = new scube_t(N, sz_p, sz_q, fp);

    // set the cut-off value for the fft size
    const int MIN_SCUBE2_FFT_EXP = 10;
    const int MIN_SCUBE3_FFT_EXP = 6;
    sfixn ex = scb->get_bounds_exp(0);
    if (N == 2) {
        if (ex <= MIN_SCUBE2_FFT_EXP) {
            delete scb;
            return NULL;
        }
    }

    if (N == 3) {
        sfixn ey = scb->get_bounds_exp(1);
        if (ex <= MIN_SCUBE3_FFT_EXP || ey <= MIN_SCUBE3_FFT_EXP) {
            delete scb;
            return NULL;   
        }
    }
    return (void *) scb;
}

#ifdef __cplusplus
extern "C" 
#endif
void free_cuda_scube(void *S) {
    scube_t *scb = (scube_t *)S;
    delete scb;
}

#ifdef __cplusplus
extern "C" 
#endif
void print_cuda_scube(const void *S) {
    scube_t *scb = (scube_t *)S;
    scb->info();
}

#ifdef __cplusplus
extern "C" 
#endif
cumodp_err 
build_cuda_scube(void *S, const sfixn *sz_p, const sfixn *P, 
    const sfixn *sz_q, const sfixn *Q) 
{
    scube_t *scb = (scube_t *)S;
    int N = scb->num_of_vars();
    bool ret;
    assert(N == 2 || N == 3); 
    if (N == 2) {
        ret = scb->build_scube_data2(sz_p[0], P, sz_q[0], Q); 
    } else {
        ret = scb->build_scube_data3(sz_p[0], sz_p[1], P, sz_q[0], sz_q[1], Q); 
    }
    // scb->info();
    return (ret == true) ? CUMODP_SUCCESS : CUMODP_FAILURE;
}

#ifdef __cplusplus
extern "C"
#endif
const sfixn *interp_subres_coeff2(sfixn *nx, void *S, sfixn i, sfixn j) 
{
    scube_t *scb = (scube_t *)S;
    sfixn w = scb->get_ldeg();
    if (j > i || i >= w || j < 0) { return NULL; }
    *nx = (sfixn(1) << scb->get_bounds_exp(0));
    // scb->info();
    return scb->subres_coeff(i, j);
}

#ifdef __cplusplus
extern "C"
#endif
const sfixn *interp_subres_coeff3(sfixn *nx, sfixn *ny, void *S, 
    sfixn i, sfixn j) 
{
    scube_t *scb = (scube_t *)S;
    sfixn w = scb->get_ldeg();
    if (j > i || i >= w || j < 0) return NULL;
    *nx = (sfixn(1) << scb->get_bounds_exp(0));
    *ny = (sfixn(1) << scb->get_bounds_exp(1));
    // scb->info();
    return scb->subres_coeff(i, j);
}

/////////////// END OF FILE ///////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//BEGIN:uni_fft_tst
///////////////////////////////////////////////////////////////////////////////
int uni_fft_tst(sfixn p, sfixn k) {
    
    sfixn i;
	//, p = 469762049, k = 24, n = (1L << k), w = 37;
	sfixn n = (1L << k);
	sfixn w = primitive_root(k,p);
    sfixn *X = (sfixn *)malloc(n*sizeof(sfixn));
    for (i = 0; i < n; ++i) X[i] = i % p;
    cumodp_err err = cumodp_fft_uni(X, n, k, w, p);
    
    if (err != CUMODP_SUCCESS) {
        //fprintf(stdout, "=====================\n");
        //fprintf(stdout, "Fail to do FFT for k=%d\n",k);
        //fprintf(stdout, "=====================\n");
        free(X); return -1;
    }

    err = cumodp_invfft_uni(X, n, k, w, p);
    if (err != CUMODP_SUCCESS) {
        //fprintf(stdout, "=======================\n");
        //fprintf(stdout, "Fail to do inverse FFT for k=%d\n",k);
        //fprintf(stdout, "=======================\n");
        free(X); return -2;
    }

    for (i = 0; i < n; ++i) {
        if (X[i] != i) {
            //fprintf(stderr, "Incorrect result for fft & inverse fft\n");
            free(X); return -3;
        }
    }
    free(X);
	return 0;
}
///////////////////////////////////////////////////////////////////////////////
//END:uni_fft_tst
///////////////////////////////////////////////////////////////////////////////
