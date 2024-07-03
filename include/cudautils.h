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



///////////////////////////////////////////////////////////////////////////////
#ifndef _CUDA_UTILS_H_
#define _CUDA_UTILS_H_

#include "defines.h"
#include <cstdio>
#include <iostream>


///////////////////////////////////////////////////////////////////////////////
//  Timing a kernel running.
//
//  Usage:
//
//  _TIMING_START_
//  kernel_run();
//  _TIMING_STOP_
//
//  cout << _elapsedTime << endl;
//
///////////////////////////////////////////////////////////////////////////////

#define _TIMING_START_ \
    cudaEvent_t _start, _stop;\
    cudaEventCreate(&_start);\
    cudaEventCreate(&_stop);\
    cudaEventRecord(_start, 0);

#define _TIMING_STOP_ \
    cudaEventRecord(_stop, 0); \
    cudaEventSynchronize(_stop); \
    float _elapsedTime; \
    cudaEventElapsedTime(&_elapsedTime, _start, _stop); \
    cudaEventDestroy(_start); \
    cudaEventDestroy(_stop);

#define start_timer(id) \
    cudaEvent_t start##id, stop##id; \
    cudaEventCreate(&start##id); \
    cudaEventCreate(&stop##id); \
    cudaEventRecord(start##id, 0)

#define stop_timer(id, elapsedTime) \
    cudaEventRecord(stop##id, 0); \
    cudaEventSynchronize(stop##id); \
    cudaEventElapsedTime(&elapsedTime, start##id, stop##id); \
    cudaEventDestroy(start##id); \
    cudaEventDestroy(stop##id)

///////////////////////////////////////////////////////////////////////////////
// Check the last error from device
static inline void checkCudaError(const char *msg) {
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (cudaSuccess != err) {
	#ifndef _mcompile_
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
	#endif
        exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////////////
static inline void printDeviceProperty() {
#ifndef _mcompile_
    int id;
    cudaError_t err = cudaGetDevice(&id);
    if (cudaSuccess != err) {

        fprintf(stderr, "Cuda error: cannot get device: %s\n",
        cudaGetErrorString(err));
    }
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, id);
    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
        fprintf(stderr, "Cuda error: no device supporting CUDA.\n");
        cudaThreadExit();
        exit(EXIT_FAILURE);     
    } else {
        fprintf(stdout, "\n++++++++++++++++++++++++++++++++++++++++++++++++\n");
        fprintf(stdout, "Device %d: \"%s\"\n", id, deviceProp.name);
        fprintf(stdout, "Capability %d.%d\n", deviceProp.major,
                        deviceProp.minor);
        fprintf(stdout, "Maximun size of each dimension of a grid [%d, %d, %d]\n",
                        deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
        fprintf(stdout, "Maximum size of each dimension of a block [%d, %d, %d]\n",
                        deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
        fprintf(stdout, "Maximum number of threads per block %d\n",
                        deviceProp.maxThreadsPerBlock);
        fprintf(stdout, "Number of multiprocessors on device %d\n",
                        deviceProp.multiProcessorCount);
        fprintf(stdout, "32-bit registers available per block %d\n",
                        deviceProp.regsPerBlock); 
        fprintf(stdout, "Shared memory available per block in bytes %d\n",
                        (int)deviceProp.sharedMemPerBlock); 
        fprintf(stdout, "Constant memory available on device in bytes %d\n", 
                        (int)deviceProp.totalConstMem);
        fprintf(stdout, "Global memory available on device in bytes %d\n",
                        (int)deviceProp.totalGlobalMem);
        fprintf(stdout, "Warp size in threads %d\n",
                        (int)deviceProp.warpSize);
        fprintf(stdout, "++++++++++++++++++++++++++++++++++++++++++++++++\n");
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////
inline sfixn get_dev_val(const sfixn * const dev_ptr) {
    sfixn t;
    cudaMemcpy(&t, dev_ptr, sizeof(sfixn), cudaMemcpyDeviceToHost);
    return t;
}

inline void set_dev_val(sfixn * const dev_ptr, sfixn t) {
    cudaMemcpy(dev_ptr, &t, sizeof(sfixn), cudaMemcpyHostToDevice);
}

///////////////////////////////////////////////////////////////////////////////
template <typename T> inline T getDeviceValue(const T* const dev_ptr) {
    T t;
    cudaMemcpy(&t, dev_ptr, sizeof(T), cudaMemcpyDeviceToHost);
    return t;
}

///////////////////////////////////////////////////////////////////////////////
template <typename T> inline void setDeviceValue(T* const dev_ptr, const T t) {
    cudaMemcpy(dev_ptr, &t, sizeof(T), cudaMemcpyHostToDevice);
}

///////////////////////////////////////////////////////////////////////////////
template <typename T> inline void printDeviceVariable(const T* const dev_ptr) {
    T t;
    cudaMemcpy(&t, dev_ptr, sizeof(T), cudaMemcpyDeviceToHost);
    //std::cout << t << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
template <typename T> inline
void printDeviceArray (const T* const arr_d, unsigned int n) {
    T *arr_h = new T[n];
    cudaMemcpy(arr_h, arr_d, sizeof(T) * n, cudaMemcpyDeviceToHost);
#ifndef _mcompile_
    for (unsigned int i = 0; i < n; ++i) {
        printf("%3d ", arr_h[i]);
    }
    printf("\n");
#endif
    delete [] arr_h;
}

template <typename T> inline
void printDeviceMatrix (const T* const arr_d, unsigned int w, unsigned int h) {
    T *arr_h = new T[w * h];
    cudaMemcpy(arr_h, arr_d, sizeof(T) * w * h, cudaMemcpyDeviceToHost);
#ifndef _mcompile_
    for (unsigned int i = 0; i < h; ++i) {
        for (unsigned int j = 0; j < w; ++j) {
            printf("%3d ", arr_h[i * w + j]);
        }
        printf("\n");
    }
    printf("\n");
#endif
    delete [] arr_h;
}

template <typename T> inline
void host_to_device(T* X_d, const T * X_h, sfixn n) {
    cudaMemcpy(X_d, X_h, n*sizeof(T), cudaMemcpyHostToDevice);
}

template <typename T> inline
void device_to_host(T* X_h, const T * X_d, sfixn n) {
    cudaMemcpy(X_h, X_d, n*sizeof(T), cudaMemcpyDeviceToHost);
}

template <typename T> inline
void device_to_device(T* X_d2, const T * X_d1, sfixn n) {
    cudaMemcpy(X_d2, X_d1, n*sizeof(T), cudaMemcpyDeviceToDevice);
}

// print the subresultant chains in the triangular encoding
static inline void print_subres_chain2(sfixn B, sfixn w, const sfixn *S) {
    sfixn sz = w * (w + 1) / 2;
    sfixn *S_h = new sfixn[B * sz]();
    cudaMemcpy(S_h, S, sizeof(sfixn)*B*sz, cudaMemcpyDeviceToHost);
#ifndef _mcompile_
    for (sfixn i = 0; i < B; ++i) {
        const sfixn *T = S_h + sz * i;
        for (sfixn k = w; k >= 1; --k) {
            for (sfixn j = 0; j < k; ++j) { printf("%4d ", T[j]); }
            T += k;
            printf("\n");
        }
        printf("\n");
    }
#endif
    delete [] S_h;
}

// print the ith subresultant chains in the triangular encoding
static inline 
void print_subres_chain2(sfixn B, sfixn w, const sfixn *S, sfixn i) {
    sfixn sz = w * (w + 1) / 2;
    sfixn *S_h = new sfixn[B * sz]();
    cudaMemcpy(S_h, S, sizeof(sfixn)*B*sz, cudaMemcpyDeviceToHost);
    const sfixn *T = S_h + sz * i;
#ifndef _mcompile_
    printf("The %d-th regular chain:\n", i);
    for (sfixn k = w; k >= 1; --k) {
        for (sfixn j = 0; j < k; ++j) { printf("%6d ", T[j]); }
        T += k;
        printf("\n");
    }
#endif
    delete [] S_h;
}

// print the subresultant chain in the triangle encoding, transposed
static inline void print_subres_chain3(sfixn B, sfixn w, const sfixn *S) {

    sfixn *S_h = new sfixn[B * w * (w+1) / 2]();
    cudaMemcpy(S_h, S, sizeof(sfixn)*B*w*(w+1)/2, cudaMemcpyDeviceToHost);

    const sfixn *T = S_h; 
#ifndef _mcompile_
    for (sfixn k = w; k >= 1; --k) {
        for (sfixn i = 0; i < B; ++i) {
            for (sfixn j = 0; j < k; ++j) { printf("%4d ", T[k*i+j]); }
            printf("\n");
        }
        T += B * k;
    }
#endif
    delete [] S_h;
}

// print the ith-subresultant chain in the triangle encoding, transposed
static inline 
void print_subres_chain3(sfixn B, sfixn w, const sfixn *S, sfixn i) {
    if (i < 0 || i >= B) return;
    sfixn *S_h = new sfixn[B * w * (w+1) / 2]();
    cudaMemcpy(S_h, S, sizeof(sfixn)*B*w*(w+1)/2, cudaMemcpyDeviceToHost);
    const sfixn *T = S_h; 
#ifndef _mcompile_
    printf("The %d-th regular chain:\n", i);
    for (sfixn k = w; k >= 1; --k) {
        for (sfixn j = 0; j < k; ++j) { printf("%4d ", T[i * k + j]); }
        printf("\n");
        T += B * k;
    }
#endif
    delete [] S_h;
}

#endif /* _CUDA_UTILS_H_ */
