#include <cassert>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "two_convolution_poly_mul.h"
#include <thrust/copy.h>
#include <thrust/scan.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "fft_aux.h"
#include "list_stockham.h"

int find_largest_bit_width_of_coefficients_naive(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{
    auto num_bits = [](const mpz_class& x) -> int {
        return mpz_sizeinbase(x.get_mpz_t(), 2);
    };
    int largest_bit_width = 0;
    for (const mpz_class& x : a)
        largest_bit_width = std::max(largest_bit_width, num_bits(x));
    for (const mpz_class& x : b)
        largest_bit_width = std::max(largest_bit_width, num_bits(x));
    return largest_bit_width;
}

// __device__ helper: compute number of bits in the limb
__device__ __forceinline__ int limb_bit_length(mp_limb_t limb) {
#if GMP_LIMB_BITS == 32
    // For 32-bit limbs use __clz
    return limb ? (32 - __clz(limb)) : 0;
#else
    // For 64-bit limbs use __clzll
    return limb ? (64 - __clzll(limb)) : 0;
#endif
}

// First kernel: compute per-element bit width and reduce within a block.
__global__ void reduce_max_bit_width_kernel(const size_t* mpz_sizes,
                                              const mp_limb_t* ms_limbs,
                                              int n,
                                              int* partial_max) {
    extern __shared__ int sdata[];
    int tid = threadIdx.x;
    // Use unrolling by a factor of 2 for improved load efficiency.
    int idx = blockIdx.x * blockDim.x * 2 + tid;
    int max_val = 0;

    while (idx < n) {
        int bit_width = 0;
        size_t limbs = mpz_sizes[idx];
        mp_limb_t limb = ms_limbs[idx];
        if (limbs > 0) {
            // If coefficient is nonzero, compute its bit width.
            bit_width = (limbs - 1) * GMP_LIMB_BITS + limb_bit_length(limb);
        }
        max_val = max(max_val, bit_width);

        // Unroll a second element per thread.
        int idx2 = idx + blockDim.x;
        if (idx2 < n) {
            int bit_width2 = 0;
            size_t limbs2 = mpz_sizes[idx2];
            mp_limb_t limb2 = ms_limbs[idx2];
            if (limbs2 > 0) {
                bit_width2 = (limbs2 - 1) * GMP_LIMB_BITS + limb_bit_length(limb2);
            }
            max_val = max(max_val, bit_width2);
        }
        idx += blockDim.x * gridDim.x * 2;
    }
    sdata[tid] = max_val;
    __syncthreads();

    // Intra-block reduction in shared memory.
    for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
        if (tid < s)
            sdata[tid] = max(sdata[tid], sdata[tid + s]);
        __syncthreads();
    }
    // Final warp-level reduction (no __syncthreads needed with volatile).
    if (tid < 32) {
        volatile int* smem = sdata;
        smem[tid] = max(smem[tid], smem[tid + 32]);
        smem[tid] = max(smem[tid], smem[tid + 16]);
        smem[tid] = max(smem[tid], smem[tid + 8]);
        smem[tid] = max(smem[tid], smem[tid + 4]);
        smem[tid] = max(smem[tid], smem[tid + 2]);
        smem[tid] = max(smem[tid], smem[tid + 1]);
    }
    if (tid == 0)
        partial_max[blockIdx.x] = sdata[0];
}

// Second kernel: reduce an array of ints (partial max values) to a single maximum.
__global__ void reduce_max_kernel(const int* d_in, int n, int* d_out) {
    extern __shared__ int sdata[];
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x * 2 + tid;
    int max_val = 0;
    if (idx < n)
        max_val = d_in[idx];
    if (idx + blockDim.x < n)
        max_val = max(max_val, d_in[idx + blockDim.x]);
    sdata[tid] = max_val;
    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
        if (tid < s)
            sdata[tid] = max(sdata[tid], sdata[tid + s]);
        __syncthreads();
    }
    if (tid < 32) {
        volatile int* smem = sdata;
        smem[tid] = max(smem[tid], smem[tid + 32]);
        smem[tid] = max(smem[tid], smem[tid + 16]);
        smem[tid] = max(smem[tid], smem[tid + 8]);
        smem[tid] = max(smem[tid], smem[tid + 4]);
        smem[tid] = max(smem[tid], smem[tid + 2]);
        smem[tid] = max(smem[tid], smem[tid + 1]);
    }
    if (tid == 0)
        d_out[blockIdx.x] = sdata[0];
}

// Host function: launch kernels to perform the complete reduction.
int find_largest_bit_width_of_coefficients_dev(const CoefficientsOnDevice& coeffs)
{
    const thrust::device_vector<size_t>& d_mpz_sizes = coeffs.mpz_sizes;
    const thrust::device_vector<mp_limb_t>& d_most_significant_mpz_limbs = coeffs.most_significant_mpz_limbs;

    int n = d_mpz_sizes.size();
    if (n == 0) return 0;

    // First-level reduction configuration.
    int threadsPerBlock = 256;
    // We unroll by a factor of 2.
    int blocks = (n + threadsPerBlock * 2 - 1) / (threadsPerBlock * 2);

    // Allocate a device vector for block–level partial maximums.
    thrust::device_vector<int> d_partial_max(blocks);
    const size_t* raw_mpz_sizes = thrust::raw_pointer_cast(d_mpz_sizes.data());
    const mp_limb_t* raw_ms_limbs = thrust::raw_pointer_cast(d_most_significant_mpz_limbs.data());
    int* raw_partial_max = thrust::raw_pointer_cast(d_partial_max.data());

    size_t sharedMemSize = threadsPerBlock * sizeof(int);
    // Launch the first kernel.
    reduce_max_bit_width_kernel<<<blocks, threadsPerBlock, sharedMemSize>>>
        (raw_mpz_sizes, raw_ms_limbs, n, raw_partial_max);
    cudaDeviceSynchronize();

    // Continue reducing the partial results until only one value remains.
    int s = blocks;
    while (s > 1) {
        int threads = (s < threadsPerBlock * 2) ? ((s + 1) / 2) : threadsPerBlock;
        int grid = (s + threads * 2 - 1) / (threads * 2);
        thrust::device_vector<int> d_out(grid);
        int* raw_in = thrust::raw_pointer_cast(d_partial_max.data());
        int* raw_out = thrust::raw_pointer_cast(d_out.data());
        reduce_max_kernel<<<grid, threads, threads * sizeof(int)>>>(raw_in, s, raw_out);
        cudaDeviceSynchronize();
        // Swap the partial results with the output for the next iteration.
        d_partial_max.swap(d_out);
        s = grid;
    }
    int result;
    thrust::copy(d_partial_max.begin(), d_partial_max.end(), &result);
    return result;
}


CoefficientsOnDevice copy_polynomial_data_to_device(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{
    thrust::host_vector<size_t> mpz_sizes(a.size() + b.size());

    size_t total_limbs = 0;
    for (size_t i = 0; i < a.size() + b.size(); ++i) {
        mpz_sizes[i] = mpz_size(i >= a.size() ? b[i - a.size()].get_mpz_t() : a[i].get_mpz_t());
        total_limbs += mpz_sizes[i];
    }
    // Would it be better instead to cudaMemcpy each coefficient individually, so that
    // each is 256-bit aligned?
    thrust::host_vector<mp_limb_t> mpz_limbs(total_limbs);
    thrust::host_vector<mp_limb_t> most_significant_mpz_limbs(a.size() + b.size());
    for (size_t i = 0, offset = 0; i < a.size() + b.size(); ++i) {
        // A lot of read-after-write hazards here, maybe could be an optimization target
        const mpz_srcptr mpz = (i >= a.size() ? b[i - a.size()].get_mpz_t() : a[i].get_mpz_t());
        std::memcpy(mpz_limbs.data() + offset, mpz->_mp_d, mpz_sizes[i] * sizeof(mp_limb_t));
        offset += mpz_sizes[i];
        most_significant_mpz_limbs[i] = mpz_sizes[i] ? mpz_limbs[offset - 1] : 0;
    }

    return CoefficientsOnDevice {
        thrust::device_vector<size_t>(mpz_sizes),
        thrust::device_vector<mp_limb_t>(mpz_limbs),
        thrust::device_vector<mp_limb_t>(most_significant_mpz_limbs)
    };
}


BivariateMPZPolynomial convert_to_modular_bivariate(const UnivariateMPZPolynomial& p, const BivariateBase& base, sfixn prime)
{
    assert(base.K * base.M == base.N);
    BivariateMPZPolynomial bi(p.size() * base.K);
    const int block_size {base.M};
    const int y_terms = p.size();

    auto convert_mpz_to_modular_univariate = [&](int y_power) {
        mpz_srcptr raw_mpz {p[y_power].get_mpz_t()};
        size_t num_limbs {mpz_size(raw_mpz)};
        
        size_t current_block_bits = 0;
        sfixn current_block = 0;
        
        size_t x_power = 0;
        // Iterate through limbs from least significant to most significant
        for (size_t limb_idx = 0; limb_idx < num_limbs; ++limb_idx) {
            mp_limb_t limb {raw_mpz->_mp_d[limb_idx]};
            
            // Process the current limb
            size_t bits_remaining {GMP_NUMB_BITS};
            while (x_power < base.K) {
                // Determine how many bits to take
                size_t bits_to_take = std::min(block_size - current_block_bits, bits_remaining);
                
                // Extract bits and add to current block
                sfixn extracted_bits {limb & ((1UL << bits_to_take) - 1)};
                current_block |= (extracted_bits << current_block_bits);
                current_block_bits += bits_to_take;
                
                // If block is full, process it
                if (current_block_bits == block_size) {
                    bi[y_power * base.K + x_power++] = current_block % prime;
                    current_block = 0;
                    current_block_bits = 0;
                }
                
                // Shift limb and update remaining bits
                limb >>= bits_to_take;
                bits_remaining -= bits_to_take;
            }
        }
        
        // Handle any remaining bits in an incomplete block
        if (current_block_bits > 0) {
            bi[y_power * base.K + x_power++] = current_block % prime;
        }
    };

    for (int y_power = 0; y_power < y_terms; ++y_power)
        convert_mpz_to_modular_univariate(y_power);
    return bi;
}

// CUDA kernel: one thread per coefficient.
__global__
void convert_kernel(const size_t* __restrict__ d_mpz_sizes,
                      const mp_limb_t* __restrict__ d_mpz_limbs,
                      const size_t* __restrict__ d_offsets,
                      sfixn* __restrict__ d_modular_bivariate,
                      int num_coeffs, int baseK, int baseM, sfixn prime)
{
    int coeff = blockIdx.x * blockDim.x + threadIdx.x;
    if (coeff >= num_coeffs) return;
    
    // Starting offset and size (number of limbs) for this coefficient.
    size_t offset = d_offsets[coeff];
    size_t size   = d_mpz_sizes[coeff];
    
    // Process each of the base.K blocks.
    for (int j = 0; j < baseK; j++) {
        // Compute the bit position of the jth block.
        int bit_offset = j * baseM;
        int limb_index = bit_offset / GMP_LIMB_BITS;
        int bit_in_limb = bit_offset % GMP_LIMB_BITS;
        
        // Initialize the block value to zero.
        mp_limb_t block_val = 0;
        
        // If the block's starting limb exists, extract from it.
        if (limb_index < size) {
            // Shift the limb so that the desired bits are in the lower part.
            block_val = d_mpz_limbs[offset + limb_index] >> bit_in_limb;
            // If the block spans a limb boundary, fetch the extra bits.
            if ((bit_in_limb + baseM) > GMP_LIMB_BITS && (limb_index + 1) < size) {
                int bits_from_next = (bit_in_limb + baseM) - GMP_LIMB_BITS;
                mp_limb_t next_part = d_mpz_limbs[offset + limb_index + 1] 
                                       & (((mp_limb_t)1 << bits_from_next) - 1);
                block_val |= next_part << (GMP_LIMB_BITS - bit_in_limb);
            }
        }
        // Mask to keep only baseM bits.
        mp_limb_t mask = (((mp_limb_t)1) << baseM) - 1;
        block_val &= mask;
        
        // Modular reduction: if the block value is at least prime, subtract prime.
        if (block_val >= (mp_limb_t)prime)
            block_val -= prime;
        
        // Store the result. The output layout is coefficient-major:
        // coefficient 0 occupies indices 0..baseK-1, coefficient 1 indices baseK..2*baseK-1, etc.
        d_modular_bivariate[coeff * baseK + j] = block_val;
    }
}

// Host function that wraps the kernel launch.
thrust::device_vector<sfixn> convert_to_modular_bivariate_dev(
    const CoefficientsOnDevice& coeffs,
    const BivariateBase& base,
    sfixn prime)
{
    const thrust::device_vector<size_t>& d_mpz_sizes = coeffs.mpz_sizes;
    const thrust::device_vector<mp_limb_t>& d_mpz_limbs = coeffs.mpz_limbs;

    int num_coeffs = d_mpz_sizes.size();
    
    // Compute per-coefficient offsets using an exclusive scan.
    thrust::device_vector<size_t> d_offsets(num_coeffs);
    thrust::exclusive_scan(d_mpz_sizes.begin(), d_mpz_sizes.end(), d_offsets.begin());
    
    // Allocate the output vector:
    // Each coefficient produces base.K blocks.
    thrust::device_vector<sfixn> d_modular_bivariate(num_coeffs * base.K);
    
    // Launch kernel with one thread per coefficient.
    int threadsPerBlock = 256;
    int blocks = (num_coeffs + threadsPerBlock - 1) / threadsPerBlock;
    convert_kernel<<<blocks, threadsPerBlock>>>(
        thrust::raw_pointer_cast(d_mpz_sizes.data()),
        thrust::raw_pointer_cast(d_mpz_limbs.data()),
        thrust::raw_pointer_cast(d_offsets.data()),
        thrust::raw_pointer_cast(d_modular_bivariate.data()),
        num_coeffs,
        base.K,
        base.M,
        prime);
        
    // Synchronize to ensure kernel completion (and check for errors if desired).
    cudaDeviceSynchronize();
    
    return d_modular_bivariate;
}

// Kernel: each thread multiplies one coefficient by theta^j mod prime,
// where j is the x–exponent given by (index % K).
__global__ void scale_kernel(sfixn* A, size_t num_elements, int K, sfixn prime, const sfixn* theta_powers) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_elements) {
        int j = idx % K; // determine the x–exponent for this coefficient
        sfixn a = A[idx];
        sfixn factor = theta_powers[j];
        A[idx] = (a * factor) % prime;
    }
}

// Host function to scale the x argument of the bivariate polynomial.
void scale_x_argument_dev(thrust::device_vector<sfixn>& A, const BivariateBase& base, sfixn theta, sfixn prime) {
    int K = base.K;
    size_t num_elements = A.size();

    // Precompute theta powers on the host.
    // h_theta_powers[j] = theta^j mod prime for j = 0,...,K-1.
    thrust::host_vector<sfixn> h_theta_powers(K);
    sfixn current = 1;
    for (int j = 0; j < K; ++j) {
        h_theta_powers[j] = current;
        current = (current * theta) % prime;
    }

    // Copy the theta powers to device memory.
    thrust::device_vector<sfixn> d_theta_powers = h_theta_powers;

    // Set up CUDA kernel execution configuration.
    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;

    // Launch the kernel.
    scale_kernel<<<gridSize, blockSize>>>(
        thrust::raw_pointer_cast(A.data()),
        num_elements,
        K,
        prime,
        thrust::raw_pointer_cast(d_theta_powers.data())
    );

    // Wait for GPU to finish before returning.
    cudaDeviceSynchronize();
}

thrust::device_vector<sfixn> two_convolution_2d_dev(const thrust::device_vector<sfixn>& A, const thrust::device_vector<sfixn>& B, const BivariateBase& base, sfixn prime)
{
    sfixn product_x_size = base.K;
    sfixn product_y_size = A.size() / base.K + B.size() / base.K - 1;

    auto num_bits = [](sfixn x) {
        return sizeof(sfixn) * 8 - __builtin_clz(x);
    };

    sfixn product_x_num_bits = num_bits(product_x_size);
    sfixn product_y_num_bits = num_bits(product_y_size);
    sfixn padded_product_x_size = 1LL << product_x_num_bits;
    sfixn padded_product_y_size = 1LL << product_y_num_bits;

    thrust::device_vector<sfixn> padded_A(padded_product_x_size * padded_product_y_size, 0);
    const sfixn* A_ptr = thrust::raw_pointer_cast(A.data());
    sfixn* padded_A_ptr = thrust::raw_pointer_cast(padded_A.data());
    expand_to_fft2_dev(product_x_num_bits, product_y_num_bits, padded_A_ptr, base.K, A.size() / base.K, A_ptr);
    
    thrust::device_vector<sfixn> padded_B(padded_product_x_size * padded_product_y_size, 0);
    const sfixn* B_ptr = thrust::raw_pointer_cast(B.data());
    sfixn* padded_B_ptr = thrust::raw_pointer_cast(padded_B.data());
    expand_to_fft2_dev(product_x_num_bits, product_y_num_bits, padded_B_ptr, base.K, B.size() / base.K, B_ptr);

    // The output is written to padded_A, so we store a copy of it for when we compute the negacyclic convolution
    thrust::device_vector<sfixn> negacyclic_convolution(padded_A);
    bi_stockham_poly_mul_dev(padded_product_y_size, product_y_num_bits, padded_product_x_size, product_x_num_bits, padded_A_ptr, padded_B_ptr, prime);


}

UnivariateMPZPolynomial two_convolution_poly_mul(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{



    
    /*
    // Find the largest bit-width of any coefficient in a or b.
    // Note that it suffices to compare the number of GMP limbs and the most significant limb in each coefficient.
    const int largest_bit_width_of_coefficients {find_largest_bit_width_of_coefficients(a, b)};

    BivariateBase base {determine_bivariate_base(largest_bit_width_of_coefficients)};
    assert(base.K * base.M == base.N);

    // Convert the univariate polynomials a and b to bivariate polynomials.

    // assuming a machine word is 64 bit
    constexpr sfixn prime_word1 = 4179340454199820289;
    constexpr sfixn prime_word2 = 2485986994308513793;

    BivariateMPZPolynomial a_bivariate1 {convert_to_modular_bivariate(a, base, prime_word1)};
    BivariateMPZPolynomial a_bivariate2 {convert_to_modular_bivariate(a, base, prime_word2)};

    BivariateMPZPolynomial b_bivariate1 {convert_to_modular_bivariate(b, base, prime_word1)};
    BivariateMPZPolynomial b_bivariate2 {convert_to_modular_bivariate(b, base, prime_word2)};

    BivariateMPZPolynomial c_minus1 {cyclic_convolution(a_bivariate1, b_bivariate1, base, prime_word1)};
    BivariateMPZPolynomial c_minus2 {cyclic_convolution(a_bivariate2, b_bivariate2, base, prime_word2)};

    BivariateMPZPolynomial c_plus1 {negacyclic_convolution(a_bivariate1, b_bivariate1, base, prime_word1)};
    BivariateMPZPolynomial c_plus2 {negacyclic_convolution(a_bivariate2, b_bivariate2, base, prime_word2)};
    */
}
