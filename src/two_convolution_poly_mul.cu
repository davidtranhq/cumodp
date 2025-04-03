#include <cuda.h>
#include <cuda_runtime.h>
#include <cassert>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "two_convolution_poly_mul.h"
#include <thrust/copy.h>
#include <thrust/scan.h>
#include "fft_aux.h"
#include "list_stockham.h"
#include "inlines.h"

#define TWO_CONV_POLY_MUL_DEBUG 0
#if TWO_CONV_POLY_MUL_DEBUG
    #define DEBUG_PRINT(x, ...) fprintf(stderr, x, ##__VA_ARGS__); fflush(stdout);
#else
    #define DEBUG_PRINT(x, ...)
#endif

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

constexpr BivariateBase determine_bivariate_base(sfixn largest_bit_width)
{
    using TwoConvolutionConstants::base_table;

    for (size_t i = 0; i < base_table.size(); ++i) {
        if (largest_bit_width <= base_table[i].N && base_table[i].K >= 32)
            return base_table[i];
    }
    assert(false && "No suitable base found, the largest bit width exceeds the maximum base size");
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
    DEBUG_PRINT("convert_to_modular_bivariate(%d (polynomial size), %d (base.K), %d (base.M), %d (prime))\n", p.size(), base.K, base.M, prime);
    assert(base.K * base.M == base.N);
    BivariateMPZPolynomial bi(p.size() * base.K);
    const int block_size {base.M};
    const int y_terms = p.size();

    auto convert_mpz_to_modular_univariate = [&](int y_power) {
        bool isNeg = 0;
        mpz_t coef {*p[y_power].get_mpz_t()};
	    int s = coef->_mp_size; //abs(s)>=1
        if (s < 0) {
            isNeg = 1;
            s = -s;
        } else if (s == 0) {
            return; // skip zero coefficients
        }

        mp_ptr limbs = coef->_mp_d;
        mp_limb_t carry = 0;
        mp_size_t carry_size = 0;

        // all the limbs except for the last one are full-size (64 bits)
        int idx = y_power * base.K;
        for (int i = 0; i < s-1; ++i) {
            mp_limb_t w = limbs[i], q;
            size_t k = 0; // current position in w
            if (!carry_size) {
                while (GMP_LIMB_BITS - k >= base.M) {
                    q = w >> base.M;
                    bi[idx] = (sfixn) ((w) ^ (q << base.M));
                    if (isNeg) {
                        bi[idx] = -bi[idx];
                        bi[idx] = bi[idx] + prime;
                    }
                    idx++;
                    w = q;
                    k += base.M;
                }
            }
            else {
                size_t l = base.M - carry_size;
                while (GMP_LIMB_BITS - k >= l) {
                    q = w >> l;
                    bi[idx] = (sfixn) ((((w) ^ (q << l)) << carry_size) + carry); //w % 2^a + carry
                    if (isNeg) {
                        bi[idx] = -bi[idx];
                        bi[idx] = bi[idx] + prime;
                    }
                    w = q;
                    idx++;
                    k += l;

                    l = base.M;
                    carry_size = 0;
                    carry = 0;
                }
            }
            carry_size = GMP_LIMB_BITS - k;
            carry = w;
        }

        // last limb, its bit size may be smaller than long
        mp_limb_t w = limbs[s-1];
        size_t b = log2(w)+1; // num of bits
        size_t k = 0; // current position in w

        if (!carry_size) {
            while (b - k >= base.M) {
                mp_limb_t q = w >> base.M;
                bi[idx] = (sfixn) ((w) ^ (q << base.M));
                if (isNeg) {
                    bi[idx] = -bi[idx];
                    bi[idx] = bi[idx] + prime;
                }
                idx++;
                w = q;
                k += base.M;
            }
            if (w) {
                bi[idx] = (sfixn) w;
                if (isNeg) {
                    bi[idx] = -bi[idx];
                    bi[idx] = bi[idx] + prime;
                }
            }
        }
        else {
            size_t l = base.M - carry_size;
            while (b - k >= l) {
                mp_limb_t q = w >> l; //quotient = w / 2^a
                bi[idx] = (sfixn) ((((w) ^ (q << l)) << carry_size) + carry); //w % 2^a + carry
                if (isNeg) {
                    bi[idx] = -bi[idx];
                    bi[idx] = bi[idx] + prime;
                }
                w = q;
                k += l;
                idx++;

                l = base.M;
                carry_size = 0;
                carry = 0;
            }
            if (w) {
                bi[idx] = (sfixn) ((w << carry_size) + carry);
                if (isNeg) {
                    bi[idx] = -bi[idx];
                    bi[idx] = bi[idx] + prime;
                }
            }
        }
    };

    for (int y_power = 0; y_power < y_terms; ++y_power)
        convert_mpz_to_modular_univariate(y_power);
    
    DEBUG_PRINT("Converted to modular bivariate polynomial:\n");
    DEBUG_PRINT("\tbi: "); for (auto x : bi) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
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
        A[idx] = ((longfixnum)a * factor) % prime;
    }
}

// Host function to scale the x argument of the bivariate polynomial.
void scale_x_argument_dev(thrust::device_vector<sfixn>& A, sfixn size, sfixn theta, sfixn prime) {
    size_t num_elements = A.size();

    // Precompute theta powers on the host.
    // h_theta_powers[j] = theta^j mod prime for j = 0,...,K-1.
    thrust::host_vector<sfixn> h_theta_powers(size);
    sfixn current = 1;
    for (int j = 0; j < size; ++j) {
        h_theta_powers[j] = current;
        current = ((longfixnum)current * theta) % prime;
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
        size,
        prime,
        thrust::raw_pointer_cast(d_theta_powers.data())
    );

    // Wait for GPU to finish before returning.
    cudaDeviceSynchronize();
}

struct ModularSumAndDifferenceResult {
    thrust::device_vector<sfixn> modular_sum;
    thrust::device_vector<sfixn> modular_difference;
};
ModularSumAndDifferenceResult modular_sum_and_difference_dev(const thrust::device_vector<sfixn>& u,
                                                            const thrust::device_vector<sfixn>& v,
                                                            sfixn prime)
{
    DEBUG_PRINT("modular_sum_and_difference_dev(%d, %d, %d)\n", u.size(), v.size(), prime);
    assert(u.size() == v.size());
    thrust::device_vector<sfixn> modular_sum(u.size());
    thrust::transform(u.begin(), u.end(),
                      v.begin(),
                      modular_sum.begin(),
                      [prime] __device__ (sfixn a, sfixn b) -> sfixn {
                          sfixn sum = a + b;
                          return (sum >= prime) ? sum - prime : sum;
                      });

    DEBUG_PRINT("\tmodular_sum: "); for (sfixn x : modular_sum) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    
    thrust::device_vector<sfixn> modular_difference(u.size());
    thrust::transform(u.begin(), u.end(),
                      v.begin(),
                      modular_difference.begin(),
                      [prime] __device__ (sfixn a, sfixn b) -> sfixn {
                          sfixn diff = a - b;
                          return (diff < 0) ? diff + prime : diff;
                      });

    DEBUG_PRINT("\tmodular_difference: "); for (sfixn x : modular_difference) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    return {
        .modular_sum = std::move(modular_sum),
        .modular_difference = std::move(modular_difference)
    };
}

/*
// From https://github.com/orcca-uwo/BPAS/blob/7a7ce86819eba3f21098b6cc71a9762ca443d871/include/FFT/src/modpn_hfiles/inlineFuncs.h#L1614
static inline sfixn MontMulModSpe_OPT3_AS_GENE_globalfunc(sfixn a, sfixn b, sfixn inv, sfixn prime)
{
    static constexpr int montgomery_base_power = 30;
    static_assert(TwoConvolutionConstants::prime1 >> montgomery_base_power, "Montgomery base must be less than prime1");
    static_assert(TwoConvolutionConstants::prime2 >> montgomery_base_power, "Montgomery base must be less than prime2");

    asm("mulq %2\n\t"
        "movq %%rax,%%rsi\n\t"
        "movq %%rdx,%%rdi\n\t"
        "imulq %3,%%rax\n\t"
        "mulq %4\n\t"
        "add %%rsi,%%rax\n\t"
        "adc %%rdi,%%rdx\n\t"
        "subq %4,%%rdx\n\t"
        "mov %%rdx,%%rax\n\t"
        "sar %%cl,%%rax\n\t"
        "andq %4,%%rax\n\t"
        "addq %%rax,%%rdx\n\t"
        : "=d"(a)
        : "a"(a), "rm"(b), "rm"(inv), "rm"(prime), "c"(montgomery_base_power)
        : "rsi", "rdi");
    return a;
}
*/

/*
// From https://github.com/orcca-uwo/BPAS/blob/7a7ce86819eba3f21098b6cc71a9762ca443d871/src/IntegerPolynomial/Multiplication/MulSSA-64_bit_arithmetic.cpp#L384
void reconstruct_mpz_with_crt(mpz_t zp, sfixn *a, sfixn *b, int size, size_t limb_bits)
{

    mpz_t zn;
    mpz_init2(zp, limb_bits);
    mpz_init2(zn, limb_bits);

    int idx = 0;
    unsigned __int128 postainer[2] = {0}, negtainer[2] = {0};
    sfixn shifter = 0;
    for (int i = 0; i < size; ++i)
    {
        // sfixn diff =(a[i] - b[i]);
        // elem = elem * P2 + b[i];
        sfixn diff = (a[i] - b[i]);
        if (diff < 0)
        {
            diff += TwoConvolutionConstants::prime1;
        }
        __int128 elem = MontMulModSpe_OPT3_AS_GENE_globalfunc(diff, TwoConvolutionConstants::u2_r1_sft, TwoConvolutionConstants::inverse_prime1, TwoConvolutionConstants::prime1);
        elem = elem * TwoConvolutionConstants::prime2 + b[i];
        if (elem > TwoConvolutionConsants::half_prime1_prime2)
        {
            elem -= TwoConvolutionConstants::prime1_prime2
        }
        else if (elem < TwoConvolutionConsants::neg_half_prime1_prime2)
        {
            elem += TwoConvolutionConstants::prime1_prime2;
        }

        if (elem < 0)
        {
            elem = -elem;
            unsigned __int128 tmp = elem << shifter;
            negtainer[0] += tmp;
            bool carry = negtainer[0] < tmp;
            if (shifter)
                tmp = elem >> (128 - shifter);
            else
            {
                tmp = 0;
            }
            negtainer[1] += tmp + carry;
        }
        else if (elem > 0)
        {
            unsigned __int128 tmp = elem << shifter;
            postainer[0] += tmp;
            bool carry = postainer[0] < tmp;
            if (shifter)
                tmp = elem >> (128 - shifter);
            else
            {
                tmp = 0;
            }
            postainer[1] += tmp + carry;
        }
        shifter += M;

        if (shifter >= 128)
        {
            if (postainer[0] > 0)
            {
                zp->_mp_d[idx] = (mp_limb_t)postainer[0];
                zp->_mp_d[idx + 1] = (mp_limb_t)(postainer[0] >> GMP_LIMB_BITS);
                zp->_mp_size = idx + 2;
            }
            else
            {
                zp->_mp_d[idx] = 0;
                zp->_mp_d[idx + 1] = 0;
            }
            postainer[0] = postainer[1];
            postainer[1] = 0;

            if (negtainer[0] > 0)
            {
                zn->_mp_d[idx] = (mp_limb_t)negtainer[0];
                zn->_mp_d[idx + 1] = (mp_limb_t)(negtainer[0] >> GMP_LIMB_BITS);
                zn->_mp_size = idx + 2;
            }
            else
            {
                zn->_mp_d[idx] = 0;
                zn->_mp_d[idx + 1] = 0;
            }
            negtainer[0] = negtainer[1];
            negtainer[1] = 0;

            shifter -= 128;
            idx += 2;
        }
    }

    if (postainer[0] > 0)
    {
        zp->_mp_d[idx] = (mp_limb_t)postainer[0];
        zp->_mp_d[idx + 1] = (mp_limb_t)(postainer[0] >> GMP_LIMB_BITS);
        zp->_mp_size = idx + 2;
    }
    if (negtainer[0] > 0)
    {
        zn->_mp_d[idx] = (mp_limb_t)negtainer[0];
        zn->_mp_d[idx + 1] = (mp_limb_t)(negtainer[0] >> GMP_LIMB_BITS);
        zn->_mp_size = idx + 2;
    }
    idx += 2;
    if (postainer[1] > 0)
    {
        zp->_mp_d[idx] = (mp_limb_t)postainer[1];
        zp->_mp_d[idx + 1] = (mp_limb_t)(postainer[1] >> GMP_LIMB_BITS);
        zp->_mp_size = idx + 2;
    }
    if (negtainer[1] > 0)
    {
        zn->_mp_d[idx] = (mp_limb_t)negtainer[1];
        zn->_mp_d[idx + 1] = (mp_limb_t)(negtainer[1] >> GMP_LIMB_BITS);
        zn->_mp_size = idx + 2;
    }

    if (zp->_mp_size && !zp->_mp_d[zp->_mp_size - 1])
        zp->_mp_size--;
    if (zn->_mp_size && !zn->_mp_d[zn->_mp_size - 1])
        zn->_mp_size--;

    if (zn > 0)
    {
        mpz_sub(zp, zp, zn);
    }
    mpz_clear(zn);
}
*/

/**
 * Convert from CRT representation to a mpz_t object
 *
 * Output:
 * @zp: Big integer
 *
 * Input:
 * @a: An array mod prime1: a[0] + a[1] * 2^M + ...
 * @b: An array mod prime2: b[0] + b[1] * 2^M + ...
 * @c: An array mod prime3: c[0] + c[1] * 2^M + ...
 * @size: Size of each array
 **/
void reconstruct_mpz_with_crt(mpz_t zp, const sfixn* a, const sfixn* b, const sfixn* c, const BivariateBase& base, int limb_bits)
{
    using namespace TwoConvolutionConstants;
	mpz_t zn;
	mpz_init2(zp, limb_bits);
	mpz_init2(zn, limb_bits);
    int size = base.K;

	int idx = 0;
	unsigned __int128 postainer[2] = {0}, negtainer[2] = {0};
	sfixn shifter = 0;
	for (int i = 0; i < size; ++i) {
		__int128 elem = (__int128) (((long int) a[i] * U23) % prime1) * P2_P3
				+ (__int128) (((long int) b[i] * U31) % prime2) * P3_P1
				+ (__int128) (((long int) c[i] * U12) % prime3) * P1_P2;
		if (elem > HALF_P1_P2_P3) { elem -= P1_P2_P3; }
		else if (elem < N_HALF_P1_P2_P3) { elem += P1_P2_P3; }

		if (elem < 0) {
			elem = -elem;
			unsigned __int128 tmp = elem << shifter;
			negtainer[0] += tmp;
			bool carry = negtainer[0] < tmp;
			if (shifter)
				tmp = elem >> (128 - shifter);
			else { tmp = 0; }
			negtainer[1] += tmp + carry;
		}
		else if (elem > 0) {
			unsigned __int128 tmp = elem << shifter;
			postainer[0] += tmp;
			bool carry = postainer[0] < tmp;
			if (shifter)
				tmp = elem >> (128 - shifter);
			else { tmp = 0; }
			postainer[1] += tmp + carry;
		}
		shifter += base.M;

		if (shifter >= 128) {
			if (postainer[0] > 0) {
				zp->_mp_d[idx] = (mp_limb_t) postainer[0];
				zp->_mp_d[idx+1] = (mp_limb_t) (postainer[0] >> GMP_LIMB_BITS);
				zp->_mp_size = idx + 2;
			}
			else { zp->_mp_d[idx] = 0; zp->_mp_d[idx+1] = 0; }
			postainer[0] = postainer[1];
			postainer[1] = 0;

			if (negtainer[0] > 0) {
				zn->_mp_d[idx] = (mp_limb_t) negtainer[0];
				zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[0] >> GMP_LIMB_BITS);
				zn->_mp_size = idx + 2;
			}
			else { zn->_mp_d[idx] = 0; zn->_mp_d[idx+1] = 0; }
			negtainer[0] = negtainer[1];
			negtainer[1] = 0;

			shifter -= 128;
			idx += 2;
		}
	}

	if (postainer[0] > 0) {
		zp->_mp_d[idx] = (mp_limb_t) postainer[0];
		zp->_mp_d[idx+1] = (mp_limb_t) (postainer[0] >> GMP_LIMB_BITS);
		zp->_mp_size = idx + 2;
	}
	if (negtainer[0] > 0) {
		zn->_mp_d[idx] = (mp_limb_t) negtainer[0];
		zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[0] >> GMP_LIMB_BITS);
		zn->_mp_size = idx + 2;
	}
	idx += 2;
	if (postainer[1] > 0) {
		zp->_mp_d[idx] = (mp_limb_t) postainer[1];
		zp->_mp_d[idx+1] = (mp_limb_t) (postainer[1] >> GMP_LIMB_BITS);
		zp->_mp_size = idx + 2;
	}
	if (negtainer[1] > 0) {
		zn->_mp_d[idx] = (mp_limb_t) negtainer[1];
		zn->_mp_d[idx+1] = (mp_limb_t) (negtainer[1] >> GMP_LIMB_BITS);
		zn->_mp_size = idx + 2;
	}

	if (zp->_mp_size && !zp->_mp_d[zp->_mp_size-1])
		zp->_mp_size--;
	if (zn->_mp_size && !zn->_mp_d[zn->_mp_size-1])
		zn->_mp_size--;

	if (mpz_cmp_ui(zn, 0) > 0) { mpz_sub (zp, zp, zn); }
	mpz_clear(zn);
}


UnivariateMPZPolynomial recover_product_dev(
    const thrust::device_vector<sfixn>& sum1,
    const thrust::device_vector<sfixn>& sum2,
    const thrust::device_vector<sfixn>& sum3,
    const thrust::device_vector<sfixn>& diff1,
    const thrust::device_vector<sfixn>& diff2,
    const thrust::device_vector<sfixn>& diff3,
    const BivariateBase& base,
    int limb_bits)
{
    DEBUG_PRINT("Recovering Product from CRT representation...\n");
    assert(sum1.size() == sum2.size());
    assert(diff1.size() == diff2.size());
    assert(sum1.size() == diff1.size());
    UnivariateMPZPolynomial product(sum1.size() / base.K);
    for (size_t i = 0; i < product.size(); ++i) {
        // Offset into the bivariate polynomial pointing to the start of y^i
        size_t offset = i * base.K;
        mpz_t uv_sum, uv_diff;

        reconstruct_mpz_with_crt(
            uv_sum,
            thrust::raw_pointer_cast(sum1.data() + offset),
            thrust::raw_pointer_cast(sum2.data() + offset),
            thrust::raw_pointer_cast(sum3.data() + offset),
            base,
            limb_bits);
        reconstruct_mpz_with_crt(
            uv_diff,
            thrust::raw_pointer_cast(diff1.data() + offset),
            thrust::raw_pointer_cast(diff2.data() + offset),
            thrust::raw_pointer_cast(diff3.data() + offset),
            base,
            limb_bits);
        mpz_mul_2exp(uv_diff, uv_diff, base.N);
        mpz_add(product[i].get_mpz_t(), uv_sum, uv_diff);
        product[i] >>= 1;
        mpz_clear(uv_sum);
        mpz_clear(uv_diff);
    }
    DEBUG_PRINT("Product recovered:\n");
    DEBUG_PRINT("\tproduct: "); for (auto x : product) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    return product;
}

// Rewrite the above function but using host_vectors instead of device_vectors
UnivariateMPZPolynomial recover_product_host(
    const thrust::host_vector<sfixn>& sum1,
    const thrust::host_vector<sfixn>& sum2,
    const thrust::host_vector<sfixn>& sum3,
    const thrust::host_vector<sfixn>& diff1,
    const thrust::host_vector<sfixn>& diff2,
    const thrust::host_vector<sfixn>& diff3,
    const BivariateBase& base,
    int limb_bits)
{
    DEBUG_PRINT("Recovering Product from CRT representation...\n");
    assert(sum1.size() == sum2.size());
    assert(diff1.size() == diff2.size());
    assert(sum1.size() == diff1.size());
    UnivariateMPZPolynomial product(sum1.size() / base.K);
    for (size_t i = 0; i < product.size(); ++i) {
        // Offset into the bivariate polynomial pointing to the start of y^i
        size_t offset = i * base.K;
        mpz_t uv_sum, uv_diff;

        reconstruct_mpz_with_crt(
            uv_sum,
            thrust::raw_pointer_cast(sum1.data() + offset),
            thrust::raw_pointer_cast(sum2.data() + offset),
            thrust::raw_pointer_cast(sum3.data() + offset),
            base,
            limb_bits);
        reconstruct_mpz_with_crt(
            uv_diff,
            thrust::raw_pointer_cast(diff1.data() + offset),
            thrust::raw_pointer_cast(diff2.data() + offset),
            thrust::raw_pointer_cast(diff3.data() + offset),
            base,
            limb_bits);
        mpz_mul_2exp(uv_diff, uv_diff, base.N);
        mpz_add(product[i].get_mpz_t(), uv_sum, uv_diff);
        product[i] >>= 1;
        mpz_clear(uv_sum);
        mpz_clear(uv_diff);
    }
    DEBUG_PRINT("Product recovered:\n");
    DEBUG_PRINT("\tproduct: "); for (auto x : product) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    return product;
}




TwoConvolutionResult two_convolution_2d_dev(const thrust::device_vector<sfixn>& A, const thrust::device_vector<sfixn>& B, const BivariateBase& base, sfixn prime)
{
    DEBUG_PRINT("Performing 2D convolution with base %d (K), %d (M), %d (N) and prime %d...\n", base.K, base.M, base.N, prime);

    sfixn product_x_size = base.K;
    sfixn product_y_size = A.size() / base.K + B.size() / base.K - 1;

    // assert product_x_size is a power of 2

    auto num_bits = [](sfixn x) -> sfixn {
        return sizeof(sfixn) * 8 - __builtin_clz(x);
    };

    auto is_power_of_two = [](sfixn x) { return (x & (x - 1)) == 0; };

    sfixn log_padded_product_x_size = num_bits(base.K) - 1;
    sfixn log_padded_product_y_size = num_bits(product_y_size) - (is_power_of_two(product_y_size) ? 1 : 0);

    log_padded_product_x_size = std::max(log_padded_product_x_size, log_padded_product_y_size);
    log_padded_product_y_size = std::max(log_padded_product_x_size, log_padded_product_y_size);

    sfixn padded_product_x_size = 1 << log_padded_product_x_size;
    sfixn padded_product_y_size =  1 << log_padded_product_y_size;

    thrust::device_vector<sfixn> padded_A(padded_product_x_size * padded_product_y_size);
    const sfixn *A_ptr = thrust::raw_pointer_cast(A.data());
    sfixn *padded_A_ptr = thrust::raw_pointer_cast(padded_A.data());
    expand_to_fft2_dev(log_padded_product_x_size, log_padded_product_y_size, padded_A_ptr, product_x_size, A.size() / base.K, A_ptr);
    DEBUG_PRINT("\tpadded_A: "); for (sfixn x : padded_A) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");

    thrust::device_vector<sfixn> padded_B(padded_product_x_size * padded_product_y_size);
    const sfixn *B_ptr = thrust::raw_pointer_cast(B.data());
    sfixn *padded_B_ptr = thrust::raw_pointer_cast(padded_B.data());
    expand_to_fft2_dev(log_padded_product_x_size, log_padded_product_y_size, padded_B_ptr, product_x_size, B.size() / base.K, B_ptr);
    DEBUG_PRINT("\tpadded_B: "); for (sfixn x : padded_B) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");

    // The output is written to padded_A, so we store a copy of it for when we compute the negacyclic convolution
    thrust::device_vector<sfixn> padded_A_copy(padded_A);
    thrust::device_vector<sfixn> padded_B_copy(padded_B);
    bi_stockham_poly_mul_dev(padded_product_y_size, log_padded_product_y_size, padded_product_x_size, log_padded_product_x_size, padded_A_ptr, padded_B_ptr, prime);
    // TODO: bi_stockham_poly_mul_dev perfomrs FFT on B here, maybe we can reuse this result for NCC 
    DEBUG_PRINT("\tFFT2(A * B): "); for (sfixn x : padded_A) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    cudaDeviceSynchronize();
    thrust::device_vector<sfixn> extracted_A(product_x_size * product_y_size);
    extract_from_fft2_dev(product_x_size, product_y_size, thrust::raw_pointer_cast(extracted_A.data()), log_padded_product_x_size, padded_A_ptr);

    // theta is a 2(padded K)th primitive root of unity
    // something to fix here, not that clean, num_bits(base.K) computed twice (once at the start of the function?)
    sfixn theta = primitive_root(log_padded_product_x_size + 1, prime);
    scale_x_argument_dev(padded_A_copy, padded_product_x_size, theta, prime);    
    scale_x_argument_dev(padded_B_copy, padded_product_x_size, theta, prime);
    cudaDeviceSynchronize();

    // padded_A_copy will store the result of the negacyclic convolution
    sfixn *padded_A_copy_ptr = thrust::raw_pointer_cast(padded_A_copy.data());
    sfixn *padded_B_copy_ptr = thrust::raw_pointer_cast(padded_B_copy.data());

    bi_stockham_poly_mul_dev(padded_product_y_size, log_padded_product_y_size, padded_product_x_size, log_padded_product_x_size, padded_A_copy_ptr, padded_B_copy_ptr, prime);

    cudaDeviceSynchronize();

    scale_x_argument_dev(padded_A_copy, padded_product_x_size, inv_mod(theta, prime), prime);
    cudaDeviceSynchronize();
    
    thrust::device_vector<sfixn> extracted_A_copy(product_x_size * product_y_size);
    extract_from_fft2_dev(product_x_size, product_y_size, thrust::raw_pointer_cast(extracted_A_copy.data()), log_padded_product_x_size, padded_A_copy_ptr);

    DEBUG_PRINT("2D convolution completed\n");
    DEBUG_PRINT("\tCyclic Convolution: "); for (sfixn x : extracted_A) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
    DEBUG_PRINT("\tNegacyclic Convolution: "); for (sfixn x : extracted_A_copy) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");

    return TwoConvolutionResult {
        .cyclic_convolution = std::move(extracted_A),
        .negacyclic_convolution = std::move(extracted_A_copy),
    };
}

UnivariateMPZPolynomial two_convolution_poly_mul(const UnivariateMPZPolynomial& a_src, const UnivariateMPZPolynomial& b_src)
{
    UnivariateMPZPolynomial a = a_src;
    UnivariateMPZPolynomial b = b_src;
    if (a.size() < 32)
        a.resize(32);
    if (b.size() < 32)
        b.resize(32);
    BivariateBase base {determine_bivariate_base(find_largest_bit_width_of_coefficients_dev(copy_polynomial_data_to_device(a, b)))};
    assert(base.K * base.M == base.N);

    // Convert the univariate polynomials a and b to bivariate polynomials.

    // assuming a machine word is 64 bit

    struct UVSumAndDifference {
        thrust::device_vector<sfixn> uv_sum;
        thrust::device_vector<sfixn> uv_diff;
    };

    auto compute_uv_sum_and_diff = [&] (sfixn prime) -> UVSumAndDifference {
        DEBUG_PRINT("Computing UV sum and difference for prime %d...\n", prime);
        BivariateMPZPolynomial a_bivariate {convert_to_modular_bivariate(a, base, prime)};
        BivariateMPZPolynomial b_bivariate {convert_to_modular_bivariate(b, base, prime)};
        const auto& [cyclic_convolution, negacyclic_convolution] = two_convolution_2d_dev(a_bivariate, b_bivariate, base, prime);
        const auto& [uv_sum, uv_diff] = modular_sum_and_difference_dev(cyclic_convolution, negacyclic_convolution, prime);
        DEBUG_PRINT("UV sum and difference computed for prime %d:\n", prime);
        DEBUG_PRINT("\tUV sum: "); for (sfixn x : uv_sum) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
        DEBUG_PRINT("\tUV diff: "); for (sfixn x : uv_diff) { DEBUG_PRINT("%d ", x); } DEBUG_PRINT("\n");
        return UVSumAndDifference {
            .uv_sum = std::move(uv_sum),
            .uv_diff = std::move(uv_diff)
        };
    };

    using namespace TwoConvolutionConstants;
    UVSumAndDifference uv_sum_and_diff1 = compute_uv_sum_and_diff(prime1);
    UVSumAndDifference uv_sum_and_diff2 = compute_uv_sum_and_diff(prime2);
    UVSumAndDifference uv_sum_and_diff3 = compute_uv_sum_and_diff(prime3);
    auto bitwidth = [](sfixn x) {
        return sizeof(sfixn) * 8 - __builtin_clz(x);
    };
    auto result = recover_product_host(
        uv_sum_and_diff1.uv_sum,
        uv_sum_and_diff2.uv_sum,
        uv_sum_and_diff3.uv_sum,
        uv_sum_and_diff1.uv_diff,
        uv_sum_and_diff2.uv_diff,
        uv_sum_and_diff3.uv_diff,
        base,
        1 + base.M + base.N + bitwidth(base.K) + bitwidth(std::max(a.size(), b.size()))
    );
    result.resize(a_src.size() + b_src.size() - 1); // resize to the correct size
    return result;
}
