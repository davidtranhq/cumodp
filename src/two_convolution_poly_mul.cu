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

constexpr BivariateBase determine_bivariate_base(size_t largest_bit_width)
{
    static constexpr std::array<BivariateBase, 432> base_table {
        {2, 2, 1},
        {4, 2, 2},
        {6, 2, 3},
        {8, 2, 4},
        {10, 2, 5},
        {12, 2, 6},
        {14, 2, 7},
        {16, 2, 8},
        {18, 2, 9},
        {20, 2, 10},
        {22, 2, 11},
        {24, 2, 12},
        {26, 2, 13},
        {28, 2, 14},
        {30, 2, 15},
        {32, 2, 16},
        {34, 2, 17},
        {36, 2, 18},
        {38, 2, 19},
        {40, 2, 20},
        {42, 2, 21},
        {44, 2, 22},
        {46, 2, 23},
        {48, 2, 24},
        {50, 2, 25},
        {52, 2, 26},
        {54, 2, 27},
        {56, 2, 28},
        {58, 2, 29},
        {60, 2, 30},
        {62, 2, 31},
        {64, 2, 32},
        {66, 2, 33},
        {68, 2, 34},
        {70, 2, 35},
        {72, 2, 36},
        {74, 2, 37},
        {76, 2, 38},
        {78, 2, 39},
        {80, 2, 40},
        {82, 2, 41},
        {84, 2, 42},
        {86, 2, 43},
        {88, 2, 44},
        {90, 2, 45},
        {92, 2, 46},
        {94, 2, 47},
        {96, 2, 48},
        {98, 2, 49},
        {100, 2, 50},
        {102, 2, 51},
        {104, 2, 52},
        {106, 2, 53},
        {108, 2, 54},
        {110, 2, 55},
        {112, 2, 56},
        {114, 2, 57},
        {116, 4, 29},
        {120, 4, 30},
        {124, 4, 31},
        {128, 4, 32},
        {132, 4, 33},
        {136, 4, 34},
        {140, 4, 35},
        {144, 4, 36},
        {148, 4, 37},
        {152, 4, 38},
        {156, 4, 39},
        {160, 4, 40},
        {164, 4, 41},
        {168, 4, 42},
        {172, 4, 43},
        {176, 4, 44},
        {180, 4, 45},
        {184, 4, 46},
        {188, 4, 47},
        {192, 4, 48},
        {196, 4, 49},
        {200, 4, 50},
        {204, 4, 51},
        {208, 4, 52},
        {212, 4, 53},
        {216, 4, 54},
        {220, 4, 55},
        {224, 4, 56},
        {232, 8, 29},
        {240, 8, 30},
        {248, 8, 31},
        {256, 8, 32},
        {264, 8, 33},
        {272, 8, 34},
        {280, 8, 35},
        {288, 8, 36},
        {296, 8, 37},
        {304, 8, 38},
        {312, 8, 39},
        {320, 8, 40},
        {328, 8, 41},
        {336, 8, 42},
        {344, 8, 43},
        {352, 8, 44},
        {360, 8, 45},
        {368, 8, 46},
        {376, 8, 47},
        {384, 8, 48},
        {392, 8, 49},
        {400, 8, 50},
        {408, 8, 51},
        {416, 8, 52},
        {424, 8, 53},
        {432, 8, 54},
        {440, 8, 55},
        {448, 16, 28},
        {464, 16, 29},
        {480, 16, 30},
        {496, 16, 31},
        {512, 16, 32},
        {528, 16, 33},
        {544, 16, 34},
        {560, 16, 35},
        {576, 16, 36},
        {592, 16, 37},
        {608, 16, 38},
        {624, 16, 39},
        {640, 16, 40},
        {656, 16, 41},
        {672, 16, 42},
        {688, 16, 43},
        {704, 16, 44},
        {720, 16, 45},
        {736, 16, 46},
        {752, 16, 47},
        {768, 16, 48},
        {784, 16, 49},
        {800, 16, 50},
        {816, 16, 51},
        {832, 16, 52},
        {848, 16, 53},
        {864, 16, 54},
        {896, 32, 28},
        {928, 32, 29},
        {960, 32, 30},
        {992, 32, 31},
        {1024, 32, 32},
        {1056, 32, 33},
        {1088, 32, 34},
        {1120, 32, 35},
        {1152, 32, 36},
        {1184, 32, 37},
        {1216, 32, 38},
        {1248, 32, 39},
        {1280, 32, 40},
        {1312, 32, 41},
        {1344, 32, 42},
        {1376, 32, 43},
        {1408, 32, 44},
        {1440, 32, 45},
        {1472, 32, 46},
        {1504, 32, 47},
        {1536, 32, 48},
        {1568, 32, 49},
        {1600, 32, 50},
        {1632, 32, 51},
        {1664, 32, 52},
        {1696, 32, 53},
        {1728, 64, 27},
        {1792, 64, 28},
        {1856, 64, 29},
        {1920, 64, 30},
        {1984, 64, 31},
        {2048, 64, 32},
        {2112, 64, 33},
        {2176, 64, 34},
        {2240, 64, 35},
        {2304, 64, 36},
        {2368, 64, 37},
        {2432, 64, 38},
        {2496, 64, 39},
        {2560, 64, 40},
        {2624, 64, 41},
        {2688, 64, 42},
        {2752, 64, 43},
        {2816, 64, 44},
        {2880, 64, 45},
        {2944, 64, 46},
        {3008, 64, 47},
        {3072, 64, 48},
        {3136, 64, 49},
        {3200, 64, 50},
        {3264, 64, 51},
        {3328, 64, 52},
        {3456, 128, 27},
        {3584, 128, 28},
        {3712, 128, 29},
        {3840, 128, 30},
        {3968, 128, 31},
        {4096, 128, 32},
        {4224, 128, 33},
        {4352, 128, 34},
        {4480, 128, 35},
        {4608, 128, 36},
        {4736, 128, 37},
        {4864, 128, 38},
        {4992, 128, 39},
        {5120, 128, 40},
        {5248, 128, 41},
        {5376, 128, 42},
        {5504, 128, 43},
        {5632, 128, 44},
        {5760, 128, 45},
        {5888, 128, 46},
        {6016, 128, 47},
        {6144, 128, 48},
        {6272, 128, 49},
        {6400, 128, 50},
        {6528, 128, 51},
        {6656, 256, 26},
        {6912, 256, 27},
        {7168, 256, 28},
        {7424, 256, 29},
        {7680, 256, 30},
        {7936, 256, 31},
        {8192, 256, 32},
        {8448, 256, 33},
        {8704, 256, 34},
        {8960, 256, 35},
        {9216, 256, 36},
        {9472, 256, 37},
        {9728, 256, 38},
        {9984, 256, 39},
        {10240, 256, 40},
        {10496, 256, 41},
        {10752, 256, 42},
        {11008, 256, 43},
        {11264, 256, 44},
        {11520, 256, 45},
        {11776, 256, 46},
        {12032, 256, 47},
        {12288, 256, 48},
        {12544, 256, 49},
        {12800, 256, 50},
        {13312, 512, 26},
        {13824, 512, 27},
        {14336, 512, 28},
        {14848, 512, 29},
        {15360, 512, 30},
        {15872, 512, 31},
        {16384, 512, 32},
        {16896, 512, 33},
        {17408, 512, 34},
        {17920, 512, 35},
        {18432, 512, 36},
        {18944, 512, 37},
        {19456, 512, 38},
        {19968, 512, 39},
        {20480, 512, 40},
        {20992, 512, 41},
        {21504, 512, 42},
        {22016, 512, 43},
        {22528, 512, 44},
        {23040, 512, 45},
        {23552, 512, 46},
        {24064, 512, 47},
        {24576, 512, 48},
        {25088, 512, 49},
        {25600, 1024, 25},
        {26624, 1024, 26},
        {27648, 1024, 27},
        {28672, 1024, 28},
        {29696, 1024, 29},
        {30720, 1024, 30},
        {31744, 1024, 31},
        {32768, 1024, 32},
        {33792, 1024, 33},
        {34816, 1024, 34},
        {35840, 1024, 35},
        {36864, 1024, 36},
        {37888, 1024, 37},
        {38912, 1024, 38},
        {39936, 1024, 39},
        {40960, 1024, 40},
        {41984, 1024, 41},
        {43008, 1024, 42},
        {44032, 1024, 43},
        {45056, 1024, 44},
        {46080, 1024, 45},
        {47104, 1024, 46},
        {48128, 1024, 47},
        {49152, 1024, 48},
        {51200, 2048, 25},
        {53248, 2048, 26},
        {55296, 2048, 27},
        {57344, 2048, 28},
        {59392, 2048, 29},
        {61440, 2048, 30},
        {63488, 2048, 31},
        {65536, 2048, 32},
        {67584, 2048, 33},
        {69632, 2048, 34},
        {71680, 2048, 35},
        {73728, 2048, 36},
        {75776, 2048, 37},
        {77824, 2048, 38},
        {79872, 2048, 39},
        {81920, 2048, 40},
        {83968, 2048, 41},
        {86016, 2048, 42},
        {88064, 2048, 43},
        {90112, 2048, 44},
        {92160, 2048, 45},
        {94208, 2048, 46},
        {96256, 2048, 47},
        {98304, 4096, 24},
        {102400, 4096, 25},
        {106496, 4096, 26},
        {110592, 4096, 27},
        {114688, 4096, 28},
        {118784, 4096, 29},
        {122880, 4096, 30},
        {126976, 4096, 31},
        {131072, 4096, 32},
        {135168, 4096, 33},
        {139264, 4096, 34},
        {143360, 4096, 35},
        {147456, 4096, 36},
        {151552, 4096, 37},
        {155648, 4096, 38},
        {159744, 4096, 39},
        {163840, 4096, 40},
        {167936, 4096, 41},
        {172032, 4096, 42},
        {176128, 4096, 43},
        {180224, 4096, 44},
        {184320, 4096, 45},
        {188416, 4096, 46},
        {196608, 8192, 24},
        {204800, 8192, 25},
        {212992, 8192, 26},
        {221184, 8192, 27},
        {229376, 8192, 28},
        {237568, 8192, 29},
        {245760, 8192, 30},
        {253952, 8192, 31},
        {262144, 8192, 32},
        {270336, 8192, 33},
        {278528, 8192, 34},
        {286720, 8192, 35},
        {294912, 8192, 36},
        {303104, 8192, 37},
        {311296, 8192, 38},
        {319488, 8192, 39},
        {327680, 8192, 40},
        {335872, 8192, 41},
        {344064, 8192, 42},
        {352256, 8192, 43},
        {360448, 8192, 44},
        {368640, 8192, 45},
        {376832, 16384, 23},
        {393216, 16384, 24},
        {409600, 16384, 25},
        {425984, 16384, 26},
        {442368, 16384, 27},
        {458752, 16384, 28},
        {475136, 16384, 29},
        {491520, 16384, 30},
        {507904, 16384, 31},
        {524288, 16384, 32},
        {540672, 16384, 33},
        {557056, 16384, 34},
        {573440, 16384, 35},
        {589824, 16384, 36},
        {606208, 16384, 37},
        {622592, 16384, 38},
        {638976, 16384, 39},
        {655360, 16384, 40},
        {671744, 16384, 41},
        {688128, 16384, 42},
        {704512, 16384, 43},
        {720896, 16384, 44},
        {753664, 32768, 23},
        {786432, 32768, 24},
        {819200, 32768, 25},
        {851968, 32768, 26},
        {884736, 32768, 27},
        {917504, 32768, 28},
        {950272, 32768, 29},
        {983040, 32768, 30},
        {1015808, 32768, 31},
        {1048576, 32768, 32},
        {1081344, 32768, 33},
        {1114112, 32768, 34},
        {1146880, 32768, 35},
        {1179648, 32768, 36},
        {1212416, 32768, 37},
        {1245184, 32768, 38},
        {1277952, 32768, 39},
        {1310720, 32768, 40},
        {1343488, 32768, 41},
        {1376256, 32768, 42},
        {1409024, 32768, 43},
        {1441792, 65536, 22},
        {1507328, 65536, 23},
        {1572864, 65536, 24},
        {1638400, 65536, 25},
        {1703936, 65536, 26},
        {1769472, 65536, 27},
        {1835008, 65536, 28},
        {1900544, 65536, 29},
        {1966080, 65536, 30},
        {2031616, 65536, 31},
        {2097152, 65536, 32},
        {2162688, 65536, 33},
        {2228224, 65536, 34},
        {2293760, 65536, 35},
        {2359296, 65536, 36},
        {2424832, 65536, 37},
        {2490368, 65536, 38},
        {2555904, 65536, 39},
        {2621440, 65536, 40},
        {2686976, 65536, 41},
        {2752512, 65536, 42},
        {2883584, 131072, 22},
        {3014656, 131072, 23},
        {3145728, 131072, 24},
        {3276800, 131072, 25},
        {3407872, 131072, 26},
        {3538944, 131072, 27},
        {3670016, 131072, 28},
        {3801088, 131072, 29},
        {3932160, 131072, 30},
        {4063232, 131072, 31},
        {4194304, 131072, 32},
    };
    constexpr auto check_base_table = [] {
        for (const auto [N, K, M] : base_table) {
            static_assert(N == K * M, "N must be the product of K and M");
            static_assert(N > 0, "N must be positive");
            static_assert(M > 0, "M must be positive");
            static_assert(K > 0 && ((K & (K - 1)) == 0), "K must be a power of two");
            static_assert(M <= sizeof(sfixn) * 8, "M must be less than or equal to the size of a machine word".);
        }
    };
    check_base_table();

    for (size_t i = 0; i < base_table.size(); ++i) {
        if (largest_bit_width <= base_table[i].N && (i + 1 >= base_table.size() || largest_bit_width > base_table[i + 1].N))
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

thrust::device_vector<sfixn> modular_sum_and_difference_dev(thrust::device_vector<sfixn>& u,
                                                            const thrust::device_vector<sfixn>& v,
                                                            sfixn prime)
{
    assert(u.size() == v.size());
    thrust::device_vector<sfixn> modular_sum(u.size());
    thrust::transform(u.begin(), u.end(),
                      v.begin(), v.end(),
                      modular_sum.begin(),
                      [prime] __device__ (sfixn a, sfixn b) {
                          sfixn sum = a + b;
                          return (sum >= prime) ? sum - prime : sum;
                      });
    thrust::transform(u.begin(), u.end(),
                      v.begin(), v.end(),
                      u.begin(),
                      [prime] __device__ (sfixn a, sfixn b) {
                          sfixn diff = a - b;
                          return (diff < 0) ? diff + prime : diff;
                      });
    return modular_sum;
}

// From https://github.com/orcca-uwo/BPAS/blob/7a7ce86819eba3f21098b6cc71a9762ca443d871/src/IntegerPolynomial/Multiplication/MulSSA-64_bit_arithmetic.cpp#L384
void reconstruct_mpz_with_crt(mpz_t zp, sfixn *a, sfixn *b, int size, int prime1, int prime2)
{
    mpz_t zn;
    mpz_init2(zp, LIMB_BITS);
    mpz_init2(zn, LIMB_BITS);

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
            diff += P1;
        }
        __int128 elem = MontMulModSpe_OPT3_AS_GENE_globalfunc(diff, U2_R1_sft, INV_PRIME1, MY_PRIME1);
        elem = elem * P2 + b[i];
        if (elem > HALF_P1_P2)
        {
            elem -= P1_P2;
        }
        else if (elem < N_HALF_P1_P2)
        {
            elem += P1_P2;
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

UnivariateMPZPolynomial recover_product_dev(
    const thrust::device_vector<sfixn>& sum1,
    const thrust::device_vector<sfixn>& sum2,
    const thrust::device_vector<sfixn>& diff1,
    const thrust::device_vector<sfixn>& diff2,
    const BivariateBase& base,
    sfixn prime1,
    sfixn prime2)
{
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
            base.K,
            prime1,
            prime2
        );
        reconstruct_mpz_with_crt(
            uv_diff,
            thrust::raw_pointer_cast(diff1.data() + offset),
            thrust::raw_pointer_cast(diff2.data() + offset),
            base.K,
            prime1,
            prime2
        );

        mpz_mul_2exp(uv_diff, uv_diff, base.N);
        mpz_add(product[i].get_mpz_t(), uv_sum, uv_diff);
        product[i] >>= 1;
        mpz_clear(uv_sum);
        mpz_clear(uv_diff);
    }
    return product;
}

struct TwoConvolutionResult {
    thrust::device_vector<sfixn> cyclic_convolution;
    thrust::device_vector<sfixn> negacyclic_convolution;
};

TwoConvolutionResult two_convolution_2d_dev(const thrust::device_vector<sfixn>& A, const thrust::device_vector<sfixn>& B, const BivariateBase& base, sfixn prime)
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
    const sfixn *A_ptr = thrust::raw_pointer_cast(A.data());
    sfixn *padded_A_ptr = thrust::raw_pointer_cast(padded_A.data());
    expand_to_fft2_dev(product_x_num_bits, product_y_num_bits, padded_A_ptr, base.K, A.size() / base.K, A_ptr);

    thrust::device_vector<sfixn> padded_B(padded_product_x_size * padded_product_y_size, 0);
    const sfixn *B_ptr = thrust::raw_pointer_cast(B.data());
    sfixn *padded_B_ptr = thrust::raw_pointer_cast(padded_B.data());
    expand_to_fft2_dev(product_x_num_bits, product_y_num_bits, padded_B_ptr, base.K, B.size() / base.K, B_ptr);

    // The output is written to padded_A, so we store a copy of it for when we compute the negacyclic convolution
    thrust::device_vector<sfixn> padded_A_copy(padded_A);
    bi_stockham_poly_mul_dev(padded_product_y_size, product_y_num_bits, padded_product_x_size, product_x_num_bits, padded_A_ptr, padded_B_ptr, prime);

    // theta is a 2Kth primitive root of unity
    sfixn theta = primitive_root(product_x_num_bits + 1, prime);
    scale_x_argument_dev(padded_A_copy, base, theta, prime);
    scale_x_argument_dev(padded_B, base, theta, prime);

    // padded_A_copy will store the result of the negacyclic convolution
    sfixn *padded_A_copy_ptr = thrust::raw_pointer_cast(padded_A_copy.data());
    bi_stockham_poly_mul_dev(padded_product_y_size, product_y_num_bits, padded_product_x_size, product_x_num_bits, padded_A_copy_ptr, padded_B_ptr, prime);

    return TwoConvolutionResult {
        .cyclic_convolution = std::move(padded_A),
        .negacyclic_convolution = std::move(padded_A_copy),
    };
}

UnivariateMPZPolynomial two_convolution_poly_mul(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{
    BivariateBase base {determine_bivariate_base(find_largest_bit_width_of_coefficients_dev(copy_polynomial_data_to_device(a, b)))};
    assert(base.K * base.M == base.N);

    // Convert the univariate polynomials a and b to bivariate polynomials.

    // assuming a machine word is 64 bit
    constexpr sfixn prime1 = 4179340454199820289;
    constexpr sfixn prime2 = 2485986994308513793;
    constexpr auto is_proth_prime = [] (sfixn n) {
        sfixn q = 1;
        while (q < n) {
            sfixn two_k = (prime - 1) / q;
            if ((two_k & (two_k - 1)) == 0) {
                return true;
            }
            ++q;
        }
        return false;

        // check if n is prime
        if (n < 2) return false;
        if (n == 2 || n == 3) return true;
        if (n % 2 == 0) return false;
        for (sfixn i = 3; i * i <= n; i += 2) {
            if (n % i == 0) return false;
        }
        return true;
    };
    static_assert(is_proth_prime(prime1), "prime1 is not a Proth prime");
    static_assert(is_proth_prime(prime2), "prime2 is not a Proth prime");

    struct UVSumAndDifference {
        thrust::device_vector<sfixn> uv_sum;
        thrust::device_vector<sfixn> uv_diff;
    };

    auto compute_uv_sum_and_diff = [&] (sfixn prime) -> UVSumAndDifference {
        BivariateMPZPolynomial a_bivariate {convert_to_modular_bivariate(a, base, prime)};
        BivariateMPZPolynomial b_bivariate {convert_to_modular_bivariate(b, base, prime)};
        TwoConvolutionResult result {two_convolution_2d_dev(a_bivariate, b_bivariate, base, prime)};
        thrust::device_vector<sfixn> uv_sum = modular_sum_and_difference_dev(result.cyclic_convolution, result.negacyclic_convolution, prime);
        thrust::device_vector<sfixn> uv_diff = modular_sum_and_difference_dev(result.cyclic_convolution, result.negacyclic_convolution, prime);
        return UVSumAndDifference {
            .uv_sum = std::move(uv_sum),
            .uv_diff = std::move(uv_diff)
        };
    };

    UVSumAndDifference uv_sum_and_diff1 = compute_uv_sum_and_diff(prime1);
    UVSumAndDifference uv_sum_and_diff2 = compute_uv_sum_and_diff(prime2);
    return recover_product_dev(
        uv_sum_and_diff1.uv_sum,
        uv_sum_and_diff2.uv_sum,
        uv_sum_and_diff1.uv_diff,
        uv_sum_and_diff2.uv_diff,
        base,
        prime_word1,
        prime_word2
    );

}
