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

int find_largest_bit_width_of_coefficients(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
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

constexpr BivariateBase determine_bivariate_base(sfixn largest_bit_width)
{
    using TwoConvolutionConstants::base_table;

    for (size_t i = 0; i < base_table.size(); ++i) {
        if (largest_bit_width <= base_table[i].N && base_table[i].K >= 32)
            return base_table[i];
    }
    assert(false && "No suitable base found, the largest bit width exceeds the maximum base size");
}

BivariateMPZPolynomial convert_to_modular_bivariate(const UnivariateMPZPolynomial& p, const BivariateBase& base, sfixn prime)
{
    assert(base.K * base.M == base.N);
    std::vector<sfixn> bi(p.size() * base.K);
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
    
    return bi;
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
    assert(u.size() == v.size());
    thrust::device_vector<sfixn> modular_sum(u.size());
    thrust::transform(u.begin(), u.end(),
                      v.begin(),
                      modular_sum.begin(),
                      [prime] __device__ (sfixn a, sfixn b) -> sfixn {
                          sfixn sum = a + b;
                          return (sum >= prime) ? sum - prime : sum;
                      });
    
    thrust::device_vector<sfixn> modular_difference(u.size());
    thrust::transform(u.begin(), u.end(),
                      v.begin(),
                      modular_difference.begin(),
                      [prime] __device__ (sfixn a, sfixn b) -> sfixn {
                          sfixn diff = a - b;
                          return (diff < 0) ? diff + prime : diff;
                      });
    return {
        .modular_sum = std::move(modular_sum),
        .modular_difference = std::move(modular_difference)
    };
}

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
    return product;
}




TwoConvolutionResult two_convolution_2d_dev(const thrust::device_vector<sfixn>& A, const thrust::device_vector<sfixn>& B, const BivariateBase& base, sfixn prime)
{
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

    thrust::device_vector<sfixn> padded_B(padded_product_x_size * padded_product_y_size);
    const sfixn *B_ptr = thrust::raw_pointer_cast(B.data());
    sfixn *padded_B_ptr = thrust::raw_pointer_cast(padded_B.data());
    expand_to_fft2_dev(log_padded_product_x_size, log_padded_product_y_size, padded_B_ptr, product_x_size, B.size() / base.K, B_ptr);

    // The output is written to padded_A, so we store a copy of it for when we compute the negacyclic convolution
    thrust::device_vector<sfixn> padded_A_copy(padded_A);
    thrust::device_vector<sfixn> padded_B_copy(padded_B);
    bi_stockham_poly_mul_dev(padded_product_y_size, log_padded_product_y_size, padded_product_x_size, log_padded_product_x_size, padded_A_ptr, padded_B_ptr, prime);
    // TODO: bi_stockham_poly_mul_dev perfomrs FFT on B here, maybe we can reuse this result for NCC 
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

    return TwoConvolutionResult {
        .cyclic_convolution = std::move(extracted_A),
        .negacyclic_convolution = std::move(extracted_A_copy),
    };
}

UnivariateMPZPolynomial two_convolution_poly_mul(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{
    BivariateBase base {determine_bivariate_base(find_largest_bit_width_of_coefficients(a, b))};
    assert(base.K * base.M == base.N);

    // Convert the univariate polynomials a and b to bivariate polynomials.

    // assuming a machine word is 64 bit

    auto compute_uv_sum_and_diff = [&] (sfixn prime) -> ModularSumAndDifferenceResult {
        BivariateMPZPolynomial a_bivariate {convert_to_modular_bivariate(a, base, prime)};
        BivariateMPZPolynomial b_bivariate {convert_to_modular_bivariate(b, base, prime)};
        const TwoConvolutionResult& two_conv = two_convolution_2d_dev(a_bivariate, b_bivariate, base, prime);
        const ModularSumAndDifferenceResult& mod_sum_and_diff = modular_sum_and_difference_dev(two_conv.cyclic_convolution, two_conv.negacyclic_convolution, prime);
        return mod_sum_and_diff;
    };

    using namespace TwoConvolutionConstants;
    ModularSumAndDifferenceResult uv_sum_and_diff1 = compute_uv_sum_and_diff(prime1);
    ModularSumAndDifferenceResult uv_sum_and_diff2 = compute_uv_sum_and_diff(prime2);
    ModularSumAndDifferenceResult uv_sum_and_diff3 = compute_uv_sum_and_diff(prime3);
    auto bitwidth = [](sfixn x) {
        return sizeof(sfixn) * 8 - __builtin_clz(x);
    };
    auto result = recover_product_host(
        uv_sum_and_diff1.modular_sum,
        uv_sum_and_diff2.modular_sum,
        uv_sum_and_diff3.modular_sum,
        uv_sum_and_diff1.modular_difference,
        uv_sum_and_diff2.modular_difference,
        uv_sum_and_diff3.modular_difference,
        base,
        1 + base.M + base.N + bitwidth(base.K) + bitwidth(std::max(a.size(), b.size()))
    );
    result.resize(a.size() + b.size() - 1); // resize to the correct size
    return result;
}
