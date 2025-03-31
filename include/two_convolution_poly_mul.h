#pragma once

#include <gmpxx.h>
#include <vector>
#include "types.h"
#include <thrust/device_vector.h>

using UnivariateMPZPolynomial = std::vector<mpz_class>;
using BivariateMPZPolynomial = std::vector<sfixn>;

// The following information is used to determine how each coefficient in the univariate polynomial P
// is decomposed as itself a univariate polynomial Q, to convert P from a univariate polynomial to a
// bivariate polynmoial.
struct BivariateBase {
    // The number of bits used to represent each coefficient in the univariate polynomial.
    int N;

    // The number of limbs used to represent the univariate coefficient. This is the partial
    // degree of the bivariate polynomial with respect to x.
    int K;

    // The bit-width of each limb in the univariate coefficient. This should be less than a machine
    // word. Note that this may be different from the size of a limb in the GMP library.
    int M;
};

struct CoefficientsOnDevice {
    thrust::device_vector<size_t> mpz_sizes;
    thrust::device_vector<mp_limb_t> mpz_limbs;
    thrust::device_vector<mp_limb_t> most_significant_mpz_limbs;
};

CoefficientsOnDevice copy_polynomial_data_to_device(const UnivariateMPZPolynomial &a, const UnivariateMPZPolynomial &b);

int find_largest_bit_width_of_coefficients_naive(const UnivariateMPZPolynomial &, const UnivariateMPZPolynomial &);
int find_largest_bit_width_of_coefficients_dev(const CoefficientsOnDevice &);

BivariateMPZPolynomial convert_to_modular_bivariate(const UnivariateMPZPolynomial&, const BivariateBase&, sfixn prime);
thrust::device_vector<sfixn> convert_to_modular_bivariate_dev(const CoefficientsOnDevice&, const BivariateBase&, sfixn prime);

void scale_x_argument_dev(thrust::device_vector<sfixn>&, const BivariateBase& base, sfixn theta, sfixn prime);
void evaluate_at_x_dev(thrust::device_vector<sfixn>&, const BivariateBase&, sfixn beta, sfixn prime);

thrust::device_vector<sfixn> two_convolution_2d_dev(const thrust::device_vector<sfixn>& A, const thrust::device_vector<sfixn>& B, const BivariateBase&, sfixn prime);