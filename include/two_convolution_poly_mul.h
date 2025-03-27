#pragma once

#include <gmpxx.h>
#include <vector>
#include "types.h"

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

int find_largest_bit_width_of_coefficients(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b);
BivariateMPZPolynomial convert_to_modular_bivariate(const UnivariateMPZPolynomial& p, const BivariateBase& base, sfixn prime);