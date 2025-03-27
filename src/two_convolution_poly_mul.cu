#include <cassert>
#include "two_convolution_poly_mul.h"

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

BivariateMPZPolynomial convert_to_modular_bivariate(const UnivariateMPZPolynomial& p, const BivariateBase& base, sfixn prime)
{
    assert(base.K * base.M == base.N);
    BivariateMPZPolynomial bi(p.size() * base.K);
    const int block_size {base.M};
    const int y_terms = p.size();

    auto convert_mpz_to_modular_univariate = [&](int y_power) {
        int x_power = 0;
        mpz_srcptr raw_mpz {p[y_power].get_mpz_t()};
        size_t num_limbs {mpz_size(raw_mpz)};
        
        size_t current_block_bits = 0;
        sfixn current_block = 0;
        
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
                    bi[y_power * base.M + x_power++] = current_block % prime;
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
            bi[y_power * base.M + x_power++] = current_block % prime;
        }
    };

    for (int y_power = 0; y_power < y_terms; ++y_power)
        convert_mpz_to_modular_univariate(y_power);
    return bi;
}

/*
UnivariateMPZPolynomial two_convolution_poly_mul(const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b)
{
    // assuming a machine word is 64 bit
    constexpr sfixn prime_word1 = 4179340454199820289;
    constexpr sfixn prime_word2 = 2485986994308513793;
    
    // Find the largest bit-width of any coefficient in a or b.
    // Note that it suffices to compare the number of GMP limbs and the most significant limb in each coefficient.
    const int largest_bit_width_of_coefficients {find_largest_bit_width_of_coefficients(a, b)};

    BivariateBase base {determine_bivariate_base(largest_bit_width_of_coefficients)};
    assert(base.K * base.M == base.N);

    // Convert the univariate polynomials a and b to bivariate polynomials.
    BivariateMPZPolynomial a_bivariate1 {convert_to_modular_bivariate(a, base, prime_word1)};
    BivariateMPZPolynomial a_bivariate2 {convert_to_modular_bivariate(a, base, prime_word2)};

    BivariateMPZPolynomial b_bivariate1 {convert_to_modular_bivariate(b, base, prime_word1)};
    BivariateMPZPolynomial b_bivariate2 {convert_to_modular_bivariate(b, base, prime_word2)};

    BivariateMPZPolynomial c_minus1 {cyclic_convolution(a_bivariate1, b_bivariate1, base, prime_word1)};
    BivariateMPZPolynomial c_minus2 {cyclic_convolution(a_bivariate2, b_bivariate2, base, prime_word2)};

    BivariateMPZPolynomial c_plus1 {negacyclic_convolution(a_bivariate1, b_bivariate1, base, prime_word1)};
    BivariateMPZPolynomial c_plus2 {negacyclic_convolution(a_bivariate2, b_bivariate2, base, prime_word2)};

}
*/
