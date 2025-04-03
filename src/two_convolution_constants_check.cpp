#include "two_convolution_poly_mul.h"

// Compute the Montgomery inverse for an odd number p modulo R = 2^k.
// That is, find INV such that p * INV ≡ -1 (mod R).
// We use Newton’s iteration which doubles the number of correct bits each time.
constexpr longfixnum montgomery_inverse(longfixnum p, sfixn k) {
    longfixnum R = (longfixnum)1 << k;
    longfixnum inv = 1; // initial guess (1 is correct modulo 2)
    // For 64-bit numbers, about 6 iterations are sufficient.
    for (int i = 0; i < 6; ++i) {
        inv = (inv * (2 - p * inv)) & (R - 1);
    }
    // We want the number such that p * ( -inv ) ≡ -1 mod R.
    return (-inv) & (R - 1);
}

// Extended Euclidean algorithm to compute the modular inverse of a modulo m.
// Assumes that a and m are coprime.
constexpr longfixnum mod_inverse(longfixnum a, longfixnum m) {
    longfixnum m0 = m;
    longfixnum x0 = 0, x1 = 1;
    while (a > 1) {
        longfixnum q = a / m;
        longfixnum t = m;
        m = a % m;
        a = t;
        longfixnum temp = x0;
        x0 = x1 - static_cast<longfixnum>(q) * x0;
        x1 = temp;
    }
    return (x1 < 0) ? static_cast<longfixnum>(x1 + m0) : static_cast<longfixnum>(x1);
}

// Convert a number x into Montgomery representation modulo p.
// That is, compute (x * R) mod p where R = 2^k.
constexpr longfixnum toMontgomery(longfixnum x, sfixn k, longfixnum p) {
    return (x * ((longfixnum)1 << k)) % p;
}

static constexpr bool is_fourier_prime(sfixn p)
{
    if (p <= 1) return false;
    if (p == 2) return true;
    if (p % 2 == 0) return false;
    
    // Check if p is a Proth number: p = k*2^n + 1 with k < 2^n
    sfixn pm1 = p - 1;
    sfixn n = 0;
    while (pm1 % 2 == 0) {
        pm1 /= 2;
        ++n;
    }
    // Now pm1 = k, and we need k < 2^n
    if (pm1 >= (1 << n)) return false;
    
    // Now verify primality (simple trial division for constexpr)
    if (pm1 == 1) return true;  // Fermat primes case
    for (sfixn i = 3; i * i <= p; i += 2) {
        if (p % i == 0) return false;
    }
    return true;
}

/* loop context isn't constepxr
static constexpr bool check_base_table()
{
    for (const auto base : TwoConvolutionConstants::base_table) {
        auto N = base.N, K = base.K, M = base.M;
        static_assert(N == K * M, "N must be the product of K and M");
        static_assert(N > 0, "N must be positive");
        static_assert(M > 0, "M must be positive");
        static_assert(K > 0 && ((K & (K - 1)) == 0), "K must be a power of two");
        static_assert(M <= sizeof(sfixn) * 8, "M must be less than or equal to the size of a machine word".);
    }
}
*/

struct TwoConvolutionConstantsCheck {
    static_assert(sizeof(sfixn) * 2 == sizeof(longfixnum), "longfixnum should fit two sfixn");


    static_assert(is_fourier_prime(TwoConvolutionConstants::prime1), "prime1 is not a Fourier prime");
    static_assert(is_fourier_prime(TwoConvolutionConstants::prime2), "prime2 is not a Fourier prime");
    static_assert(is_fourier_prime(TwoConvolutionConstants::prime3), "prime3 is not a Fourier prime");
    /*
    static_assert(check_base_table(), "Base table check failed");

    // A simple static assertion to check INV_PRIME1:
    // We want: MY_PRIME1 * INV_PRIME1 ≡ (2^K - 1) mod 2^K, i.e. ≡ -1 mod 2^K.
    static_assert((prime1 * inv_prime1) & (1LL << montgomery_power - 1) == (powerOfTwo(K) - 1),
              "inv_prime1 does not satisfy the required congruence.");

    */
    
};