#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include "two_convolution_poly_mul.h"

#include <vector>
#include <gmpxx.h>
#include <algorithm>
#include <cassert>

// Multiply two univariate integer polynomials A and B represented as std::vector<mpz_class>
// using Kronecker substitution. Returns the product polynomial.
UnivariateMPZPolynomial kronecker_poly_mul(const UnivariateMPZPolynomial& A,
                                     const UnivariateMPZPolynomial& B) {
    // Sizes and result length.
    size_t nA = A.size(), nB = B.size();
    size_t resSize = nA + nB - 1;

    // STEP 1. Compute shifts to make coefficients nonnegative.
    mpz_class shiftA = 0, shiftB = 0;
    for (const auto& a : A) {
        if (a < 0 && (-a) > shiftA)
            shiftA = -a;
    }
    for (const auto& b : B) {
        if (b < 0 && (-b) > shiftB)
            shiftB = -b;
    }
    // Create shifted copies: Aprime[i] = A[i] + shiftA, and similarly for B.
    std::vector<mpz_class> Aprime(nA), Bprime(nB);
    for (size_t i = 0; i < nA; i++)
        Aprime[i] = A[i] + shiftA;
    for (size_t j = 0; j < nB; j++)
        Bprime[j] = B[j] + shiftB;
    
    // STEP 2. Find maximum coefficients in the shifted polynomials.
    mpz_class maxA_prime = 0, maxB_prime = 0;
    for (size_t i = 0; i < nA; i++) {
        if (Aprime[i] > maxA_prime)
            maxA_prime = Aprime[i];
    }
    for (size_t j = 0; j < nB; j++) {
        if (Bprime[j] > maxB_prime)
            maxB_prime = Bprime[j];
    }

    // Choose base r such that every product digit is less than r.
    // A safe choice is:
    //    r = 1 + nA * nB * (maxA_prime * maxB_prime)
    mpz_class r = 1 + (mpz_class(nA) * nB * maxA_prime * maxB_prime);

    // STEP 3. Pack Aprime and Bprime into GMP integers by evaluating at x = r.
    mpz_class valA = 0, power = 1;
    for (size_t i = 0; i < nA; i++) {
        valA += Aprime[i] * power;
        power *= r;
    }
    mpz_class valB = 0;
    power = 1;
    for (size_t j = 0; j < nB; j++) {
        valB += Bprime[j] * power;
        power *= r;
    }

    // STEP 4. Multiply the two large integers.
    mpz_class valC = valA * valB;

    // STEP 5. Unpack the coefficients for Cprime(x) = A'(x)*B'(x)
    std::vector<mpz_class> Cprime(resSize, 0);
    for (size_t k = 0; k < resSize; k++) {
        mpz_class coeff = valC % r; // The coefficient of x^k.
        Cprime[k] = coeff;
        valC /= r;
    }

    // STEP 6. Correct for the shifts.
    // We need to compute the convolution with "all-ones" vectors.
    // Let F_A(x) = 1 + x + ... + x^(nA-1) and F_B(x) = 1 + x + ... + x^(nB-1).
    // Compute conv1 = (F_A * Bprime)(x) and conv2 = (F_B * Aprime)(x)
    // and conv3 = (F_A * F_B)(x). Since F_A and F_B are vectors of ones,
    // these can be computed by simple loops.
    std::vector<mpz_class> conv1(resSize, 0), conv2(resSize, 0), conv3(resSize, 0);

    for (size_t k = 0; k < resSize; k++) {
        // conv1[k] = sum_{i=0}^{nA-1} (Bprime[k-i]) for valid indices.
        for (size_t i = 0; i < nA; i++) {
            if (k >= i && (k - i) < nB)
                conv1[k] += Bprime[k - i];
        }
        // conv2[k] = sum_{j=0}^{nB-1} (Aprime[k-j]) for valid indices.
        for (size_t j = 0; j < nB; j++) {
            if (k >= j && (k - j) < nA)
                conv2[k] += Aprime[k - j];
        }
        // conv3[k] = number of pairs (i, j) with i in [0, nA-1], j in [0, nB-1] and i+j == k.
        // This is equivalent to: max(0, min(k, nA-1) - max(0, k - (nB-1)) + 1).
        size_t lower = (k < nB) ? 0 : k - (nB - 1);
        size_t upper = (k < nA) ? k : nA - 1;
        if (upper >= lower)
            conv3[k] = upper - lower + 1;
        else
            conv3[k] = 0;
    }

    // Now, recover the actual product polynomial.
    // Note that:
    //   A(x) = A'(x) - shiftA * F_A(x)
    //   B(x) = B'(x) - shiftB * F_B(x)
    // Hence,
    //   A(x)*B(x) = A'(x)*B'(x)
    //               - shiftA*(F_A*B')(x)
    //               - shiftB*(F_B*A')(x)
    //               + shiftA*shiftB*(F_A*F_B)(x)
    std::vector<mpz_class> result(resSize, 0);
    for (size_t k = 0; k < resSize; k++) {
        result[k] = Cprime[k]
                    - shiftA * conv1[k]
                    - shiftB * conv2[k]
                    + shiftA * shiftB * conv3[k];
    }
    return result;
}

template<typename PolyMulFunc>
auto time_poly_mul(
    const std::string& input_file,
    PolyMulFunc poly_mul)
{
    std::ifstream input(input_file);
    if (!input)
        throw std::runtime_error("Failed to open input file: " + input_file);
    
    auto read_polynomial_from_file = [](std::ifstream& input) {
        std::string line;
        if (!std::getline(input, line))
            throw std::runtime_error("Failed to read polynomial from input file: the input file should have 3 lines, one for each polynomial. The 3rd line is the expected product of the first 2 lines.");
        UnivariateMPZPolynomial poly;
        std::istringstream iss(line);
        mpz_class coeff;
        while (iss >> coeff)
            poly.push_back(coeff);
        return poly;
    };
    UnivariateMPZPolynomial a, b, c;
    a = read_polynomial_from_file(input);
    b = read_polynomial_from_file(input);
    c = read_polynomial_from_file(input);
    auto start = std::chrono::high_resolution_clock::now();
    UnivariateMPZPolynomial result = poly_mul(a, b);
    auto end = std::chrono::high_resolution_clock::now();
    if (result != c) {
        std::ofstream output (input_file + ".actual");
        for (const auto& coeff : result)
            output << coeff << " ";
        output.flush();
        throw std::runtime_error("The result of polynomial multiplication does not match the expected product.");
    }
    auto duration = end - start;
    return duration;
}

int main(int argc, char* argv[])
{
    auto print_summary = [](const std::string& input_name, const std::string& algo_name, auto duration) {
        std::cout << input_name << ',' 
                  << algo_name << ','
                  << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " ms" << std::endl;
    };

    if (argc > 1) {
        auto duration = time_poly_mul(argv[1], kronecker_poly_mul);
        print_summary(argv[1], "kronecker_poly_mul", duration);
        duration = time_poly_mul(argv[1], two_convolution_poly_mul);
        print_summary(argv[1], "two_convolution_poly_mul", duration);
        return 0;
    }

    std::vector<std::string> input_files = {
        "n128_m128_p10000000000.test",
        "n256_m256_p10000000000.test",
        "n512_m512_p10000000000.test",
        "n1024_m1024_p10000000000.test",
        "n2048_m2048_p10000000000.test",
        "n4096_m4096_p10000000000.test",
    };

    for (const auto& input_file : input_files) {
        auto duration = time_poly_mul(input_file, kronecker_poly_mul);
        print_summary(input_file, "kronecker_poly_mul", duration);
        duration = time_poly_mul(input_file, two_convolution_poly_mul);
        print_summary(input_file, "two_convolution_poly_mul", duration);
    }
}