#include <gtest/gtest.h>
#include "two_convolution_poly_mul.h"
#include <iostream>
#include <exception>

TEST(TwoConvolutionPolyMul, FindLargestBitWidthOfCoefficientsNaive)
{
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_naive(a, a), 2);
    }
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_naive(a, a), 16);
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_naive(a, a), 3);
    }
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        UnivariateMPZPolynomial b = {4, 5, 6};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_naive(a, b), 3);
    }
}

TEST(TwoConvolutionPolyMul, FindLargestBitWidthOfCoefficientsHost)
{
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        UnivariateMPZPolynomial b = {4, 5, 6};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 3);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 3);
    }
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        UnivariateMPZPolynomial b = {0x9abc, 0xdef0, 0x1234};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 16);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 16);
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        UnivariateMPZPolynomial b = {7, 8, 9};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 4);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 4);
    }
    {
        UnivariateMPZPolynomial a = {-1, -2, -3};
        UnivariateMPZPolynomial b = {-4, -5, -6};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 3);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 3);
    }
    {
        UnivariateMPZPolynomial a = {0};
        UnivariateMPZPolynomial b = {0};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 0);
    }
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        UnivariateMPZPolynomial b = {0};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 2);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 2);
    }
    {
        UnivariateMPZPolynomial a = {1, 0};
        UnivariateMPZPolynomial b = {0xffffffff, 0};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 32);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 32);
    }
    {
        UnivariateMPZPolynomial a = {mpz_class("ffffffffffffffffffffffffffffffffffffffff", 16), 0};
        UnivariateMPZPolynomial b = {1, 0};
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(a, b), 160);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_host(b, a), 160);
    }
}

TEST(TwoConvolutionPolyMul, ConvertUnivariateToBivariate)
{
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({1, 0, 2, 0, 3, 0}));
    }/*
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        BivariateBase base = {.N = 16, .K = 2, .M = 8};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({0xcd, 0xab, 0x34, 0x12, 0x78, 0x56}));
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({0, 1, 1, 1, 2, 1}));
    }*/
}

/*
TEST(TwoConvolutionPolyMul, ComputeRecoveryPrime)
{
    auto primePropertiesHold = [](int64_t p, size_t d, int K) -> bool {
        return p > 0 && (p - 1) % K == 0 && (p - 1) % (2 * d - 1) == 0;
    };
    {
        const int d = 3;
        const int K = 2;
        const int M = 2;
        const int64_t p = computeRecoveryPrime(d, K, M);
        EXPECT_TRUE(primePropertiesHold(p, d, K));
        EXPECT_EQ(p, 401);
    }
    {
        const int d = 3;
        const int K = 2;
        const int M = 8;
        const int64_t p = computeRecoveryPrime(d, K, M);
        EXPECT_TRUE(primePropertiesHold(p, d, K));
        EXPECT_EQ(p, 1573081);
    }
}

TEST(TwoConvolutionPolyMul, BinaryModularPower)
{
    EXPECT_EQ(binaryModularPower(2, 3, 5), 3);
    EXPECT_EQ(binaryModularPower(2, 3, 7), 1);
    EXPECT_EQ(binaryModularPower(2, 3, 11), 8);
}

TEST(TwoConvolutionPolyMul, ModularInverse)
{
    EXPECT_EQ(modularInverse(3, 5), 2);
    EXPECT_EQ(modularInverse(3, 7), 5);
    EXPECT_EQ(modularInverse(3, 11), 4);
}

TEST(TwoConvolutionPolyMul, FindPrimitiveRootOfUnity)
{
    EXPECT_EQ(findPrimitiveRootOfUnity(16, 17), 3);
    EXPECT_EQ(findPrimitiveRootOfUnity(5, 401), 39);
    EXPECT_EQ(findPrimitiveRootOfUnity(5, 1573081), 8621);
}

TEST(TwoConvolutionPolyMul, Dft)
{
    {
        const int64_t p = 401;
        const std::vector<int64_t> a = {1, 2, 3};
        const std::vector<int64_t> rootPowers = {1, 39, 318, 372, 72};
        const std::vector<int64_t> expected = {6, 231, 51, 60, 58};
        EXPECT_EQ(dft(a, p, rootPowers), expected);
    }
}

TEST(TwoConvolutionPolyMul, Dft2d)
{
    {
        const int64_t p = 401;
        const std::vector<int64_t> A = {1, 0, 2, 0, 3, 0};
        const std::vector<int64_t> xPowers = {1, 400};
        const std::vector<int64_t> yPowers = {1, 39, 318, 372, 72};
        const std::vector<int64_t> expected = {6, 6, 231, 231, 51, 51, 60, 60, 58, 58};
        EXPECT_EQ(dft2d(A, 2, p, xPowers, yPowers), expected);
    }
}

TEST(TwoConvolutionPolyMul, CyclicConvolution)
{
    {
        const int64_t p = 401;
        const std::vector<int64_t> A = {1, 0, 2, 0, 3, 0};
        const std::vector<int64_t> B = {1, 0, 2, 0, 3, 0};
        const std::vector<int64_t> xPowers = {1, 400};
        const std::vector<int64_t> yPowers = {1, 39, 318, 372, 72};
        const std::vector<int64_t> expected = {1, 0, 4, 0, 10, 0, 12, 0, 9, 0};
        EXPECT_EQ(cyclicConvolution(A, B, 2, p, 400, 39), expected);
    }
}

TEST(TwoConvolutionPolyMul, NegacyclicConvolution)
{
    {
        const int64_t p = 401;
        const std::vector<int64_t> A = {1, 0, 2, 0, 3, 0};
        const std::vector<int64_t> B = {1, 0, 2, 0, 3, 0};
        const std::vector<int64_t> xPowers = {1, 400};
        const std::vector<int64_t> yPowers = {1, 39, 318, 372, 72};
        const std::vector<int64_t> expected = {1, 0, 4, 0, 10, 0, 12, 0, 9, 0};
        EXPECT_EQ(negacyclicConvolution(A, B, 2, p, 400, 39), expected);
    }
}

TEST(TwoConvolutionPolyMul, Fast)
{
    std::vector<int64_t> a = {1, 2, 3};
    std::vector<int64_t> b = {1, 2, 3};
    std::vector<int64_t> expected = {1, 4, 10, 12, 9};
    EXPECT_EQ(multiplyPolynomials(a, b, 2), expected);

    std::vector<int64_t> c = {1, 2, 3};
    std::vector<int64_t> d = {4, 5, 6};
    std::vector<int64_t> expected2 = {4, 13, 28, 27, 18};
    EXPECT_EQ(multiplyPolynomials(c, d, 2), expected2);

    std::vector<int64_t> e = {1, 2, 3, 4, 5};
    std::vector<int64_t> f = {6, 7, 8, 9, 10};
    std::vector<int64_t> expected3 = {6, 19, 40, 70, 110, 114, 106, 85, 50};
    EXPECT_EQ(multiplyPolynomials(e, f, 2), expected3);
}
*/