#include <gtest/gtest.h>
#include "two_convolution_poly_mul.h"
#include <thrust/host_vector.h>

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

TEST(TwoConvolutionPolyMul, FindLargestBitWidthOfCoefficientsDev)
{
    auto run_test = [](const UnivariateMPZPolynomial& a, const UnivariateMPZPolynomial& b, int expected_bit_width) {
        CoefficientsOnDevice coeffs = copy_polynomial_data_to_device(a, b);
        EXPECT_EQ(find_largest_bit_width_of_coefficients_dev(coeffs), expected_bit_width);
    };
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        UnivariateMPZPolynomial b = {4, 5, 6};
        run_test(a, b, 3);
    }
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        UnivariateMPZPolynomial b = {0x9abc, 0xdef0, 0x1234};
        run_test(a, b, 16);
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        UnivariateMPZPolynomial b = {7, 8, 9};
        run_test(a, b, 4);
    }
    {
        UnivariateMPZPolynomial a = {-1, -2, -3};
        UnivariateMPZPolynomial b = {-4, -5, -6};
        run_test(a, b, 3);
    }
    {
        UnivariateMPZPolynomial a = {0};
        UnivariateMPZPolynomial b = {0};
        run_test(a, b, 0);
    }
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        UnivariateMPZPolynomial b = {0};
        run_test(a, b, 2);
    }
    {
        UnivariateMPZPolynomial a = {1, 0};
        UnivariateMPZPolynomial b = {0xffffffff, 0};
        run_test(a, b, 32);
    }
    {
        UnivariateMPZPolynomial a = {mpz_class("ffffffffffffffffffffffffffffffffffffffff", 16), 0};
        UnivariateMPZPolynomial b = {1, 0};
        run_test(a, b, 160);
    }
}

TEST(TwoConvolutionPolyMul, ConvertToModularBivariateNaive)
{
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({1, 0, 2, 0, 3, 0}));
    }
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        BivariateBase base = {.N = 16, .K = 2, .M = 8};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({0xcd, 0xab, 0x34, 0x12, 0x78, 0x56}));
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        EXPECT_EQ(convert_to_modular_bivariate(a, base, 401), BivariateMPZPolynomial({0, 1, 1, 1, 2, 1}));
    }
}

TEST(TwoConvolutionPolyMul, ConvertToModularBivariateDev)
{
    auto run_test = [](const UnivariateMPZPolynomial& a, const BivariateBase& base, int prime, const std::vector<sfixn>& expected) {
        CoefficientsOnDevice coeffs = copy_polynomial_data_to_device(a, {});
        thrust::host_vector<sfixn> result = convert_to_modular_bivariate_dev(coeffs, base, prime);
        thrust::host_vector<sfixn> expected_result = expected;
        EXPECT_EQ(result, expected_result);
    };
    {
        UnivariateMPZPolynomial a = {1, 2, 3};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        run_test(a, base, 401, {1, 0, 2, 0, 3, 0});
    }
    {
        UnivariateMPZPolynomial a = {0xabcd, 0x1234, 0x5678};
        BivariateBase base = {.N = 16, .K = 2, .M = 8};
        run_test(a, base, 401, {0xcd, 0xab, 0x34, 0x12, 0x78, 0x56});
    }
    {
        UnivariateMPZPolynomial a = {4, 5, 6};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        run_test(a, base, 401, {0, 1, 1, 1, 2, 1});
    }
}

TEST(TwoConvolutionPolyMul, ScaleXArgumentDev)
{
    auto run_test = [](const BivariateMPZPolynomial& A, const BivariateBase& base, sfixn theta, sfixn prime) {
        thrust::device_vector<sfixn> A_dev(A.begin(), A.end());
        scale_x_argument_dev(A_dev, base, theta, prime);
        thrust::host_vector<sfixn> result = A_dev;
        return result;
    };
    {
        BivariateMPZPolynomial A = {1, 0, 2, 0, 3, 0};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        sfixn theta = 400;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, theta, prime), (BivariateMPZPolynomial{1, 0, 2, 0, 3, 0}));
    }
    {
        BivariateMPZPolynomial A = {1, 1, 1, 2, 2, 2, 4, 4, 4};
        BivariateBase base = {.N = 6, .K = 3, .M = 2};
        sfixn theta = 2;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, theta, prime), (BivariateMPZPolynomial{1, 2, 4, 2, 4, 8, 4, 8, 16}));
    }
    {
        BivariateMPZPolynomial A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        BivariateBase base = {.N = 6, .K = 3, .M = 2};
        sfixn theta = 2;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, theta, prime), (BivariateMPZPolynomial{1, 4, 12, 4, 10, 24, 7, 16, 36}));
    }
}

TEST(TwoConvolutionPolyMul, EvaluateAtXDev)
{
    auto run_test = [](const BivariateMPZPolynomial& A, const BivariateBase& base, sfixn beta, sfixn prime) {
        thrust::device_vector<sfixn> A_dev(A.begin(), A.end());
        evaluate_at_x_dev(A_dev, base, beta, prime);
        thrust::host_vector<sfixn> result = A_dev;
        return result;
    };
    {
        BivariateMPZPolynomial A = {1, 0, 2, 0, 3, 0};
        BivariateBase base = {.N = 4, .K = 2, .M = 2};
        sfixn beta = 1;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, beta, prime), (BivariateMPZPolynomial{1, 2, 3}));
    }
    {
        BivariateMPZPolynomial A = {1, 1, 1, 2, 2, 2};
        BivariateBase base = {.N = 6, .K = 3, .M = 2};
        sfixn beta = 2;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, beta, prime), (BivariateMPZPolynomial{7, 14}));
    }
    {
        BivariateMPZPolynomial A = {1, 2, 3};
        BivariateBase base = {.N = 6, .K = 3, .M = 2};
        sfixn beta = -1;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, beta, prime), (BivariateMPZPolynomial{2}));
    }
    {
        BivariateMPZPolynomial A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        BivariateBase base = {.N = 12, .K = 3, .M = 4};
        sfixn beta = 2;
        sfixn prime = 401;
        EXPECT_EQ(run_test(A, base, beta, prime), (BivariateMPZPolynomial{17, 38, 59}));
    }
}

TEST(TwoConvolutionPolyMul, ReconstructBivariatePolynomialWithCRTDev)
{
    auto run_test = [](const BivariateMPZPolynomial& p1, const BivariateMPZPolynomial& p2, const BivariateBase& base, sfixn prime1, sfixn prime2) {
        thrust::device_vector<sfixn> result(p1);
        reconstruct_bivariate_polynomial_with_crt_dev(result, p2, base, prime1, prime2);
        return result;
    };
    {
        BivariateMPZPolynomial p1 = {1, 4, 1, 5, 2, 3, 5, 2, 3, 4, 3, 3, 2, 0, 5, 4};
        BivariateMPZPolynomial p2 = {4, 5, 2, 9, 9, 6, 6, 3, 10, 5, 6, 7, 9, 3, 6, 8};
        BivariateBase base = {.N = 16, .K = 4, .M = 4};
        sfixn prime1 = 7;
        sfixn prime2 = 11
        thrust::device_vector<sfixn> expected = {15, 60, 57, 75, 9, 17, 61, 58, 10, 60, 17, 73, 9, 14, 61, 74};
        EXPECT_EQ(run_test(p1, p2, base, prime1, prime2) == ;
    }
}

TEST(TwoConvolutionPolyMul, RecoverProductDev)
{
    auto run_test = [](const UnivariateMPZPolynomial& cyclic_convolution, const UnivariateMPZPolynomial& negacyclic_convolution, sfixn largest_bit_width_coefficient) {
        thrust::device_vector<sfixn> d_cyclic_convolution(cyclic_convolution.begin(), cyclic_convolution.end());
        const thrust::device_vector<sfixn> d_negacyclic_convolution(negacyclic_convolution.begin(), negacyclic_convolution.end());
        recover_product_dev(d_cyclic_convolution, d_negacyclic_convolution, largest_bit_width_coefficient);
        thrust::host_vector<sfixn> result = d_cyclic_convolution;
        return result;
    };
    {
        UnivariateMPZPolynomial cyclic_convolution = {1, 2, 3};
        UnivariateMPZPolynomial negacyclic_convolution = {4, 5, 6};
        sfixn largest_bit_width_coefficient = 3;
        EXPECT_EQ(run_test(cyclic_convolution, negacyclic_convolution, largest_bit_width_coefficient), (UnivariateMPZPolynomial{10, 11, 12}));
    }
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