#include <gtest/gtest.h>
#include <duzp.hpp>

TEST(Duzp, Naive)
{
    std::vector<int64_t> a = {1, 2, 3};
    std::vector<int64_t> b = {4, 5, 6};
    std::vector<int64_t> expected = {4, 13, 28, 27, 18};
    EXPECT_EQ(naivelyMultiplyPolynomials(a, b), expected);
}

TEST(Duzp, LargestBitSize)
{
    {
        std::vector<int64_t> a = {1, 2, 3};
        EXPECT_EQ(largestBitSize(a), 2);
    }
    {
        std::vector<int64_t> a = {0xabcd, 0x1234, 0x5678};
        EXPECT_EQ(largestBitSize(a), 16);
    }
    {
        std::vector<int64_t> a = {4, 5, 6};
        EXPECT_EQ(largestBitSize(a), 3);
    }
}

TEST(Duzp, ConvertUnivariateToBivariate)
{
    {
        std::vector<int64_t> a = {1, 2, 3};
        std::vector<int64_t> expected = {1, 0, 2, 0, 3, 0};
        EXPECT_EQ(convertUnivariateToBivariate(a, 2, 2), expected);
    }
    {
        std::vector<int64_t> a = {0xabcd, 0x1234, 0x5678};
        std::vector<int64_t> expected = {0xcd, 0xab, 0x34, 0x12, 0x78, 0x56};
        EXPECT_EQ(convertUnivariateToBivariate(a, 2, 8), expected);
    }
    {
        std::vector<int64_t> a = {4, 5, 6};
        std::vector<int64_t> expected = {0, 1, 1, 1, 2, 1};
        EXPECT_EQ(convertUnivariateToBivariate(a, 2, 2), expected);
    }
}

TEST(Duzp, ComputeRecoveryPrime)
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

TEST(Duzp, BinaryModularPower)
{
    EXPECT_EQ(binaryModularPower(2, 3, 5), 3);
    EXPECT_EQ(binaryModularPower(2, 3, 7), 1);
    EXPECT_EQ(binaryModularPower(2, 3, 11), 8);
}

TEST(Duzp, ModularInverse)
{
    EXPECT_EQ(modularInverse(3, 5), 2);
    EXPECT_EQ(modularInverse(3, 7), 5);
    EXPECT_EQ(modularInverse(3, 11), 4);
}

TEST(Duzp, FindPrimitiveRootOfUnity)
{
    EXPECT_EQ(findPrimitiveRootOfUnity(16, 17), 3);
    EXPECT_EQ(findPrimitiveRootOfUnity(5, 401), 39);
    EXPECT_EQ(findPrimitiveRootOfUnity(5, 1573081), 8621);
}

TEST(Duzp, Dft)
{
    {
        const int64_t p = 401;
        const std::vector<int64_t> a = {1, 2, 3};
        const std::vector<int64_t> rootPowers = {1, 39, 318, 372, 72};
        const std::vector<int64_t> expected = {6, 231, 51, 60, 58};
        EXPECT_EQ(dft(a, p, rootPowers), expected);
    }
}

TEST(Duzp, Dft2d)
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

TEST(Duzp, CyclicConvolution)
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

TEST(Duzp, NegacyclicConvolution)
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

TEST(Duzp, Fast)
{
    std::vector<int64_t> a = {1, 2, 3, 4};
    std::vector<int64_t> b = {1, 2, 3, 4};
    std::vector<int64_t> expected = {1, 4, 10, 20, 30, 34, 30, 20, 10, 4};
    EXPECT_EQ(multiplyPolynomials(a, b, 2), expected);

    std::vector<int64_t> c = {1, 2, 3, 4};
    std::vector<int64_t> d = {4, 5, 6};
    std::vector<int64_t> expected2 = {4, 13, 28, 27, 18};
    EXPECT_EQ(multiplyPolynomials(c, d, 2), expected2);

    /*std::vector<int64_t> e = {1, 2, 3, 4, 5};*/
    /*std::vector<int64_t> f = {6, 7, 8, 9, 10};*/
    /*std::vector<int64_t> expected3 = {6, 19, 40, 70, 110, 114, 106, 85, 50};*/
    /*EXPECT_EQ(multiplyPolynomials(e, f, 2), expected3);*/
}
