#include "duzp.hpp"
#include <cassert>
#include <iostream>

std::vector<Coefficient> naivelyMultiplyPolynomials(std::span<Coefficient> a, std::span<Coefficient> b)
{
    std::vector<Coefficient> result(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j)
            result[i + j] += a[i] * b[j];
    }
    return result;
}

int largestBitSize(std::span<const Coefficient> a)
{
    int N = 0;
    for (auto x : a)
        N = std::max(N, static_cast<int>(sizeof(Coefficient) * 8 - __builtin_clzll(x)));
    return N;
}

std::vector<Coefficient> convertUnivariateToBivariate(std::span<Coefficient> a, int K, int M)
{
    std::vector<Coefficient> bivariate(a.size() * K);
    for (size_t i = 0; i < a.size(); ++i) {
        Coefficient ai = a[i];
        for (int j = 0; j < K; ++j) {
            // Assign the coefficient a_{i, j} as the j-th M-bit chunk of a_i.
            bivariate[i * K + j] = ai & ((1 << M) - 1);
ai >>= M;
        }
    }
    return bivariate;
}

int64_t computeRecoveryPrime(size_t d, int K, int M)
{
    // Check that 2^k >= 2d - 1, where k is such that p = q*2^k + 1.
    auto containsLargeEnoughTwoPower = [](int64_t p, int64_t d) -> bool {
        assert(p > 1);
        assert(d > 0);
        int64_t twoPower = 1;
        // Find the largest power of 2 that divides p - 1.
        // We also check that twoPower is not 0 to avoid overflow.
        while ((p - 1) % (twoPower << 1) == 0 && (twoPower << 1) != 0)
            twoPower <<= 1;
        return twoPower >= 2 * d - 1;
    };

    auto isPrime = [](int64_t p) -> bool {
        if (p < 2)
            return false;
for (int64_t divisor = 2; divisor < p; ++divisor) {
            if (p % divisor == 0)
                return false;
        }
        return true;
    };
    // See the paper for why p >= 4*d*K*(2^(2*M)) + 1.
    for (int64_t p = 4 * d * K * (1 << (2 * M)) + 1; p < 0xffff'ffff'ffff'ffff; ++p) {
        if (isPrime(p) &&
            (p - 1) % (2 * d - 1) == 0 &&
            (p - 1) % K == 0 &&
            containsLargeEnoughTwoPower(p, d)) {
            return p;
        }
    };
    assert(false && "Failed to find a valid prime p under 0xffff'ffff'ffff'ffff.");
}

int64_t binaryModularPower(int64_t base, int64_t exponent, int64_t p)
{
    int64_t result = 1;
    while (exponent) {
        if (exponent & 1)
            result = (result * base) % p;
        base = (base * base) % p;
        exponent >>= 1;
    }
    return result;
}

int64_t modularInverse(int64_t x, int64_t p)
{
    return binaryModularPower(x, p - 2, p);
}

std::vector<Coefficient> dft(std::span<const Coefficient> a, int64_t p, std::span<const Coefficient> rootPowers)
{
    const int numPowers = rootPowers.size();
    // evaluated[i] corresponds to evaluating a(rootPowers[i]) in Z_p.
    std::vector<Coefficient> evaluated(numPowers);
    for (int i = 0; i < numPowers; ++i) {
        evaluated[i] = 0;
        // We have a(rootPowers[i]) = sum_{j=0}^{a.size()-1} (a[j] * rootPowers[i]^j).
        // Since rootPowers[i] = rootPowers[0]^i, we have rootPowers[i]^j = rootPowers[0]^((i*j) % rootPowers.size()).
        for (int j = 0; j < a.size(); ++j)
            evaluated[i] = (evaluated[i] + a[j] * rootPowers[i * j % numPowers] % p) % p;
    }
    return evaluated;
}

std::vector<Coefficient> dft2d(std::span<const Coefficient> A, int K, int64_t p, std::span<const Coefficient> xPowers, std::span<const Coefficient> yPowers)
{
    // We use the row-column algorithm to compute the 2-D DFT.
    // First, we compute the 1-D DFT on the rows of the matrix A.
    // Then, we compute the 1-D DFT on the columns of the matrix A.
    // Note that this is possible due to the separability of the DFT over the 2 variables.
    
    // A (not completely sufficient) check that the representation of the bivariate polynomial is correct, that is, each term in y is represented as a polynomial of degree K-1 in x.
    assert(A.size() % K == 0);

    // The matrix representation of the bivariate polynomial A(x, y) is such that the coefficient a_{i, j} of x^j * y^i is at (i, j).
    const int numRows = A.size() / K; // also, the degree in y + 1
    const int numCols = K; // also, the degree in x + 1

    std::vector<Coefficient> evaluatedRows(xPowers.size() * numRows);
    // Perform a 1-D DFT on the rows of the matrix A.
    for (int row = 0; row < numRows; ++row) {
        auto evaluatedRow = dft(A.subspan(row * numCols, numCols), p, xPowers);
        std::copy(evaluatedRow.begin(), evaluatedRow.end(), evaluatedRows.begin() + row * xPowers.size());
    }

    // Perform a 1-D DFT on the cols of the result of the row-wise DFT above.
    std::vector<Coefficient> evaluated(xPowers.size() * yPowers.size());
    for (int col = 0; col < numCols; ++col) {
        std::vector<Coefficient> column(numRows);
        for (int row = 0; row < numRows; ++row)
            column[row] = evaluatedRows[row * numCols + col];
        auto transformedColumn = dft(column, p, yPowers);
        for (int row = 0; row < yPowers.size(); ++row)
            evaluated[row * numCols + col] = transformedColumn[row];
    }
    return evaluated;
}

Coefficient findPrimitiveRootOfUnity(int n, int64_t p)
{
    if (n == 1)
        return 1;
    // Find a primitive root of unity of order K in Z_p.
    // That is, find a number g such that g^K = 1 (mod p) and
    // g^k != 1 (mod p) for all k < K.
    for (int64_t g = 2; g < p; ++g) {
        bool isPrimitive = true;
        int64_t gk = g;
        for (int k = 1; k < n; ++k) {
            if (gk == 1) {
                isPrimitive = false;
                break;
            }
            gk = (gk * g) % p;
        }
        if (isPrimitive && gk == 1)
            return g;
    }
    return -1;
} 

std::vector<Coefficient> cyclicConvolution(std::span<const Coefficient> A, std::span<const Coefficient> B, int K, int64_t p, Coefficient xPrimitiveRoot, Coefficient yPrimitiveRoot)
{
    // The degree of the polynomials A(x, y) and B(x, y) in x and y should be the same.
    assert(A.size() == B.size());
    // The degree of the polynomials A(x, y) and B(x, y) in x should be K.
    assert(A.size() % K == 0);
    const int yDegree = A.size() / K;
    const int productYDegree = yDegree * 2 - 1;

    auto precomputePowers = [](int64_t primitiveRoot, int64_t p, int numPowers) -> std::vector<Coefficient> {
        std::vector<Coefficient> powers(numPowers);
        powers[0] = 1;
        powers[1] = primitiveRoot;
        for (int i = 2; i < numPowers; ++i)
            powers[i] = (powers[i - 1] * primitiveRoot) % p;
        return powers;
    };

    std::vector<Coefficient> xPowers = precomputePowers(xPrimitiveRoot, p, K);
    std::vector<Coefficient> yPowers = precomputePowers(yPrimitiveRoot, p, productYDegree);

    const std::vector<Coefficient> transformedA = dft2d(A, K, p, xPowers, yPowers);
    const std::vector<Coefficient> transformedB = dft2d(B, K, p, xPowers, yPowers);
    std::vector<Coefficient> transformedProduct(K * productYDegree);
    assert(transformedA.size() == transformedB.size());
    assert(transformedA.size() == transformedProduct.size());

    for (int i = 0; i < transformedA.size(); ++i)
        transformedProduct[i] = transformedA[i] * transformedB[i] % p;

    // We use the fact that r^(-1) = r^(n - 1) (mod p) for any nth primitive root of unity in Z_p.
    std::vector<Coefficient> inverseXPowers = xPowers;
    std::reverse(inverseXPowers.begin() + 1, inverseXPowers.end());
    std::vector<Coefficient> inverseYPowers = yPowers;
    std::reverse(inverseYPowers.begin() + 1, inverseYPowers.end());

    std::vector<Coefficient> product = dft2d(transformedProduct, K, p, inverseXPowers, inverseYPowers);

    const int normalizationFactor = modularInverse(K * productYDegree, p);
    for (Coefficient& a : product)
        a = a * normalizationFactor % p;
    return product;
}

std::vector<Coefficient> negacyclicConvolution(std::span<const Coefficient> A, std::span<const Coefficient> B, int K, int64_t p, Coefficient xPrimitiveRoot, Coefficient yPrimitiveRoot)
{
    // Given the polynomial A(x, y) = sum_{i=0}^{d-1} (A(x)_i * y^i), return A'(x, y) = A(theta*x, y).
    auto transformByTheta = [](std::span<const Coefficient> A, int K, int64_t p, Coefficient theta) -> std::vector<Coefficient> {
        std::vector<Coefficient> transformed(A.size());
        for (int i = 0; i < A.size(); ++i)
            transformed[i] = (i % K ? A[i] * theta % p : A[i]);
        return transformed;
    };

    // TODO: theta is just the square root of xPrimitiveRoot, we should be able to compute it more efficiently?
    const Coefficient theta = findPrimitiveRootOfUnity(2 * K, p);
    const std::vector<Coefficient> APrime = transformByTheta(A, K, p, theta);
    const std::vector<Coefficient> BPrime = transformByTheta(B, K, p, theta);
    const std::vector<Coefficient> CPrime = cyclicConvolution(APrime, BPrime, K, p, xPrimitiveRoot, yPrimitiveRoot);
    const Coefficient thetaInverse = modularInverse(theta, p);
    const std::vector<Coefficient> CPlus = transformByTheta(CPrime, K, p, thetaInverse);

    return CPlus;
}

std::vector<Coefficient> evaluateAtBeta(std::span<const Coefficient> C, int K, int64_t p, Coefficient Beta)
{
    const int yDegree = C.size() / K;
    std::vector<Coefficient> c(yDegree);
    for (int i = 0; i < yDegree; ++i) {
        c[i] = 0;
        for (int j = 0; j < K; ++j)
            c[i] = (c[i] + C[i * K + j] * binaryModularPower(Beta, j, p) % p) % p;
    }
    return c;
}

std::vector<Coefficient> multiplyPolynomials(std::span<Coefficient> a, std::span<Coefficient> b, int M)
{
    assert(a.size() == b.size());
    const int N = std::max(largestBitSize(a), largestBitSize(b));
    assert(N <= 64);

    const int K = N / M + (N % M != 0);

    const std::vector<Coefficient> A = convertUnivariateToBivariate(a, K, M);
    const std::vector<Coefficient> B = convertUnivariateToBivariate(b, K, M);
    assert(A.size() == B.size());

    const int64_t p = computeRecoveryPrime(a.size(), K, M);

    const Coefficient xPrimitiveRoot = findPrimitiveRootOfUnity(K, p);
    assert(xPrimitiveRoot != -1);
    const Coefficient yPrimitiveRoot = findPrimitiveRootOfUnity(a.size() + b.size() - 1, p);
    assert(yPrimitiveRoot != -1);

    const std::vector<Coefficient> CMinus = cyclicConvolution(A, B, K, p, xPrimitiveRoot, yPrimitiveRoot);
    const std::vector<Coefficient> CPlus = negacyclicConvolution(A, B, K, p, xPrimitiveRoot, yPrimitiveRoot);

    const std::vector<Coefficient> u = evaluateAtBeta(CPlus, K, p, 1 << M);
    const std::vector<Coefficient> v = evaluateAtBeta(CMinus, K, p, 1 << M);
    std::vector<Coefficient> c(a.size() + b.size() - 1);

    assert(u.size() == v.size());
    assert(c.size() == u.size());

    for (size_t i = 0; i < c.size(); ++i)
        c[i] = (u[i] + v[i]) / 2 + (((v[i] - u[i]) / 2) << N);

    return c;
}
