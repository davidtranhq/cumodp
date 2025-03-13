#include <cstdint>
#include <span>
#include <vector>

using Coefficient = int64_t;

std::vector<Coefficient> naivelyMultiplyPolynomials(std::span<Coefficient> a, std::span<Coefficient> b);

// Return the largest bit size of any of the elements in a.
int largestBitSize(std::span<const Coefficient> a);

// Convert a univariate polynomial of y with coefficients a_i of degree d: 
// a(y) = sum_{i=0}^{d-1} (a_i * y^i)
// to a polynomial of degree d in y with coefficients A_i(x) of degree K in x:
// A(x, y) = sum_{i=0}^{d-1} (A_i * y^i)
// where A_i = sum_{j=0}^{K-1} (a_{i,j} * x^j).
//
// In other words, we represent the coefficients a_i of the univariate polynomial
// as themselves polynomials of degree K-1 in x, whose coefficients a_{i,j} are
// the representation of a_i in base 2^M. That is, A(2^M, y) = a(y).
//
// The elements in the resultant vector representing the bivariate polynomial are
// such that the coefficient a_{i, j} of x^j * y^i is at index i * K + j.
//
// In essence, each coefficient a_i of the univariate polynomial is represented as
// K blocks of M bits each, where the j-th block is the j-th least-significant M-bit chunk of a_i.
std::vector<Coefficient> convertUnivariateToBivariate(std::span<Coefficient> a, int K, int M);

// Return a prime p such that all the coefficients of C^+(x, y) and C^-(x, y)
// may be uniquely identified in the field Z_p, and such that
// Z_p supports the FFT (admits primitive roots of unity for computing the polynomials
// C^+(x, y) and C^-(x, y), that is, both 2d - 1 and K divide p - 1, and 2^k >= 2d - 1, where
// k is such that p = q*2^k + 1.
int64_t computeRecoveryPrime(size_t d, int K, int M);

// Compute the modular power of a number in Z_p.
int64_t binaryModularPower(int64_t base, int64_t exponent, int64_t p);

// Compute the modular inverse of x in Z_p, using
// Fermat's little theorem which states K^(p - 2) = K^{-1} (mod p) for all K in Z_p.
int64_t modularInverse(int64_t x, int64_t p);

// Perform a DFT on the polynomial a in Z_p at the points rootPowers.
// That is, return a vector `DFT` of size rootPowers.size(), where 
// `DFT[i] = a(rootPowers[i]) = sum_{j=0}^{a.size()-1} (a[j] * rootPowers[i]^j)`.
std::vector<Coefficient> dft(std::span<const Coefficient> a, int64_t p, std::span<const Coefficient> rootPowers);

std::vector<Coefficient> fastDft(std::span<const Coefficient> a, int64_t p, std::span<const Coefficient> rootPowers);

// Perform a 2-D DFT on the polynomial a in Z_p at the points (xPowers[i], yPowers[j]).
// That is, return a vector `DFT` of size xPowers.size() * yPowers.size(), where
// `DFT[i * xPowers.size() + j] = a(xPowers[i], yPowers[j]) = sum_{i=0}^{a.size() / K} sum_{j=0}^{K} a[i * K + j] * xPowers[0]^j * yPowers[0]^i`.
//
// The polynomial a is expected to be represented as a bivariate polynomial of degree K in x and y
// as described in the function `convertUnivariateToBivariate`.
//
std::vector<Coefficient> dft2d(std::span<const Coefficient> A, int K, int64_t p, std::span<const Coefficient> xPowers, std::span<const Coefficient> yPowers);

// Find an nth primitive root of unity in Z_p. If no such root exists, return -1.
Coefficient findPrimitiveRootOfUnity(int n, int64_t p);

// Compute the convolution of two polynomials A(x, y) and B(x, y) in Z_p.
// That is, compute A(x,y) * B(x,y) mod (x^k - 1, p).
// The polynomials of A and B are of degree K in x and should have the same degree in y.
// 
// xPrimitiveRoot is a primitive 2K-th root of unity in Z_p.
// yPrimitiveRoot is a primitive 2d-th root of unity in Z_p.
std::vector<Coefficient> cyclicConvolution(std::span<const Coefficient> A, std::span<const Coefficient> B, int K, int64_t p, Coefficient xPrimitiveRoot, Coefficient yPrimitiveRoot);

// Given a bivariate polynomial C(x, y) of degree K in x, return a univariate polynomial c(y) = C(Beta, y) in Z_p.
std::vector<Coefficient> evaluateAtBeta(std::span<const Coefficient> C, int K, int64_t p, Coefficient Beta);

// Compute the negacyclic convolution of two polynomials A(x, y) and B(x, y) in Z_p.
// That is, compute A(x,y) * B(x,y) mod (x^k + 1, p).
// theta is a primitive 2K-th root of unity in Z_p.
std::vector<Coefficient> negacyclicConvolution(std::span<const Coefficient> A, std::span<const Coefficient> B, int K, int64_t p, Coefficient xPrimitiveRoot, Coefficient yPrimitiveRoot);

// Multiply two polynomials using the two-convolution method.
std::vector<Coefficient> multiplyPolynomials(std::span<Coefficient> a, std::span<Coefficient> b, int M);
