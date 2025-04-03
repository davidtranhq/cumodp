# Two Convolution Polynomial Multiplication

In the below, we describe the implementation of the two convolution polynomial multiplication algorithm on the GPU.

## DetermineBase

**Input**: The univariate polynomials $a(y)$ and $b(y)$ represented as `std::vector<mpz_class>`'s, where `a[i]` is the coefficient of `y^i` in $a(y)$, and symmetrically for $b$.

**Output**: The `BivariateBase` for the bivariate expansion:
```cpp
BivariateBase {
    int N; // the largest number of bits required to represent any coefficient of a or b.

    int K; // the number of blocks by which to partition each coefficient. The partial degree with respect to x after bivariate expansion.

    int M; // the bit-width of each block. So, each coefficient of the bivariate expansion can be represented with M bits.
};
```

Note that $N = KM$. We consider $K$ and $M$ as functions of $N$, which can be found by lookup in a table. 