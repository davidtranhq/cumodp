UnivariatePolynomial multiply(const UnivariatePolynomial& a, const UnivariatePolynomial& b)
{
    
    const auto aBivariate = a.toModularBivariate(prime1, prime2);
    const auto bBivariate = b.toModularBivariate(prime1, prime2);
    const auto outputDegree = a.degree() + b.degree();
    const auto exponentBits = logceiling(outputDegree);
    const auto xDegree = aBivariate.xDegree();
    if ((1 << exponentBits) == outputDegree) {
        // degree is an exact power of 2
        const auto cMinus = cyclicConvolution(aBivariate, bBivariate, prime1, xDegree);
        const auto cPlus = negacyclicConvolution(aBivariate, bBivariate, prime2, xDegree);
    }
    const auto c = recoverProduct(cMinus, cPlus);
    return c;
}
