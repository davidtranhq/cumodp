#include <vector>

template<typename CoefficientType>
class BivariatePolynomial
{
public:
    explicit BivariatePolynomial(std::vector<std::vector<CoefficientType>> coefficients);
private:
    std::vector<std::vector<CoefficientType>> m_coefficients;
};

