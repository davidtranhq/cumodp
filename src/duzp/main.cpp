#include "duzp.hpp"
#include <iostream>

int main()
{
    int n;
    std::cin >> n;
    std::vector<Coefficient> a(n), b(n);
    for (auto& x : a)
        std::cin >> x;
    for (auto& x : b)
        std::cin >> x;
    std::vector<Coefficient> c = multiplyPolynomials(a, b);
    for (auto x : c)
        std::cout << x << ' ';
    return 0;
}
