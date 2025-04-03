import random
import json

def mod_inverse(a, m):
    """
    Compute the modular inverse of a modulo m using the extended Euclidean algorithm.
    Raises an exception if the inverse does not exist.
    """
    def egcd(a, b):
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = egcd(b % a, a)
            return (g, x - (b // a) * y, y)
    
    g, x, _ = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist for a = {} and m = {}'.format(a, m))
    return x % m

def generate_test_case(degree_x, degree_y, prime1, prime2, coeff_range=(-20, 20)):
    """
    Generates a test case for the CRT reconstruction.

    The bivariate polynomial is assumed to be stored in row-major order:
      coefficient of x^i y^j is at index i * baseK + j,
    where baseK = degree_y + 1.

    Returns a dictionary containing:
      - baseK: number of y coefficients per x degree,
      - prime1, prime2: the two primes,
      - original: the original integer coefficients,
      - product1: coefficients reduced modulo prime1,
      - product2: coefficients reduced modulo prime2,
      - expected: the expected reconstruction computed via CRT.
    """
    baseK = degree_y + 1
    original = []
    # Generate coefficients for polynomial with degrees 0..degree_x in x and 0..degree_y in y.
    for i in range(degree_x + 1):
        for j in range(degree_y + 1):
            coeff = random.randint(coeff_range[0], coeff_range[1])
            original.append(coeff)
    
    # Compute the modular representations.
    product1 = [x % prime1 for x in original]
    product2 = [x % prime2 for x in original]
    
    # Precompute the modular inverse of prime1 modulo prime2.
    inv_prime1 = mod_inverse(prime1, prime2)
    
    # Compute the expected result using the CRT formula.
    expected = []
    for a, b in zip(product1, product2):
        # Compute the difference modulo prime2.
        diff = (b - a) % prime2
        # Multiply by the modular inverse and reduce modulo prime2.
        factor = (diff * inv_prime1) % prime2
        # Reconstruct the full coefficient.
        expected_val = a + factor * prime1
        expected.append(expected_val)
    
    return {
        "baseK": baseK,
        "prime1": prime1,
        "prime2": prime2,
        "original": original,
        "product1": product1,
        "product2": product2,
        "expected": expected
    }

def main():
    # Set parameters for the test case.
    degree_x = 3  # Maximum degree in x (polynomial has degree 0..3 in x)
    degree_y = 3  # Maximum degree in y (polynomial has degree 0..3 in y)
    prime1 = 7
    prime2 = 11
    
    # Generate the test case.
    test_case = generate_test_case(degree_x, degree_y, prime1, prime2)
    
    # Output the test case in JSON format.
    print(json.dumps(test_case, indent=2))

if __name__ == "__main__":
    main()