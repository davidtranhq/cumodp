import numpy as np
import random

def generate_poly(degree_x, degree_y):
    """
    Generate a bivariate polynomial represented as a 2D numpy array.
    The array shape is (degree_y+1, degree_x+1) so that the coefficient
    of x^i * y^j is at poly[j, i]. Coefficients are random integers between -10000 and 10000.
    """
    width = degree_x + 1
    height = degree_y + 1
    total_coeffs = width * height
    coeffs = [random.randint(-10000, 10000) for _ in range(total_coeffs)]
    poly = np.array(coeffs, dtype=np.int64).reshape((height, width))
    return poly

def poly_mult_naive(poly1, poly2):
    """
    Multiply two bivariate polynomials represented as 2D arrays using the naive method.
    If poly1 has shape (h1, w1) and poly2 has shape (h2, w2), then the product has shape
    (h1+h2-1, w1+w2-1) with:
        product[j, i] = sum_{j1+j2=j, i1+i2=i} poly1[j1, i1]*poly2[j2, i2]
    """
    h1, w1 = poly1.shape
    h2, w2 = poly2.shape
    out_h = h1 + h2 - 1
    out_w = w1 + w2 - 1
    product = np.zeros((out_h, out_w), dtype=np.int64)
    
    for j1 in range(h1):
        for i1 in range(w1):
            coeff1 = poly1[j1, i1]
            # Only proceed if the coefficient is nonzero to save time
            if coeff1 == 0:
                continue
            for j2 in range(h2):
                for i2 in range(w2):
                    product[j1 + j2, i1 + i2] += coeff1 * poly2[j2, i2]
    return product

def flatten_poly(poly):
    """
    Flatten the 2D polynomial into a single line of space-separated integers.
    The flattening is done in row-major order, so that the coefficient of x^i * y^j
    appears at index j * (width) + i.
    """
    flat = poly.flatten()
    return " ".join(map(str, flat.tolist()))

def main():
    # Number of test cases to generate.
    num_test_cases = 5

    # Open the output file for writing.
    with open("bi_poly_mul_tests_cases.txt", "w") as f:
        for _ in range(num_test_cases):
            # Randomly choose degrees for the two input polynomials.
            dx1 = random.randint(0, 8)
            dy1 = random.randint(0, 8)
            dx2 = random.randint(0, 8)
            dy2 = random.randint(0, 8)

            # Generate the two polynomials.
            poly1 = generate_poly(dx1, dy1)
            poly2 = generate_poly(dx2, dy2)

            # Compute their product using the naive multiplication algorithm.
            product = poly_mult_naive(poly1, poly2)

            # Write the three lines (polynomial1, polynomial2, product) to the file.
            f.write(flatten_poly(poly1) + "\n")
            f.write(flatten_poly(poly2) + "\n")
            f.write(flatten_poly(product) + "\n")

if __name__ == "__main__":
    main()