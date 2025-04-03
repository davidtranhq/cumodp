import random

# Function to generate a polynomial with random coefficients
def generate_polynomial(degree, coef_range):
    return [random.randint(*coef_range) for _ in range(degree + 1)]

# Function to multiply two polynomials
def multiply_polynomials(poly1, poly2):
    result_degree = len(poly1) + len(poly2) - 2
    result = [0] * (result_degree + 1)

    for i in range(len(poly1)):
        for j in range(len(poly2)):
            result[i + j] += poly1[i] * poly2[j]
    
    return result

def main():
    # User input for degree and coefficient range
    degree1 = int(input("Enter the degree of the first polynomial: "))
    degree2 = int(input("Enter the degree of the second polynomial: "))
    min_coef = int(input("Enter the minimum coefficient value: "))
    max_coef = int(input("Enter the maximum coefficient value: "))
    
    # Generate two random polynomials
    poly1 = generate_polynomial(degree1, (min_coef, max_coef))
    poly2 = generate_polynomial(degree2, (min_coef, max_coef))
    
    # Multiply the polynomials
    product = multiply_polynomials(poly1, poly2)
    
    # Print the coefficients as comma-separated lists
    print(f"{{{", ".join(map(str, poly1))}}}")
    print(f"{{{", ".join(map(str, poly2))}}}")
    print(f"{{{", ".join(map(str, product))}}}")

if __name__ == "__main__":
    main()
