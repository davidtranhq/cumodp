# Procedure to print the coefficients of a polynomial to a file.
printPolyFile := proc(p, x, d, fd)
    local j:
    for j from 0 by 1 to d do 
        fprintf(fd, "%d ", coeff(p, x, j)):
    od:
end proc:
#--------------------------------------------------
# Parameter lists.
# n and m are the numbers of coefficients for the two polynomials.
# p (here) is simply the upper limit for the coefficients; coefficients are chosen in the range 0..p-1.
nList := [1024]:
mList := [1024]:
pList := [99999]:  # Examples: p=13 will generate coefficients from 0 to 12.

# Loop over each combination of parameters.
for n in nList do
    for m in mList do
        for p in pList do
            # Generate random polynomials M and N with coefficients in 0 .. p-1.
            # M is of degree n-1 and N is of degree m-1.
            M := RandomTools[Generate](
                     polynom(integer(range = 0..(p-1)), x, degree = n-1)
                 ):
            N := RandomTools[Generate](
                     polynom(integer(range = 0..(p-1)), x, degree = m-1)
                 ):
                 
            # Multiply M and N and expand the result.
            C := M * N:
            
            # Build an output filename that encodes the parameters.
            # For example: "n256_m256_p13.dat"
            filename := cat("n", convert(n, string), "_m", convert(m, string),
                            "_p", convert(p, string), ".test"):
            fd := fopen(filename, WRITE):
            
            # Print the three polynomials on separate lines.
            printPolyFile(M, x, n-1, fd):
            fprintf(fd, "\n"):
            printPolyFile(N, x, m-1, fd):
            fprintf(fd, "\n"):
            printPolyFile(C, x, (n + m - 2), fd):
            
            fclose(fd):
        od:
    od:
od: