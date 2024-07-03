####################################################################################
## Implementation of Newton Iteration for comuting the modular inverse of a polynomial.
## (From Evaluation, written by Sardar Haque)
#####################################
newtonIteration := proc(f, i, p)
	local g, j:
	g := 1:
	for j from 1 by 1 to i do 		
		g := Rem( (2*g - f*g*g) , x^(2^j), x) mod p:
	od:
	return g:
end proc:

####################################################################################
## Prints polynomial's coefficients in to a file.
## (From Evaluation, written by Sardar Haque)
#####################################
printPolyFile := proc(p, x, d, fd)
	local j:
	for j from 0 by 1 to d do 		
		fprintf(fd, "%d ", coeff(p, x, j)  ):
	od:
end proc:

####################################################################################
## This procedure evaluates polynomial with maple the function.
## (From Evaluation, written by Sardar Haque)
#####################################
simpleEvaluate := proc(k,p,fd)
	local j:
	global Value:
	
	for j from 1 by 1 to 2^(k-1) do 
		Value[j] :=  Eval(F[1][1], x= Points[j] )mod p:
		fprintf(fd, "%d ", Value[j]):
	od:
end proc:

####################################################################################
## This procedure construct the subproduct-tree & subInverse-tree for polynomial evaluation.
## (From Evaluation, written by Sardar Haque)
#####################################
buildSubproductTree := proc(k,p,fdM,fdMinv)
	local j,l:
	global M,Inv:

	for j from 1 by 1 to 2^(k-1) do 
		M[1][j] := x - Points[j]:
		Inv[1][j] :=  1:
		printPolyFile(M[1][j], x, 1  , fdM):
		fprintf(fdM, "\n"):
	od:

	for j from 2 by 1 to (k-1) do
		for l from 1 by 1 to 2^(k-j) do  
			M[j][l] := Expand( M[j-1][2*l-1]* M[j-1][2*l] ) mod p:
			Inv[j][l] := eval(M[j][l], x = 1/x):
	
			Inv[j][l] := expand(Inv[j][l]*x^(2^(j-1))):
			Inv[j][l] := newtonIteration(Inv[j][l], j-1 ,p):
			
			printPolyFile(M[j][l], x, 2^(j-1), fdM):
			if (j-1) >= plainMulLimit then 
				printPolyFile(Inv[j][l], x, 2^(j-1)-1  , fdMinv):
				fprintf(fdMinv, "\n"):
			end if:
	
			fprintf(fdM, "\n"): 
	
		od:
	od:
end proc:

####################################################################################
## This procedure evaluates polynomial.
## (From Evaluation, written by Sardar Haque)
#####################################
evaluate := proc(k,p,fd)
	local j,l:
	global F, Points:
	printPolyFile(F[1][1], x, 2^(k-1)-1  , fd):
	for j from 2 by 1 to k do
		for l from 1 by 1 to 2^(j-1) do
			F[j][l] := Rem(F[j-1][ceil(l*0.5)], M[k-j+1][l], x) mod p:
			fprintf(fd, "\n"):
			printPolyFile(F[j][l], x, 2^(k-j)-1  , fd):
		od:
	od:
end proc:

####################################################################################
## This procedure compute the modular interpolation of the polynomial over a finite field. For doing this,
## we construct a interpolation-tree down-to-top in which the final result will be the root of this tree.
## This tree contains polynomials with 1 degree less than the corresponding member in the sobproduct-tree.  
## first we initialize the leaves with C_{i}s; at every step we multiply every polynomial with adjacent polynomial
## of the correspoding one in the subproduct-tree and sum the result with the other polynomial in order to 
## generate the polynomial at the upper layer.
#####################################
fastInterpolation := proc(k,p,fdI)
	local j,l, tmp, fMaple:
	global INTER:
	
	for j from 1 by 1 to 2^(k-1) do
		INTER[1][j] := c[j]:
	od:

	for j from 2 by 1 to (k-1) do
		for l from 1 by 1 to 2^(k-j) do  
			INTER[j][l] := (Expand( M[j-1][2*l-1]* INTER[j-1][2*l])+Expand(M[j-1][2*l]* INTER[j-1][2*l-1])) mod p:
			printPolyFile(INTER[j][l], x, 2^(j-1)-1, fdI):
			fprintf(fdI, "\n"):
		od:
	od:
	INTER[k][1] := (Expand( M[k-1][1]* INTER[k-1][2])+Expand(M[k-1][2]* INTER[k-1][1])) mod p:
	fMaple := fopen("finalMaple.dat", WRITE):
	printPolyFile(INTER[k][1], x, 2^(k-1)-1, fMaple):
	fclose(fMaple):
end proc:

####################################################################################
## This procedure generates c_{i} = v_{i}*s_{i} in which v_{i} is the evaluated value of the polynomial for u_{i} (i-th point); 
## And s_{i}=1/(u_{i}-u_{0})...(u_{i}-u_{i-1})(u_{i}-u_{i+1})...(u_{i}-u_{n}). 
## For computing s_{i} we compute the derivation of (x-u_{0})(x-u_{1})...(x-u_{n}) and evaluate it on point u_{i}.
#####################################
generateInterpolationCoefficients := proc(k,p,fdC)
	global c:
	local s, sInv, j, mPrime, MM:
	
	MM := Expand(M[k-1][2]*M[k-1][1]) mod p:

        mPrime := diff(MM, x) mod p:
	
	for j from 1 by 1 to 2^(k-1) do
		sInv := Eval(mPrime,x= Points[j])mod p:
		s[j] := sInv^(-1) mod p:
		c[j] := s[j]*Value[j] mod p:
		fprintf(fdC, "%d ", c[j]):
	od:

end proc:

################################################################################################################################
################################################################################################################################
############################################ Here's the start of the program ###################################################

##########################################################
## Initialization...
##########################################################

##############
#### Read value of "k" & "p" from a file.
fd := fopen("KP.dat", READ):
X := readdata(fd, [integer, integer]):
fclose(fd):
w := X[1][1]:
p := X[1][2]:

k := w+1:
plainMulLimit := 8:

##############
#### Generate different random points and write them to a file.

points:={}:
while(nops(points)<2^(k-1) )do
        points:=points union {modp(rand(),p)}:
end do:

Points:=[op(points)]:
fd := fopen("Points.dat", WRITE):
for j from 1 by 1 to 2^(k-1) do
        fprintf(fd, "%d ", Points[j]):
od:
fclose(fd):

##############
#### Generate random values for polynomial's coefficients.
F[1][1] :=  RandomTools[Generate](polynom( integer(range = -0..(p-1))  , x, degree =  2^(k-1)-1)):

##########################################################
## First we evaluate polynomial.
##########################################################

##############
#### Evaluate F on Points and store the results to a file with maple fuction.
fd := fopen("value.dat", WRITE):
simpleEvaluate(k,p,fd):
fclose(fd):

##############
#### Build subproduct & subInverse trees and write them in files.
fdM := fopen("PolyM.dat", WRITE):
fdMinv := fopen("PolyMinv.dat", WRITE):
buildSubproductTree(k,p,fdM,fdMinv):
fclose(fdM):
fclose(fdMinv):

##############
#### Evaluate F on Points and write them in files.
fd := fopen("PolyF.dat", WRITE):
evaluate(k,p,fd):
fclose(fd):

##############
#### Check if the results from the algorithm & maple are the same or not.
for j from 1 by 1 to 2^(k-1) do 
	if Value[j] <> F[k][j] then print(false, j)
	end if;
od:

##########################################################
## Here we start polynomial interpolation.
##########################################################

##############
## Compute c_{i}=v_{i}*s_{i} which are leaves of the interpolation-tree.
fdC := fopen("c.dat", WRITE):
generateInterpolationCoefficients(k,p,fdC):
fclose(fdC):

##############
## Compute polynomial interpolation by constructins interpolation-tree.
fdI := fopen("PolyI.dat", WRITE):
fastInterpolation(k,p,fdI):
fclose(fdI):

##############
#### Check if the results of the interpolation & F are equal or not.
#for j from 1 by 1 to 2^(k-1) do 
#	if coeff(INTER[k][1], x, j) <> coeff(F[1][1], x, j) then print(false, j, coeff(F[1][1], x, j), coeff(INTER[k][1], x, j))
#	end if:
#od:
