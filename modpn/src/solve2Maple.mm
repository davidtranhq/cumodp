printPolyFile := proc(A, x, y, pd1, pd2, fd)
	local j:
	local i:
	for j from 0 by 1 to pd2 do 	
		for i from 0 by 1 to pd1 do	
			fprintf(fd, "%d ", coeff(coeff(A, x, j), y, i)  ):
		od:
	od:
end proc:




fd := fopen("PDs.dat", READ):
X := readdata(fd, [integer, integer, integer, integer, integer]):
fclose(fd):


p := X[1][1]:
pd11 := X[1][2]:
pd12 := X[1][3]:
pd21 := X[1][4]:
pd22 := X[1][5]:

A := randpoly([x], coeffs = proc() randpoly([y], dense, degree = X[1][2]) end proc, degree = X[1][3]):
#A := (x+y)^pd11- 10:
A := A mod X[1][1]:

fd := fopen("Poly1.dat", WRITE):
printPolyFile(A, x, y, pd11, pd12, fd):

fclose(fd):

B := randpoly([x], coeffs = proc() randpoly([y], dense, degree = X[1][4]) end proc, degree = X[1][5]):
#B := (x+y)^pd21- 20:
B := B mod X[1][1]:

fd := fopen("Poly2.dat", WRITE):
printPolyFile(B, x, y, pd21, pd22, fd):

fclose(fd):

A;
B;
#save A, B, "mapleSession.txt":


#R := RegularChains:-PolynomialRing([x,y], p);
#RegularChains:-Triangularize([A, B], R);
#cs2 := RegularChains:-ConstructibleSetTools:-GeneralConstruct([A,B], [], R);
#RegularChains:-Display(cs2, R);


