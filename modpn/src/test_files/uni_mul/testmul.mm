
printPolyFile := proc(p, x, d, fd)
	local j:
	for j from 0 by 1 to d do 		
		fprintf(fd, "%d ", coeff(p, x, j)  ):
	od:
end proc:


fd := fopen("MNP.dat", READ):
X := readdata(fd, [integer, integer, integer]):
fclose(fd):


n := X[1][1]:
m := X[1][2]:
p := X[1][3]:

M:
M :=  RandomTools[Generate](polynom( integer(range = -0..(p-1))  , x, degree =  (m-1))):
N:
N :=  RandomTools[Generate](polynom( integer(range = -0..(p-1))  , x, degree =  (n-1))):

C:
C := Expand( M*N ) mod p:

fdM := fopen("PolyM.dat", WRITE):
printPolyFile(M, x, m-1  , fdM):
fclose(fdM):

fdN := fopen("PolyN.dat", WRITE):
printPolyFile(N, x, n-1  , fdN):
fclose(fdN):

fdC := fopen("PolyC.dat", WRITE):
printPolyFile(C, x, (m  + n -2)  , fdC):
fclose(fdC):

