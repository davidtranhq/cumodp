
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

M:=Primpart(M, x) mod p:

fdM := fopen("PolyM.dat", WRITE):
printPolyFile(M, x, m-1  , fdM):
fclose(fdM):

fdN := fopen("PolyN.dat", WRITE):
printPolyFile(N, x, n-1  , fdN):
fclose(fdN):



fdC := fopen("PolyQuo.dat", WRITE):
C:
C := Quo(N, M, x) mod p:
printPolyFile(C, x, (n  - m )  , fdC):
fclose(fdC):

fdRem := fopen("PolyRem.dat", WRITE):
E:
E := Rem(N, M, x) mod p:
printPolyFile(E, x, (m  - 2 )  , fdRem):
fclose(fdRem):

