/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/


printPolyFile := proc(A, x, y, pd1, pd2, p, fd)
	local j:
	local i:
	for j from 0 by 1 to pd2 do 	
		for i from 0 by 1 to pd1 do
			#print(i,j, coeff(coeff(A, x, j), y, i) );	
			fprintf(fd, "%d ", (coeff(coeff(A, x, i), y, j)) mod p):
		od:
	od:
end proc:

fd5 := fopen("prime.txt", READ):
RR := readdata(fd5, [integer]);#, integer]):
fclose(fd5):

p := RR[1]:



fd := fopen("inputFile.txt", READ):
RRR := readdata(fd, [string, string]):
fclose(fd):

s1 := RRR[1][1]:
s2 := RRR[1][2]:

read s1:
read s2:

#vars[2] > vars1

#name should be vars, n_vars, eqs

#vars := ListTools:-Reverse(vars);
vars := ListTools:-Reverse([op(vars)]);

if n_vars = 2 then
eqs1 := eqs[1]:# mod p
pd11 := degree(eqs1,vars[1]):
pd12 := degree(eqs1,vars[2]):

eqs2 := eqs[2]:# mod p;
pd21 := degree(eqs2,vars[1]):
pd22 := degree(eqs2,vars[2]):



fd2 := fopen("outputFile.txt", WRITE):
fprintf(fd2,"1 "):
fprintf(fd2,"%d ",pd11):
fprintf(fd2,"%d ",pd12):
fprintf(fd2,"%d ",pd21):
fprintf(fd2,"%d ",pd22):
fclose(fd2):


fd3 := fopen("Poly1.dat", WRITE):
printPolyFile(eqs1, vars[1], vars[2], pd11, pd12,p,fd3):
fclose(fd3):

fd7 := fopen("Poly2.dat", WRITE):
printPolyFile(eqs2, vars[1], vars[2], pd21, pd22, p,fd7):
fclose(fd7):



else 
print("not a bivariate system.");
fd1 := fopen("outputFile.txt", WRITE):
fprintf(fd1,"-1"):
fclose(fd1):

end if;


