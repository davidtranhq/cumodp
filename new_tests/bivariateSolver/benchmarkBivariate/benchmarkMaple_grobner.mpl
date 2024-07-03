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


with(RegularChains): 
with(FastArithmeticTools):
with(ChainTools):
with(PolynomialIdeals):
with(Groebner): 
with(FileTools):

writeFile:=fopen("results_maple_grobner.txt",WRITE);

if Exists("bivariateSystem.mpl")=true then
	read("bivariateSystem.mpl");
else
	p := 469762049: 
	nElements:=5:
	stepSize:=15:
	a_0:=40:
	b_0:=20:
end if;

vars := [x, y]: 
R := PolynomialRing(vars, p):
a:=Array(1..nElements, i->(stepSize*i)+a_0):
b:=Array(1..nElements, i->(stepSize*i)+b_0):

fprintf(writeFile, "Bivariate system: [x^a + y - 1, x + y^b - 1]\n"):
fprintf(writeFile, "--------------------------------------------\n"):
fprintf(writeFile, "a \t b \t elapsed \t algorithm\n"):
fprintf(writeFile, "--------------------------------------------\n"):

for i from 1 to nElements do
    f1 := x^a[i] + y - 1:
    f2 := x + y^b[i] - 1:

    J := < f1, f2, characteristic = p >:
    startTime := time():
    Basis(J, plex(x,y)):
    elapsedTime:=time() - startTime;
    fprintf(writeFile,"%d \t %d \t %f \t (Groebner Basis) \n" ,a[i], b[i], elapsedTime):
    print("Groebner", a[i], b[i], elapsedTime):
end do:

fprintf(writeFile, "\n--------------------------------------------\n\n"):