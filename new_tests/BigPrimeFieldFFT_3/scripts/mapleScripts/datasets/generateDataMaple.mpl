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


restart;
with(ArrayTools):
with(Ordinals):
kernelopts(printbytes=false):

r:=2^63+2^34;
p:=r^8+1;
k:=8;
#####################################
#convert array[1..8] to n, bigPrimeNumber element
toBigInt:=proc(v,r,p)
  n:=0:
  for i from 1 to k do
    n:=n+v[i]*r^(i-1):
  end do:
  return n:
  end proc:
#######################################
#recursive function for converting big int n to array [1..8]
toBigFieldElementRecursive:=proc(n,i,r,p) 
  R:=r^(2^(i-1)):
  q,s:= Div(n,R);
  v:=Array([s,q]):
  if i>1 then
    v:= Concatenate(1,toBigFieldElementRecursive(s,i-1,r,p),toBigFieldElementRecursive(q,i-1,r,p)):
  end;
  return v:
end proc:  

#######################################
#function for converting big int n to array [1..8]
toBigFieldElement:=proc(n,r,p,k)
  i:=3; #=log[2](k)
  # i:= log(k,2)[1];
  x:=toBigFieldElementRecursive(n,i,r,p):
  x:=convert(x,vector):
  x:=convert(x,Array):
  return x:
end proc;
#######################################

vecSize:=1024;

#xN=-1
system("cp data_0_11_22_1k/xData data_5_11_22_1k/xData");
writeFile:=fopen("data_5_11_22_1k/yData",WRITE,TEXT):
for i from 1 to vecSize do
	xN:=modp(-1,p);
	# print(modp(2^i*xN,p));
	xArray:=toBigFieldElement(xN,r,p,k);
	v:=convert(xArray,list):
	# print(v);
	writedata(writeFile,v,integer):
	writeline(writeFile,""):
end do:
fclose(writeFile);


#inverse of powers of 2 (1/2^n)
system("cp data_0_11_22_1k/xData data_6_11_22_1k/xData");
writeFile:=fopen("data_6_11_22_1k/yData",WRITE,TEXT):
for i from 1 to vecSize do
	xN:=modp(1/2^i,p);
	# print(modp(2^i*xN,p));
	xArray:=toBigFieldElement(xN,r,p,k);
	v:=convert(xArray,list):
	writedata(writeFile,v,integer):
	writeline(writeFile,""):
end do:
fclose(writeFile);


system("cp data_0_11_22_1k/xData data_7_11_22_1k/xData");
# system("cp $CUMODP_HOME/new_tests/BigPrimeFieldFFT_3/src/data_0_11_22_1k/xData data_6_11_22_1k/xData");


#xN=-1
xArray:=Array(1..k,i->r-1);
v:=convert(xArray,list):
writeFile:=fopen("data_8_11_22_1k/xData",WRITE,TEXT):
for i from 1 to vecSize do
	# print(v);
	writedata(writeFile,v,integer):
	writeline(writeFile,""):
end do:
fclose(writeFile);

system("cp data_8_11_22_1k/xData data_8_11_22_1k/yData");