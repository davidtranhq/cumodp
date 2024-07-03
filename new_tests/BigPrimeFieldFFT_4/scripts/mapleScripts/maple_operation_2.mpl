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
r:=2^62+2^36;
k:=16;
p:=r^k+1:
ULMAX:= 2^64-1;

vecSize:=1024:
#read data files for x and y
#stored in form of vectors of k=8 coefficients
xFile:=fopen("../../test/tmp/xData",READ,TEXT):
yFile:=fopen("../../test/tmp/yData",READ,TEXT):
writeFile:=fopen("../../test/tmp/uData_maple_2",WRITE):

#convert array[1..8] to n, bigPrimeNumber element
toBigInt:=proc(v,r,p)
  n:=0:
  for i from 1 to k do
    n:=n+v[i]*r^(i-1):
  end do:
  return n:
  end proc:

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

#function for converting big int n to array [1..8]
toBigFieldElement:=proc(n,r,p,k)
  i:=4; #=log[2](k)
  x:=toBigFieldElementRecursive(n,i,r,p):
  x:=convert(x,vector):
  #special case when x[k]==r
  # if x[k]=r then
    # x[k]:=0;
    # for i from 1 to 7 do
    # x[1]:=1;
    # end do;
  # end if;
  # if x[k]=r then
  #   x[k]:=0;
  #   for i from 1 to 7 do
  #     x[i]:=1;
  #   end do;
  # end if;
  return x:
end proc;

#computing bigSubtraction using Maple arithmetic
bigSubMaple:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs,r,p):
  yN:=toBigInt(ys,r,p):

  #addition modulo big prime "p"
  result:=(xN-yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,r,p,k):
  #write back result of addition to "uData_add_maple"
  resultArray:=convert(resultVector,array):
  writedata(writeFile,resultArray,integer):
  writeline(writeFile,""):
end proc:

xData:=readdata(xFile,integer):
yData:=readdata(yFile,integer):

for i from 1 to vecSize do
  xs:=convert(xData[(i-1)*k+1..i*k],Array):
  ys:=convert(yData[(i-1)*k+1..i*k],Array):
  bigSubMaple(xs,ys,r,p);
end do:

fclose(xFile):
fclose(yFile):
fclose(writeFile):