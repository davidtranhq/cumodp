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


########################################
# maple code for revised multiplication 
# of two elements of big prime field
#########################################
restart;
with(ArrayTools):
with(Ordinals):
kernelopts(printbytes=false);

#r: radix of big prime 
r:=2^63+2^34;
#k: power of radix in big prime
k:=8;
#p: big prime in z/pz
p:=r^k+1:
#vecSize: size of input vector
vecSize:=16:

#read data files for x and y
#stored in form of vectors of k=8 coefficients
# make sure you have executed "make generateData" in "BigPrimeFieldFFT_3"
xFile:=fopen("datasets/tmp/xData",READ,TEXT):
yFile:=fopen("datasets/tmp/yData",READ,TEXT):
#######################################
# v is a vector of k machine words, 
# representing one element of z/pz
# toBigInt converts v to an integer n 
toBigInt:=proc(v)
  n:=0:
  for i from 1 to k do
    n:=n+v[i]*r^(i-1):
  end do:
  return n:
  end proc:
#######################################
# n is an integer of size k machine-words
# toBigFieldElementRecursive converts n to a vector v
# of k machine-words, representing one element of z/pz

toBigFieldElementRecursive:=proc(n,i) 
  R:=r^(2^(i-1)):
  q,s:= Div(n,R);
  v:=Array([s,q]):
  if i>1 then
    v:= Concatenate(1,toBigFieldElementRecursive(s,i-1),toBigFieldElementRecursive(q,i-1)):
  end;
  return v:
end proc:  

#######################################
# n is an integer of size of k machine-words
# toBigFieldElement converts n to a vector x of k
# machine-words which will represent an element of z/pz

toBigFieldElement:=proc(n,k)
  local i;
  i:=3; #=log[2](k)
  x:=toBigFieldElementRecursive(n,i):
  x:=convert(x,vector):
  return x:
end proc;

#######################################
#xs and ys are elements of z/pz, each storing k digits 
# equalBigFieldElements is used for checking equality 
# of two elements of z/pz

equalBigFieldElements:=proc(xs,ys)
  local i;
  for i from 1 to k do
    if xs[i] <> ys[i] then
      print(" FAILED, inequal at index =",i);
      quit;
      return -1;
    fi;
  end do;
  return 1;
end proc;

#######################################
#xs and ys are elements of z/pz, each storing k digits 
#bigAddMaple computes addition in z/pz using Maple arithmetic

bigAddMaple:=proc(xs,ys)
  local xN, yN,result, resultVector;
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs):
  yN:=toBigInt(ys):

  #addition modulo big prime "p"
  result:=(xN+yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,k):
  return resultVector
end proc:

#######################################
# xs and ys are elements of z/pz, each storing k digits 
# bigAddMaple computes subtraction in z/pz using Maple arithmetic

bigSubMaple:=proc(xs,ys)
  local xN, yN,result, resultVector;
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs):
  yN:=toBigInt(ys):

  #addition modulo big prime "p"
  result:=(xN-yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,8):
  return resultVector
end proc:
#######################################
# xs is an element of z/pz, storing k digits 
# simpleRotateRight rotates digits of one element of z/pz 
# one unit to the right

simpleRotateRight:=proc(xs)
  local xShifted, i, t;
  xShifted:=Array(1..k,i->0);
  t:=xs[k];
  for i from k to 2 by -1 do
    xShifted[i]:=xs[i-1];
  end do;
  xShifted[1]:=t*(1);
  return xShifted;
end proc;
#######################################
# xs and ys are elements of z/pz, each storing k digits 
# bigMult0 computes multiplication in z/pz using Maple arithmetic

bigMult0:=proc(xs,ys)
  local xN, yN, result, resultVector;
  
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs):
  yN:=toBigInt(ys):

  #addition modulo big prime "p"
  result:=(xN*yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,k):
  if resultVector[k]=r then
    resultVector[k]:=0;
  	resultVector[1]:=r-1;
  fi;
  return resultVector
end proc:
#######################################
# s is sum of k products of the form x*y
# toLHC converts s to (l,h,c) such that
# s:=l+h*r+c*r^2

toLHC:=proc (s)
  local sN, l, h, c;
  sN:=s;
  c:=iquo(sN,r^2);
  sN:=irem(sN,r^2);
  h:=iquo(sN,r);
  sN:=irem(sN,r);
  l:=sN;
  return l,h,c;
end proc;
####################################
# xs is an element of z/pz, storing k digits 
# bigModulo converts xs to three vectors 
# lVector,hVector, and cVector such that
# xs = lVector + hVector * r + cVector*r^2;

bigModulo:=proc(xs)
  local lVector, hVector, cVector, res;
  local s,i;

	lVector:=Array(1..k,i->0):
	hVector:=Array(1..k,i->0):
	cVector:=Array(1..k,i->0):

  for i from 1 to k do
    s:=xs[i];
    res:=toLHC(s);
    cVector[i]:=res[3];;
    hVector[i]:=res[2];
    lVector[i]:=res[1];
  end do;
  return [lVector,hVector,cVector];
end proc;

#######################################
# xs and ys are elements of z/pz, each storing k digits 
# bigMult1 computes multiplication in z/pz using bigModulo 

bigMult1:=proc(xs,ys)
  local i,j,s;
  local ly, u, result, lVector, hVector, cVector;
  ly:=ys;
  result:=Array(1..k,i->0);
  u:=Array(1..k,i->0);
  
  for i from 1 to k do
  	if i>1 then
  		ly:=simpleRotateRight(ly);
      ly[1]:=(-1)*ly[1];
  	fi;

    s:=0;
    for j from 1 to k do 
      s:=s+xs[j]*ly[k+1-j];
    end do;
    u[k+1-i]:=s; 
  end do;
  ###################
  result:=bigModulo(u);# u-> l,h,c -> result:=l+hr+cr^2
  lVector:=result[1]:   
  hVector:=result[2]:  
  cVector:=result[3]:

  hVector:=simpleRotateRight(hVector);   
  cVector:=simpleRotateRight(cVector);    cVector:=simpleRotateRight(cVector);
  hVector[1]:=(-1)*hVector[1];
  cVector[1]:=(-1)*cVector[1];    cVector[2]:=(-1)*cVector[2];

  #result:=l+h*r
  result:=bigAddMaple(lVector,hVector);
  #result:=result+c*r^2
  result:=bigAddMaple(result,cVector);
  return result;
end proc:
#######################################
#reading data from file
xData:=readdata(xFile,integer):
yData:=readdata(yFile,integer):
###################
nSuccess:=0:
for i from 1 to vecSize do
  xs:=xData[(i-1)*k+1..i*k]:
  ys:=yData[(i-1)*k+1..i*k]:
  xs:=xs[1..k]:
  ys:=ys[1..k]:
  m0:=bigMult0(xs,ys);
  m1:=bigMult1(xs,ys);
  # print("idx",i);
  nSuccess:=nSuccess+ equalBigFieldElements(m0,m1);
  # print("Maple Mult:",m0):
  # print("revised Mult:",m1):
end do:
###################
if nSuccess = vecSize then
  print("Passed!");
  # `quit`(128):
else
  print("Failed!");
  # `quit`(255):
fi;
#######################################