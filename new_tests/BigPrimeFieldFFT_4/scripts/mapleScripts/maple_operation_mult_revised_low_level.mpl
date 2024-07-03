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


###############################
restart:
with(ArrayTools):
with(Ordinals):
r:=2^63+2^34:
k:=8:
p:=r^k+1:
ULMAX:= 2^64-1:
# size of input vector
vecSize:=1024:

#read data files for x and y
#stored in form of vectors of k=8 coefficients
xFile:=fopen("../../../test/data/data_4_11_22_1k/xData",READ,TEXT):
yFile:=fopen("../../../test/data/data_4_11_22_1k/yData",READ,TEXT):
writeFile:=fopen("../../../test/tmp/uData_maple_0",WRITE, TEXT):
#######################################
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
  q,s:= Div(n,R):
  v:=Array([s,q]):
  if i>1 then
    v:= Concatenate(1,toBigFieldElementRecursive(s,i-1,r,p),toBigFieldElementRecursive(q,i-1,r,p)):
  end:
  return v:
end proc:  

#######################################
#function for converting big int n to array [1..8]
toBigFieldElement:=proc(n,r,p,k)
  i:=3: #=log[2](k)
  # i:= log(k,2)[1]:
  x:=toBigFieldElementRecursive(n,i,r,p):
  x:=convert(x,vector):
  return x:
end proc:

#######################################
equalBigFieldElements:=proc(xs,ys)
for i from 1 to k do
  if xs[i] <> ys[i] then
  print("FAIL FAIL FAIL FAIL FAIL inequal at index =",i):
  quit:
  return:
  fi:
  end do:
  print("equal"):
end proc:

#######################################
# #computing bigAddition using Maple arithmetic
bigAddMaple:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs,r,p):
  yN:=toBigInt(ys,r,p):

  #addition modulo big prime "p"
  result:=(xN+yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,r,p,k):
  return resultVector
end proc:

#######################################
# #computing bigSubtraction using Maple arithmetic
bigSubMaple:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs,r,p):
  yN:=toBigInt(ys,r,p):

  #addition modulo big prime "p"
  result:=(xN-yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,r,p,8):
  return resultVector
end proc:

#######################################
# #computing bigSubtraction using Maple arithmetic
bigSubMaple3_v1:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  # xN:=xs[1]+xs[2]*r+xs[3]*r^2;
  # yN:=ys[1]+ys[2]*r+ys[3]*r^2;

resultVector:=Array(1..3,i->0);
s:=ys[1]:
if xs[1]>=s then
  s:=xs[1]-s:
  c:=0:
else:
  s:=r+xs[1]-s:
  c:=1:
fi:
if s>=r then
  s:=s-r:
  c:=1:
fi:
resultVector[1]:=s;

s:=ys[2]+c:
if xs[2]>=s then
  s:=xs[2]-s:
  c:=0:
else:
  s:=r+xs[2]-s:
  c:=1:
fi:
if s>=r then
  s:=s-r:
  c:=1:
fi:
resultVector[2]:=s;

s:=ys[3]+c:
if xs[3]>=s then
  s:=xs[3]-s:
  c:=0:
else:
  s:=r+xs[3]-s:
  c:=1:
fi:
if s>=r then
  s:=s-r:
  c:=1:
fi:
resultVector[3]:=s;

  # #addition modulo big prime "p"
  # result:=(xN-yN) mod p;
  # s2:=iquo(result,r^2);
  # result:=irem(result,r^2);
  # s1:=iquo(result,r);
  # s0:=irem(result,r);
  # #convert result of addition to vector of k=8 coefficients
  # # resultVector:=toBigFieldElement(result,r,p,8):
  # resultVector:=Array(1..3,i->0);
  # resultVector[1]:=s0;
  # resultVector[2]:=s1;
  # resultVector[3]:=s2;
  return resultVector
end proc:

#######################################
bigSubMaple3:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=xs[1]+xs[2]*r+xs[3]*r^2;
  yN:=ys[1]+ys[2]*r+ys[3]*r^2;
  #addition modulo big prime "p"
  result:=(xN-yN) mod p;
  # s2:=iquo(result,r^2);
  # result:=irem(result,r^2);
  # s1:=iquo(result,r);
  # s0:=irem(result,r);
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,r,p,8):
  # resultVector:=Array(1..3,i->0);
  # resultVector[1]:=s0;
  # resultVector[2]:=s1;
  # resultVector[3]:=s2;
  print(resultVector);
  return resultVector
end proc:

#######################################
#multiplication by a power of radix
multPowR:=proc(xs,s,r,p)
  if s<=0 or s>=k
  then
  print("in multPowR: s is either leq 0 or geq k=",k):
  return xs:
  fi:
  xN:=toBigInt(xs,r,p):
  xN:=(xN*r^s)mod (p):
  result:=toBigFieldElement(xN,r,p,8):
  if result[8]=r then
    result[8]:=0:
    result[1]:=r-1:
    print("special case in multPowR"):
  fi:
  return result:
end proc:
#######################################
simpleRightShift:=proc(xs,r,p)
  xShifted:=Array(1..k,i->0):
  t:=xs[k]:
  for i from k to 2 by -1 do
    xShifted[i]:=xs[i-1]:
  end do:
  xShifted[1]:=t*(1):
  return xShifted:
end proc:
#######################################
#computing bigMult using Maple arithmetic
bigMult0:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  xN:=toBigInt(xs,r,p):
  yN:=toBigInt(ys,r,p):

  #addition modulo big prime "p"
  result:=(xN*yN) mod p:
  #convert result of addition to vector of k=8 coefficients
  resultVector:=toBigFieldElement(result,r,p,k):
  #write back result of addition to "uData_add_maple"
  # resultArray:=convert(resultVector,array):
  if resultVector[k]=r then
	resultVector[k]:=0:
  	resultVector[1]:=r-1:
  	# print("special case in mult0"):
  	fi:
  return resultVector
end proc:
#######################################
# subtraction modulo r
subR:=proc(x,y,c0,r)
s:=y:
if x>=s then
  s:=x-s:
  c:=0:
else:
  s:=r+x-s:
  c:=1:
fi:
if s>=r then
  s:=s-r:
  c:=1:
fi:
if c=1 then c:=r-1: fi:
c:=c+c0:
return s,c:
end proc:
#######################################
# addition modulo r
addR:=proc(x,y,c0,r)
  s:=irem(x+y,2^64):
  c:=0:
  if s<x or s<y then
    s:=s+2^64-r:
    c:=1:
  fi:
  if s>=r then:
    # print("here"):
    s:=s-r:
    c:=1:
  fi:
  c:=c+c0:
  c:=irem(c,r):
  return s,c:
end proc:
#######################################
divR_v4:=proc(xs,r)
  r0:=2^34:
  s:=xs:
  q:=0:
  m:=0:
  #s can be 128 bit, heavy arithmetic?
  f:=1:
  # while s>=r do
  if s>=2^64 then print ("s is huge"): fi:
  if s>=r then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    q0:=iquo(s,r-r0):
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0)):
      f2:=-1:
    else
      s:=irem(s,(r-r0))-q0*r0:
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
  fi:
  if s>=r then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    q0:=iquo(s,r-r0):
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0)):
      f2:=-1:
    else
      s:=irem(s,(r-r0))-q0*r0:
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
  fi:

  if s>=r then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    q0:=iquo(s,r-r0):
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0)):
      f2:=-1:
    else
      s:=irem(s,(r-r0))-q0*r0:
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
  fi:

  # print(s,q,f):
  if f=-1 and s<>0 then
    # print("f = -1"):
    s:=r-s:
    q:=q-1:
  fi:
  m:=s:
  return q,m:
end proc:
#######################################
divR_v5:=proc(xs,r)
  r0:=2^34:
  s:=xs:
  q:=0:
  m:=0:
  f:=1:
  s1:=iquo(s,2^64): #iquo(s,r-r0)
  s0:=irem(s,2^64): #irem(s,r-r0)

  # if s>=r then
  if s1>=2^63 then print ("s is huge"): fi:
  t1:=2*s1:
  # t11:=iquo(t1,2^64):
  # t10:=irem(t1,2^64):
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  if q0>0 then
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)
  # print("s0,s1",s0,s1):
  t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    # t1:=2*s1:
    # t0:=iquo(s0,2^63):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)
 t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    # t1:=2*s1:
    # t0:=iquo(s0,2^63):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)

  if f=-1 and s0<>0 then
    # print("f = -1"):
    print("s0,s1",s0,s1):
    # s:=r-s:
    s0:=r-s0:
    q:=q-1:
  fi:
  m:=s0:
  return q,m:
end proc:
#######################################
divR_v6:=proc(xs,r)
  r0:=2^34:
  s:=xs:
  q:=0:
  m:=0:
  f:=1:
  s1:=iquo(s,2^64): #iquo(s,r-r0)
  s0:=irem(s,2^64): #irem(s,r-r0)

  # if s>=r then
  # if s1>=2^63 then print ("s is huge"): fi:
  t1:=2*s1:
  if t1>=2^64 then print ("t1 is huge"): quit: fi:
  # never happens, because xs<2^128 -> q<2^64

  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  if q0>=2^64 then print ("q0 is huge"): quit: fi:
  if q0>0 then
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    # if m0 < q0*r0 then
    if m0<s0 or s1>0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)
  # print("s0,s1",s0,s1):
  t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    # t1:=2*s1:
    # t0:=iquo(s0,2^63):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    # if m0 < q0*r0 then
    if m0<s0 or s1>0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)
 t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print(s0,s1):
    # t1:=2*s1:
    # t0:=iquo(s0,2^63):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    # if m0 < q0*r0 then
    if m0<s0 or s1>0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    s0:=sl[1]:
    s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  # s:=s0+s1*2^64:
  # s1:=iquo(s,2^64): #iquo(s,r-r0)
  # s0:=irem(s,2^64): #irem(s,r-r0)

  if f=-1 and s0<>0 then
    # print("f = -1"):
    print("s0,s1",s0,s1):
    # s:=r-s:
    s0:=r-s0:
    q:=q-1:
  fi:
  m:=s0:
  if q>=2^64 then print ("q is huge"): quit: fi:
  return q,m:
end proc:
#######################################
divR_v7:=proc(xs,r)
  r0:=2^34:
  q:=0:
  m:=0:
  f:=1:
  s:=xs:
  s1:=iquo(s,2^64): #iquo(s,r-r0)
  s0:=irem(s,2^64): #irem(s,r-r0)

  #round 1
  # t1:=2*s1:
  # t0:=iquo(s0,2^63):
  # q0:=t0+t1:
  # q:=q0:
  # if q0>0 then
  #   # s1:=iquo(s,2^64):
  #   # s0:=irem(s,2^64):
  #   m0:=irem(s0,2^63):
  #   # q0:=iquo(s,r-r0):
  #   # q0:=t1+t0:
  #   # q0:=iquo(s0+s1*2^64,r-r0):
  #   # q:=q+q0*f:
  #   print("m0,q,q0,f",m0,q,q0,f):
  #   # if irem(s,(r-r0))<q0*r0 then
  #   s1:=iquo(q0*r0,2^64):
  #   s0:=irem(q0*r0,2^64):
  #   if m0 < q0*r0 then
  #     print("m0<q0r0"):
  #     # if m0<s0 or s1>0 then
  #     # s:=q0*r0-irem(s,(r-r0)):
  #     # s:=q0*r0-m0:
  #     sl:=sub_128(s0,s1,m0,0):
  #     # s:=s0+s1*2^64:
  #     # s:=s-m0:
  #     # s1:=iquo(s,2^64):
  #     # s0:=irem(s,2^64):
  #     # print("sub-maple:",s0,s1):
  #     s0:=sl[1]:
  #     s1:=sl[2]:
  #     # print("sub-arithmetic:",s0,s1):

  #     print("========================"):
  #     f2:=-1:
  #     # q:=q-q0:
  #   else
  #   print("m0 >> q0r0"):
  #     # s:=irem(s,(r-r0))-q0*r0:
  #     # s:=m0-q0*r0:
  #     sl:=sub_128(m0,0,s0,s1):
  #     # s:=s0+s1*2^64:
  #     # s:=m0-s:
  #     # s1:=iquo(s,2^64):
  #     # s0:=irem(s,2^64):
  #     # print("sub-maple:",s0,s1):
  #     s0:=sl[1]:
  #     s1:=sl[2]:
  #     # print("sub-arithmetic:",s0,s1):
  #     print("#########################"):
  #     f2:=1:
  #     # q:=q+q0:
  #   fi:
  #   # if s<0 then f:=-1: else f:=1: fi:
  #   # s:=s*f:
  #   # f:=f*f2:
  #   f:=f2:
  #   # print("s",s):
  #   # s0:=sl[1]:
  #   # s1:=sl[2]:
  #   # print("s0,s1",s0,s1):
  # fi:
q:=0:
for it from 1 to 3 do
  t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  if q0>0 then
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:
    print("m0,q,q0,f",m0,q,q0,f):
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    if m0 < q0*r0 then
      print("m0<q0r0"):
    # if m0<s0 or s1>0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      # s:=s0+s1*2^64:
      # s:=s-m0:
      # s1:=iquo(s,2^64):
      # s0:=irem(s,2^64):
      # print("sub-maple:",s0,s1):
      s0:=sl[1]:
      s1:=sl[2]:
      # print("sub-arithmetic:",s0,s1):

      print("========================"):
      f2:=-1:
      # q:=q-q0:
    else
    print("m0 >> q0r0"):
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      # s:=s0+s1*2^64:
      # s:=m0-s:
      # s1:=iquo(s,2^64):
      # s0:=irem(s,2^64):
      # print("sub-maple:",s0,s1):
      s0:=sl[1]:
      s1:=sl[2]:
      # print("sub-arithmetic:",s0,s1):
      print("#########################"):
      f2:=1:
      # q:=q+q0:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    # s0:=sl[1]:
    # s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  end do:

  if f=-1 and s0<>0 then
    # print("s0,s1",s0,s1):
    # s:=r-s:
    s0:=r-s0:
    q:=q-1:
  fi:
  m:=s0:
  if q>=2^64 then print ("q is huge"): quit: fi:
  print("divR,q,m,",q,m):
  print("divR,maple",iquo(xs,r),irem(xs,r)):
  return q,m:
end proc:
#######################################
divR_v8:=proc(xs,r)
  r0:=2^34:
  s:=xs:
  q:=0:
  m:=0:
  #s can be 128 bit, heavy arithmetic?
  f:=1:
  # while s>=r do

  for it from 1 to 3 do
  if s>=r then
    # print(s):
    # q0 is a monomial, can be computed by shifting s 
    s1:=iquo(s,2^64):
    s0:=irem(s,2^64):
    # print(s0,s1):
    q0:=iquo(s,r-r0):
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:

    m0:=irem(s0,2^63):
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    if m0 <q0*r0 then
      # s:=q0*r0-irem(s,(r-r0)):  
      sl:=sub_128(s0,s1,m0,0):
      s:=sl[1]+sl[2]*2^64:
      f2:=-1:
    else
      # s:=irem(s,(r-r0))-q0*r0:
      print("m0,s0,s1",m0,s0,s1):
      sl:=sub_128(m0,0,s0,s1):
      s:=sl[1]+sl[2]*2^64:
      print("s",s):
      f2:=1:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
  fi:
  end do:

  # print(s,q,f):
  if f=-1 and s<>0 then
    # print("f = -1"):
    s:=r-s:
    q:=q-1:
  fi:
  m:=s:
  return q,m:
end proc:
#########################################
divR:=proc(xs,r)
  r0:=2^34:
  q:=0:
  m:=0:
  f:=1:
  s:=xs:
  s1:=iquo(s,2^64): #iquo(s,r-r0)
  s0:=irem(s,2^64): #irem(s,r-r0)

q:=0:
for it from 1 to 3 do
  t1:=2*s1:
  t0:=iquo(s0,2^63):
  q0:=t0+t1:
  if q0>0 then
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    m0:=irem(s0,2^63):
    # q0:=iquo(s,r-r0):
    # q0:=t1+t0:
    # q0:=iquo(s0+s1*2^64,r-r0):
    q:=q+q0*f:
    # print("m0,q,q0,f",m0,q,q0,f):
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64):
    s0:=irem(q0*r0,2^64):
    # print("s0,s1",s0,s1):
    # if m0 < q0*r0 then
      # print("m0<q0r0"):
    if m0<s0 or s1>0 then
      # s:=q0*r0-irem(s,(r-r0)):
      # s:=q0*r0-m0:
      sl:=sub_128(s0,s1,m0,0):
      # s:=s0+s1*2^64:
      # s:=s-m0:
      # s1:=iquo(s,2^64):
      # s0:=irem(s,2^64):
      # print("sub-maple:",s0,s1):
      s0:=sl[1]:
      s1:=sl[2]:
      # print("sub-arithmetic:",s0,s1):

      # print("========================"):
      f2:=-1:
      # q:=q-q0:
    else
    # print("m0 >> q0r0"):
      # s:=irem(s,(r-r0))-q0*r0:
      # s:=m0-q0*r0:
      sl:=sub_128(m0,0,s0,s1):
      # s:=s0+s1*2^64:
      # s:=m0-s:
      # s1:=iquo(s,2^64):
      # s0:=irem(s,2^64):
      # print("sub-maple:",s0,s1):
      s0:=sl[1]:
      s1:=sl[2]:
      # print("sub-arithmetic:",s0,s1):
      # print("#########################"):
      f2:=1:
      # q:=q+q0:
    fi:
    # if s<0 then f:=-1: else f:=1: fi:
    # s:=s*f:
    f:=f*f2:
    # print("s",s):
    # s0:=sl[1]:
    # s1:=sl[2]:
    # print("s0,s1",s0,s1):
  fi:
  end do:

  if f=-1 and s0<>0 then
    # print("s0,s1",s0,s1):
    # s:=r-s:
    s0:=r-s0:
    q:=q-1:
  fi:
  # print(s0,r):
  if s0>=r then
    s0:=s0-r:
    q:=q+1:
  fi:

  m:=s0:
  # if q>=2^64 then print ("q is huge"): quit: fi:
  # print("divR,q,m,",m,q):
  # print("divR,maple",irem(xs,r),iquo(xs,r)):
  return q,m:
end proc:
#######################################
# s (three machine-words) -> (l,h,c) modulo r
toLHC:=proc (s,r)
  s0:=s[1]:
  s1:=s[2]:
  s2:=s[3]:
  r0:=2^34:
  q:=iquo(2^64,r-r0):
  # p2:=s2*q^2:
  p2:=4*s2:
  # p1:=s1*q:
  # p1:=2*s1:
  if s1<2^63 then 
    p1:=2*s1:
  else
    p1:=2*s1:
    p2:=p2+iquo(p1,r):
    p1:=irem(p1,r):
  fi:

  p0:=s0:
  t0:=(2*s1*r0):
  t1:=(4*s2*r0^2):

  # print("t00,t01,t10,t11",irem(t0,2^64),iquo(t0,2^64),irem(t1,2^64),iquo(t1,2^64)):
  t2:=8*r0*s2:
  # if t0> 2^64 then print("t0 is huge"): fi:
  # if t1>2^64 then print("t1 is huge"): fi:
  t0q:=iquo(t0,r):
  t0r:=irem(t0,r):
  d:=divR(t0,r):
  t0q:=d[1]:
  t0r:=d[2]:

  t1q:=iquo(t1,r):
  t1r:=irem(t1,r):
  d:=divR(t1,r):
  t1q:=d[1]:
  t1r:=d[2]:

  # print("tvalues", t0q,t0r,t1q,t1r):
  # p1:=p1-iquo(2*s1*r0,r)+iquo(4*s2*r0^2,r)-8*r0*s2:
  # p0:=p0-rem(2*s1*r0,r)+rem(4*s2*r0^2,r):
  # p1:=p1-iquo(t0,r)+iquo(t1,r)-t2:
  # p0:=p0-rem(t0,r)+rem(t1,r):
  # p1:=p1-t0q+t1q-t2:
  # p1:=p1-t0q:
  
  # print("ps before final",p0,p1,p2):  
  
  # print("before",p1,t1q):
  c:=0:
  tmp:=0:
  p1,tmp:=addR(p1,t1q,tmp,r):
  # print("after",p1,tmp):
  p2:=p2+tmp:
  # p1:=p1-t0q-t2:

  # tmp:=0:
  # p1,tmp:=subR(p1,t0q,tmp,r):
  # p2:=p2+tmp:

  

  a1:=Array(1..8,i->0):
  a2:=Array(1..8,i->0):
  a1[1]:=p1:
  a1[2]:=p2:
  a2[1]:=t0q:
  a1:=bigSubMaple(a1,a2,r,p):
  p1:=a1[1]:
  p2:=a1[2]:



  tmp:=0:
  # p1:=p1-iquo(t0,r)-t2:
  p1,tmp:=subR(p1,t2,tmp,r):
  p2:=p2+tmp:
  tmp:=0:

  # # print(p1,"after"):
  # # p1:=2*s1-t0q-t2+t1q:
  # # print(p1-r,"before"):
  # # p1:=p1-iquo(t0,r)-t2:
  # p0:=p0-t0r+t1r:

  a1:=Array(1..8,i->0):
  a2:=Array(1..8,i->0):
  a1[1]:=p0:
  a1[2]:=p1:
  a1[3]:=p2:
  a2[1]:=t1r:
  # tmp:=0:
  # p0,tmp:=addR(p0,t1r,tmp,r):
  # p1,tmp:=addR(p1,tmp,0,r): 
  # p2:=p2+tmp:
  a1:=bigAddMaple(a1,a2,r,p):
  a2[1]:=t0r:
  a1:=bigSubMaple(a1,a2,r,p):
  p0:=a1[1]:
  p1:=a1[2]:
  p2:=a1[3]: 
  

  # p1:=p1+tmp:
  # p0:=p0-t0r:  
  # tmp:=0:
  # p0,tmp:=subR(p0,t0r,tmp,r):
  # # p1,tmp:=addR(p1,tmp,0,r):
  # p1:=p1+tmp:

  # p1:=p1+iquo(p0,r):
  # p0:=irem(p0,r):
  # p2:=p2+iquo(p1,r):
  # p1:=irem(p1,r):
  pN:=p0+p1*r+p2*(r^2):
  if p1<0 then
    p2:=p2-1:
    p1:=r+p1:
  fi:
  if p0>=r then
    p0:=p0-r:
    p1:=p1+1:
  fi:

  if p1>=r then
    p1:=p1-r:
    p2:=p2+1:
  fi:
  if p2>=r then
    p2:=p2-r:
    # p2:=p2+1:
  fi:

  # print("second",p0,p1,p2):

  sN:=s0+s1*(2^64)+s2*(2^128):
  c:=iquo(sN,r^2):
  sN:=irem(sN,r^2):
  h:=iquo(sN,r):
  sN:=irem(sN,r):
  l:=sN:
  # print("first,",l,h,c):

  # if p0=l then print("l equal"): fi:
  # if p1=h then print("h equal"): fi:
  # if p2=c then print("c equal"): fi:
  # sN:=s0+s1*(2^64)+s2*(2^128):
  # if p0=l and p1=h and p2=c then
  #   print("verified"):
  # fi:
  return p0,p1,p2:

end proc:
####################################

bigModulo:=proc(xs,r,p)
	lVector:=Array(1..k,i->0):
	hVector:=Array(1..k,i->0):
	cVector:=Array(1..k,i->0):

  for i from 1 to k do
    # l:=0:
    # h:=0:
    # c:=0:
    # #####################
    # s:=xs[i]:
    # c:=iquo(s,r^2):
    # s:=irem(s,r^2):
    # # print(s):
    # # print(i,"s",s):
    # h:=iquo(s,r):
    # s:=irem(s,r):
    # l:=s:

    # print(l,h,c):
    # if l<0 then 
    #   l:=l+r:
    #   h:=h-1:
    # fi:
    # if h<0 then
    #   h:=h+r:
    #   c:=c-1:
    # fi:
    # if c<0 then
    #   c:=r+c:
    # fi:
    #####################
    s:=xs[i]:
    s2:=iquo(s,2^128):
    s:=irem(s,2^128):
    s0:=irem(s,2^64):
    s1:=iquo(s,2^64):
    # print("u0,u1,u2",s0,s1,s2): 
    # print(s):
    # print(i,"s",s):
    # h:=iquo(s,r):
    # s:=irem(s,r):
    # l:=s:
    sArray:=Array(1..3,i->0):
    sArray:=[s0,s1,s2]:
    # print("third",l,h,c):
    res:=toLHC(sArray,r):
    # print("res",res):
    l:=res[1]:
    h:=res[2]:
    c:=res[3]:
    ####################
    # s:=xs[i]:
    # s2:=iquo(s,((2^64)^2)):
    # s:=irem(s,((2^64)^2)):
    # s1:=iquo(s,((2^64)^1)):
    # s:=irem(s,((2^64)^1)):
    # s0:=s:
    # printf("s0=%d, s1=%d, s2=%d \n",s0,s1,s2):
    # l:=irem(s0,r):
    # t:=iquo(s0,r)+irem(s1,r^2):
    # h:=irem(t,r^2):
    # c:=iquo(t,r^2)+s2:
    # print(l,h,c,"second-step"):
    # if l<0 then
    # l:=l+r:
    # h:=h-1:
    # fi:
    # if h<0 then
    # h:=h+r:
    # c:=c-1:
    # fi:
    cVector[i]:=c:
    hVector[i]:=h:
    lVector[i]:=l:
  end do:
  print("lvector ", lVector): 
  print("hvector ", hVector): 
  print("cvector ", cVector): 
  hVector:=simpleRightShift(hVector,r,p):
  
  h_tmp:=Array(1..k,i->0):
  h_tmp[1]:=hVector[1]:
  hVector[1]:=0:

  cVector:=simpleRightShift(cVector,r,p):
  cVector:=simpleRightShift(cVector,r,p):
  c_tmp:=Array(1..k,i->0):
  c_tmp[1]:=cVector[1]:
  c_tmp[2]:=cVector[2]:
  cVector[1]:=0:
  cVector[2]:=0:
  
  #l:=l+hr
	lVector:=bigAddMaple(lVector,hVector,r,p):
  #l:=l+cr^2
	lVector:=bigAddMaple(lVector,cVector,r,p):
  #l:=l+(-h)r
  lVector:=bigSubMaple(lVector,h_tmp,r,p):
  #l:=l+(-c)r^2
  lVector:=bigSubMaple(lVector,c_tmp,r,p):
	# print(lVector):
	return lVector:
end proc:

#######################################
bigModulo_v1:=proc(xs,r,p)
  lVector:=Array(1..k,i->0):
  hVector:=Array(1..k,i->0):
  cVector:=Array(1..k,i->0):

  for i from 1 to k do
    #####################
    s:=xs[i]:
    s2:=iquo(s,2^128):
    s:=irem(s,2^128):
    s0:=irem(s,2^64):
    s1:=iquo(s,2^64):
  
    sArray:=Array(1..3,i->0):
    sArray:=[s0,s1,s2]:
    # print("third",l,h,c):
    res:=toLHC(sArray,r):
    # print("res",res):
    l:=res[1]:
    h:=res[2]:
    c:=res[3]:
    ####################
    cVector[i]:=c:
    hVector[i]:=h:
    lVector[i]:=l:
  end do:
  return lVector,hVector,cVector;
  # print("lvector ", lVector): 
  # print("hvector ", hVector): 
  # print("cvector ", cVector): 
  # hVector:=simpleRightShift(hVector,r,p):
  
  # h_tmp:=Array(1..k,i->0):
  # h_tmp[1]:=hVector[1]:
  # hVector[1]:=0:

  # cVector:=simpleRightShift(cVector,r,p):
  # cVector:=simpleRightShift(cVector,r,p):
  # c_tmp:=Array(1..k,i->0):
  # c_tmp[1]:=cVector[1]:
  # c_tmp[2]:=cVector[2]:
  # cVector[1]:=0:
  # cVector[2]:=0:
  
  # #l:=l+hr
  # lVector:=bigAddMaple(lVector,hVector,r,p):
  # #l:=l+cr^2
  # lVector:=bigAddMaple(lVector,cVector,r,p):
  # #l:=l+(-h)r
  # lVector:=bigSubMaple(lVector,h_tmp,r,p):
  # #l:=l+(-c)r^2
  # lVector:=bigSubMaple(lVector,c_tmp,r,p):
  # # print(lVector):
  # return lVector:
end proc:

#######################################
add_128:=proc(x0,x1,y0,y1,x2)
  s:=0:
  c:=0:
  s:=x0+y0:
  s:=irem(s,2^64):
  if s<x0 or s<y0 then
    c:=1:
  fi:
  u0:=s:
  u1:=y1+c:
  s:=x1+u1:
  s:=irem(s,2^64):
  c:=0:
  if s<x1 or s<u1 then
    c:=1:
  fi:
  u1:=s:
  u2:=0:
  u2:=x2+c:
  return u0,u1,u2:
end proc:
#######################################
sub_128:=proc(x0,x1,y0,y1)
  s:=0:
  c:=0:
  if x0 >=y0 then 
    s:=x0-y0:
  else
    s:=2^64+x0-y0:
    c:=1:
  fi:
  u0:=s:
  c:=y1+c:
  if x1>c then
    s:=x1-c:
  else
    s:=2^64+x1-c:
  fi:
  u1:=s:
  if u1=2^64 then u1:=0: fi:
  return u0,u1:
end proc:
#######################################
sub_192_v0:=proc(x0,x1,x2,y0,y1,y2)
  # s:=0:
  # c:=0:
  # if x0 >=y0 then 
  #   s:=x0-y0:
  # else
  #   s:=2^64+x0-y0:
  #   c:=1:
  # fi:
  # u0:=s:
  # c:=y1+c:
  # if x1>=c then
  #   s:=x1-c:
  # else
  #   s:=2^64+x1-c:
  # fi:
  # u1:=s:

  # c:=y2+c:
  # if x2>=c then
  #   s:=x2-c:
  # else
  #   s:=2^64+x2-c:
  # fi:
  # u2:=s:
  # if u1=2^64 then u1:=0: u2:=1: fi:
  sX:=x0+x1*2^64+x2*2^128:
  sY:=y0+y1*2^64+y2*2^128:
  if sX>= sY then 
    s:=sX-sY:
    sign:=1;
  else
    s:=sY-sX:
    sign:=-1;
  fi;

  s2:=iquo(s,2^128):
  s:=irem(s,2^128):
  s1:=iquo(s,2^64):
  s0:=irem(s,2^64):
  # if s2<0 then s2:=2^64-s2; fi;
  # if s1<0 then s1:=2^64-s1; fi;
  # if s0<0 then s0:=2^64-s0; fi;
  # return s0,s1,s2:
  return s0,s1,s2,sign:
end proc:
#######################################
sub_192:=proc(x0,x1,x2,y0,y1,y2)
  sX:=x0+x1*2^64+x2*2^128:
  sY:=y0+y1*2^64+y2*2^128:
  if sX>= sY then 
    # s:=sX-sY:
    sign:=1;
  else
    # s:=sY-sX:
    sign:=-1;
  fi;

  s:=0:
  c:=0:
if sign =1 then
  if x0 >=y0 then 
    s:=x0-y0:
  else
    s:=2^64+x0-y0:
    c:=1:
  fi:
  u0:=s:
  c:=y1+c:
  if x1>=c then
    s:=x1-c:
    c:=0;
  else
    s:=2^64+x1-c:
    c:=1;
  fi:
  u1:=s:

  c:=y2+c:
  if x2>=c then
    s:=x2-c:
    c:=0;
  else
    s:=2^64+x2-c:
    c:=1;
  fi:
  u2:=s:
  if u1=2^64 then u1:=0: u2:=1: fi:
fi;

if sign =-1 then
  if y0 >=x0 then 
    s:=y0-x0:
  else
    s:=2^64+y0-x0:
    c:=1:
  fi:
  u0:=s:
  c:=x1+c:
  if y1>=c then
    s:=y1-c:
    c:=0;
  else
    s:=2^64+y1-c:
    c:=1;
  fi:
  u1:=s:

  c:=x2+c:
  if y2>=c then
    s:=y2-c:
    c:=0;
  else
    s:=2^64+y2-c:
    c:=1;
  fi:
  u2:=s:
  if u1=2^64 then u1:=0: u2:=1: fi:
fi;



  # s2:=iquo(s,2^128):
  # s:=irem(s,2^128):
  # s1:=iquo(s,2^64):
  # s0:=irem(s,2^64):
  # if s2<0 then s2:=2^64-s2; fi;
  # if s1<0 then s1:=2^64-s1; fi;
  # if s0<0 then s0:=2^64-s0; fi;
  # return s0,s1,s2:
  return u0,u1,u2,sign:
end proc:
#######################################
#computing bigAddition using Maple arithmetic
bigMult1_v0:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  # xN:=toBigInt(xs,r,p):
  # yN:=toBigInt(ys,r,p):
  ly:=ys:
  result:=Array(1..k,i->0):
  u:=Array(1..k,i->0):
  uSub:=Array(1..k,i->0):
  s:=0:
  for i from 1 to k do
  	if i>1 then
  		# ly:=multPowR(ly,1,r,p):
  		ly:=simpleRightShift(ly,r,p):
      # ly[1]:=ly[1]*(-1):
  		# ly:=normalize(ly,r,p):
  		# print("normalized",ly):
  	fi:
  	# for j from 1 to k do 
  	# 	s:=s+xs[j]*ly[k+1-j]:
  	# end do:

    #sum of negative values
    s:=0:
    # for j from k+2-i to k do 
    #   s:=s+xs[j]*ly[k+1-j]:
    # end do:
    # s2:=iquo(s,2^128):
    # s:=irem(s,2^128):
    # s1:=iquo(s,2^64):
    # s0:=irem(s,2^64):
    # print("before s0,s1,s2",s0,s1,s2):

    s0:=0:
    s1:=0:
    s2:=0:
    sArray:=Array(1..3,i->0):
    for j from k+2-i to k do 
      sArray:=[0,0,0]:
      m:=xs[j]*ly[k+1-j]:
      m1:=iquo(m,2^64):
      m0:=irem(m,2^64):
      # print("s2 before",s2):
      sArray:=add_128(s0,s1,m0,m1,s2):
      s0:=sArray[1]:
      s1:=sArray[2]:
      s2:=sArray[3]:
    end do:
    # print("afterr s0,s1,s2",s0,s1,s2):
    s:=s0+s1*(2^64)+s2*(2^128):
  	uSub[k+1-i]:=s:
    #sum of positive values
    s:=0:
    # for j from 1 to k+1-i do 
    #   s:=s+xs[j]*ly[k+1-j]:
    # end do:

    s0:=0:
    s1:=0:
    s2:=0:
    for j from 1 to k+1-i do 
      sArray:=[0,0,0]:
      m:=xs[j]*ly[k+1-j]:
      m1:=iquo(m,2^64):
      m0:=irem(m,2^64):
      # print("s2 before",s2):
      # print("m0,m1,j",m0,m1,j):
      sArray:=add_128(s0,s1,m0,m1,s2):
      s0:=sArray[1]:
      s1:=sArray[2]:
      s2:=sArray[3]:
    end do:
    print("i,s0,s1,s2,",i,s0,s1,s2):
    s:=s0+s1*(2^64)+s2*(2^128):
    u[k+1-i]:=s: 
  end do:

  resultSub:=bigModulo(uSub,r,p):# uSub-> l,h,c -> resultSub:=l+hr+cr^2
  result:=bigModulo(u,r,p):# u-> l,h,c -> result:=l+hr+cr^2
  result:=bigSubMaple(result,resultSub,r,p):
  return result:
end proc:
#######################################
bigMult1:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  # xN:=toBigInt(xs,r,p):
  # yN:=toBigInt(ys,r,p):
  ly:=ys:
  result:=Array(1..k,i->0):
  u:=Array(1..k,i->0):
  uSub:=Array(1..k,i->0):
  s:=0:

  lVector:=Array(1..k,i->0):
  hVector:=Array(1..k,i->0):
  cVector:=Array(1..k,i->0):
  lVectorSub:=Array(1..k,i->0):
  hVectorSub:=Array(1..k,i->0):
  cVectorSub:=Array(1..k,i->0):
  signVector:=Array(1..k,i->0):
  
  for i from 1 to k do
    if i>1 then
      # ly:=multPowR(ly,1,r,p):
      ly:=simpleRightShift(ly,r,p):
      # ly[1]:=ly[1]*(-1):
      # ly:=normalize(ly,r,p):
      # print("normalized",ly):
    fi:
    # for j from 1 to k do 
    #   s:=s+xs[j]*ly[k+1-j]:
    # end do:

    #sum of negative values
    s:=0:
    s0:=0:
    s1:=0:
    s2:=0:
    sArray:=Array(1..3,i->0):
    for j from k+2-i to k do 
      sArray:=[0,0,0]:
      m:=xs[j]*ly[k+1-j]:
      m1:=iquo(m,2^64):
      m0:=irem(m,2^64):
      # print("s2 before",s2):
      sArray:=add_128(s0,s1,m0,m1,s2):
      s0:=sArray[1]:
      s1:=sArray[2]:
      s2:=sArray[3]:
    end do:

    s0Sub:=s0;
    s1Sub:=s1:
    s2Sub:=s2;

    # print("afterr s0,s1,s2",s0,s1,s2):
    s:=s0+s1*(2^64)+s2*(2^128):
    uSub[k+1-i]:=s:
    #sum of positive values
    s:=0:
    # for j from 1 to k+1-i do 
    #   s:=s+xs[j]*ly[k+1-j]:
    # end do:
    # print("iSub,s0,s1,s2,",i,s0,s1,s2):
    s0:=0:
    s1:=0:
    s2:=0:
    for j from 1 to k+1-i do 
      sArray:=[0,0,0]:
      m:=xs[j]*ly[k+1-j]:
      m1:=iquo(m,2^64):
      m0:=irem(m,2^64):
      # print("s2 before",s2):
      # print("m0,m1,j",m0,m1,j):
      sArray:=add_128(s0,s1,m0,m1,s2):
      s0:=sArray[1]:
      s1:=sArray[2]:
      s2:=sArray[3]:
    end do:
    # print("i,s0,s1,s2,",i,s0,s1,s2):

    u[k+1-i]:=s:
    u4:=sub_192_v0(s0,s1,s2,s0Sub,s1Sub,s2Sub);
    # t:=u4;
    t:=Array(1..3,i->0);
    t[1]:=u4[1]:
    t[2]:=u4[2]:
    t[3]:=u4[3]:
    # print(t);
    s:=t[1]+t[2]*(2^64)+t[3]*(2^128):
    # print("s to lhc in mult1",s);
    l3:=toLHC(t,r):
    l8:=Array(1..k,i->0);
    l8[1]:=l3[1]; l8[2]:=l3[2]; l8[3]:=l3[3];
    if u4[4]=-1 then
      zArray:=Array(1..k,i->0):
      # print("sign is negative");
      l8:=bigSubMaple(zArray,l8,r,p);
    fi;
    lVector[k+1-i]:=l8[1]:
    hVector[k+1-i]:=l8[2]:
    cVector[k+1-i]:=l8[3]:
    if u4[4] = -1 then
      signVector[k+1-i]:=1:
    else
      signVector[k+1-i]:=0;
    fi;
    # print("------------------")
  end do:
  # print("lvector",lVector);
  # print("hvector",hVector);
  # print("cvector",cVector);
  hVector:=simpleRightShift(hVector,r,p):
  complementN:=r^8-r^3;
  # complementN:= complementN*r^3 mod p;
  # complement:=toBigFieldElement(complementN,r,p);
  # print("complement",complement);
  # for i from 1 to k do
  #   if signVector[i]=1 then
  #     tN:=complementN*r^(i-1) mod p;
  #     t:=toBigFieldElement(tN,r,p);
  #     print(t);
  #   fi;
  # end do;

  h_tmp:=Array(1..k,i->0):
  h_tmp[1]:=hVector[1]:
  hVector[1]:=0:
  # signVector:=bigSubMaple(zArray, signVector,r,p);
  subVector:=Array(1..k,i->0);
  # print("sign0",signVector);
  for i from 1 to k do
    if signVector[i]=1 then
      tN:=(complementN*r^(i-1)) mod p;
      t:=toBigFieldElement(tN,r,p);
      subVector:=bigAddMaple(subVector,t,r,p);
    fi;
  end do;
  # subVector:=bigSubMaple(zArray,subVector,r,p);
  lVector:=bigAddMaple(lVector,subVector,r,p);
  # print("sub",subVector);
  # signVector:=simpleeRightShift(signVector,r,p);
  # signVector:=simpleRightShift(signVector,r,p);
  # signVector:=simpleRightShift(signVector,r,p);
  # for i from 1 to 3 do
  #   signVector[i]:=r-signVector[i];
  # end do;
  # signN:=toBigInt(signVector,r,p);
  # signN:=signN*r^3 mod p;
  # signVector:=toBigFieldElement(signN,r,p);
  # print("sign2",signVector);
  # lVector:=bigSubMaple(lVector,signVector,r,p);

  cVector:=simpleRightShift(cVector,r,p):
  cVector:=simpleRightShift(cVector,r,p):
  # print("cvetor 2 shift",cVector);
  c_tmp:=Array(1..k,i->0):
  c_tmp[1]:=cVector[1]:
  c_tmp[2]:=cVector[2]:
  cVector[1]:=0:
  cVector[2]:=0:
  
  #l:=l+hr
  lVector:=bigAddMaple(lVector,hVector,r,p):
  # print("lvecto1",lVector);
  #l:=l+cr^2
  lVector:=bigAddMaple(lVector,cVector,r,p):
  # print("lvecto2",lVector);
  #l:=l+(-h)r
  lVector:=bigSubMaple(lVector,h_tmp,r,p):
  # print("lvecto3",lVector);
  #l:=l+(-c)r^2
  lVector:=bigSubMaple(lVector,c_tmp,r,p):
  # print("lvecto4",lVector);
  

  # print(lVector):
  ########################
  # result:=bigSubMaple(lVector,lVectorSub,r,p):
  return lVector:
end proc:
#####################################
bigMult2:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  # xN:=toBigInt(xs,r,p):
  # yN:=toBigInt(ys,r,p):
  ly:=ys:
  result:=Array(1..k,i->0):
  u:=Array(1..k,i->0):
  uSub:=Array(1..k,i->0):
  s:=0:

  lVector:=Array(1..k,i->0):
  hVector:=Array(1..k,i->0):
  cVector:=Array(1..k,i->0):
  lVectorSub:=Array(1..k,i->0):
  hVectorSub:=Array(1..k,i->0):
  cVectorSub:=Array(1..k,i->0):
  zArray:=Array(1..k,i->0):
  for i from 1 to k do
    if i>1 then
      # ly:=multPowR(ly,1,r,p):
      ly:=simpleRightShift(ly,r,p):
      ly[1]:=ly[1]*(-1):
      # ly:=normalize(ly,r,p):
      # print("normalized",ly):
    fi:
    # for j from 1 to k do 
    #   s:=s+xs[j]*ly[k+1-j]:
    # end do:

    #sum of negative values
    s:=0:
    s0:=0:
    s1:=0:
    s2:=0:
    sArray:=Array(1..3,i->0):
    for j from 1 to k do 
      sArray:=[0,0,0]:
      m:=xs[j]*ly[k+1-j]:
      # m1:=iquo(m,2^64):
      # m0:=irem(m,2^64):
      # # print("s2 before",s2):
      # # sArray:=add_128(s0,s1,m0,m1,s2):
      # # s0:=sArray[1]:
      # # s1:=sArray[2]:
      # # s2:=sArray[3]:
      s:=s+m;
    end do:
    
    u[k+1-i]:=s:
    # s:=s mod p:
    # print("s to lhc in mult2",s);
    # sList:=toBigFieldElement(s,r,p);
    s2:=iquo(s,r^2);
    s:=irem(s,r^2):
    s1:=iquo(s,r):
    s0:=irem(s,r):
    # s0:=sList[1];
    # s1:=sList[2];
    # s2:=sList[3];
    l3:=Array(1..3,i->0):
    l3[1]:=s0;
    l3[2]:=s1;
    l3[3]:=s2:
    # l3:=toLHC(t,r):

    l8:=Array(1..k,i->0);
    l8[1]:=l3[1]; l8[2]:=l3[2]; l8[3]:=l3[3];
    print("lhc in mult2",l8[1],l8[2],l8[3]);
    lVector[k+1-i]:=l8[1]:
    hVector[k+1-i]:=l8[2]:
    cVector[k+1-i]:=l8[3]:
    print("------------------")
  end do:

  # lVector:=result[1];
  # hVector:=result[2];
  # cVector:=result[3];
  print("lvector",lVector);
  print("hvector",hVector);
  print("cvector",cVector);
  hVector:=simpleRightShift(hVector,r,p):
  
  h_tmp:=Array(1..k,i->0):
  h_tmp[1]:=hVector[1]:
  hVector[1]:=0:

  cVector:=simpleRightShift(cVector,r,p):
  cVector:=simpleRightShift(cVector,r,p):
  c_tmp:=Array(1..k,i->0):
  c_tmp[1]:=cVector[1]:
  c_tmp[2]:=cVector[2]:
  cVector[1]:=0:
  cVector[2]:=0:
  
  #l:=l+hr
  lVector:=bigAddMaple(lVector,hVector,r,p):
  #l:=l+cr^2
  lVector:=bigAddMaple(lVector,cVector,r,p):
  #l:=l+(-h)r
  lVector:=bigSubMaple(lVector,h_tmp,r,p):
  #l:=l+(-c)r^2
  lVector:=bigSubMaple(lVector,c_tmp,r,p):
  # print(lVector):
  ########################
  # result:=bigSubMaple(lVector,lVectorSub,r,p):
  return lVector:
end proc:
#####################################
xData:=readdata(xFile,integer):
yData:=readdata(yFile,integer):

for i from 1 to 10 do
  # xs:=convert(xData[(i-1)*8+1..i*8],Array):
  # ys:=convert(yData[(i-1)*8+1..i*8],Array):
  xs:=xData[(i-1)*k+1..i*k]:
  ys:=yData[(i-1)*k+1..i*k]:
  xs:=xs[1..k]:
  ys:=ys[1..k]:
  m0:=bigMult0(xs,ys,r,p):
  m1:=bigMult1(xs,ys,r,p):
  equalBigFieldElements(m0,m1):
  # print("m1-m0",bigSubMaple(m1,m0,r,p));
  print(m0);
  print(m1);
  # print(xs):
  # quit;
end do:

# xs:=Array(1..k,i->r-2*i):
# ys:=Array(1..k,i->3*i):
# xs:=Array(1..k,i->0):
# ys:=Array(1..k,i->0):
# xs[8]:=r:
# ys[8]:=r:
# xs:=[1,1,1,1,2,0,0,0]:
# xs:=
# [4122069494891546112,
# 416448630046793856,
# 3768218334765421568,
# 4312036060060413440,
# 2597503149375243776,
# 6581614685471126528,
# 3507585816910166528,
# 7624981294075398144]:

# ys:=[1,1,1,1,2,0,0,0]:
# ys:=[
# 8244138994078059520,
# 832897264388555008,
# 7536436669530843136,
# 8624072120120826880,
# 5195006303045454848,
# 3939857316907610112,
# 7015171638115300352,
# 6026590529821184000]:

# m0:=bigMult0(xs,ys,r,p):
# m1:=bigMult1(xs,ys,r,p):
# m2:=bigMult2(xs,ys,r,p):
# print(m0):
# print(m1):
# print(m2):
# # equalBigFieldElements(m0,m1):
# print(bigSubMaple(m1,m0,r,p));


# subArray:=Array(1..8,i->0);
# subArray[1]:=1;
# zArray:=Array(1..8,i->0);
# print("zeroArray",zArray);
# res:=bigSubMaple(zArray,subArray,r,p);
# res[8]:=0;
# res[1]:=r-1;
# res[2]:=r-1;
# res[3]:=r-1;
# print(bigAddMaple(res,subArray,r,p));
# a:=(r-2)^2:
# d:=divR(a,r):
# m:=irem(a,r):
# q:=iquo(a,r):
# print(d,q,m):
# print(divR(rr-2))
# x[1]:=1:
# y[1]:=1:
# y[2]:=1:
# x:=multPowR(x,8,r,p):

# s:=Array(1..3,i->0):
# # ## for 8(r-2)^2 we have:
# # # s[3]:=2:
# # # s[2]:=137438953583:
# # # s[1]:=18446743523953737760:

# # s[3]:=2:
# # s[2]:=137438953583:
# # s[1]:=18446743523953737760:


# a:=Array(1..8,i->0);

# a[3]:= 10;
# a[2]:= 10;
# a[1]:= 10;

# s:=Array(1..8,i->0);
# s[3]:=11:
# s[2]:=11:
# s[1]:=11:
# zArray:=Array(1..8,i->0):
# sub3:=bigSubMaple(zArray,s,r,p):
# sub3:=bigSubMaple3_v1(zArray,s,r,p):
# print(sub3);
# print(sub3_v1);
# for i from 4 to 8 do
#   sub3[i]:=0;
# end do;
# # sub3[1]:=r-1;
# # sub3[2]:=r-2;
# # sub3[3]:=r-2;
# add3:=bigAddMaple(sub3,s,r,p);
# print(add3);

# print(sub_192_v0(a[1],a[2],a[3],s[1],s[2],s[3]));
# print(sub_192(a[1],a[2],a[3],s[1],s[2],s[3]));

# # s[1]:=18446743523953737760: s[2]:=137438953583: s[3]:=2 :
# s[1]:=18446743592673214492: s[2]:=13835058175541248097: s[3]:=1 :
# # s[1]:=18446743661392691224: s[2]:=9223372139933990995: s[3]:=1 :
# # s[1]:=18446743730112167956: s[2]:=4611686104326733893: s[3]:=1 :
# # s[1]:=18446743798831644688: s[2]:=68719476791: s[3]:=1 :
# # s[1]:=18446743867551121420: s[2]:=13835058106821771305: s[3]:=0 :
# # s[1]:=18446743936270598152: s[2]:=9223372071214514203: s[3]:=0 :
# # s[1]:=18446744004990074884: s[2]:=4611686035607257101: s[3]:=0 :

# d:=toLHC(s,r):
# print(d):
# sN:=s[1]+s[2]*2^64:
# print(divR(sN,r)):

# fclose(xFile):
# fclose(yFile):
# fclose(writeFile):


#need to do simple bigSubMaple on 8 digits, but the problem 
#is that with the current order, it takes 8*8 digits to store.

