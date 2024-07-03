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
# toLHC does not verify!
###############################
restart;
with(ArrayTools):
with(Ordinals):
r:=2^63+2^34;
k:=8;
p:=r^k+1:
ULMAX:= 2^64-1;
# size of input vector
vecSize:=1024:

#read data files for x and y
#stored in form of vectors of k=8 coefficients
xFile:=fopen("../../test/data/data_0_11_22_1k/xData",READ,TEXT):
yFile:=fopen("../../test/data/data_0_11_22_1k/yData",READ,TEXT):
writeFile:=fopen("../../test/tmp/uData_maple_0",WRITE, TEXT):
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
  return x:
end proc;

#######################################
equalBigFieldElements:=proc(xs,ys)
for i from 1 to k do
  if xs[i] <> ys[i] then
  print("FAIL FAIL FAIL FAIL FAIL inequal at index =",i);
  quit;
  return;
  fi;
  end do;
  print("equal");
end proc;

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
#multiplication by a power of radix
multPowR:=proc(xs,s,r,p)
  if s<=0 or s>=k
  then
  print("in multPowR: s is either leq 0 or geq k=",k);
  return xs;
  fi;
  xN:=toBigInt(xs,r,p):
  xN:=(xN*r^s)mod (p):
  result:=toBigFieldElement(xN,r,p,8):
  if result[8]=r then
    result[8]:=0;
    result[1]:=r-1;
    print("special case in multPowR");
  fi;
  return result:
end proc:
#######################################
simpleRightShift:=proc(xs,r,p)
  xShifted:=Array(1..k,i->0);
  t:=xs[k];
  for i from k to 2 by -1 do
    xShifted[i]:=xs[i-1];
  end do;
  xShifted[1]:=t*(1);
  return xShifted;
end proc;
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
	resultVector[k]:=0;
  	resultVector[1]:=r-1;
  	# print("special case in mult0");
  	fi;
  return resultVector
end proc:
#######################################
# subtraction modulo r
subR:=proc(x,y,c0,r)
s:=y;
if x>=s then
  s:=x-s;
  c:=0;
else:
  s:=r+x-s;
  c:=1;
fi;
if s>=r then
  s:=s-r;
  c:=1;
fi;
if c=1 then c:=r-1; fi;
c:=c+c0;
return s,c;
end proc;
#######################################
# addition modulo r
addR:=proc(x,y,c0,r)
  s:=irem(x+y,2^64);
  c:=0;
  if s<x or s<y then
    s:=s+2^64-r;
    c:=1;
  fi;
  if s>=r then:
    print("here");
    s:=s-r;
    c:=1;
  fi;
  c:=c+c0;
  c:=irem(c,r);
  return s,c;
end proc;
#######################################
divR_v4:=proc(xs,r)
  r0:=2^34;
  s:=xs;
  q:=0;
  m:=0;
  #s can be 128 bit, heavy arithmetic?
  f:=1;
  # while s>=r do
  if s>=2^64 then print ("s is huge"); fi;
  if s>=r then
    # print(s);
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    # print(s0,s1);
    q0:=iquo(s,r-r0);
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0));
      f2:=-1;
    else
      s:=irem(s,(r-r0))-q0*r0;
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
  fi;
  if s>=r then
    # print(s);
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    # print(s0,s1);
    q0:=iquo(s,r-r0);
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0));
      f2:=-1;
    else
      s:=irem(s,(r-r0))-q0*r0;
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
  fi;

  if s>=r then
    # print(s);
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    # print(s0,s1);
    q0:=iquo(s,r-r0);
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;

    if irem(s,(r-r0))<q0*r0 then
      s:=q0*r0-irem(s,(r-r0));
      f2:=-1;
    else
      s:=irem(s,(r-r0))-q0*r0;
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
  fi;

  # print(s,q,f);
  if f=-1 and s<>0 then
    # print("f = -1");
    s:=r-s;
    q:=q-1;
  fi;
  m:=s;
  return q,m;
end proc;
#######################################
divR:=proc(xs,r)
  r0:=2^34;
  s:=xs;
  q:=0;
  m:=0;
  f:=1;
  s1:=iquo(s,2^64); #iquo(s,r-r0)
  s0:=irem(s,2^64); #irem(s,r-r0)

  # if s>=r then
  if s1>0 then print ("s is huge"); fi;
  t1:=2*s1;
  t0:=iquo(s0,2^63);
  q0:=t0+t1;
  if q0>0 then
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    m0:=irem(s0,2^63);
    # q0:=iquo(s,r-r0);
    # q0:=t1+t0;
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;
    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64);
    s0:=irem(q0*r0,2^64);
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0));
      # s:=q0*r0-m0;
      sl:=sub_128(s0,s1,m0,0);
      f2:=-1;
    else
      # s:=irem(s,(r-r0))-q0*r0;
      # s:=m0-q0*r0;
      sl:=sub_128(m0,0,s0,s1);
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
    s0:=sl[1];
    s1:=sl[2];
    # print("s0,s1",s0,s1);
  fi;
  # s:=s0+s1*2^64;
  # s1:=iquo(s,2^64); #iquo(s,r-r0)
  # s0:=irem(s,2^64); #irem(s,r-r0)
  # print("s0,s1",s0,s1);
  t1:=2*s1;
  t0:=iquo(s0,2^63);
  q0:=t0+t1;
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s);
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    # print(s0,s1);
    # t1:=2*s1;
    # t0:=iquo(s0,2^63);
    m0:=irem(s0,2^63);
    # q0:=iquo(s,r-r0);
    # q0:=t1+t0;
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64);
    s0:=irem(q0*r0,2^64);
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0));
      # s:=q0*r0-m0;
      sl:=sub_128(s0,s1,m0,0);
      f2:=-1;
    else
      # s:=irem(s,(r-r0))-q0*r0;
      # s:=m0-q0*r0;
      sl:=sub_128(m0,0,s0,s1);
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
    s0:=sl[1];
    s1:=sl[2];
    # print("s0,s1",s0,s1);
  fi;
  # s:=s0+s1*2^64;
  # s1:=iquo(s,2^64); #iquo(s,r-r0)
  # s0:=irem(s,2^64); #irem(s,r-r0)
 t1:=2*s1;
  t0:=iquo(s0,2^63);
  q0:=t0+t1;
  # if 2*s1+iquo(s0,2^63)>0 then
  if q0>0 then
    # print(s);
    # q0 is a monomial, can be computed by shifting s 
    # s1:=iquo(s,2^64);
    # s0:=irem(s,2^64);
    # print(s0,s1);
    # t1:=2*s1;
    # t0:=iquo(s0,2^63);
    m0:=irem(s0,2^63);
    # q0:=iquo(s,r-r0);
    # q0:=t1+t0;
    # q0:=iquo(s0+s1*2^64,r-r0);
    q:=q+q0*f;

    # if irem(s,(r-r0))<q0*r0 then
    s1:=iquo(q0*r0,2^64);
    s0:=irem(q0*r0,2^64);
    if m0 < q0*r0 then
      # s:=q0*r0-irem(s,(r-r0));
      # s:=q0*r0-m0;
      sl:=sub_128(s0,s1,m0,0);
      f2:=-1;
    else
      # s:=irem(s,(r-r0))-q0*r0;
      # s:=m0-q0*r0;
      sl:=sub_128(m0,0,s0,s1);
      f2:=1;
    fi;
    # if s<0 then f:=-1; else f:=1; fi;
    # s:=s*f;
    f:=f*f2;
    # print("s",s);
    s0:=sl[1];
    s1:=sl[2];
    # print("s0,s1",s0,s1);
  fi;
  # s:=s0+s1*2^64;
  # s1:=iquo(s,2^64); #iquo(s,r-r0)
  # s0:=irem(s,2^64); #irem(s,r-r0)

  if f=-1 and s0<>0 then
    # print("f = -1");
    print("s0,s1",s0,s1);
    # s:=r-s;
    s0:=r-s0;
    q:=q-1;
  fi;
  m:=s0;
  return q,m;
end proc;
#######################################
# s (three machine-words) -> (l,h,c) modulo r
toLHC:=proc (s,r)
  s0:=s[1]:
  s1:=s[2];
  s2:=s[3];
  r0:=2^34;
  q:=iquo(2^64,r-r0);
  # p2:=s2*q^2;
  p2:=4*s2;
  # p1:=s1*q;
  # p1:=2*s1;
  if s1<2^63 then 
    p1:=2*s1;
  else
    p1:=2*s1;
    p2:=p2+iquo(p1,r);
    p1:=irem(p1,r);
  fi;

  p0:=s0;
  t0:=(2*s1*r0);
  t1:=(4*s2*r0^2);
  t2:=8*r0*s2;
  if t0> 2^64 then print("t0 is huge"); fi;
  if t1>2^64 then print("t1 is huge"); fi;
  t0q:=iquo(t0,r):
  t0r:=irem(t0,r);
  d:=divR(t0,r);
  t0q:=d[1];
  t0r:=d[2];

  t1q:=iquo(t1,r);
  t1r:=irem(t1,r);
  d:=divR(t1,r);
  t1q:=d[1];
  t1r:=d[2];

  # p1:=p1-iquo(2*s1*r0,r)+iquo(4*s2*r0^2,r)-8*r0*s2;
  # p0:=p0-rem(2*s1*r0,r)+rem(4*s2*r0^2,r);
  # p1:=p1-iquo(t0,r)+iquo(t1,r)-t2;
  # p0:=p0-rem(t0,r)+rem(t1,r);
  # p1:=p1-t0q+t1q-t2;
  # p1:=p1-t0q;

  
  print("before",p1,t1q);
  c:=0;
  tmp:=0;
  p1,tmp:=addR(p1,t1q,tmp,r);
  print("after",p1,tmp);
  p2:=p2+tmp;
  # p1:=p1-t0q-t2;
  tmp:=0;
  p1,tmp:=subR(p1,t0q,tmp,r);
  p2:=p2+tmp;
  
  tmp:=0;
  # p1:=p1-iquo(t0,r)-t2;
  p1,tmp:=subR(p1,t2,tmp,r);
  p2:=p2+tmp;
  tmp:=0;
  # # print(p1,"after");
  # # p1:=2*s1-t0q-t2+t1q;
  # # print(p1-r,"before");
  # # p1:=p1-iquo(t0,r)-t2;
  # p0:=p0-t0r+t1r;

  a1:=Array(1..8,i->0);
  a2:=Array(1..8,i->0);
  a1[1]:=p0;
  a1[2]:=p1;
  a1[3]:=p2;
  a2[1]:=t1r;
  # tmp:=0;
  # p0,tmp:=addR(p0,t1r,tmp,r);
  # p1,tmp:=addR(p1,tmp,0,r); 
  # p2:=p2+tmp;
  a1:=bigAddMaple(a1,a2,r,p);
  a2[1]:=t0r;
  a1:=bigSubMaple(a1,a2,r,p);
  p0:=a1[1];
  p1:=a1[2];
  p2:=a1[3];
  # p1:=p1+tmp;
  # p0:=p0-t0r;  
  # tmp:=0;
  # p0,tmp:=subR(p0,t0r,tmp,r);
  # # p1,tmp:=addR(p1,tmp,0,r);
  # p1:=p1+tmp;

  # p1:=p1+iquo(p0,r);
  # p0:=irem(p0,r);
  # p2:=p2+iquo(p1,r);
  # p1:=irem(p1,r);
  pN:=p0+p1*r+p2*(r^2);
  if p1<0 then
    p2:=p2-1;
    p1:=r+p1;
  fi;
  if p0>=r then
    p0:=p0-r;
    p1:=p1+1;
  fi;

  if p1>=r then
    p1:=p1-r;
    p2:=p2+1;
  fi;
  if p2>=r then
    p2:=p2-r;
    # p2:=p2+1;
  fi;

  print("second",p0,p1,p2);

  sN:=s0+s1*(2^64)+s2*(2^128);
  c:=iquo(sN,r^2);
  sN:=irem(sN,r^2);
  h:=iquo(sN,r);
  sN:=irem(sN,r);
  l:=sN;
  print("first,",l,h,c);

  # if p0=l then print("l equal"); fi;
  # if p1=h then print("h equal"); fi;
  # if p2=c then print("c equal"); fi;
  # sN:=s0+s1*(2^64)+s2*(2^128);
  if p0=l and p1=h and p2=c then
    print("verified");
  fi;
  return p0,p1,p2;

end proc;
####################################
bigModulo:=proc(xs,r,p)
	lVector:=Array(1..k,i->0):
	hVector:=Array(1..k,i->0):
	cVector:=Array(1..k,i->0):

  for i from 1 to k do
    l:=0;
    h:=0;
    c:=0;
    # #####################
    s:=xs[i];
    c:=iquo(s,r^2);
    s:=irem(s,r^2);
    # print(s);
    # print(i,"s",s);
    h:=iquo(s,r):
    s:=irem(s,r);
    l:=s;

    # print(l,h,c);
    # if l<0 then 
    #   l:=l+r;
    #   h:=h-1;
    # fi;
    # if h<0 then
    #   h:=h+r;
    #   c:=c-1;
    # fi;
    # if c<0 then
    #   c:=r+c;
    # fi;
    #####################
    s:=xs[i];
    s2:=iquo(s,2^128);
    s:=irem(s,2^128);
    s0:=irem(s,2^64);
    s1:=iquo(s,2^64); 
    # print(s);
    # print(i,"s",s);
    # h:=iquo(s,r):
    # s:=irem(s,r);
    # l:=s;
    sArray:=Array(1..3,i->0);
    sArray:=[s0,s1,s2];    
    print("third",l,h,c);
    res:=toLHC(sArray,r);
    # print("res",res);
    l:=res[1];
    h:=res[2];
    c:=res[3];
    ####################
    # s:=xs[i];
    # s2:=iquo(s,((2^64)^2));
    # s:=irem(s,((2^64)^2));
    # s1:=iquo(s,((2^64)^1));
    # s:=irem(s,((2^64)^1));
    # s0:=s;
    # printf("s0=%d, s1=%d, s2=%d \n",s0,s1,s2);
    # l:=irem(s0,r);
    # t:=iquo(s0,r)+irem(s1,r^2);
    # h:=irem(t,r^2);
    # c:=iquo(t,r^2)+s2;
    # print(l,h,c,"second-step");
    # if l<0 then
    # l:=l+r;
    # h:=h-1;
    # fi;
    # if h<0 then
    # h:=h+r;
    # c:=c-1;
    # fi;
    cVector[i]:=c;
    hVector[i]:=h;
    lVector[i]:=l;
  end do;
  # lVector:=normalize(lVector,r,p);
  print("lvector",lVector);
  print("hvector",hVector);
  print("cvector",cVector);
	# print(lVector,hVector,cVector);
	# print(hVector);
	# hVector:=multPowR(hVector,1,r,p);
  # print(hVector);
  hVector:=simpleRightShift(hVector,r,p);
  # print(hVector);
  h_tmp:=Array(1..k,i->0);
  h_tmp[1]:=hVector[1];
  hVector[1]:=0;

  cVector:=simpleRightShift(cVector,r,p);
  cVector:=simpleRightShift(cVector,r,p);
  c_tmp:=Array(1..k,i->0);
  c_tmp[1]:=cVector[1];
  c_tmp[2]:=cVector[2];
  cVector[1]:=0;
  cVector[2]:=0;
	# print(hVector);
	# cVector:=multPowR(cVector,2,r,p);
  
  #l:=l+hr
	lVector:=bigAddMaple(lVector,hVector,r,p);
  #l:=l+cr^2
	lVector:=bigAddMaple(lVector,cVector,r,p);
  #l:=l+(-h)r
  lVector:=bigSubMaple(lVector,h_tmp,r,p);
  #l:=l+(-c)r^2
  lVector:=bigSubMaple(lVector,c_tmp,r,p);
	# print(lVector);
	return lVector;
end proc;

#######################################
add_128:=proc(x0,x1,y0,y1,c0)
  s:=0;
  c:=0;
  s:=x0+y0;
  s:=irem(s,2^64);
  if s<x0 or s<y0 then
    c:=1;
  fi;
  u0:=s;
  u1:=y1+c;
  s:=x1+u1;
  s:=irem(s,2^64);
  c:=0;
  if s<x1 or s<u1 then
    c:=1;
  fi;
  u1:=s;
  u2:=0;
  # print("c0",c0,"c",c,"u2",u2);
  u2:=c0+c;
  return u0,u1,u2;
end proc;
#######################################
sub_128:=proc(x0,x1,y0,y1)
  s:=0;
  c:=0;
  if x0 >y0 then 
    s:=x0-y0;
  else
    s:=2^64+x0-y0;
    c:=1;
  fi;
  u0:=s;
  c:=y1+c;
  if x1>c then
    s:=x1-c;
  else
    s:=2^64+x1-c;
  fi;
  u1:=s;
  if u1=2^64 then u1:=0; fi;
  return u0,u1;
end proc;
#######################################
#computing bigAddition using Maple arithmetic
bigMult1:=proc(xs,ys,r,p)
  #convert xs and ys to bigInt xN and yN, respectively
  # xN:=toBigInt(xs,r,p):
  # yN:=toBigInt(ys,r,p):
  ly:=ys;
  result:=Array(1..k,i->0);
  u:=Array(1..k,i->0);
  uSub:=Array(1..k,i->0);
  s:=0;
  for i from 1 to k do
  	if i>1 then
  		# ly:=multPowR(ly,1,r,p);
  		ly:=simpleRightShift(ly,r,p);
      # ly[1]:=ly[1]*(-1);
  		# ly:=normalize(ly,r,p);
  		# print("normalized",ly);
  	fi;
  	# for j from 1 to k do 
  	# 	s:=s+xs[j]*ly[k+1-j];
  	# end do;

    #sum of negative values
    s:=0;
    for j from k+2-i to k do 
      s:=s+xs[j]*ly[k+1-j];
    end do;
    s2:=iquo(s,2^128);
    s:=irem(s,2^128);
    s1:=iquo(s,2^64);
    s0:=irem(s,2^64);
    # print("before s0,s1,s2",s0,s1,s2);

    s0:=0;
    s1:=0;
    s2:=0;
    sArray:=Array(1..3,i->0);
    for j from k+2-i to k do 
      sArray:=[0,0,0];
      m:=xs[j]*ly[k+1-j];
      m1:=iquo(m,2^64);
      m0:=irem(m,2^64);
      # print("s2 before",s2);
      sArray:=add_128(s0,s1,m0,m1,s2);
      s0:=sArray[1];
      s1:=sArray[2];
      s2:=sArray[3];
    end do;
    # print("afterr s0,s1,s2",s0,s1,s2);
    s:=s0+s1*(2^64)+s2*(2^128);
  	uSub[k+1-i]:=s;
    #sum of positive values
    s:=0;
    # for j from 1 to k+1-i do 
    #   s:=s+xs[j]*ly[k+1-j];
    # end do;

    s0:=0;
    s1:=0;
    s2:=0;
    for j from 1 to k+1-i do 
      sArray:=[0,0,0];
      m:=xs[j]*ly[k+1-j];
      m1:=iquo(m,2^64);
      m0:=irem(m,2^64);
      # print("s2 before",s2);
      sArray:=add_128(s0,s1,m0,m1,s2);
      s0:=sArray[1];
      s1:=sArray[2];
      s2:=sArray[3];
    end do;
    s:=s0+s1*(2^64)+s2*(2^128);
    u[k+1-i]:=s; 
  end do;

  result:=bigModulo(u,r,p);# u-> l,h,c -> result:=l+hr+cr^2
  resultSub:=bigModulo(uSub,r,p);# uSub-> l,h,c -> resultSub:=l+hr+cr^2
  result:=bigSubMaple(result,resultSub,r,p);
  return result;
end proc:

#######################################
xData:=readdata(xFile,integer):
yData:=readdata(yFile,integer):

# for i from 1 to 10 do
#   # xs:=convert(xData[(i-1)*8+1..i*8],Array):
#   # ys:=convert(yData[(i-1)*8+1..i*8],Array):
#   xs:=xData[(i-1)*k+1..i*k]:
#   ys:=yData[(i-1)*k+1..i*k]:
#   xs:=xs[1..k]:
#   ys:=ys[1..k]:
#   m0:=bigMult0(xs,ys,r,p):
#   m1:=bigMult1(xs,ys,r,p):
#   equalBigFieldElements(m0,m1):
#   print(m0):
#   print(m1):
#   # print(xs);
# end do:

xs:=Array(1..k,i->r-2);
ys:=Array(1..k,i->r-2);
# xs[8]:=1;
# xs:=[1,1,1,1,2,0,0,0];
# ys:=[1,1,1,1,2,0,0,0];
m0:=bigMult0(xs,ys,r,p):
m1:=bigMult1(xs,ys,r,p):
print(m0):
print(m1):
equalBigFieldElements(m0,m1):

# a:=(r-2)^2;
# d:=divR(a,r);
# m:=irem(a,r);
# q:=iquo(a,r);
# print(d,q,m);
# print(divR(rr-2))
# x[1]:=1;
# y[1]:=1;
# y[2]:=1;
# x:=multPowR(x,8,r,p);

# fclose(xFile):
# fclose(yFile):
# fclose(writeFile):