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
r:=2^62+2^36;
#k: power of radix in big prime
k:=16;
#p: big prime in z/pz
p:=r^k+1:
#e in N=K^e
e:=3;
s:=e-2;
#vecSize: size of input vector
vecSize:=(2*k)^e:
K:=2*k;
#read data files for x and y
#stored in form of vectors of k=8 coefficients
# make sure you have executed "make generateData" in "BigPrimeFieldFFT_3"
libraryPath:=getenv("CUMODP_HOME");
projectPath:=sprintf("%s/new_tests/BigPrimeFieldFFT_4",libraryPath);
powersOfOmegaPath:=sprintf("%s/src/powers_of_omega_K%d",projectPath,e);
resultPath:=sprintf("DFT_K%d_result", e);
print(dataPath);
xDataPath:=sprintf("%s/test/data/data_0_11_22_1m/xData",projectPath);
# yDataPath:=sprintf("%s/test/data/data_0_11_22_1k/yData",projectPath);
# print(xDataPath);
# quit;
xFile:=fopen(xDataPath,READ,TEXT):
yFile:=fopen(powersOfOmegaPath,READ,TEXT):
# xFile:=fopen("datasets/tmp/xData",READ,TEXT):
# yFile:=fopen("datasets/tmp/yData",READ,TEXT):
# yFile:=fopen(dataPath,READ,TEXT):
uFile:=fopen(resultPath,WRITE,TEXT):
#######################################
# v is a vector of k machine words, 
# representing one element of z/pz
# toBigInt converts v to an integer n 


computeR_array:=proc(r_array)
  for i from 1 to k do 
    r_array[i]:=r^(i-1);
  end do;
  return r_array;
end proc;

toBigInt:=proc(v)
  n:=0:
  for i from 1 to k do
    n:=n+v[i]*r^(i-1):
  end do:
  return n:
  end proc:
#######################################
toBigIntPrecomputedR:=proc(v,r_array)
  n:=0:
  for i from 1 to k do
    n:=n+v[i]*r_array[i]:
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
  i:=4; #=log[2](k)
  x:=toBigFieldElementRecursive(n,i):
  x:=convert(x,vector):
  return x:
end proc;

# toBigFieldElement_inPlace:=proc(n,k,r_array)
#   local i;
#   x:=Vector(1..k,0);
#   m:=n;
#   for i from k-1 to 1 by -1 do 
#     x[i]:=iquo(m,r_array[i+1]);
#     m:=irem(m,r_array[i+1]);
#   end do;
#   x[1]:=m;
#   x:=convert(x,vector):
#   return x:
# end proc;

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
  cVector:=simpleRotateRight(cVector);    
  cVector:=simpleRotateRight(cVector);
  hVector[1]:=(-1)*hVector[1];
  cVector[1]:=(-1)*cVector[1];    cVector[2]:=(-1)*cVector[2];

  #result:=l+h*r
  result:=bigAddMaple(lVector,hVector);
  #result:=result+c*r^2
  result:=bigAddMaple(result,cVector);
  return result;
end proc:
#######################################
twiddle:=proc(xN,omega,i)
  resultN:=modp(xN*omega^(i),p);
  result:=toBigFieldElement(resultN,k);
  return result;
end proc:
#######################################
twiddle_inplace:=proc(xN,omega_in,i_in,s_in)
  K:=2*k;
  stride:=vecSize;
  # for s from 1 to s_in1 do
  #   stride:=stride/K;
  # end do;
  stride:=stride/ (K^(e-s_in));
  print("stide=",stride);
  omega_root:=omega_in;
  omega_root:=modp(omega_root^(K*(e-s_in-1)),p);
  i:=i_in;
  i:=irem(i,stride*K);
  # print("i",i);
  print(stride*K);
  i0:=irem(i,stride);
  j0:= iquo(i,stride);
  # print("i= i0,j0",i, i0,j0);
  resultN:=xN:
  if i0*j0 >0  then
    # print("i0*j0",i0*j0);
    h:=iquo((i0*j0),stride);
    w:=irem((i0*j0),stride);
    # h:=15;
    # print("w,h",w,h);
    resultN:=modp(xN*r^(h),p);
      if w>0 then
        resultN:=modp(resultN*omega_root^w,p);
      fi;
      # w:=1;
    # resultN:=modp(omega_root^(w),p);
  fi;

  result:=toBigFieldElement(resultN,k);
  # if result[1]=1898003780834354255 then print ("special",i); quit; fi;
  return result;
end proc:
#######################################
#L_{k}^{m}
permutation:=proc(xs,k,m)
tmp:=Array(1..k*m,i->0);
  for i from 1 to m do
    for j from 1 to k do
      tmp[j*m+i]:=xs[i*k+j];
    end do;
  end do;
  return tmp;
end proc:
#######################################
DFT_n:=proc(xs,omega_in,n)
  resultArray:=Array(1..n,i->0);
  omega_array:=Array(1..K^(e-1),i->0);
  omega_array_complete:=Array(1..n,i->0);
  op:=1;
  stride:=iquo(n,K);
  print("stride",stride);

  # for i from 1 to K^(e-1) do
  for i from 1 to stride do
      if irem(i,1000)=1 then 
        print("omega pow i",i-1);
      fi;
      # omega_array[i]:=modp(omega_in^((i-1)),p);
      omega_array[i]:=modp(op,p);
      op:=op*omega_in;
      # print(i, omega_array[i]);
  end do;

  xNN_array:=Array(1..n,i->0);
  # r_array:=Array(1..k,i->0);
  # r_array:=computeR_array(r_array);
  r_array:=Array(1..K,i->0);
  op:=1;
  for i from 1 to K do
    print("computing power of radix = ",i-1);
    # r_array[i]:=op;
    # op:=modp(op*r,p);
    r_array[i]:=modp(r^(i-1),p);
    print(r_array[i]);
    if r_array[i]=p then 
    print ("pow of r == p");
    quit;
    fi;
    print("===================");
  end do;

  for i from 1 to n do
    w:=irem(i-1, stride);
    h:=iquo(i-1, stride);
    omega_array_complete[i]:=r_array[h+1]*omega_array[w+1];
  end do;

  for j from 1 to n do
    if irem (j,1000)=1 then print("converting element no. to bigInt", j); fi;
    xN:=xs[(j-1)*k+1..j*k];
    # xNN:=toBigInt(xN);
    xNN:=toBigIntPrecomputedR(xN,r_array);
    # xNN:=1;
    # print(xNN);
    xNN_array[j]:=xNN;
    # if j>1 then 
    #   if(xNN_array[j]=xNN_array[j-1]) then
    #     print("Error! two consecutive are the same");
    #     quit;
    #     fi;
    # fi;
  end do;

  # result_array:=Array(1..n,i->0);
  # for i from 1 to n by stride do # to make it faster and have a have look at 
  #results of all powers of omega
  for i from 1 to n do
  result:=0;
  # print(i);
    for j from 1 to n do
      # h:=irem((i-1)*(j-1),n)+1;
      h:=irem((i-1)*(j-1),n);
      # w:=irem(h, stride);
      # h:=iquo(h, stride);

      # if i=2 then
      # print(i,j,h,w);
      # fi;
      # xN:=xs[(j-1)*k+1..j*k];
      # xNN:=toBigInt(xN);
      # result:=result+xNN_array[j]*omega_array[h]*r^w;
      
      # result:=modp(result+xNN_array[j]*modp(omega_array[h]*r_array[w+1],p),p);
      # if w=0 then 
      # mult:=xNN_array[j]*r_array[h+1]*omega_array[w+1];
      mult:=xNN_array[j]*omega_array_complete[h+1];
      # else
      
      # mult:=r_array[h+1];
      # print(mult);

      # if i=2 then
      #   print(toBigFieldElement(mult));
      # # print(r_array[h+1]);
      # # print(mult);
      # fi;
      # fi;
      result:=mult+result;
      result:=modp(result,p);  
      # result_array[i]:=modp(result_array[i]+xNN_array[j]*modp(omega_array[h]*r_array[w+1],p),p);
      # result:=result+xNN*omega_in^((i-1)*(j-1));
      # omega_in^((i-1)*(j-1));
      # result[i]:=result[i]+((j-1)*(i-1));
      # print("res[i],i,",result[i],i);
    end do;
    result:=modp(result,p);
    
    m0:=toBigFieldElement(result,k);
    # m0:=toBigFieldElement_inPlace(result,k,r_array);
    # print("i,j,h,w",i,j,h,w);
    print("i",m0);
    # if irem (i,10)=1 then print(i); fi;
    writedata(uFile,m0,integer);
    writeline(uFile,""); 
  end do;
  # return resultArray;
end proc:
#######################################

DFT_CT_n:=proc(xs,omega_in,n)
  resultArray:=Array(1..n,i->0);
  omega_array:=Array(1..K^(e-1),i->0);
  omega_array_complete:=Array(1..n,i->0);
  op:=1;
  stride:=iquo(n,K);
  print("stride",stride);

  # for i from 1 to K^(e-1) do
  for i from 1 to stride do
      if irem(i,1000)=1 then 
        print("omega pow i",i-1);
      fi;
      # omega_array[i]:=modp(omega_in^((i-1)),p);
      omega_array[i]:=modp(op,p);
      op:=op*omega_in;
      # print(i, omega_array[i]);
  end do;

  xNN_array:=Array(1..n,i->0);
  # r_array:=Array(1..k,i->0);
  # r_array:=computeR_array(r_array);
  r_array:=Array(1..K,i->0);
  op:=1;
  for i from 1 to K do
    print("computing power of radix = ",i-1);
    # r_array[i]:=op;
    # op:=modp(op*r,p);
    r_array[i]:=modp(r^(i-1),p);
    print(r_array[i]);
    if r_array[i]=p then 
    print ("pow of r == p");
    quit;
    fi;
    print("===================");
  end do;

  for i from 1 to n do
    w:=irem(i-1, stride);
    h:=iquo(i-1, stride);
    omega_array_complete[i]:=r_array[h+1]*omega_array[w+1];
  end do;

  for j from 1 to n do
    if irem (j,1000)=1 then print("converting element no. to bigInt", j); fi;
    xN:=xs[(j-1)*k+1..j*k];
    # xNN:=toBigInt(xN);
    xNN:=toBigIntPrecomputedR(xN,r_array);
    # xNN:=1;
    # print(xNN);
    xNN_array[j]:=xNN;
    if j>1 then 
      if(xNN_array[j]=xNN_array[j-1]) then
        print("Error! two consecutive are the same");
        quit;
        fi;
    fi;
  end do;

  # result_array:=Array(1..n,i->0);
  # for i from 1 to n by stride do # to make it faster and have a have look at 
  #results of all powers of omega

  # xNN_array:=permutation(xNN_array,K,stride);
  # for i from 1 to k do
  #   for j from 1 to stride do

  #   end do;
  # end do;
    
  #   m0:=toBigFieldElement(result,k);
  #   # m0:=toBigFieldElement_inPlace(result,k,r_array);
  #   # print("i,j,h,w",i,j,h,w);
  #   print("i",m0);
  #   # if irem (i,10)=1 then print(i); fi;
  #   writedata(uFile,m0,integer);
  #   writeline(uFile,""); 
  # end do;
  # return resultArray;
end proc:
#######################################
#reading data from file
xData:=readdata(xFile,integer):
yData:=readdata(yFile,integer):
###################
nSuccess:=0:
ys:=yData[1..k]:
yN:=toBigInt(ys);
# print(modp(yN^256,p));
# for i from 1 to vecSize do 
#   if modp(yN^i,p) =1 then
#   print (i,modp(yN^i,p));
#   fi;
# end do; 
# quit;

xs:=xData:
DFT_n(xs,yN,vecSize):
# print(result):

# for i from 1 to vecSize do
#   xs:=xData[(i-1)*k+1..i*k]:
#   # ys:=yData[(i-1)*k+1..i*k]:
#   # xs:=xs[1..k]:
#   # ys:=ys[1..k]:
#   xN:=toBigInt(xs);
#   # xN_twiddle:=modp(xN*yN^(i-1),p);
#   # m0:=toBigFieldElement(xN_twiddle,k);
#   # m0:=twiddle(xN,yN,i-1);
#   m1:=twiddle_inplace(xN,yN,i-1,s);
#   # m0:=bigMult0(xs,ys);
#   # print(m0);
#   print(m1);
#   m0:=m1;
#   # if irem(i,100)=1 then
#     print("i/vecSize",i, vecSize);
#     print("========================");
#   # fi;
#   # quit;
#   # m1:=bigMult1(xs,ys);
#   # print("idx",i);
#   # nSuccess:=nSuccess+ equalBigFieldElements(m0,m1);
#   # print("Maple Mult:",m0):
#   # print("revised Mult:",m1):
#   writedata(uFile,m0, integer);
#   writeline(uFile,"");
# end do:
###################
# if nSuccess = vecSize then
#   print("Passed!");
#   # `quit`(128):
# else
#   print("Failed!");
#   # `quit`(255):
# fi;
#######################################d cd 