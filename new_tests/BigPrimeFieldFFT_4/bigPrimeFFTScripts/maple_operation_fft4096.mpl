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

r:=2^63+2^34;
p:=r^8+1:
ULMAX:= 2^64-1;

#read data files for x and y
#stored in form of vectors of k=8 coefficients
xFile:=fopen("../../test/tmp/xData",READ,TEXT):
# yFile:=fopen("../../test/tmp/yData",READ,TEXT):
writeFile:=fopen("../../test/tmp/uData_maple_fft4096",WRITE,TEXT):

# xFile:=fopen("/home/dave/Desktop/DESK_SEPTEMBER/FFT/repository/cumodp/new_tests/BigPrimeFieldFFT_3/test/data/data_2_11_22_1k/xData",READ,TEXT);

# xFile:=fopen("../test/data/data_2_11_22_1k/xData",READ,TEXT):

# writeFile:=fopen("uData_maple_fft4096",WRITE,TEXT):

############################################################
#convert array[1..8] to n, bigPrimeNumber element
toBigInt:=proc(v,r,p)
  n:=0:
  for i from 1 to 8 do
    n:=n+v[i]*r^(i-1):
  end do:
  return n:
  end proc:

############################################################
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

############################################################
#function for converting big int n to array [1..8]
toBigFieldElement:=proc(n,r,p,k)
  i:=3; #=log[2](k)
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

##################################################
#This function is to calculate FFT with size 16
FFT16:=proc(xx,N, w,p)
  #input x array  w array  
  #in-place
  local i, num1,num2,x,ws;
  ws:=Array(1..8):
  for i from 1 to 8 do 
    ws[i]:=w^(i-1) mod p:
  end do:
  x:=xx;
  #first round
  num1:=x[1]:
  num2:=x[9]:
  x[1]:=(num1+num2) mod p:
  x[9]:=(num1-num2) mod p:
  num1:=x[5]:
  num2:=x[13]:
  x[5]:=(num1+num2) mod p:
  x[13]:=(num1-num2) mod p:
  num1:=x[3]:
  num2:=x[11]:
  x[3]:=(num1+num2) mod p:
  x[11]:=(num1-num2) mod p:
  num1:=x[7]:
  num2:=x[15]:
  x[7]:=(num1+num2) mod p:
  x[15]:=(num1-num2) mod p:
  num1:=x[2]:
  num2:=x[10]:
  x[2]:=(num1+num2) mod p:
  x[10]:=(num1-num2) mod p:
  num1:=x[6]:
  num2:=x[14]:
  x[6]:=(num1+num2) mod p:
  x[14]:=(num1-num2) mod p:
  num1:=x[4]:
  num2:=x[12]:
  x[4]:=(num1+num2) mod p:
  x[12]:=(num1-num2) mod p:
  num1:=x[8]:
  num2:=x[16]:
  x[8]:=(num1+num2) mod p:
  x[16]:=(num1-num2) mod p:
  
  #second round
  num1:=x[1]:
  num2:=x[5]:
  x[1]:=(num1+num2) mod p:
  x[5]:=(num1-num2) mod p:
  num1:=x[9]:
  num2:=x[13]:
  x[9]:=(num1+ws[5]*num2) mod p:
  x[13]:=(num1-ws[5]*num2) mod p:
  num1:=x[3]:
  num2:=x[7]:
  x[3]:=(num1+num2) mod p:
  x[7]:=(num1-num2) mod p:
  num1:=x[11]:
  num2:=x[15]:
  x[11]:=(num1+ws[5]*num2) mod p:
  x[15]:=(num1-ws[5]*num2) mod p:
  num1:=x[2]:
  num2:=x[6]:
  x[2]:=(num1+num2) mod p:
  x[6]:=(num1-num2) mod p:
  num1:=x[10]:
  num2:=x[14]:
  x[10]:=(num1+ws[5]*num2) mod p:
  x[14]:=(num1-ws[5]*num2) mod p:  
  num1:=x[4]:
  num2:=x[8]:
  x[4]:=(num1+num2) mod p:
  x[8]:=(num1-num2) mod p:
  num1:=x[12]:
  num2:=x[16]:
  x[12]:=(num1+ws[5]*num2) mod p:
  x[16]:=(num1-ws[5]*num2) mod p:
  
  #third round
  num1:=x[1]:
  num2:=x[3]:
  x[1]:=(num1+num2) mod p:
  x[3]:=(num1-num2) mod p:
  
  num1:=x[9]:
  num2:=x[11]:
  x[9]:=(num1+ws[3]*num2) mod p:
  x[11]:=(num1-ws[3]*num2) mod p:
  
  num1:=x[5]:
  num2:=x[7]:
  x[5]:=(num1+ws[5]*num2) mod p:
  x[7]:=(num1-ws[5]*num2) mod p:
  
  num1:=x[13]:
  num2:=x[15]:
  x[13]:=(num1+ws[7]*num2) mod p:
  x[15]:=(num1-ws[7]*num2) mod p:
  
  num1:=x[2]:
  num2:=x[4]:
  x[2]:=(num1+num2) mod p:
  x[4]:=(num1-num2) mod p:
  
  num1:=x[10]:
  num2:=x[12]:
  x[10]:=(num1+ws[3]*num2) mod p:
  x[12]:=(num1-ws[3]*num2) mod p:
  
  num1:=x[6]:
  num2:=x[8]:
  x[6]:=(num1+ws[5]*num2) mod p:
  x[8]:=(num1-ws[5]*num2) mod p:
  
  num1:=x[14]:
  num2:=x[16]:
  x[14]:=(num1+ws[7]*num2) mod p:
  x[16]:=(num1-ws[7]*num2) mod p:
  
  #fourth
  num1:=x[1]:
  num2:=x[2]:
  x[1]:=(num1+num2) mod p:
  x[2]:=(num1-num2) mod p:
  num1:=x[9]:
  num2:=x[10]:
  x[9]:=(num1+ws[2]*num2) mod p:
  x[10]:=(num1-ws[2]*num2) mod p:
  num1:=x[5]:
  num2:=x[6]:
  x[5]:=(num1+ws[3]*num2) mod p:
  x[6]:=(num1-ws[3]*num2) mod p:
  num1:=x[13]:
  num2:=x[14]:
  x[13]:=(num1+ws[4]*num2) mod p:
  x[14]:=(num1-ws[4]*num2) mod p:
  num1:=x[3]:
  num2:=x[4]:
  x[3]:=(num1+ws[5]*num2) mod p:
  x[4]:=(num1-ws[5]*num2) mod p:
  num1:=x[11]:
  num2:=x[12]:
  x[11]:=(num1+ws[6]*num2) mod p:
  x[12]:=(num1-ws[6]*num2) mod p:
  num1:=x[7]:
  num2:=x[8]:
  x[7]:=(num1+ws[7]*num2) mod p:
  x[8]:=(num1-ws[7]*num2) mod p:
  num1:=x[15]:
  num2:=x[16]:
  x[15]:=(num1+ws[8]*num2) mod p:
  x[16]:=(num1-ws[8]*num2) mod p:
  
  #2 9
  num1:=x[2]:
  x[2]:=x[9]:
  x[9]:=num1:
  #3 5
  num1:=x[3]:
  x[3]:=x[5]:
  x[5]:=num1:
  #4 13
  num1:=x[4]:
  x[4]:=x[13]:
  x[13]:=num1:
  #6 11
  num1:=x[6]:
  x[6]:=x[11]:
  x[11]:=num1:
  #8 15
  num1:=x[8]:
  x[8]:=x[15]:
  x[15]:=num1:
  #12 14
  num1:=x[12]:
  x[12]:=x[14]:
  x[14]:=num1:
  return x;
end proc:

##################################################
#This function is to calculate FFT with size 256. It calls FFT16.
FFT256:=proc(a,N,w,p)
  local J,K,rs,cs,ds,es,b,j,k:
  
  J:=16:
  K:=16:
  rs:=Array(1..16,1..16):
  cs:=Array(1..16,1..16):
  ds:=Array(1..16,1..16):
  es:=Array(1..16,1..16):
  b:=Array(1..256): #store the final result
  
  for j from 1 to J do
    for k from 1 to K do
      rs[j,k]:=a[(k-1)*J+j]:
    end do:
  end do:
  
  for j from 1 to J do
    cs[j]:=FFT16(rs[j], 16, (w^16) mod p, p); 
  end do:
  
  for k from 1 to K do
    for j from 1 to J do
      ds[k,j]:=(cs[j,k]*w^((j-1)*(k-1))) mod p:
    end do:
  end do:
  
  for k from 1 to K do
    es[k]:=FFT16(ds[k],16, (w^16) mod p, p);
  end do:  
  
  for k from 1 to K do 
    for j from 1 to J do
      b[(j-1)*K+k]:=es[k,j]:
    end do:
  end do:
  return b;
end proc:

#Test for FFT256 
# b := FFT256(convert([seq(i, i = 1 .. 256)], Array), 256, 3, 257):
# convert(b, list);

##################################################
#This function is to calculate FFT with size 256. It calls FFT16.
FFT4096:=proc(a,N,w,p)
  local J,K,rs,cs,ds,es,b,j,k:
  
  J:=256:
  K:=16:
  rs:=Array(1..J,1..K):
  cs:=Array(1..J,1..K):
  ds:=Array(1..K,1..J):
  es:=Array(1..K,1..J):
  b:=Array(1..N): #store the final result
  
  for j from 1 to J do
    for k from 1 to K do
      rs[j,k]:=a[(k-1)*J+j]:
    end do:
  end do:
  
  for j from 1 to J do
    cs[j]:=FFT16(rs[j], 16, (w^J) mod p, p); 
  end do:
  
  for k from 1 to K do
    for j from 1 to J do
      ds[k,j]:=(cs[j,k]*w^((j-1)*(k-1))) mod p:
    end do:
  end do:
  
  for k from 1 to K do
    es[k]:=FFT256(ds[k], 256, (w^K) mod p, p);
  end do:  
  
  for k from 1 to K do 
    for j from 1 to J do
      b[(j-1)*K+k]:=es[k,j]:
    end do:
  end do:
  return b;
end proc:

#Test for FFT4096
# b := FFT4096(convert([seq(i, i = 1 .. 4096)], Array), 4096, 41, 12289):

omega=41:

vecSize:=4096:

# read vectors of size 8 from file, each vector represents a big prime field element
xData:=readdata(xFile,integer):
xN:=Array(1..vecSize):

# convert each vector to a big integer
for i from 1 to vecSize do
  xs:=convert(xData[(i-1)*8+1..i*8],Array):
  xN[i]:=toBigInt(xs,r,p):
  # print(xN[i]);
end do:

# compute FFT256 over vector of big integers, of size 256
result:=FFT4096(xN, 4096, omega, p):

# convert result to big prime field elements, to vectors of size 8, write result to uData file
for i from 1 to vecSize do
  ys:=toBigFieldElement(result[i], r, p, 3):
  # print(ys);
  writedata(writeFile, ys, integer):
  writeline(writeFile, ""):
end do: