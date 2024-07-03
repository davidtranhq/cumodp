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


with(ArrayTools):
#This function is to calculate FFT with size 16
FFT16:=proc(xx,w,p)
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

#This function is to calculate FFT with size 256. It calls FFT16.
with(ArrayTools):
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
    cs[j]:=FFT16(rs[j], (w^16) mod p, p); 
  end do:
  
  for k from 1 to K do
    for j from 1 to J do
      ds[k,j]:=(cs[j,k]*w^((j-1)*(k-1))) mod p:
    end do:
  end do:
  
  for k from 1 to K do
    es[k]:=FFT16(ds[k],(w^16) mod p, p);
  end do:  
  
  for k from 1 to K do 
    for j from 1 to J do
      b[(j-1)*K+k]:=es[k,j]:
    end do:
  end do:
  return b;
end proc:

#Test for FFT256 
b := FFT256(convert([seq(i, i = 1 .. 256)], Array), 256, 3, 257):
convert(b, list);



#This function is to calculate FFT with size 256. It calls FFT16.
with(ArrayTools):
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
    cs[j]:=FFT16(rs[j], (w^J) mod p, p); 
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
b := FFT4096(convert([seq(i, i = 1 .. 4096)], Array), 4096, 41, 12289):
convert(b, list);
