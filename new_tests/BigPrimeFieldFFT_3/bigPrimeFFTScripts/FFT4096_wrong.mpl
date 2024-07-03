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


with(LinearAlgebra):
with(ArrayTools):

N:=4096:
decompose := proc (n) 
  local i, result, r, p, nn; 
  result := Array(1 .. 8); 
  r := 2^63+2^34; 
  p := r^8+1; 
  if p <= n or n < 0 then 
    print("the input number must be in [0,r^8)."); 
    return;  
  end if; 
  nn := n; 
  for i to 8 do 
    result[i] := `mod`(nn, r); 
    nn := (nn-result[i])/r: 
  end do; 
  return convert(result, list): 
end proc:

FFT16:=proc(xx,ww,p)
  #input x array  ww omega  
  #in-place
  local num1,num2,x,w,i;
  x:=Array(1..16):
  ArrayTools[Copy](xx,x):
  w:=Array(1..8):
  for i from 1 to 8 do
    w[i]:=ww^(i-1) mod p:
  end do:
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
  x[9]:=(num1+w[5]*num2) mod p:
  x[13]:=(num1-w[5]*num2) mod p:
  num1:=x[3]:
  num2:=x[7]:
  x[3]:=(num1+num2) mod p:
  x[7]:=(num1-num2) mod p:
  num1:=x[11]:
  num2:=x[15]:
  x[11]:=(num1+w[5]*num2) mod p:
  x[15]:=(num1-w[5]*num2) mod p:
  num1:=x[2]:
  num2:=x[6]:
  x[2]:=(num1+num2) mod p:
  x[6]:=(num1-num2) mod p:
  num1:=x[10]:
  num2:=x[14]:
  x[10]:=(num1+w[5]*num2) mod p:
  x[14]:=(num1-w[5]*num2) mod p:  
  num1:=x[4]:
  num2:=x[8]:
  x[4]:=(num1+num2) mod p:
  x[8]:=(num1-num2) mod p:
  num1:=x[12]:
  num2:=x[16]:
  x[12]:=(num1+w[5]*num2) mod p:
  x[16]:=(num1-w[5]*num2) mod p:
  
  #third round
  num1:=x[1]:
  num2:=x[3]:
  x[1]:=(num1+num2) mod p:
  x[3]:=(num1-num2) mod p:
  
  num1:=x[9]:
  num2:=x[11]:
  x[9]:=(num1+w[3]*num2) mod p:
  x[11]:=(num1-w[3]*num2) mod p:
  
  num1:=x[5]:
  num2:=x[7]:
  x[5]:=(num1+w[5]*num2) mod p:
  x[7]:=(num1-w[5]*num2) mod p:
  
  num1:=x[13]:
  num2:=x[15]:
  x[13]:=(num1+w[7]*num2) mod p:
  x[15]:=(num1-w[7]*num2) mod p:
  
  num1:=x[2]:
  num2:=x[4]:
  x[2]:=(num1+num2) mod p:
  x[4]:=(num1-num2) mod p:
  
  num1:=x[10]:
  num2:=x[12]:
  x[10]:=(num1+w[3]*num2) mod p:
  x[12]:=(num1-w[3]*num2) mod p:
  
  num1:=x[6]:
  num2:=x[8]:
  x[6]:=(num1+w[5]*num2) mod p:
  x[8]:=(num1-w[5]*num2) mod p:
  
  num1:=x[14]:
  num2:=x[16]:
  x[14]:=(num1+w[7]*num2) mod p:
  x[16]:=(num1-w[7]*num2) mod p:
  
  #fourth
  num1:=x[1]:
  num2:=x[2]:
  x[1]:=(num1+num2) mod p:
  x[2]:=(num1-num2) mod p:
  num1:=x[9]:
  num2:=x[10]:
  x[9]:=(num1+w[2]*num2) mod p:
  x[10]:=(num1-w[2]*num2) mod p:
  num1:=x[5]:
  num2:=x[6]:
  x[5]:=(num1+w[3]*num2) mod p:
  x[6]:=(num1-w[3]*num2) mod p:
  num1:=x[13]:
  num2:=x[14]:
  x[13]:=(num1+w[4]*num2) mod p:
  x[14]:=(num1-w[4]*num2) mod p:
  num1:=x[3]:
  num2:=x[4]:
  x[3]:=(num1+w[5]*num2) mod p:
  x[4]:=(num1-w[5]*num2) mod p:
  num1:=x[11]:
  num2:=x[12]:
  x[11]:=(num1+w[6]*num2) mod p:
  x[12]:=(num1-w[6]*num2) mod p:
  num1:=x[7]:
  num2:=x[8]:
  x[7]:=(num1+w[7]*num2) mod p:
  x[8]:=(num1-w[7]*num2) mod p:
  num1:=x[15]:
  num2:=x[16]:
  x[15]:=(num1+w[8]*num2) mod p:
  x[16]:=(num1-w[8]*num2) mod p:
  
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

uniTranspose:=proc(a,s,dnm,he,wi)
  #a matrix  s number of submatrix  dnm elements of submatrix
  #he height of submatrix  wi width of submatrix
  local i,x,x1,x2;
  x:=Array(1..dnm):
  for i from 1 to s do 
    ArrayTools[Copy](dnm,a,(i-1)*dnm,1,x,0,1):
    x1:=Matrix(he,wi,convert(x,list),order=C_order):
    x2:=Transpose(x1):
    x:=convert(convert(x2,list),Array):
    ArrayTools[Copy](dnm,x,0,1,a,(i-1)*dnm,1):
  end do:
end proc:


uniFFT16:=proc(a,s,w,p)
  #a matrix  s number of submatrix  dnm elements of submatrix
  #he height of submatrix  wi width of submatrix
  local i,x,y:
  x:=Array(1..16):
  for i from 1 to s do 
    ArrayTools[Copy](16,a,(i-1)*16,1,x,0,1):
    y:=FFT16(x,w,p);
    ArrayTools[Copy](16,y,0,1,a,(i-1)*16,1):
  end do:
end proc:

uniTwiddle:=proc(a,s,dnm,he,wi,w)
  local x,i,j,pos,nh,nw;
  x:=Array(1..dnm):
  for i from 1 to s do 
    ArrayTools[Copy](dnm,a,(i-1)*dnm,1,x,0,1):
    for j from 1 to dnm do
      pos:=j-1:
      nh:=floor(pos/wi):
      nw:=pos-nh*wi:
      x[j]:=x[j]*(w^(nh*nw) mod p) mod p;
    end do:
    ArrayTools[Copy](dnm,x,0,1,a,(i-1)*dnm,1):
  end do:
end proc:



r:=2^63+2^34:
p:=r^8+1:

#d1:=FileTools[Binary][Read](cat("/home/lychen/svn/cumodp/new_tests/BigPrimeFieldFFT/FFT_4096_input.dat"), integer[8], 4096*8, byteorder=native):
d1:=FileTools[Binary][Read](cat("C:/chenly/ScienceSource/uwosvn/cumodp/new_tests/BigPrimeFieldFFT/FFT_4096_input.dat"), integer[8], 4096*8, byteorder=native):
d1:=convert(d1,Array):
e1:=Array(1..4096):
for i from 1 to 4096 do
  pos:=(i-1)*8+1:
  e1[i]:=d1[pos]+d1[pos+1]*r^1+d1[pos+2]*r^2+d1[pos+3]*r^3+d1[pos+4]*r^4+d1[pos+5]*r^5+d1[pos+6]*r^6+d1[pos+7]*r^7:
end do:

uniTranspose(e1,1,4096,256,16);
uniTranspose(e1,16,256,16,16);
uniFFT16(e1,N/16,r,p);
w4096 := 36850024757647431031555179090680053993233791031683935463357793094299055666867852024663701514624880435798127524087306906649662755441422342494620897280641:
uniTwiddle(e1,16,256,16,16,w4096);
uniTranspose(e1,16,256,16,16);
uniFFT16(e1,N/16,r,p);
uniTranspose(e1,16,256,16,16);
uniTwiddle(e1,1,4096,16,256,w4096);
uniTranspose(e1,1,4096,16,256);
uniFFT16(e1,N/16,r,p);
uniTranspose(e1,1,4096,256,16);

