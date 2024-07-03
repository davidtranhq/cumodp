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
decompose := proc (n) 
  local i, result, r, p, nn; 
  result := Array(1 .. 8); 
  r := 9223372054034644992; 
  p := r^8+1; 
  if p <= n or n < 0 then 
    print("the input number must be in [0,r^8)."); 
    return;
  end if; 
  nn := n; for i to 8 do 
  result[i] := `mod`(nn, r); 
  nn := (nn-result[i])/r end do; 
  return convert(result, list): 
end proc:

FFTeval := proc (a::Vector, N::integer, W::integer, p::integer) 
  local v, aeven, aodd, veven, vodd, halfN, k, Wsqr, zk; 
  v := Vector(N); 
  if N = 1 then 
    v[1] := a[1]: 
  else 
    aeven := Vector((1/2)*N); 
    aodd := Vector((1/2)*N); 
    veven := Vector((1/2)*N); 
    vodd := Vector((1/2)*N); 
    halfN := (1/2)*N; 
    aeven := Vector([seq(a[2*k+1], k = 0 .. halfN-1)]); 
    aodd := Vector([seq(a[2*k+2], k = 0 .. halfN-1)]); 
    Wsqr := `mod`(W^2, p); 
    veven := FFTeval(aeven, halfN, Wsqr, p); 
    vodd := FFTeval(aodd, halfN, Wsqr, p); 
    zk := 1; 
    for k from 0 to halfN-1 do 
      v[k+1] := `mod`(veven[k+1]+`mod`(zk*vodd[k+1], p), p); 
      v[halfN+k+1] := `mod`(veven[k+1]-(`mod`(zk*vodd[k+1], p)), p); 
      zk := `mod`(zk*W, p): 
    end do: 
  end if; 
  return v: 
end proc:

forwardFFT := proc (f::Vector, w::integer, p::integer) 
  local N, Winv; N := LinearAlgebra:-Dimension(f); 
  Winv := w; 
  return `mod`(FFTeval(f, N, Winv, p), p): 
end proc:

w:=36850024757647431031555179090680053993233791031683935463357793094299055666867852024663701514624880435798127524087306906649662755441422342494620897280641:
N:=4096:
a := FileTools[Binary][Read](cat("/home/lychen/svn/cumodp/new_tests/BigPrimeFieldFFT/FFT_4096_input.dat"), integer[8], N*8, byteorder = native): 
b := Array(1 .. N):
r := 2^63+2^34;
p := r^8+1;
for i to N do 
  d1 := 0; 
  pos := (i-1)*8; 
  for j to 8 do 
    d1 := d1+a[pos+j]*r^(j-1): 
  end do; 
  b[i] := d1: 
end do:

c := forwardFFT(convert(b, Vector), w, p):


result:=Array(1..N*8):
pos:=1:
for i from 1 to N do 
  d2:=convert(decompose(c[i]),Array):
  for j from 1 to 8 do
    result[pos]:=d2[j]:
    pos:=pos+1:
  end do:
end do:

#result is an array to store the final result.