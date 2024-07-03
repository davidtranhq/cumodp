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
r:=2^63+2^34;
p:=r^8+1:
xs:=convert([2,3,4,0,0,1,2,0],Array);
ys:=convert([9223372054034644991,9223372054034644990,9223372054034644989,9223372054034644988,9223372054034644987,9223372054034644986,9223372054034644985,9223372054034644991],Array);


bigSub:=proc(xs,ys,r,p)
  local c,us,ui,i,j,nonzero,x,y,u1,u2,pos;
  c:=0:
  us:=Array(1..8):
  for i from 1 to 8 do
    if xs[i] < ys[i]+c then      
      ui := r-ys[i]+xs[i]-c;
      c := 1 ;
    else
      ui := xs[i]-ys[i]-c;
      c := 0;
    end if:
    us[i]:=ui;
  end do:
  
  if c=1 then
    pos:=0:
    for i from 1 to 8 do 
      if us[i] <> r-1 then
        pos:=i:
        break;
      end if:
    end do:
    if pos>0 then
      for j from 1 to pos-1 do
        us[j]:=0:
      end do:
      us[pos]:=us[pos]+1:
    else
      us[1]:=2^64-1;
      for i from 2 to 8 do
        us[i]:=0:
      end do:
    end if:
  end if:
  
  #verify
  x:=0:
  y:=0:
  u1:=0:
  for i from 1 to 8 do
    x := x + xs[i]*r^(i-1):
    y := y + ys[i]*r^(i-1):
    u1 := u1 + us[i]*r^(i-1):
  end do:
  u2:=(x-y) mod p;
  if u1-u2=0 then
    print("verify right");
  else
    print("verify wrong");
    print(cat("u1=",u1,",u2=",u2));
  end if:

  uss:=convert(us,array):
  f:=fopen("result_BigSubtractionMaple.txt",WRITE):
  writedata(f,uss):
  #print(us[]):
  #return us:
end proc:

bigSub(xs,ys,r,p);
