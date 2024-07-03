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


#restart;
with(ArrayTools):
r:=2^63+2^34;
p:=r^8+1:
us:=convert([2^131-8,2^131-7,2^131-6,2^131-5,2^131-4,2^131-3,2^131-2,2^131-1],Array):

calRQ := proc (num, p) 
  local r0, q0; 
  r0 := `mod`(num, p); 
  q0 := (num-r0)/p; 
  return [r0, q0]: 
end proc:

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
  return us;  
end proc:


bigModulo:=proc(us,r,p)
  local i,us9,rq,ri,qi,ss1,ss2,ss3;
  #DEBUG();
  for i from 1 to 7 do
    while us[i] >= 2^64 or us[i] < 0 do
      rq:=calRQ(us[i],2^63);
      ri:=op(1,rq):
      qi:=op(2,rq):
      us[i+1]:=us[i+1]+qi;
      us[i]:=ri-qi*2^34:      
    end do:
    if us[i]>=r then
      us[i]:=us[i]-r;
    end if:
  end do:
  us9:=0:
  while us[8] >= 2^64 or us[8] < 0 do
    #[ri,qi]:=calRQ(us[8]);
    rq:=calRQ(us[8],2^63);
    ri:=op(1,rq):
    qi:=op(2,rq):
    us9:=us9+qi;
    us[8]:=ri-qi*2^34:      
  end do:
  if us[8]>=r then
    us[8]:=us[8]-r;
  end if:

  if us9 <>0 then
    #us - us9. 
    #since us9 > 2^64, we need modulo us9
    #ss1:=new Array(1..8):
    ss1:= Array(1..8):
    ss1[1]:=us9:
    ss2:=bigModulo(ss1,r,p);
    ss3:=bigSub(us,ss2,r,p);
    return ss3;
  else
    return us;
  end if:
end proc:

ans:=convert(bigModulo(us,r,p),array);
fid:=fopen("result_bigModuleMaple.txt",WRITE);
writedata(fid,ans,integer);