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


with(Ordinals);
r:=2^63+2^34;

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

#######################################
divR_v0:=proc(s,r)
	q:=iquo(s,r);
	m:=irem(s,r);
	return q,m;
end proc;
#######################################
divR_v1:=proc(xs,r)
	r0:=2^34;
	s:=xs;
	q:=0;
	m:=0;
	while abs(s)>=r do
		q0:=iquo(s,r-r0);
		q:=q+q0;
		s:=irem(s,(r-r0))-q0*r0;
	end do;
	m:=s;
	return q,m;
end proc;
#######################################
divR_v2:=proc(xs,r)
	r0:=2^34;
	s:=xs;
	q:=0;
	m:=0;
	#s can be 128 bit, heavy arithmetic?
	f:=1;
	while s>=r do
		# print(s);
		#q0 is a monomial, can be computed by shifting s 
		q0:=iquo(s,r-r0);
		q:=q+q0*f;
		# print("f",f);
		#q0*r0 is a monomial, can be computed by shifting q0 
		# if irem(s,(r-r0)) > q0*r0 then
		# # s:=irem(s,(r-r0))-q0*r0;
		# f:=f*(-1);
		# print("set to 1");
		# else
		# # s:=-irem(s,(r-r0))+q0*r0;
		# f:=-1;
		# fi;
		# s:=irem(s,(r-r0))-q0*r0;
		# s:=s*f;
		# if s<0 then f:=-1; else f:=1; fi;
		# s:=s*f;
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
	end do;
	# print(s,q,f);
	if f=-1 and s<>0 then
		print("f = -1");
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
	#s can be 128 bit, heavy arithmetic?
	f:=1;
	while s>=r do
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
	end do;
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
divR_128:=proc(xs0,xs1,r)
	r0:=2^34;
	# s:=xs;

		s1:=xs1;
		s0:=xs0;
	
	s:=s0+s1*2^64;
	q:=0;
	m:=0;
	#s can be 128 bit, heavy arithmetic?
	f:=1;
	while s>=r do
		# print(s);
		#q0 is a monomial, can be computed by shifting s 

		s1:=iquo(s,2^64);
		s0:=irem(s,2^64);
		# s1:=xs1;
		# s0:=xs0;
		# print(s0,s1);
		# q0:=iquo(s,r-r0);
		q0:=iquo(s0+s1*2^64,r-r0);
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
	end do;
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
# how to implement this function on 64 bit arithmetic?
#######################################
divR2:=proc(s,r)
	# q:=iquo(s,r^2);
	# m:=irem(s,r^2);
	# q0:=iquo(s,r);
	# m0:=irem(s,r);
	q0,m0:=divR(s,r);
	# q1:=iquo(q0,r);
	# m1:=irem(q0,r);
	q1,m1:=divR(q0,r);
	m:=m0+m1*r;
	q:=q1;
	# q:=0;
	# m:=s;
	return q,m;
end proc;
#######################################
add_128:=proc(x0,x1,y0,y1)
	s:=0;
	c:=0;
	s:=x0+y0;
	c:=iquo(s,2^64);
	u0:=irem(s,2^64);
	u1:=x1+y1+c;
	if u1=2^64 then u1:=0; fi;
	return u0,u1;
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
toLHC_v0:=proc (s,r)
	s0:=s[1]:
	s1:=s[2];
	s2:=s[3];
	sN:=s0+s1*(2^64)+s2*(2^128);
	c:=iquo(sN,r^2);
	sN:=irem(sN,r^2);
	h:=iquo(sN,r);
	sN:=irem(sN,r);
	l:=sN;
	print("first,",l,h,c);
	sN:=s0+s1*(2^64)+s2*(2^128);
	r0:=2^34;
	q:=iquo(2^64,r-r0);
	p2:=s2*q^2;
	p1:=s1*q;
	p0:=s0;
	t0:=s2*(q^2)*(r0^2);
	t1:=s2*(2)*(q^2)*(r*r0);
	t2:=s1*q*r0;
	#can avoid this step by having t0<r^2 and t1<r^2
	t0R2:=divR2(t0,r);
	t1R2:=divR2(t1,r);
	# p2:=p2+iquo(t0,r^2)-iquo(t1,r^2);
	# t0:=irem(t0,r^2);
	# t1:irem(t1,r^2);
	p2:=p2+t0R2[1]-t1R2[1];
	t0:=t0R2[2];
	t1:=t1R2[2];
	# tmp:=iquo(t0,r)-iquo(t1,r)-iquo(t2,r);
	print("t0");
	t0R:=divR(t0,r);
	print(t0R);
	print(divR_v0(t0,r));
	print("t1");
	t1R:=divR(t1,r);
	print(t1R);
	print(divR_v0(t1,r));
	print("t2");	
	t2R:=divR(t2,r);
	print(t2R);
	print(divR_v0(t2,r));
	tmp:=t0R[1]-t1R[1]-t2R[1];
	p1:=p1+tmp;
	# if abs(tmp) > 2^64 then print ("larg etmp"); fi;
	print("p1 is ",p1);
	p2:=p2+iquo(p1,r);
	p1:=irem(p1,r);
	if p1>=r then
		p1:=p1-r;
		p2:=p2+1;
	fi;
	if p1<0 then
		p1:=p1+r;
		p2:=p2-1;
	fi;
	# t0:=irem(t0,r);
	# t1:=irem(t1,r);
	# t2:=irem(t2,r);
	# tmp:=t0-t1-t2;
	tmp:=t0R[2]-t1R[2]-t2R[2];
	p0:=p0+tmp;
	if p0>=r then
		p0:=p0-r;
		p1:=p1+1;
	fi;
	if p0<0 then
		p0:=p0+r;
		p1:=p1-1;
	fi;
		# print(p0,p1,p2);
	pN:=p0+p1*r+p2*(r^2);
	print("second",p0,p1,p2);
	if p0=l and p1=h and p2=c then
		print("verified");
	fi;
	#	l:=irem(s0,r);
	#	t:=iquo(s0,r)+irem(s1,r^2);
	#	h:=irem(t,r^2);
	#	c:=iquo(t,r^2)+s2;
    # print("second-step",l,h,c);
	# l:=0;
	# h:=0;
	# c:=0;
	# return p0,p1,p2;
	# print("second,",l,h,c);
end proc;

#######################################
toLHC_v1:=proc (s,r)
	s0:=s[1]:
	s1:=s[2];
	s2:=s[3];
	sN:=s0+s1*(2^64)+s2*(2^128);
	c:=iquo(sN,r^2);
	sN:=irem(sN,r^2);
	h:=iquo(sN,r);
	sN:=irem(sN,r);
	l:=sN;
	print("first,",l,h,c);
	sN:=s0+s1*(2^64)+s2*(2^128);
	r0:=2^34;
	q:=iquo(2^64,r-r0);
	p2:=s2*q^2;
	p1:=s1*q;
	p0:=s0;
	t0:=s2*(q^2)*(r0^2);
	t1:=s2*(2)*(q^2)*(r*r0);
	t2:=s1*q*r0;
	#can avoid this step by having t0<r^2 and t1<r^2
	t0R2:=divR2(t0,r);
	t1R2:=divR2(t1,r);
	# p2:=p2+iquo(t0,r^2)-iquo(t1,r^2);
	# t0:=irem(t0,r^2);
	# t1:irem(t1,r^2);
	p2:=p2+t0R2[1]-t1R2[1];
	p20:=p2;
	p21:=0;
	p20,p21:= add_128(p20,p21,irem(t0R2[1],2^64),iquo(t0R2[1],2^64));
	p20,p21:= sub_128(p20,p21,irem(t1R2[1],2^64),iquo(t1R2[1],2^64));
	# p2:=p20+p21*2^64;
	t0:=t0R2[2];
	t1:=t1R2[2];
	# tmp:=iquo(t0,r)-iquo(t1,r)-iquo(t2,r);
	print("t0");
	t0R:=divR(t0,r);
	print(t0R);
	print(divR_v0(t0,r));
	print("t1");
	t1R:=divR(t1,r);
	print(t1R);
	print(divR_v0(t1,r));
	print("t2");	
	t2R:=divR(t2,r);
	print(t2R);
	print(divR_v0(t2,r));
	tmp:=t0R[1]-t1R[1]-t2R[1];
	p10:=irem(p1,2^64);
	p11:=iquo(p1,2^64);
	# p10:=0;
	# p11:=0;
	p10,p11:=add_128(p10,p11,irem(t0R[1],2^64),iquo(t0R[1],2^64));
	p10,p11:=sub_128(p10,p11,irem(t1R[1],2^64),iquo(t1R[1],2^64));
	p10,p11:=sub_128(p10,p11,irem(t2R[1],2^64),iquo(t2R[1],2^64));
	# print("i,j",irem(t2R[1],2^64),iquo(t2R[1],2^64));
	# print("p0,p1,tmp0,tmp1",p10,p11,irem(tmp,2^64),iquo(tmp,2^64));
	# sum0,sum1:=add_128(p10,p11,abs(irem(tmp,2^64)),iquo(tmp,2^64));
	# print(sum0,sum1);
	# tmp:=irem(p10+0*2^64,2^64);
	# print(tmp);
	p1:=p1+tmp;

	# p10,p11:=add_128(p10,p11,irem(tmp,2^64),iquo(tmp,2^64));
	# p1:=p10+p11*2^64;
	print(p1, p10+p11*2^64);
	# if abs(tmp) > 2^64 then print ("large tmp"); fi;
	print("p1 is ",p1);
	# if p1<0 then
	# 	p1:=(2^64-1)*2^64+2^64+p1;
	# fi;
	pR:=divR(p1,r);
	# print("t",pR[1],pR[2],iquo(p1,r),irem(p1,r));
	# p2:=p2+iquo(p1,r);
	p2:=p2+pR[1];
	# p1:=irem(p1,r);
	p1:=pR[2];
	if p1>=r then
		p1:=p1-r;
		p2:=p2+1;
	fi;
	if p1<0 then
		p1:=p1+r;
		p2:=p2-1;
	fi;
	# t0:=irem(t0,r);
	# t1:=irem(t1,r);
	# t2:=irem(t2,r);
	# tmp:=t0-t1-t2;
	tmp:=t0R[2]-t1R[2]-t2R[2];
	p0:=p0+tmp;
	if p0>=r then
		p0:=p0-r;
		p1:=p1+1;
	fi;
	if p0<0 then
		p0:=p0+r;
		p1:=p1-1;
	fi;
		# print(p0,p1,p2);
	pN:=p0+p1*r+p2*(r^2);
	print("second",p0,p1,p2);
	if p0=l and p1=h and p2=c then
		print("verified");
	fi;
	#	l:=irem(s0,r);
	#	t:=iquo(s0,r)+irem(s1,r^2);
	#	h:=irem(t,r^2);
	#	c:=iquo(t,r^2)+s2;
    # print("second-step",l,h,c);
	# l:=0;
	# h:=0;
	# c:=0;
	# return p0,p1,p2;
	# print("second,",l,h,c);
end proc;

#######################################
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
toLHC:=proc (s,r)
	s0:=s[1]:
	s1:=s[2];
	s2:=s[3];
	sN:=s0+s1*(2^64)+s2*(2^128);
	c:=iquo(sN,r^2);
	sN:=irem(sN,r^2);
	h:=iquo(sN,r);
	sN:=irem(sN,r);
	l:=sN;
	print("first,",l,h,c);
	sN:=s0+s1*(2^64)+s2*(2^128);
	r0:=2^34;
	q:=iquo(2^64,r-r0);
	# p2:=s2*q^2;
	p2:=4*s2;
	# p1:=s1*q;
	p1:=2*s1;
	p0:=s0;
	t0:=(2*s1*r0);
	t1:=(4*s2*r0^2);
	t2:=8*r0*s2;
	# if t0>2^64-1 then 
	# 	print("t0 is large");
	# fi;

	# if t1 > 2^64-1 then 
	# print("t1 is large");
	# fi;
	# if t2 > 2^64-1 then 
	# print("t2 is large");
	# fi;
	t0q:=iquo(t0,r):
	t0r:=irem(t0,r);
	t1q:=iquo(t1,r);
	t1r:=irem(t1,r);

	# p1:=p1-iquo(2*s1*r0,r)+iquo(4*s2*r0^2,r)-8*r0*s2;
	# p0:=p0-rem(2*s1*r0,r)+rem(4*s2*r0^2,r);
	# p1:=p1-iquo(t0,r)+iquo(t1,r)-t2;
	p0:=p0-rem(t0,r)+rem(t1,r);
	# p1:=p1-t0q+t1q-t2;
	# c:=0;
	tmp:=0;
	p1,tmp:=addR(p1,t1q,tmp,r);
	p2:=p2+tmp;
	tmp:=0;
	p1,tmp:=subR(p1,t0q,tmp,r);
	p2:=p2+tmp;
	tmp:=0;
	# p1:=p1-iquo(t0,r)-t2;
	p1,tmp:=subR(p1,t2,tmp,r);
	p2:=p2+tmp;
	tmp:=0;
	# print(p1,"after");
	# p1:=2*s1-t0q-t2+t1q;
	# print(p1-r,"before");
	# p1:=p1-iquo(t0,r)-t2;
	# p0:=p0-t0r+t1r;
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

	if p0=l then print("l equal"); fi;
	if p1=h then print("h equal"); fi;
	if p2=c then print("c equal"); fi;
	print("second",p0,p1,p2);
	if p0=l and p1=h and p2=c then
		print("verified");
	fi;
end proc;

#######################################
s:=Array(1..3,i->4*i+1);
s:=[r-2,1,7];
m0:=toLHC(s,r):

# x:=r-1;
# y:=r-3;
# c:=0;
# s1:=0;
# s0,s1:=subR(x,y,c,r);
# s2,s3:=addR(y,s0,s1,r);
# print(x,0);
# s:= (x-y) mod r:
# print(subR(x,y,c,r));
# x0:=1;
# x1:=1;
# y0:=1;
# y1:=1;
# z0:=1;
# z1:=1;
# s0,s1:=add_128(x0,x1,y0,y1);
# s0,s1:=sub_128(z0,z1,s0,s1);

# print(x0,x1);
# print(y0,y1);
# print(s0,s1);
# b0,b1:=sub_128(s0,s1,y0,y1);
# print(b0,b1);
# b2,b3:=sub_128(s0,s1,x0,x1);
# print(b2,b3);

# if x0=b0 and x1=b1 and y0=b2 and y1=b3 then
# 	print("sub and add verified");
# fi;

# x:=2^126+2^125+2^124;
# q,m:=divR_v0(x,r);
# q2,m2:=divR_128(irem(x,2^64),iquo(x,2^64),r);
# print(m,q);
# print(m2,q2);
# if m=m2 and q=q2 then print("div verified"); fi;

# division by r is a bug and do the principle from the paper,
# add comments and tests and specifications