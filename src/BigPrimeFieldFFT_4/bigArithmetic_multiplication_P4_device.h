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


// for 16-word prime
#ifndef BIG_ARITHMETIC_MULTIPLICATION_DEVICE_H_
#define BIG_ARITHMETIC_MULTIPLICATION_DEVICE_H_

#include "../../include/BigPrimeFieldFFT_4/bigPrimeField_P4.h"

#include "bigArithmetic_addition_P4_device.h"
#include "bigArithmetic_subtraction_P4_device.h"

//__global__ void
//mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
//												 const usfixn64* __restrict__ ys, usfixn64* parameters,
//												 usfixn64* __restrict__ lVector,
//												 usfixn64* __restrict__ hVector,
//												 usfixn64* __restrict__ cVector,
//												 usfixn32* __restrict__ signVector);
//__global__ void
//mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
//												 usfixn64* __restrict__ lVector,
//												 usfixn64* __restrict__ hVector,
//												 usfixn64* __restrict__ cVector,
//												 const usfixn32* __restrict__ signVector);
//
//__global__ void
//mult_revised_8lhc_step3 (usfixn64* __restrict__ xs,
//												 const usfixn64* __restrict__ lVector,
//												 const usfixn64 * __restrict__ hVector,
//												 const usfixn64* __restrict__ cVector,
//												 const usfixn32* __restrict__ signVector,
//												 const usfixn64 * __restrict__ parameters);
//
//__global__ void
//kernel_one_shift_right (usfixn64 *__restrict__ xs,
//												const usfixn64 * __restrict__ parameters);

//[l,h,c]
//#include "bigArithmetic_addition.h"

//#include "bigArithmetic_subtraction.h"
//#include "bigArithmetic_cyclicShift.h"
//#include "bigArithmetic_fft.h"

/**********************************************/
//# subtraction modulo r
//subR:=proc(x,y,c0,r)
//s:=y;
//if x>=s then
//  s:=x-s;
//  c:=0;
//else:
//  s:=r+x-s;
//  c:=1;
//fi;
//if s>=r then
//  s:=s-r;
//  c:=1;
//fi;
//if c=1 then c:=r-1; fi;
//c:=c+c0;
//return s,c;
//end proc;


__device__ __inline__ void
device_p4_subR (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym)
{
//	usfixn64 um[8];
	short i, pos;
	unsigned short c = 0;
	usfixn64 num1;
	num1 = 0;
//	num2 = R - 1;

	for (i = 0; i < 3; i++)
	{
		num1 = ym[i] + c;
		if (xm[i] < num1) //there is not enough to do subtraction
		{
			c = 1;
			xm[i] = R - num1 + xm[i];
		}
		else
		{
			c = 0;
			xm[i] = xm[i] - num1;
		}
	}
//	if (c > 0)
//	{
//		pos = -1;
//		for (i = 0; i < 3; i++)
//		{
//			if (xm[i] <= R_MINUS_ONE)
//			{
//				pos = i;
//				break;
//			}
////			printf("pos = %d \n",pos);
//		}
//		if (pos >= 0)
//		{
//			for (i = 0; i < pos; i++)
//			{
//				xm[i] = 0;
//			}
//			xm[pos]++;
//		}
//		else
//		{
//			xm[0] = ULMAX;
//			for (i = 1; i < 3; i++)
//			{
//				xm[i] = 0;
//			}
//		}
//	}
//	xm[0] = um[0];
//	xm[1] = um[1];
//	xm[2] = um[2];
//	xm[3] = um[3];
//	xm[4] = um[4];
//	xm[5] = um[5];
//	xm[6] = um[6];
//	xm[7] = um[7];
}
/**********************************************/
//# addition modulo r
//addR:=proc(x,y,c0,r)
//  s:=irem(x+y,2^64);
//  c:=0;
//  if s<x or s<y then
//    s:=s+2^64-r;
//    c:=1;
//  fi;
//  if s>=r then:
//    print("here");
//    s:=s-r;
//    c:=1;
//  fi;
//  c:=c+c0;
//  c:=irem(c,r);
//  return s,c;
//end proc;
__device__ __inline__ void
device_p4_addR (usfixn64 *__restrict__ xm, usfixn64 *__restrict__ ym)
{
	unsigned short c = 0;
//	short pos1;
	usfixn64 num1, num2;

//	usfixn64 um[8];

	num1 = 0;
//	num2 = R - 1;
	short i;
//	usfixn64 tmp = 0;

	for (i = 0; i < 3; i++)
	{
		num1 = xm[i] + ym[i] + c;
		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
		{
//			xm[i] = num1 + RC;
			num1 = num1 + RC;
			c = 1;
		}
		else if (num1 >= R)
		{
			c = 1;
//			xm[i] = num1 - R;
			num1 = num1 - R;
		}
		else
		{
//			xm[i] = num1;
			c = 0;
		}
		xm[i] = num1;
	}
//
//	if (c > 0)
//	{
//		pos1 = -1;
//		for (i = 0; i < 8; i++)
//		{
//			if (xm[i] != 0)
//			{
//				pos1 = i;
//				break;
//			}
//		}
//		if (pos1 >= 0)
//		{
//			for (i = 0; i < pos1; i++)
//			{
//				xm[i] = num2;
//			}
//			xm[pos1]--;
//		}
//		else
//		{
//			xm[0] = ULMAX;
//			for (i = 1; i < 8; i++)
//			{
//				xm[i] = 0;
//			}
//		}
//	}
}
/**********************************************/
//add_128:=proc(x0,x1,y0,y1,c0)
//  s:=0;
//  c:=0;
//  s:=x0+y0;
//  s:=irem(s,2^64);
//  if s<x0 or s<y0 then
//    c:=1;
//  fi;
//  u0:=s;
//  u1:=y1+c;
//  s:=x1+u1;
//  s:=irem(s,2^64);
//  c:=0;
//  if s<x1 or s<u1 then
//    c:=1;
//  fi;
//  u1:=s;
//  u2:=0;
//  # print("c0",c0,"c",c,"u2",u2);
//  u2:=c0+c;
//  return u0,u1,u2;
//end proc;
__device__ __inline__ void
device_p4_add128 (usfixn64 x0, usfixn64 x1, usfixn64 y0, usfixn64 y1, usfixn64 c0,
				usfixn64 & u0, usfixn64 & u1, usfixn64 &u2)
{
	usfixn64 s = 0;
	usfixn64 c = 0;
	usfixn64 c1 = 0;

//	s = x0 + y0;
//	if (s < x0 || s < y0)
//	{
//		c = 1;
//	}

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.cc.u64 %1,0,0;\n\t"
			"}"
			:"=l"(s),"=l"(c):"l"(x0),"l"(y0));

	u0 = s;

//	s = y1 + c;
//	if (s < y1 || s < c)
//	{
//		c = 1;
//	}
//	else
//	{
//		c = 0;
//	}
	asm("{\n\t"
			"add.cc.u64 %0,%1,%2;\n\t"
			"addc.cc.u64 %1,0,0;\n\t"
			"}"
			:"=l"(s),"+l"(c):"l"(y1));

	u1 = s;

//	s = x1 + u1;
//
//	if (s < x1 || s < u1)
//	{
//		c++;
//	}

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.cc.u64 %1,%1,0;\n\t"
			"}"
			:"=l"(s),"+l"(c):"l"(x1),"l"(u1));

	u1 = s;
	u2 = 0;
	u2 = c0 + c;
}

/**********************************************/
__device__ __forceinline__ void
device_p4_sub128 (usfixn64 x0, usfixn64 x1, usfixn64 y0, usfixn64 y1, usfixn64 & u0,
				usfixn64 & u1)
{
//	s:=0;
//	c:=0;
//	if x0 >y0 then
//	s:=x0-y0;
//	else
//	s:=2^64+x0-y0;
//	c:=1;
//	fi;
//	u0:=s;
//	c:=y1+c;
//	if x1>c then
//	s:=x1-c;
//	else
//	s:=2^64+x1-c;
//	fi;
//	u1:=s;
//	if u1=2^64 then u1:=0; fi;
//	return u0, u1;
//}
	usfixn64 s = 0;
	usfixn64 c = 0;
	if (x0 >= y0)
		s = x0 - y0;
	else
	{
		s = x0 - y0;
		c = 1;
	}
	u0 = s;
	c = y1 + c;

	if (x1 > c)
		s = x1 - c;
	else
	{
		s = x1 - c;
		c = 1;
	}
	u1 = s;
}
/**********************************************/
__device__ __inline__ void //not tested, sub of 3 machine-words based on 64-bit arithmetic
device_p4_sub192_v0 (usfixn64 &x0, usfixn64 &x1, usfixn64 &x2, usfixn64 &y0, usfixn64 &y1,
					 usfixn64 &y2, usfixn64 & sign)
{
//	s:=0;
//	c:=0;
//	if x0 >y0 then
//	s:=x0-y0;
//	else
//	s:=2^64+x0-y0;
//	c:=1;
//	fi;
//	u0:=s;
//	c:=y1+c;
//	if x1>c then
//	s:=x1-c;
//	else
//	s:=2^64+x1-c;
//	fi;
//	u1:=s;
//	if u1=2^64 then u1:=0; fi;
//	return u0, u1;
//}

	usfixn64 u0 = 0, u1 = 0, u2 = 0;

	sign = 0;
	if (x2 < y2)
		sign = -1;
	else if (x2 > y2)
		sign = 1;

	if ((x1 < y1) && sign == 0)
		sign = -1;
	else if ((x1 > y1) && sign == 0)
		sign = 1;

	if ((x0 < y0) && sign == 0)
		sign = -1;
	else if ((x0 >= y0) && sign == 0)
		sign = 1;
	if (sign == 1)
		sign = 0;
	else
		sign = 1;

	usfixn64 s = 0, c = 0;

	if (sign == 0)
	{
		if (x0 >= y0)
		{
			s = x0 - y0;
		}
		else
		{
			s = x0 - y0; //s=2^64+x0-y0;
			c = 1;
		}

		u0 = s;
		c = y1 + c;
		if (x1 >= c)
		{
			s = x1 - c;
			c = 0;
		}
		else
		{
			s = x1 - c; //s=2^64+x1-c;
			c = 1;
		}

		u1 = s;
		c = y2 + c;
		if (x2 >= c)
		{
			s = x2 - c;
			c = 0;
		}
		else
		{
			s = x2 - c; //2^64+x2-c:
			c = 1;
		}

		u2 = s;
//		if (u1 == 0)
//			u2++;

	}
	else
	{
		if (y0 >= x0)
		{
			s = y0 - x0;
		}
		else
		{
			s = y0 - x0; //s=2^64+y0-x0;
			c = 1;
		}

		u0 = s;

		c = x1 + c;
		if (y1 >= c)
		{
			s = y1 - c;
			c = 0;
		}
		else
		{
			s = y1 - c; //s=2^64+y1-c;
			c = 1;
		}

		u1 = s;
		c = x2 + c;
		if (y2 >= c)
		{
			s = y2 - c;
			c = 0;
		}
		else
		{
			s = y2 - c; //2^64+y2-c:
			c = 1;
		}

		u2 = s;
//		if (u1 == 0)
//			u2++;
	}
	x0 = u0;
	x1 = u1;
	x2 = u2;
}
/**********************************************/
__device__ __inline__ void //not tested, sub of 3 machine-words based on 64-bit arithmetic
device_p4_sub192_v1 (usfixn64 &x0, usfixn64 &x1, usfixn64 &x2, usfixn64 &y0, usfixn64 &y1,
					 usfixn64 &y2, usfixn64 & sign)
{
//	s:=0;
//	c:=0;
//	if x0 >y0 then
//	s:=x0-y0;
//	else
//	s:=2^64+x0-y0;
//	c:=1;
//	fi;
//	u0:=s;
//	c:=y1+c;
//	if x1>c then
//	s:=x1-c;
//	else
//	s:=2^64+x1-c;
//	fi;
//	u1:=s;
//	if u1=2^64 then u1:=0; fi;
//	return u0, u1;
//}

	usfixn64 u0 = 0, u1 = 0, u2 = 0;

	sign = 0;
	if (x2 < y2)
		sign = -1;
	else if (x2 > y2)
		sign = 1;

	if ((x1 < y1) && sign == 0)
		sign = -1;
	else if ((x1 > y1) && sign == 0)
		sign = 1;

	if ((x0 < y0) && sign == 0)
		sign = -1;
	else if ((x0 >= y0) && sign == 0)
		sign = 1;
	if (sign == 1)
		sign = 0;
	else
		sign = 1;

	usfixn64 s = 0, c = 0;

	if (sign == 0)
	{
		if (x0 >= y0)
		{
			s = x0 - y0;
		}
		else
		{
			s = x0 - y0; //s=2^64+x0-y0;
			c = 1;
		}

		u0 = s;
		c = y1 + c;
		if (x1 >= c)
		{
			s = x1 - c;
			c = 0;
		}
		else
		{
			s = x1 - c; //s=2^64+x1-c;
			c = 1;
		}

		u1 = s;
		c = y2 + c;
		if (x2 >= c)
		{
			s = x2 - c;
			c = 0;
		}
		else
		{
			s = x2 - c; //2^64+x2-c:
			c = 1;
		}

		u2 = s;
//		if (u1 == 0)
//			u2++;

	}
	else
	{
		if (y0 >= x0)
		{
			s = y0 - x0;
		}
		else
		{
			s = y0 - x0; //s=2^64+y0-x0;
			c = 1;
		}

		u0 = s;

		c = x1 + c;
		if (y1 >= c)
		{
			s = y1 - c;
			c = 0;
		}
		else
		{
			s = y1 - c; //s=2^64+y1-c;
			c = 1;
		}

		u1 = s;
		c = x2 + c;
		if (y2 >= c)
		{
			s = y2 - c;
			c = 0;
		}
		else
		{
			s = y2 - c; //2^64+y2-c:
			c = 1;
		}

		u2 = s;
//		if (u1 == 0)
//			u2++;
	}
	x0 = u0;
	x1 = u1;
	x2 = u2;
}
/**********************************************/
__device__ __inline__ void //not tested, sub of 3 machine-words based on 64-bit arithmetic
device_p4_sub192 (usfixn64 &x0, usfixn64 &x1, usfixn64 &x2, usfixn64 &y0, usfixn64 &y1,
				usfixn64 &y2, usfixn64 & sign)
{
//	s:=0;
//	c:=0;
//	if x0 >y0 then
//	s:=x0-y0;
//	else
//	s:=2^64+x0-y0;
//	c:=1;
//	fi;
//	u0:=s;
//	c:=y1+c;
//	if x1>c then
//	s:=x1-c;
//	else
//	s:=2^64+x1-c;
//	fi;
//	u1:=s;
//	if u1=2^64 then u1:=0; fi;
//	return u0, u1;
//}

//	usfixn64 u0 = 0, u1 = 0, u2 = 0;

	usfixn64 s = 0, c = 0;

	sign = 0;
//	usfixn16 b0=0,b1=0,b2=0;
//	b2= (x2>y2);
//	b1= (x1>y1);
//	b0= (x0>y0);
//	b0 = b0+2*b1+4*b2;
//
//	usfixn16 b0 = (x0>y0)+((x1>y1)<<1)+((x2>y2)<<2);
	s = (x2 > y2) << 1;
	s += (x1 > y1);
	s <<= 1;
	s += (x0 > y0);
//	usfixn16 m0=0,m1=0,m2=0;
//	m2= (x2<y2);
//	m1= (x1<y1);
//	m0= (x0<y0);
//	m0 = m0+2*m1+4*m2;
//
//	usfixn16 m0 = (x0<y0)+((x1<y1)<<1)+((x2<y2)<<2);
	c = (x2 < y2) << 1;
	c += (x1 < y1);
	c <<= 1;
	c += (x0 < y0);

//	if (b0>=m0)
//		sign=0;
//	else
//		sign =1;
	sign = (s < c);
//	if(threadIdx.x==0 && blockIdx.x==0)
//	{
//		printf("sign-lex=%lu \n",sign);
//	}
//	sign=0;
//	if (x2 < y2)
//		sign = -1;
//	else if (x2 > y2)
//		sign = 1;
//
//	if ((x1 < y1) && sign == 0)
//		sign = -1;
//	else if ((x1 > y1) && sign == 0)
//		sign = 1;
//
//	if ((x0 < y0) && sign == 0)
//		sign = -1;
//	else if ((x0 >= y0) && sign == 0)
//		sign = 1;
//	if (sign == 1)
//		sign = 0;
//	else
//		sign = 1;
//
//	if(threadIdx.x==0 && blockIdx.x==0)
//		{
//			printf("sign-long=%lu \n ======================== \n",sign);
//		}

	if (sign == 0)
	{
		asm("{\n\t"
				"sub.cc.u64 %0,%2,%3;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(x0),"l"(y0)
		);
//		c=~c;
		c = c ^ 1;

//		if(threadIdx.x==0 && blockIdx.x==0)
//			{
//			printf("c_asm=%lu \n",c);
//			}
//
//
//		if (x0 >= y0)
//		{
//			s = x0 - y0;
//			c=0;
//		}
//		else
//		{
//			s = x0 - y0; //s=2^64+x0-y0;
//			c = 1;
//		}
//		if(threadIdx.x==0 && blockIdx.x==0)
//			{
//			printf("c_normal=%lu \n ====================== \n",c);
//			}
//		u0 = s;
		x0 = s;
		c = y1 + c;
//		if (x1 >= c)
//		{
//			s = x1 - c;
//			c = 0;
//		}
//		else
//		{
//			s = x1 - c; //s=2^64+x1-c;
//			c = 1;
//		}
//
		asm("{\n\t"
				"sub.cc.u64 %0,%2,%1;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(x1)
		);
		c = c ^ 1;

//		u1 = s;
		x1 = s;
		c = y2 + c;
//		if (x2 >= c)
//		{
//			s = x2 - c;
//			c = 0;
//		}
//		else
//		{
//			s = x2 - c; //2^64+x2-c:
//			c = 1;
//		}

		asm("{\n\t"
				"sub.cc.u64 %0,%2,%1;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(x2)
		);
		c = c ^ 1;

//		u2 = s;
		x2 = s;
//		if (u1 == 0)
//			u2++;

	}
	else
	{
//		if (y0 >= x0)
//		{
//			s = y0 - x0;
//		}
//		else
//		{
//			s = y0 - x0; //s=2^64+y0-x0;
//			c = 1;
//		}

		asm("{\n\t"
				"sub.cc.u64 %0,%2,%3;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(y0),"l"(x0)
		);
//		c=~c;
		c = c ^ 1;

//		u0 = s;
		x0 = s;

		c = x1 + c;
//		if (y1 >= c)
//		{
//			s = y1 - c;
//			c = 0;
//		}
//		else
//		{
//			s = y1 - c; //s=2^64+y1-c;
//			c = 1;
//		}
		asm("{\n\t"
				"sub.cc.u64 %0,%2,%1;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(y1)
		);
		c = c ^ 1;

//		u1 = s;
		x1 = s;
		c = x2 + c;
//		if (y2 >= c)
//		{
//			s = y2 - c;
//			c = 0;
//		}
//		else
//		{
//			s = y2 - c; //2^64+y2-c:
//			c = 1;
//		}
		asm("{\n\t"
				"sub.cc.u64 %0,%2,%1;\n\t"
				"addc.u64 %1,0,0;\n\t"
				"}"
				:"=l"(s),"+l"(c)
				:"l"(y2)
		);
		c = c ^ 1;

//		u2 = s;
		x2 = s;
//		if (u1 == 0)
//			u2++;
	}
//	x0 = u0;
//	x1 = u1;
//	x2 = u2;
}
/**********************************************/
__device__ __inline__ void
device_p4_divR (usfixn64 &xs0, usfixn64 &xs1, usfixn64 &q, usfixn64 &m)
{
	q = 0;
	m = 0;
	usfixn64 p63 = 0x7FFFFFFFFFFFFFFFL; //(1 << 63) - 1;
	usfixn64 r0 = 1;
	r0 = r0 << 34;
	short f = 1, f2 = 1;
	usfixn64 s0 = xs0, s1 = xs1;
	usfixn64 m0;
//	  s1:=iquo(s,2^64); #iquo(s,r-r0)
//	  s0:=irem(s,2^64); #irem(s,r-r0)

//	  # if s>=r then
	usfixn64 t0, t1, q0;

// round 1

	for (short iteration = 0; iteration < 3; iteration++)
	{
//		t1 = 2 * s1;
		t1 = s1 << 1;
//for this r, t1 will always be less 2^64
//	  t0 =iquo(s0,2^63);
		t0 = s0;
		t0 = t0 >> 63;
		q0 = t0 + t1;
//		printf ("round=%d, q0 =%d \n", iteration, q0);
		if (q0 > 0)
		{
//	    # s1:=iquo(s,2^64);
//	    # s0:=irem(s,2^64);
//	    m0:=irem(s0,2^63);
			m0 = (s0 & p63);
//	    # q0:=iquo(s,r-r0);
//	    # q0:=t1+t0;
//	    # q0:=iquo(s0+s1*2^64,r-r0);

//			q = q + q0 * f;
			if (f == 1)
				q += q0;
			else
				q -= q0;

//	    # if irem(s,(r-r0))<q0*r0 then
//	    s1:=iquo(q0*r0,2^64);
//	    s0:=irem(q0*r0,2^64);
			asm("{\n\t"
					"mul.lo.u64 %0,%2,%3;\n\t"
					"mul.hi.u64 %1,%2,%3;\n\t"
					"}"
					:"=l"(s0),"=l"(s1)
					:"l"(q0),"l"(r0)
			);
//						s0=q0;
//						s1=q0;
//						s0=s0 & (0x1FFFFFFF);
//						s0=s0<<34;
////						s1=s0>>29;

//			printf ("s0=%d, s1=%d \n", s0, s1);

			if (m0 < s0 || s1 > 0)
			{
//	      # s:=q0*r0-irem(s,(r-r0));
//	      # s:=q0*r0-m0;
				device_p4_sub128 (s0, s1, m0, 0, s0, s1);
				f2 = -1;
			}
			else
			{
//	      # s:=irem(s,(r-r0))-q0*r0;
//	      # s:=m0-q0*r0;
				device_p4_sub128 (m0, 0, s0, s1, s0, s1);
				f2 = 1;
//	    fi;
			}
			f = f * f2;
		}
	}

	if (f == -1 && s0 != 0)
	{
		s0 = R - s0;
		q = q - 1;
	}
	m = s0;
}
/**********************************************/
/**********************************************/
//
//device_p4_toLHC:=proc (s,r)
//  s0:=s[1]:
//  s1:=s[2];
//  s2:=s[3];
//  r0:=2^34;
//  q:=iquo(2^64,r-r0);
//  # p2:=s2*q^2;
//  p2:=4*s2;
//  # p1:=s1*q;
//  p1:=2*s1;
//  p0:=s0;
//  t0:=(2*s1*r0);
//  t1:=(4*s2*r0^2);
//  t2:=8*r0*s2;
//  if t0> 2^64 then print("t0 is huge"); fi;
//  if t1>2^64 then print("t1 is huge"); fi;
//  # t0q:=iquo(t0,r):
//  # t0r:=irem(t0,r);
//  d:=device_p4_divR(t0,r);
//  t0q:=d[1];
//  t0r:=d[2];
//
//  # t1q:=iquo(t1,r);
//  # t1r:=irem(t1,r);
//  d:=device_p4_divR(t1,r);
//  t1q:=d[1];
//  t1r:=d[2];
//
//  # p1:=p1-iquo(2*s1*r0,r)+iquo(4*s2*r0^2,r)-8*r0*s2;
//  # p0:=p0-rem(2*s1*r0,r)+rem(4*s2*r0^2,r);
//  # p1:=p1-iquo(t0,r)+iquo(t1,r)-t2;
//  # p0:=p0-rem(t0,r)+rem(t1,r);
//  # p1:=p1-t0q+t1q-t2;
//  # p1:=p1-t0q;
//
//  c:=0;
//  tmp:=0;
//  p1,tmp:=addR(p1,t1q,tmp,r);
//  p2:=p2+tmp;
//
//  tmp:=0;
//  p1,tmp:=subR(p1,t0q,tmp,r);
//  p2:=p2+tmp;
//
//  tmp:=0;
//  # p1:=p1-iquo(t0,r)-t2;
//  p1,tmp:=subR(p1,t2,tmp,r);
//  p2:=p2+tmp;
//  tmp:=0;
//  # print(p1,"after");
//  # p1:=2*s1-t0q-t2+t1q;
//  # print(p1-r,"before");
//  # p1:=p1-iquo(t0,r)-t2;
//  # p0:=p0-t0r+t1r;
//
//  tmp:=0;
//  p0,tmp:=addR(p0,t1r,tmp,r);
//  # p1,tmp:=addR(p1,tmp,0,r);
//  # p2:=p2+tmp;
//  p1:=p1+tmp;
//
//  tmp:=0;
//  p0,tmp:=subR(p0,t0r,tmp,r);
//  # p1,tmp:=addR(p1,tmp,0,r);
//  p1:=p1+tmp;
//
//  # p1:=p1+iquo(p0,r);
//  # p0:=irem(p0,r);
//  # p2:=p2+iquo(p1,r);
//  # p1:=irem(p1,r);
//  pN:=p0+p1*r+p2*(r^2);
//  if p1<0 then
//    p2:=p2-1;
//    p1:=r+p1;
//  fi;
//  if p0>=r then
//    p0:=p0-r;
//    p1:=p1+1;
//  fi;
//
//  if p1>=r then
//    p1:=p1-r;
//    p2:=p2+1;
//  fi;
//  if p2>=r then
//    p2:=p2-r;
//    # p2:=p2+1;
//  fi;
//
//  return p0,p1,p2;
//end proc;
__device__ __forceinline__ void
device_p4_toLHC_v0 (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
					usfixn64&p2)
{
//	usfixn64 p2, p1, p0;
	usfixn64 t0, t1, t2;
	usfixn64 t00, t01;
	usfixn64 t10, t11;
	usfixn64 l = 0, h = 0, c = 0, tmp = 0;
//	usfixn64 p0,p1,p2;
	p2 = 4 * s2;
	p1 = 2 * s1;
	p0 = s0;
	usfixn64 t0q = 0, t0r = 0, t1q = 0, t1r = 0;
	usfixn64 a1[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 a2[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 r1 = 0x8000000000000000;
	if (s1 < r1)
	{
		p1 = 2 * s1;
	}
	else
	{
		tmp = 2;
		asm("{\n\t"
				"mul.lo.u64 %0,%2,%3;\n\t"
				"mul.hi.u64 %1,%2,%3;\n\t"
				"}"
				:"=l"(t00),"=l"(t01)
				:"l"(s1),"l"(tmp)
		);
		device_p4_divR (t00, t01, t0q, t0r);

		p2 = p2 + t0q;
		p1 = t0r;
		t0q = 0;
		t0r = 0;
		t00 = 0;
		t01 = 0;
	}
	usfixn64 r0;
	r0 = 1;
	r0 <<= 34;
//	t0 = (2 * s1 * r0);
	usfixn64 r02 = 2 * r0;

	asm("{\n\t"
			"mul.lo.u64 %0,%2,%3;\n\t"
			"mul.hi.u64 %1,%2,%3;\n\t"
			"}"
			:"=l"(t00),"=l"(t01)
			:"l"(s1),"l"(r02)
	);

//	t1=(4*s2*r0^2); //to be corrected
//	t1 = 0;
	t1 = 4 * s2;
	t10 = 0;
	t11 = t1 << 4; //2^(68-64);

//	printf("tvalues, %lu, %lu, %lu, %lu\n",t00,t01,t10,t11); //correct up to this point

	t2 = 8 * r0 * s2;

//	device_p4_divR (t0, t0q, t0r);
//	device_p4_divR (t1, t1q, t1r);
	device_p4_divR (t00, t01, t0q, t0r);
	device_p4_divR (t10, t11, t1q, t1r);

//	printf("tqs , %lu, %lu, %lu, %lu\n",t0q,t0r,t1q,t1r); //correct up to this point

//	tmp = 0;
//	addR (p1, t1q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;

//	printf("ps before final %lu, %lu, %lu \n",p0,p1,p2);
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t1q;
	device_bigPrimeAdd_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	subR (p1, t0q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t0q;
	device_p4_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	subR (p1, t2, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t2;
	device_p4_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

	a1[1] = p0;
	a1[2] = p1;
	a1[3] = p2;
	a2[1] = t1r;

	device_bigPrimeAdd_plain (a1, a2);
	a2[1] = t0r;
	device_p4_bigSub_plain (a1, a2);
	p0 = a1[1];
	p1 = a1[2];
	p2 = a1[3];

//		return ;

	if (p1 < 0)
	{
		p2--;
		p1 = R + p1;
	}
	if (p0 >= R)
	{
		p0 = p0 - R;
		p1++;
	}

	if (p1 >= R)
	{
		p1 = p1 - R;
		p2 = p2 + 1;
	}

	if (p2 >= R)
	{
		p2 = p2 - R;
	}

}

/**********************************************/
__device__ __forceinline__ void
device_p4_toLHC_v1 (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
					usfixn64&p2)
{
//	usfixn64 p2, p1, p0;
	usfixn64 t0, t1, t2;
	usfixn64 t00, t01;
	usfixn64 t10, t11;
	usfixn64 l = 0, h = 0, c = 0, tmp = 0;
//	usfixn64 p0,p1,p2;
	p2 = 4 * s2;
	p1 = 2 * s1;
	p0 = s0;
	usfixn64 t0q = 0, t0r = 0, t1q = 0, t1r = 0;
	usfixn64 a1[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 a2[COEFFICIENT_SIZE] =
		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 r1 = 0x8000000000000000;
	if (s1 < r1)
	{
//		p1=2*s1;
		p1 = s1 << 1;
	}
	else
	{
		tmp = 2;
		asm("{\n\t"
				"mul.lo.u64 %0,%2,%3;\n\t"
				"mul.hi.u64 %1,%2,%3;\n\t"
				"}"
				:"=l"(t00),"=l"(t01)
				:"l"(s1),"l"(tmp)
		);

		device_p4_divR (t00, t01, t0q, t0r);

		p2 = p2 + t0q;
		p1 = t0r;
		t0q = 0;
		t0r = 0;
		t00 = 0;
		t01 = 0;
	}
	usfixn64 r0;
	r0 = 1;
	r0 <<= 34;
//	t0 = (2 * s1 * r0);
	usfixn64 r02 = 2 * r0;

//	asm("{\n\t"
//			"mul.lo.u64 %0,%2,%3;\n\t"
//			"mul.hi.u64 %1,%2,%3;\n\t"
//			"}"
//			:"=l"(t00),"=l"(t01)
//			:"l"(s1),"l"(r02)
//	);
	t00 = s1 & ((1 << 29) - 1);
	t00 = t00 << 35;
	t01 = s1;
	t01 >>= 29;

//	t1=(4*s2*r0^2); //to be corrected
//	t1 = 0;
//	t1 = 4 * s2;
	t10 = 0;
//	t11 = t1 << 4;		//2^(68-64);

	t11 = s2 << 6;
//	printf("tvalues, %lu, %lu, %lu, %lu\n",t00,t01,t10,t11); //correct up to this point

//	t2 = 8 * r0 * s2;
	t2 = s2 << 37;

//	device_p4_divR (t0, t0q, t0r);
//	device_p4_divR (t1, t1q, t1r);
	device_p4_divR (t00, t01, t0q, t0r);
	device_p4_divR (t10, t11, t1q, t1r);

//	printf("tqs , %lu, %lu, %lu, %lu\n",t0q,t0r,t1q,t1r); //correct up to this point

//	tmp = 0;
//	addR (p1, t1q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;

//	printf("ps before final %lu, %lu, %lu \n",p0,p1,p2);
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t1q;
	device_bigPrimeAdd_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	subR (p1, t0q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t0q;
	device_p4_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	subR (p1, t2, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t2;
	device_p4_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

	a1[1] = p0;
	a1[2] = p1;
	a1[3] = p2;
	a2[1] = t1r;

	device_bigPrimeAdd_plain (a1, a2);
	a2[1] = t0r;
	device_p4_bigSub_plain (a1, a2);
	p0 = a1[1];
	p1 = a1[2];
	p2 = a1[3];

//		return ;

	if (p1 < 0)
	{
		p2--;
		p1 = R + p1;
	}
	if (p0 >= R)
	{
		p0 = p0 - R;
		p1++;
	}

	if (p1 >= R)
	{
		p1 = p1 - R;
		p2 = p2 + 1;
	}

	if (p2 >= R)
	{
		p2 = p2 - R;
	}

}

/**********************************************/
__device__ __inline__ void
device_p4_toLHC_v2 (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
					usfixn64&p2)
{
//	s2:=iquo(s,2^128);
//		s1:=irem(s,2^128);
//		s0:=irem(s1,2^64);
//		s1:=iquo(s1,2^64);
//
//		# print("s210,",s2,s1,s0);

//		r0:=2^36;

	usfixn64 q;
//		q:=iquo(2^64,r-r0);
	q = 4;
//		Q2:=s2*q^2;
	usfixn64 Q2 = s2 << 4;

//		q1:=iquo(s1,r-r0);
//	usfixn64 q1 = s1 >> 62;
	usfixn64 q1 = s1;
	q1 >>= 62;

//		R1:=irem(s1,r-r0);
	usfixn64 hex_2_62 = 0x3FFFFFFFFFFFFFFF;
	usfixn64 R11 = s1 & hex_2_62;

//		Q1:=q*q1;
	usfixn64 Q1 = q1 << 2;

//		# print("Q2,Q1",Q2,Q1);
//		q0:=iquo(s0,r-r0); #shift
//		m0:=irem(s0,r-r0); #shift
	usfixn64 q0 = s0;
	q0 >>= 62;
	usfixn64 m0 = s0 & hex_2_62;
	usfixn64 l = 0, h = 0;
	usfixn64 c = Q2 + Q1;

//		##################
//		#computing Q_12 only by doing shifts
//			Q_12:=(Q1+Q2)*2^10;
	usfixn64 Q_12 = (Q1 + Q2) << 10;

//		# R_12:=irem((Q1+Q2)*r0^2,r);
//		# print("R_12",R_12);
//		R_12:=r-Q_12*r0;
//		Q_12:=Q_12-1;
	usfixn64 R_12 = R - (Q_12 << 36);
	Q_12--;

//		##################
//		h0:=Q_12+q0; #always less than a mword -> okay
	usfixn64 h0 = Q_12 + q0;
//		# R1:=2^62-1;
//		#######################
//		# I think we can specifically optimize this part

//		h1:=q*R1; #at most will take "mword-1"
//		for i from 1 to 4 do
//			if h1>=r then
//				h1:=h1-r;
//				c:=c+1;
//			fi;
//		end do;
	usfixn64 h1 = R11 << 2;

	for (short i = 0; i < 4; i++)
	{
		if (h1 >= R)
		{
			h1 = h1 - R;
			c++;
		}
	}

//		#######################
//		h2:=2*r0*(Q2+Q1); # always less than a mword -> okay
//		h3:=iquo(q*r0*R1,r); # two machine-words
//	usfixn64 h2 = (Q2 + Q1) << 37;
	usfixn64 h2 = (Q2 + Q1);
	h2 <<= 37;
//		# print("h3:", Log(h3,2));
//		# print("h3",h3);
//		# h3:=iquo(R1, 2^24); # R1*2^38 /2^62
//		h3:=iquo(q*r0*R1, r-r0); # R1*2^38 /2^62
//		# print("h3:",Log(h3,2));
//		# h3:=iquo(h3,2^62);
//		# print("h3:",Log(h3,2));
//		####################################
//		#device_p4_divR (q*r0*R1) ->h3:=q, h3r:=m
//		####################################
//		# h31:=iquo(R1*2^38,2^64);
//		# h31:=h31*4; #2^64 >> 62 = 2^2
//		# h30:=irem(R1*2^38,2^64);
//		# h30:=iquo(h30,2^62);
//		###################

//		h3q:=iquo(R1,2^24); #<<38 then >> 62 : >>24
//		h3:=h3q;

	usfixn64 h3q = R11 >> 24;
	usfixn64 h3 = h3q;

//		h3r:=irem(R1*2^38,2^64);
//		h3r:=irem(h3r,2^62);
	usfixn64 h3r = R11;
//	 h3r = R11 << 38;
	h3r <<= 38;
	h3r = h3r & hex_2_62;

//		R0:=h3r;
//		R0:=R0-h3q*r0;
	usfixn64 R00 = h3r;
	R00 >>= 36; // /r0
//		sign:=1;
//		if R0 <0 then
//			sign:=-1;
//			R0:=R0*(-1);
//		fi;
	short sign = 1;
	if (R00 < h3q)
	{
		sign = -1;
		h3q--;
	}
	else
	{
		sign = 1;
	}
//		if R0>=r then
//			print("negative");
//			h3q:=iquo(R0,r-r0);
//			h3r:=irem(R0,r-r0);
//			h3:=h3+sign*h3q-1;
//			h3r:=h3r-sign*r0*h3q;
//		fi;
//	h3q >>=36;
//	printf("R00 =%lu \n",R00);
	usfixn64 R_one_part = (1 << 26);
	if (h3q >= (R_one_part + 1))
	{
//		h3q = (R0 >> 62);

		h3q >>= 26;
//			h3r=(R0&hex_2_62);
		if (sign == 1)
			h3 = h3 + h3q - 1;
		else
			h3 = h3 - h3q - 1;
//			h3r=h3r-sign*h3q;
	}

//		printf("h3=%lu \n",h3);
	usfixn64 a1 = q * R11;

//	printf("r1 =%lu \n",R_one_part);
//			a_rem:=a1-(h3*r_r)-h3;
//	printf("a1=%lu\n",a1);
	usfixn64 a_rem = a1 - (h3 * (R_one_part + 1));
//	usfixn64 a_rem = (h3 * (R_one_part));
//	printf("arem = %lu \n",a_rem);
	h3r = a_rem << 36; //#+n_rem;
	if (h3r >= R)
	{
		h3r = h3r - R;
		h3++;
	}
//	printf("h3=%lu \t h3q=%lu \n", h3,h3r);
//		# print("h3",h3);
//		# print("h3:", Log(h3,2));
//		# h3r:=iquo(R1*2^38,2^64);
//		# h3r:=irem(h3r,r-r0);
//		# print("h3r",h3r);
//		# if h3r>=q*r0 then
//		# 	h3r:=h3r - q*r0;
//		# else
//		# 	h3:=h3-1;
//		# 	h3r:=r-q*r0+h3r;
//		# fi;
//
//		# h:=h0+h1+h2+h3;
//		h0:=h0+h1; # h0 <2*(r-r0) << 2^64 : always safe
//		h2:=h2+h3; # h2 <2*(r-r0) << 2^64 : always safe
//		if h0>2^64 then print("h0 is large"); fi;
//		if h2>2^64 then print("h2 is large"); fi;
//		if h0 >= h2 then
//			h:=h0-h2;
//		else
//			h:=r-h2+h0;
//			c:=c-1;
//		fi;

//	printf ("h0=%lu,\t h1=%lu,\t h2=%lu,\t h3=%lu \n", h0, h1, h2, h3);
	h0 = h0 + h1;
	h2 = h2 + h3;
	if (h0 >= h2)
	{
		h = h0 - h2;
	}
	else
	{
		h = R - h2 + h0;
		c--;
	}
//
//		# l:=irem(Q_12*r0^2,r)-q0*r0-irem(q*r0*R1,r)+m0;
//		# l:=R_12-q0*r0-irem(q*r0*R1,r)+m0;
//		l0:=R_12;
//		l1:=m0;
//		l0:= l0+l1;
//
	usfixn64 l0 = R_12;
	usfixn64 l1 = m0;

//		l2:=q0*r0;
//		l3:=h3r;
//		l2:=l2+l3;
//		if l0>=l2 then
//			l:=l0-l2;
//		else
//			l:=r-l2+l0;
//			h:=h-1;
//		fi;
//
//	usfixn64 l2 = (q0 << 36);
	usfixn64 l2 = (q0);
	l2 <<= 36;
	usfixn64 l3 = h3r;

//	printf ("l0=%lu,\t l1=%lu,\t l2=%lu,\t l3=%lu \n", l0, l1, l2, l3);
	l0 = l0 + l1;
	l2 += l3;
	if (l0 >= l2)
		l = l0 - l2;
	else
	{
		l = R - l2 + l0;
		h--;
	}

//		if l>=r then
//			l:=l-r;
//			h:=h+1;
//		fi;
//
//		if h>=r then
//			h:=h-r;
//			c:=c+1;
//		fi;
//
	if (l >= R)
	{
		l -= R;
		h++;
	}
	if (h >= R)
	{
		h -= R;
		c++;
	}

	p0 = l;
	p1 = h;
	p2 = c;
//	  return l,h,c;
}

/**********************************************/
__device__ __inline__ void
device_p4_toLHC (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
			 usfixn64&p2)
{
//	s2:=iquo(s,2^128);
//		s1:=irem(s,2^128);
//		s0:=irem(s1,2^64);
//		s1:=iquo(s1,2^64);
//
//		# print("s210,",s2,s1,s0);

//		r0:=2^36;

	usfixn64 q;
//		q:=iquo(2^64,r-r0);
	q = 4;
//		Q2:=s2*q^2;
	usfixn64 Q2 = s2 << 4;

//		q1:=iquo(s1,r-r0);
//	usfixn64 q1 = s1 >> 62;
	usfixn64 q1 = s1;
	q1 >>= 62;

//		R1:=irem(s1,r-r0);
	usfixn64 hex_2_62 = 0x3FFFFFFFFFFFFFFF;
	usfixn64 R11 = s1 & hex_2_62;

//		Q1:=q*q1;
	usfixn64 Q1 = q1 << 2;

//		# print("Q2,Q1",Q2,Q1);
//		q0:=iquo(s0,r-r0); #shift
//		m0:=irem(s0,r-r0); #shift
	usfixn64 q0 = s0;
	q0 >>= 62;
	usfixn64 m0 = s0 & hex_2_62;
	usfixn64 l = 0, h = 0;
	usfixn64 c = Q2 + Q1;

//		##################
//		#computing Q_12 only by doing shifts
//			Q_12:=(Q1+Q2)*2^10;
	usfixn64 Q_12 = (Q1 + Q2) << 10;

//		# R_12:=irem((Q1+Q2)*r0^2,r);
//		# print("R_12",R_12);
//		R_12:=r-Q_12*r0;
//		Q_12:=Q_12-1;
	usfixn64 R_12 = R - (Q_12 << 36);
	Q_12--;

//		##################
//		h0:=Q_12+q0; #always less than a mword -> okay
	usfixn64 h0 = Q_12 + q0;
//		# R1:=2^62-1;
//		#######################
//		# I think we can specifically optimize this part

//		h1:=q*R1; #at most will take "mword-1"
//		for i from 1 to 4 do
//			if h1>=r then
//				h1:=h1-r;
//				c:=c+1;
//			fi;
//		end do;
	usfixn64 h1 = R11 << 2;

//	for (short i = 0; i < 4; i++)
//	{
//		if (h1 >= R)
//		{
//			h1 = h1 - R;
//			c++;
//		}
//	}
	if (h1 >= R)
	{
		h1 = h1 - R;
		c++;
	}
	if (h1 >= R)
	{
		h1 = h1 - R;
		c++;
	}
	if (h1 >= R)
	{
		h1 = h1 - R;
		c++;
	}
	if (h1 >= R)
	{
		h1 = h1 - R;
		c++;
	}

//		#######################
//		h2:=2*r0*(Q2+Q1); # always less than a mword -> okay
//		h3:=iquo(q*r0*R1,r); # two machine-words
//	usfixn64 h2 = (Q2 + Q1) << 37;
	usfixn64 h2 = (Q2 + Q1);
	h2 <<= 37;
//		# print("h3:", Log(h3,2));
//		# print("h3",h3);
//		# h3:=iquo(R1, 2^24); # R1*2^38 /2^62
//		h3:=iquo(q*r0*R1, r-r0); # R1*2^38 /2^62
//		# print("h3:",Log(h3,2));
//		# h3:=iquo(h3,2^62);
//		# print("h3:",Log(h3,2));
//		####################################
//		#device_p4_divR (q*r0*R1) ->h3:=q, h3r:=m
//		####################################
//		# h31:=iquo(R1*2^38,2^64);
//		# h31:=h31*4; #2^64 >> 62 = 2^2
//		# h30:=irem(R1*2^38,2^64);
//		# h30:=iquo(h30,2^62);
//		###################

//		h3q:=iquo(R1,2^24); #<<38 then >> 62 : >>24
//		h3:=h3q;

	usfixn64 h3q = R11 >> 24;
	usfixn64 h3 = h3q;

//		h3r:=irem(R1*2^38,2^64);
//		h3r:=irem(h3r,2^62);
	usfixn64 h3r = R11;
//	 h3r = R11 << 38;
	h3r <<= 38;
	h3r = h3r & hex_2_62;

//		R0:=h3r;
//		R0:=R0-h3q*r0;

////////////////////////////////////
//	R0:=h3r;
//		R0:=R0-h3q*r0;
//		sign:=1;
//		if R0 <0 then
//			sign:=-1;
//			R0:=R0*(-1);
//			print("negative");
//		fi;
//		print("R00",R0);
//		# if R0>=2^64 then print("R0 is large"); fi;
//		if R0>=r then
//			h3q:=iquo(R0,r-r0); # div 2^62
//			# print(h3q);
//			# h3r:=irem(R0,r-r0);
//			h3:=h3+sign*h3q-1;
//			# print("h3r",h3r);
//			# h3r:=h3r-r0*h3q;
//			# h3r:=irem(q*r0*R1,r);
//		fi;
//	R00= h3r;
//	if (R))

	////////////////////////////////
	usfixn64 R00 = h3r;
	R00 >>= 36; // /r0
//		sign:=1;
//		if R0 <0 then
//			sign:=-1;
//			R0:=R0*(-1);
//		fi;
	short sign = 1;
	if (R00 < h3q)
	{
		sign = -1;
		h3q--;
//		h3q-=R00;
//		printf("h3q= %lu negative !!! \n", h3q);
	}
	else
	{
		sign = 1;
	}
//		if R0>=r then
//			print("negative");
//			h3q:=iquo(R0,r-r0);
//			h3r:=irem(R0,r-r0);
//			h3:=h3+sign*h3q-1;
//			h3r:=h3r-sign*r0*h3q;
//		fi;
//	h3q >>=36;
//	printf("R00 =%lu \n",R00);
	usfixn64 R_one_part = (1 << 26);
	if (h3q >= (R_one_part + 1))
	{
//		h3q = (R0 >> 62);

		h3q >>= 26;
//			h3r=(R0&hex_2_62);
		if (sign == 1)
			h3 = h3 + h3q - 1;
		else
			h3 = h3 - h3q - 1;
//			h3r=h3r-sign*h3q;
	}

//		printf("h3=%lu \n",h3);
//	usfixn64 a1 = q * R11;
	usfixn64 a1 = R11 << 2;
//	printf("r1 =%lu \n",R_one_part);
//			a_rem:=a1-(h3*r_r)-h3;
//	printf("a1=%lu\n",a1);
//	usfixn64 a_rem = a1 - (h3 * (R_one_part + 1));

	/**********************************************/
//	usfixn64 a_rem = a1 - h3;
//	a_rem -= (h3 << 26);
	//a_rem = (h3*(R_one_part+1));
	usfixn64 a_rem = (h3 << 26);
	a_rem += h3;

	sign = 1;
	if (a1 < a_rem)
	{
		sign = -1;
	}
	a_rem = a1 - a_rem;

	if (sign == -1)
	{
		a_rem += (R_one_part + 1);
		a_rem = a_rem & (0x3FFFFFF);
		h3--;
	}
	/**********************************************/

//	printf("arem = %lu \n",a_rem);
	h3r = a_rem << 36; //#+n_rem;

	if (h3r >= R) //unroll 3 times
	{
		h3r = h3r - R;
		h3++;
	}

//	printf("h3=%lu \t h3r=%lu \n", h3,h3r);
//		# print("h3",h3);
//		# print("h3:", Log(h3,2));
//		# h3r:=iquo(R1*2^38,2^64);
//		# h3r:=irem(h3r,r-r0);
//		# print("h3r",h3r);
//		# if h3r>=q*r0 then
//		# 	h3r:=h3r - q*r0;
//		# else
//		# 	h3:=h3-1;
//		# 	h3r:=r-q*r0+h3r;
//		# fi;
//
//		# h:=h0+h1+h2+h3;
//		h0:=h0+h1; # h0 <2*(r-r0) << 2^64 : always safe
//		h2:=h2+h3; # h2 <2*(r-r0) << 2^64 : always safe
//		if h0>2^64 then print("h0 is large"); fi;
//		if h2>2^64 then print("h2 is large"); fi;
//		if h0 >= h2 then
//			h:=h0-h2;
//		else
//			h:=r-h2+h0;
//			c:=c-1;
//		fi;

//	printf ("h0=%lu,\t h1=%lu,\t h2=%lu,\t h3=%lu \n", h0, h1, h2, h3);
	h0 = h0 + h1;
	h2 = h2 + h3;
	if (h0 >= h2)
	{
		h = h0 - h2;
	}
	else
	{
		h = R - h2 + h0;
		c--;
	}
//
//		# l:=irem(Q_12*r0^2,r)-q0*r0-irem(q*r0*R1,r)+m0;
//		# l:=R_12-q0*r0-irem(q*r0*R1,r)+m0;
//		l0:=R_12;
//		l1:=m0;
//		l0:= l0+l1;
//
	usfixn64 l0 = R_12;
	usfixn64 l1 = m0;

//		l2:=q0*r0;
//		l3:=h3r;
//		l2:=l2+l3;
//		if l0>=l2 then
//			l:=l0-l2;
//		else
//			l:=r-l2+l0;
//			h:=h-1;
//		fi;
//
//	usfixn64 l2 = (q0 << 36);
	usfixn64 l2 = (q0);
	l2 <<= 36;
	usfixn64 l3 = h3r;

//	printf ("l0=%lu,\t l1=%lu,\t l2=%lu,\t l3=%lu \n", l0, l1, l2, l3);
	l0 = l0 + l1;
	l2 += l3;
	//ptx here

	if (l0 >= l2)
		l = l0 - l2;
	else
	{
		l = R - l2 + l0;
		h--;
	}

//		if l>=r then
//			l:=l-r;
//			h:=h+1;
//		fi;
//
//		if h>=r then
//			h:=h-r;
//			c:=c+1;
//		fi;
//
	//round 1 of checking remainders
	if (l >= R)
	{
		l -= R;
		h++;
	}
	if (h >= R)
	{
		h -= R;
		c++;
	}
	//round 2 of checking remainders
	if (l >= R)
	{
		l -= R;
		h++;
	}
	if (h >= R)
	{
		h -= R;
		c++;
	}

	//round 3 of checking remainders
	if (l >= R)
	{
		l -= R;
		h++;
	}
	if (h >= R)
	{
		h -= R;
		c++;
	}

	p0 = l;
	p1 = h;
	p2 = c;
//	  return l,h,c;
}

/**********************************************/
#endif
