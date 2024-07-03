#ifndef BIG_ARITHMETIC_MULTIPLICATION_H_
#define BIG_ARITHMETIC_MULTIPLICATION_H_

#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"

#include "BigPrimeFieldFFT_3/bigArithmetic_addition_P3.h"
#include "BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.h"

/**********************************************/
__device__ void
device_p3_mulLong (usfixn64 x, usfixn64 y, usfixn64 *s)
{
	short x1, y1;
	usfixn64 l, h, c, x0, y0, v2, v5, v9, v10, v11, v14, v15, v16, v17, q, t;
	usfixn64 a0, a1, b0, b1, c0, c1, c1prime, d0, d1, d2, e0, e1;

	if (x <= SQRTR && y <= SQRTR)
	{
		s[0] = x * y;
		s[1] = 0;
		s[2] = 0;
		return;
	}

	x1 = (x >= R ? 1 : 0);
	x0 = (x1 > 0 ? x - R : x);
	y1 = (y >= R ? 1 : 0);
	y0 = (y1 > 0 ? y - R : y);

	v2 = x0 * y1; //[0,v2,0];
	v5 = x1 * y0; //[0,v5,0];
	v9 = x1 * y1; //[0,0,1];

	c = v9;
	l = 0;
	h = v5 + v2;
	h < v5 || h < v2 ? (c = c + 1) : (c = c);
	c > v9 ? (h = h + RC) : (h = h);

	if (x0 <= SQRTR && y0 <= SQRTR)
	{
		s[0] = x0 * y0;
		s[1] = h;
		s[2] = c;
		return;
	}

//lhc
//x0*y0
	a1 = x0 >> 32;
	a0 = x0 - (a1 << 32);
	b1 = y0 >> 32;
	b0 = y0 - (b1 << 32);

	c0 = 0;
	c1 = a1 * b1;

	t = a0 * b1;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	c0 += t;  //safe

	t = a1 * b0;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	q = c0 + t;               //here, is not related to r.
	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
	c0 = q;

	t = a0 * b0;
	q = c0 + t;
	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=x0*y0
	c0 = q;

	c1prime = c1 << 1;

	c0 >= R ? (v11 = 1) : (v11 = 0);
	v11 > 0 ? (v10 = c0 - R) : (v10 = c0);
//v12=0;

	q = l + v10;  //[l,h,c] + [v10,v11,0]
	q < l || q < v10 ? (v11 = v11 + 1) : (v11 = v11);
	q < l || q < v10 ? (l = q + RC) : (l = q);
	if (l >= R)
	{
		l = l - R;
		v11++;
	}
	q = h + v11;
	q < h || q < v11 ? (c = c + 1) : (c = c);
	q < h || q < v11 ? (h = q + RC) : (h = q);
	if (h >= R)
	{
		h = h - R;
		c++;
	}
//v13=0;
	c1prime >= R ? (v15 = 1) : (v15 = 0);
	v15 > 0 ? (v14 = c1prime - R) : (v14 = c1prime); //v13=0;

	q = h + v14;  //[l,h,c]+[0,v14,v15]
	q < h || q < v14 ? (c = c + v15 + 1) : (c = c + v15);
	q < h || q < v14 ? (h = q + RC) : (h = q);
	if (h >= R)
	{
		h = h - R;
		c++;
	}
//[l,h,c]

	d1 = c1prime >> 29;
	d0 = c1prime - (d1 << 29);
	if (d0 >= d1)
	{
		d2 = d0 - d1;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v16 = (e0 - e1) << 34) : (v16 = R - (e1 << 34) + (e0 << 34));
		e0 >= e1 ? (v17 = e1 + d1) : (v17 = e1 + d1 - 1);
		/*
		 if(e0>=e1)
		 {
		 v16=(e0-e1)<<34;
		 v17=e1+d1;
		 }
		 else
		 {
		 v17=e1+d1-1;
		 v16=R-(e1<<34)+(e0<<34);
		 }
		 */
	}
	else
	{
		//d1>d0
		d2 = d1 - d0;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v16 = R - ((e0 - e1) << 34)) : (v16 = (e1 - e0) << 34);
		e0 >= e1 ? (v17 = d1 - e1 - 1) : (v17 = d1 - e1);
		/*
		 if(e0>=e1)
		 {
		 v16=R-((e0-e1)<<34);
		 v17=d1-e1-1;
		 }
		 else
		 {
		 v16=(e1-e0)<<34;
		 v17=d1-e1;
		 }
		 */
	}
//[l,h,c]-[v16,v17,0]
//q
	q = 0;
	if (l >= v16)
	{
		l = l - v16;
	}
	else
	{
		l = R - v16 + l;
		q = 1;
	}
//t
	if (h < q + v17)
	{
		c = c - 1;
		h = R - q - v17 + h;
	}
	else
	{
		h = h - q - v17;
	}
	s[0] = l;
	s[1] = h;
	s[2] = c;
}

/**********************************************/
__device__ void
device_p3_mulLong_2_plain (usfixn64 const &x, usfixn64 const &y, usfixn64 &s0,
								 usfixn64 &s1, usfixn64 &s2)
{
	usfixn64 x1, y1;
	usfixn64 l, h, c, x0, y0, v2, v5, v9, v10, v11, v14, v15, v16, v17, q, t;
	usfixn64 a0, a1, b0, b1, c0, c1, c1prime, d0, d1, d2, e0, e1;

	if (x <= SQRTR && y <= SQRTR)
	{
		s0 = x * y;
		s1 = 0;
		s2 = 0;
		return;
	}

	x1 = (x >= R ? 1 : 0);
	x0 = (x1 > 0 ? x - R : x);
	y1 = (y >= R ? 1 : 0);
	y0 = (y1 > 0 ? y - R : y);

	v2 = x0 * y1; //[0,v2,0];
	v5 = x1 * y0; //[0,v5,0];
	v9 = x1 * y1; //[0,0,1];

	c = v9;
	l = 0;
	h = v5 + v2;
	h < v5 || h < v2 ? (c = c + 1) : (c = c);
	c > v9 ? (h = h + RC) : (h = h);

	if (x0 <= SQRTR && y0 <= SQRTR)
	{
		s0 = x0 * y0;
		s1 = h;
		s2 = c;
		return;
	}

//lhc
//x0*y0
	a1 = x0 >> 32;
	a0 = x0 - (a1 << 32);
	b1 = y0 >> 32;
	b0 = y0 - (b1 << 32);

	c0 = 0;
	c1 = a1 * b1;

	t = a0 * b1;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	c0 += t;  //safe

	t = a1 * b0;
	q = t >> 32;
	t = (t - (q << 32)) << 32;
	c1 += q;
	q = c0 + t;               //here, is not related to r.
	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
	c0 = q;

	t = a0 * b0;
	q = c0 + t;
	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=x0*y0
	c0 = q;

	c1prime = c1 << 1;

	c0 >= R ? (v11 = 1) : (v11 = 0);
	v11 > 0 ? (v10 = c0 - R) : (v10 = c0);
//v12=0;

	q = l + v10;  //[l,h,c] + [v10,v11,0]
	q < l || q < v10 ? (v11 = v11 + 1) : (v11 = v11);
	q < l || q < v10 ? (l = q + RC) : (l = q);
	if (l >= R)
	{
		l = l - R;
		v11++;
	}
	q = h + v11;
	q < h || q < v11 ? (c = c + 1) : (c = c);
	q < h || q < v11 ? (h = q + RC) : (h = q);
	if (h >= R)
	{
		h = h - R;
		c++;
	}
//v13=0;
	c1prime >= R ? (v15 = 1) : (v15 = 0);
	v15 > 0 ? (v14 = c1prime - R) : (v14 = c1prime); //v13=0;

	q = h + v14;  //[l,h,c]+[0,v14,v15]
	q < h || q < v14 ? (c = c + v15 + 1) : (c = c + v15);
	q < h || q < v14 ? (h = q + RC) : (h = q);
	if (h >= R)
	{
		h = h - R;
		c++;
	}
//[l,h,c]

	d1 = c1prime >> 29;
	d0 = c1prime - (d1 << 29);
	if (d0 >= d1)
	{
		d2 = d0 - d1;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v16 = (e0 - e1) << 34) : (v16 = R - (e1 << 34) + (e0 << 34));
		e0 >= e1 ? (v17 = e1 + d1) : (v17 = e1 + d1 - 1);
		/*
		 if(e0>=e1)
		 {
		 v16=(e0-e1)<<34;
		 v17=e1+d1;
		 }
		 else
		 {
		 v17=e1+d1-1;
		 v16=R-(e1<<34)+(e0<<34);
		 }
		 */
	}
	else
	{
		//d1>d0
		d2 = d1 - d0;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v16 = R - ((e0 - e1) << 34)) : (v16 = (e1 - e0) << 34);
		e0 >= e1 ? (v17 = d1 - e1 - 1) : (v17 = d1 - e1);
		/*
		 if(e0>=e1)
		 {
		 v16=R-((e0-e1)<<34);
		 v17=d1-e1-1;
		 }
		 else
		 {
		 v16=(e1-e0)<<34;
		 v17=d1-e1;
		 }
		 */
	}
//[l,h,c]-[v16,v17,0]
//q
	q = 0;
	if (l >= v16)
	{
		l = l - v16;
	}
	else
	{
		l = R - v16 + l;
		q = 1;
	}
//t
	if (h < q + v17)
	{
		c = c - 1;
		h = R - q - v17 + h;
	}
	else
	{
		h = h - q - v17;
	}
	s0 = l;
	s1 = h;
	s2 = c;
}

/**********************************************/
__device__ __inline__ void
device_p3_mult_32bit_with_carry (usfixn32 & __restrict__ x,
															usfixn32 & __restrict__ y,
															usfixn64 & __restrict__ u,
															usfixn64 & __restrict__ c)
{

	usfixn32 carry;
	usfixn32 mult;

	//l=(x*y).low
	//h=(x*y).high
	asm("mul.lo.u32 %0, %1, %2;":"=r"(mult):"r"(x),"r"(y));
	asm("mul.hi.u32 %0, %1, %2;":"=r"(carry):"r"(x),"r"(y));

	u = mult;
	c = carry;

//	c = 0;
//	usfixn32 a0, a1, b0, b1;
//	usfixn32 c0, c1, c2, c3;
//	usfixn32 s;
//
//	a0 = x & 0xFFFF;
//	a1 = x >> 16;
//	b0 = y & 0xFFFF;
//	b1 = y >> 16;
//
//	c0 = a0 * b0;
//	c1 = a0 * b1;
//	c2 = a1 * b0;
//	c3 = a1 * b1;
//
//	c1 <<= 16;
//	c2 <<= 16;
//	s = c0 + c1;
//	if (s < c0 || s < c1)
//		c3++;
//	c0 = s;
//	s = c0 + c2;
//	if (s < c0 || s < c2)
//		c3++;
//	c1 = a0 * b1;
//	c2 = a1 * b0;
//	c1 >>= 16;
//	c2 >>= 16;
//	c3 += (c1 + c2);
//	u = s;
//	c = c3;

}

/**********************************************/
__device__ __inline__ void
device_p3_mult_32bit_no_carry (usfixn32 & __restrict__ x,
														usfixn32 & __restrict__ y, usfixn64 &__restrict__ u) //u should be sfixn
{
	usfixn32 carry;
	usfixn32 mult;

	//l=(x*y).low
	//h=(x*y).high
	asm("mul.lo.u32 %0, %1, %2;":"=r"(mult):"r"(x),"r"(y));
	asm("mul.hi.u32 %0, %1, %2;":"=r"(carry):"r"(x),"r"(y));
	u = carry;
	u <<= 32;
//	u = (carry<<32);
//	u+=mult;
	u |= (mult);
//	c = carry;

////	c=0;
//	usfixn32 a0, a1, b0, b1;
//	usfixn32 c0, c1, c2, c3;
//	usfixn32 s;
//
//	a0 = x & 0xFFFF;
//	a1 = x >> 16;
//	b0 = y & 0xFFFF;
//	b1 = y >> 16;
//
//	c0 = a0 * b0;
//	c1 = a0 * b1;
//	c2 = a1 * b0;
//	c3 = a1 * b1;
//
////	asm( "shl.u32 %0,%1,0x100;":"=r(c1)":"r"(c1));
//	c1 <<= 16;
//	c2 <<= 16;
//
////	s = c0 + c1;
//		asm("add.u32 %0,%1,%2;":"=r"(s) : "r" (c0), "r" (c1));
//
//	if (s < c0 || s < c1)
//		c3++;
//	c0 = s;
//	s = c0 + c2;
//	if (s < c0 || s < c2)
//		c3++;
//	c1 = a0 * b1;
//	c2 = a1 * b0;
//	c1 >>= 16;
//	c2 >>= 16;
//	c3 += (c1 + c2);
//	//c3 can't handle 32-bits of shift, first laod u=c3; then do the shift
//	u = (c3);
//	u = (u << 32) + s;
////	c=c3;

}

/**********************************************/
__device__ __inline__ void
device_p3_mult_32bit_no_carry_0_working (usfixn32 & __restrict__ x,
																			usfixn32 & __restrict__ y,
																			usfixn64 &__restrict__ u) //u should be sfixn
{
//	c=0;
	usfixn32 a0, a1, b0, b1;
	usfixn32 c0, c1, c2, c3;
	usfixn32 s;
	a0 = x & 0xFFFF;
	a1 = x >> 16;
	b0 = y & 0xFFFF;
	b1 = y >> 16;

	c0 = a0 * b0;
	c1 = a0 * b1;
	c2 = a1 * b0;
	c3 = a1 * b1;

	c1 <<= 16;
	c2 <<= 16;
	s = c0 + c1;
	if (s < c0 || s < c1)
		c3++;
	c0 = s;
	s = c0 + c2;
	if (s < c0 || s < c2)
		c3++;
	c1 = a0 * b1;
	c2 = a1 * b0;
	c1 >>= 16;
	c2 >>= 16;
	c3 += (c1 + c2);
	//c3 can't handle 32-bits of shift, first laod u=c3; then do the shift
	u = (c3);
	u = (u << 32) + s;
//	c=c3;
}

/**********************************************/
__device__ __inline__ void
device_p3_mult_32bit_with_carry_working (usfixn32 & __restrict__ x,
																			usfixn32 & __restrict__ y,
																			usfixn64 & __restrict__ u,
																			usfixn64 & __restrict__ c)
{
	c = 0;
	usfixn32 a0, a1, b0, b1;
	usfixn32 c0, c1, c2, c3;
	usfixn32 s;

	a0 = x & 0xFFFF;
	a1 = x >> 16;
	b0 = y & 0xFFFF;
	b1 = y >> 16;

	c0 = a0 * b0;
	c1 = a0 * b1;
	c2 = a1 * b0;
	c3 = a1 * b1;

	c1 <<= 16;
	c2 <<= 16;
	s = c0 + c1;
	if (s < c0 || s < c1)
		c3++;
	c0 = s;
	s = c0 + c2;
	if (s < c0 || s < c2)
		c3++;
	c1 = a0 * b1;
	c2 = a1 * b0;
	c1 >>= 16;
	c2 >>= 16;
	c3 += (c1 + c2);
	u = s;
	c = c3;

}

/**********************************************/
__device__ __inline__ void
device_p3_mult_64bit_with_carry (usfixn64 & __restrict__ x,
															usfixn64 & __restrict__ y,
															usfixn64 & __restrict__ u,
															usfixn64 & __restrict__ c)
{

	usfixn64 carry;
	usfixn64 mult;

	//l=(x*y).low
	//h=(x*y).high
	asm("mul.lo.u64 %0, %1, %2;":"=l"(mult):"l"(x),"l"(y));
	asm("mul.hi.u64 %0, %1, %2;":"=l"(carry):"l"(x),"l"(y));

	u = mult;
	c = carry;
}

/**********************************************/
//__device__ __inline__ void mulLong_2_revised_1
__device__ __inline__ void
device_p3_mulLong_2_plain_2 (usfixn64 const &x, usfixn64 const &y, usfixn64 &s0_out,
									 usfixn64 &s1_out, usfixn64 &s2_out)
{
	usfixn64 u0, u1, u2, u3;
//	usfixn64 r0 = R;
//	usfixn64 r1 = R;
//	usfixn64 rc0 = RC;
//	usfixn64 rc1 = RC;

//	r0 = r0 & (0xFFFFFFFF);
//	r1 >>= (32);
//	r0=R0;
//	r1=R1;

//	rc0 = rc0 & (0xFFFFFFFF);
//	rc1 >>= (32);

//		rc0=RC0;
//		rc1=RC1;

	usfixn32 x0 = (x & 0xFFFFFFFF);
	usfixn32 y0 = (y & 0xFFFFFFFF);

	bool y2 = (y >= R);
	bool x2 = (x >= R);
//	usfixn64 x1 = (x - x2 * R);
//	x1 = (x2 * R);
//	x1 = x -  x1;
//	x1 >>= 32;
//	usfixn32 x1 = (x - x2 * R) >> 32;
	usfixn32 x1;
	usfixn32 y1;

	x1 = (x >> 32);
	y1 = (y >> 32);
	if (x2)
		x1 -= R1;
	//	y1 = (y - y2 * R) >> 32;

	if (y2)
		y1 -= R1;

//	x1=(x-x0-r*x2)>>32

//	y0 = y % (1 << 32)
//	y2 = int(y >= r)
//	y1 = (y - y0 - r * y2) >> 32

//	usfixn64 y1 = (y - y0 - y2 * R);
//	y1 = (y2 * R);
//	y1 = y - y0 - y1;
//	y1 >>= 32;

	usfixn64 s;
	usfixn64 c;
//	usfixn64 t=x1;
//	usfixn64 t1;
//	t<<=32;
////	s=(x0 + (x1 << 32));
//	t= x0 + t;
////	* (y0 + (y1 << 32));
//	t1=y1;
//	t1<<=32;
//	t1+=y0;
//	s=t1*t;
//	c=0;
//	if(s<t1||s<t)
//		c=1;

	/////////////////////////////////////////////// problem is here, in mutliplication at the beginnging,
//	usfixn64 a0, a1, a2, a3;
//	c = 0;
//	a0 = (x0 * y0) & (0xFFFFFFFF);
//	a1 = (x0 * y1) & (0xFFFFFFFF);
//	a2 = (x1 * y0) & (0xFFFFFFFF);
//	a3 = (x1 * y1) & (0xFFFFFFFF);
//
//	a0 = (a0) & (0xFFFFFFFF);
//	a1 = (a1) & (0xFFFFFFFF);
//	a2 = (a2) & (0xFFFFFFFF);
//	a3 = (a3) & (0xFFFFFFFF);
//
//	a1 <<= 32;
//	s = a0 + (a1);
//	if (s < a0 || s < a1)
//		c++;
//	a0 = s;
//	a2 <<= 32;
//	s = a0 + a2;
//	if (s < a0 || s < a2)
//		c++;
//
//	a1 = (x0 * y1);
//	a2 = (x1 * y0);
//	a1 = (a1) & (0xFFFFFFFF);
//	a2 = (a2) & (0xFFFFFFFF);
//	a1 >>= 32;
//	a2 >>= 32;
//	c = c + a3 + a1 + a2;

//	usfixn64 d0, d1, e0, e1;
	usfixn64 v0, v1;
//	v2;

//	s1_out=x0;
//	usfixn64 X, Y;
	v0 = x1;
	v1 = y1;
	v0 <<= 32;
	v1 <<= 32;
	v0 = v0 + x0;
	v1 = v1 + y0;

	usfixn64 c2 = 0;

//	device_p3_mult_64bit_with_carry(X, Y, s, c);
//	device_p3_mult_64bit_with_carry(X, Y, s, c2);
//	device_p3_mult_64bit_with_carry(v0, v1, s, c2);
//	device_p3_mult_64bit_with_carry(v0, v1, u0, c2);

	asm("mul.lo.u64 %0, %1, %2;":"=l"(u0):"l"(v0),"l"(v1));
	asm("mul.hi.u64 %0, %1, %2;":"=l"(c2):"l"(v0),"l"(v1));

	v0 = 0;
	v1 = 0;
//	v2 = 0;

//	v0, v1;
//	c = c2 * 2;
	c = c2 << 1;
//	d1 = c >> 29;
////	d0 = c - (d1 << 29);
//	d0 = (d1 << 29);
//	d0 = c - (d0);
//
//	if (d0 >= d1)
//	{
//		v2 = d0 - d1;
//		e1 = v2 >> 29;
//		e0 = v2 - (e1 << 29);
//		e0 >= e1 ? (v0 = (e0 - e1) << 34) : (v0 = R - (e1 << 34) + (e0 << 34));
//		e0 >= e1 ? (v1 = e1 + d1) : (v1 = e1 + d1 - 1);
//		/*
//		 if(e0>=e1)
//		 {
//		 v16=(e0-e1)<<34;
//		 v17=e1+d1;
//		 }
//		 else
//		 {
//		 v17=e1+d1-1;
//		 v16=R-(e1<<34)+(e0<<34);
//		 }
//		 */
//	}
//	else
//	{
//		//d1>d0
//		v2 = d1 - d0;
//		e1 = v2 >> 29;
//		e0 = v2 - (e1 << 29);
//		e0 >= e1 ? (v0 = R - ((e0 - e1) << 34)) : (v0 = (e1 - e0) << 34);
//		e0 >= e1 ? (v1 = d1 - e1 - 1) : (v1 = d1 - e1);
//	}

//	usfixn64 u0_aux = 0;
	u1 = 0;
	u2 = 0;
	u3 = 0;
	u1 = c >> 29;
	usfixn32 u1_p = u1 & (0x1FFFFFFF);
//	u0_aux = c - (u1 << 29);
//	u0_aux = (u1 << 29);
//	u0_aux = c - (u0_aux);
//	u0_aux = (c & 0x1FFFFFFF);  //c mod 2^29 = c & (2^29-1)

//
//	if (u0_aux >= u1)
//	{
//		v2 = u0_aux - u1;
//		u3 = v2 >> 29;
//		u2 = v2 - (u3 << 29);
//		u2 >= u3 ? (v0 = (u2 - u3) << 34) : (v0 = R - (u3 << 34) + (u2 << 34));
//		u2 >= u3 ? (v1 = u3 + u1) : (v1 = u3 + u1 - 1);
//		/*
//		 if(u2>=u3)
//		 {
//		 v16=(u2-u3)<<34;
//		 v17=u3+u1;
//		 }
//		 else
//		 {
//		 v17=u3+u1-1;
//		 v16=R-(u3<<34)+(u2<<34);
//		 }
//		 */
//	}
//	else
//	{
//		//u1>u0_aux
//		v2 = u1 - u0_aux;
//		u3 = v2 >> 29;
//		u2 = v2 - (u3 << 29);
//		u2 >= u3 ? (v0 = R - ((u2 - u3) << 34)) : (v0 = (u3 - u2) << 34);
//		u2 >= u3 ? (v1 = u1 - u3 - 1) : (v1 = u1 - u3);
//	}

	/******************************************************/
	usfixn32 r0, r1, r2, r3;
//	usfixn64 r = R;
//	usfixn64 d0;
//	usfixn64 d1;

	r0 = c & (0xFFFFFFFF);
	r2 = c >> 32;

	r1 = r0 >> 29;
	r3 = r2 >> 29;
	r0 = r0 & (0x1FFFFFFF);
	r2 = r2 & (0x1FFFFFFF);
//	usfixn64 e0;
//	usfixn32 e1;

//	r0 = c&(0xFFFFFFFF);
//	r2 = c>>32;
//
//	r1 = r0 >>29;
//	r0 = r0 & (0x1FFFFFFF);

//	r3 = r2 >>29;
//	r2 = r2 & (0x1FFFFFFF);
//	d0 = r0;
//	d1 = (c - d0) >> 29;
//	d1 = u1;
//	d0 = u0_aux;
//
////# d0 always less than d1
//	if (r3 != 0)
//	{
//
//		if (r0 < (r1 + 8 * r2))
//		{
////# print(d0,d1)
////# e1=(d1-d0)>>29
//			e0 = (r1 + 8 * r2 - r0);
////# print(e0,e1)
//			e1 = (e0 >> 29) + 8 * r3;
////# e0 = e0&(0xFFFFFFFF)
////# print(e1<<29 < (1<<63))
////# e0=(d1-d0-(e1<<29))
//			e0 = (d1 - r0 - (e1 << 29)) & (0x1FFFFFFF);
////# e0=(d1-r0-(e1<<29))
////# print(((e0-e1)<<34))
//			if (e0 >= e1)
//			{
////# v0=(r-((e0-e1)<<34))
////				v0 = (((r >> 32) - (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1 - (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v1 = (d1 - e1 - 1);
//			}
//			else
//			{
//				v0 = ((e1 - e0) << 2) << 32;
//				v1 = (d1 - e1);
//			}
//		}
//		else
//		{
////# print(d0,d1)
////# e1=(d1-d0)>>29
////# print((0xFFFFFFFF-r0+1))
////# print(r1+8*r2)
//			e0 = (r1 + 8 * r2) + (0xFFFFFFFF - r0 + 1);
////# e0=e0&(0xFFFFFFFF)
////# e0 =(r1+8*r2-r0)
////# print("st",e0)
////# print(1<<32)
////# +((-(1<<32)+r0)&0xFFFFFFFF))
////# e0 = e0&(0xFFFFFFFF)
//			e1 = (e0 >> 29) + 8 * (r3 - 1);
////			e0 = (d1 - d0) & (0x1FFFFFFF);
//			e0 = (d1 - r0) & (0x1FFFFFFF);
////		# print(e0,e1)
////		# print(e1<<29 < (1<<63))
////		# e0=(d1-d0-(e1<<29))
////		# e0=(d1-r0-(e1<<29))
////		# print(((e0-e1)<<34))
//			if (e0 >= e1)
//			{
////			# v0=(r-((e0-e1)<<34))
////				v0 = (((r >> 32) - (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1- (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v1 = (d1 - e1 - 1);
//			}
//			else
//			{
//				v0 = ((e1 - e0) << 2) << 32;
//				v1 = (d1 - e1);
//			}
//		}
//	}
//	if (r3 == 0)
//	{
////	#d0 less than d1
//		if (r0 <= (r1 + 8 * r2))
//		{
////# print(d0,d1)
////# e1=(d1-d0)>>29
//			e0 = (r1 + 8 * r2 - r0);
//			e1 = (e0 >> 29) + 8 * r3;
////# e0 = e0&(0xFFFFFFFF)
////# print(e1<<29 < (1<<63))
////# e0=(d1-d0-(e1<<29))
//			e0 = (d1 - r0) & (0x1FFFFFFF);
////# print(((e0-e1)<<34))
//			if (e0 >= e1)
//			{
////# v0=(r-((e0-e1)<<34))
////				v0 = (((r >> 32) - (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1 - (4 * (e0 - e1))) & 0xFFFFFFFF) << 32;
//				v1 = (d1 - e1 - 1);
//			}
//			else
//			{
//				v0 = ((e1 - e0) << 2) << 32;
//				v1 = (d1 - e1);
//			}
//		}
//
////########################
////# d0>=d1
//		if (r0 > (r1 + 8 * r2))
//		{
////# print(d0,d1)
////			e1 = (d1 - d0) >> 29;
//			e1=0;
//			e0=(r0-r1-8*r2);
////# e0=(d0-d1-(e1<<29))
////# e0 =r0-(r1+8*r2)
////# e1=(e0>>29)-8*r3
////# e0 = e0&(0xFFFFFFFF)
////# e0=(d1-d0-(e1<<29))
////			e0 = (d1 - d0) & (0x1FFFFFFF);
//			if (e0>=e1)
//			{
//				v0=(e0-e1)<<34;
//				v1=(e1+d1);
//			}
//			else
//			{
//				v0 =(R-(e1<<34)+(e0<<34));
//				v1 = (e1+d1);
//			}
//		}
//	}
	/******************************************************************/

//	r1=r1+(r2*8);
	r1 = r1 + (r2 << 3);
	r3 = r3 << 3;

	bool r1Stat = (r0 <= r1);
	bool r3Stat = (r3 > 0);

	usfixn32 e0;
	//# d0 always less than d1
	//	if (r3 != 0)
	{

		//		if (r3 && (r0 <= (r1)))
		if (r3Stat && r1Stat)
		{
			//# print(d0,u1)
			//# r2=(u1-d0)>>29
//			e0 = (r1 - r0) >> 29;
//			//# print(e0,r2)
//			//			r2 = (e0 >> 29) + 8 * r3;
//			//			r2 = (e0 >> 29) + (r3<<3);
//			//			r2 = (e0) + (r3 << 3);
//			r2 = (e0) + (r3);
			r1 -= r0;
			r1 >>= 29;
			r2 = r1 + r3;
			//# print(e0,r2)
			//			r2 = (e0 >> 29) + 8 * r3;
			//			r2 = (e0 >> 29) + (r3<<3);
			//			r2 = (e0) + (r3 << 3);
//						r2 = (e0) + (r3);
			//# e0 = e0&(0xFFFFFFFF)
			//# print(r2<<29 < (1<<63))
			//# e0=(u1-d0-(r2<<29))

//			e0 = (u1_p - r0 - (r2 << 29)) & (0x1FFFFFFF);
			r1 = r2 << 29;
			r1 = u1_p - r1;
			r1 -= r0;
			e0 = r1 & (0x1FFFFFFF);
			//# e0=(u1-r0-(r2<<29))
			//# print(((e0-r2)<<34))
			v1 = u1 - r2;
			if (e0 >= r2)
			{
				//# v0=(r-((e0-r2)<<34))
				//				v0 = (((r >> 32) - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
				//				v0 = ((R1 - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1 - ((e0 - r2) << 2)) & 0xFFFFFFFF) << 32;
				e0 -= r2;
				e0 <<= 2;
				e0 = R1 - e0;
				e0 &= (0xFFFFFFFF);
				v0 = e0;
				v0 <<= 32;
//				v0 = ((R1 - ((e0 - r2) << 2)) & 0xFFFFFFFF) << 32;
//				v1 = (u1 - r2 - 1);
				v1--;
			}
			else
			{
				//				v0 = ((r2 - e0) << 2) << 32;
//				v0 = ((r2 - e0) << 34);
				r1 = (r2 - e0);
				v0 = r1;
				v0 <<= 34;
//				v1 = (u1 - r2);
			}
		}
		//		if (r3 && (r0 > (r1)))
		if (r3Stat && (!r1Stat))
		{
			//# print(d0,u1)
			//# r2=(u1-d0)>>29
			//# print((0xFFFFFFFF-r0+1))
			//# print(r1+8*r2)
//			e0 = ((r1) + (0xFFFFFFFF - r0 + 1)) >> 29;
//
//			//# e0=e0&(0xFFFFFFFF)
//			//# e0 =(r1+8*r2-r0)
//			//# print("st",e0)
//			//# print(1<<32)
//			//# +((-(1<<32)+r0)&0xFFFFFFFF))
//			//# e0 = e0&(0xFFFFFFFF)
//			//			r2 = (e0 >> 29) + 8 * (r3 - 1);
//			//			r2 = (e0 >> 29) + ((r3 - 1)<<3);
//			//			r2 = (e0) + ((r3 - 1) << 3);
//			r2 = (e0) + ((r3 - 8));
			r1 -= r0;
			r1 >>= 29;
			r1 += (r3);
			r2 = r1 - 8;

			//			e0 = (u1 - d0) & (0x1FFFFFFF);
//			e0 = (u1_p - r0) & (0x1FFFFFFF);
			e0 = (u1_p - r0) & (0x1FFFFFFF);
			//		# print(e0,r2)
			//		# print(r2<<29 < (1<<63))
			//		# e0=(u1-d0-(r2<<29))
			//		# e0=(u1-r0-(r2<<29))
			//		# print(((e0-r2)<<34))
			v1 = u1 - r2;
			if (e0 >= r2)
			{
				//			# v0=(r-((e0-r2)<<34))
				//				v0 = (((r >> 32) - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
				//				v0 = ((R1 - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1 - ((e0 - r2) << 2)) & 0xFFFFFFFF) << 32;
				r1 = (e0 - r2);
				r1 <<= 2;
				r1 = R1 - (r1);
				r1 &= (0xFFFFFFFF);
				v0 = r1;
				v0 <<= 32;
//				v0 = ((R1 - ((e0 - r2) << 2)) & 0xFFFFFFFF) << 32;
//				v1 = (u1 - r2 - 1);
				v1--;
			}
			else
			{
//				v0 = ((r2 - e0) << 2) << 32;
				r1 = r2 - e0;
				v0 = r1;
				v0 <<= 34;
//				v1 = (u1 - r2);
			}
		}
	}
	//	if (r3 == 0)

	{
		//	#d0 less than u1
		//		if (r0 <= (r1))
		//		if (!r3 && (r0 <= (r1)))
		if ((!r3Stat) && r1Stat)
		{
			//# print(d0,u1)
			//# r2=(u1-d0)>>29
//			e0 = (r1 - r0) >> 29;
////			//			r2 = (e0 >> 29) + (8 * r3);
////			//			r2 = (e0 >> 29) + (r3<<3);
////			//			r2 = (e0) + (r3 << 3);
//			r2 = (e0) + (r3);
			r1 -= r0;
			r1 >>= 29;
			r1 += r3;
			r2 = r1;

//			asm("{\n\t"
//					".reg .u32 t;\n\t"
//					"sub.u32 %1,%1,%2;\n\t"
//					"shl.b32 %1,%1,0x1D;\n\t"
//					"add.u32 %0,%1,%3;\n\t"
//					"}"
//					:"=r"(r2),"+r"(r1)
//					:"r"(r0),"r"(r3));
			//# e0 = e0&(0xFFFFFFFF)
			//# print(r2<<29 < (1<<63))
			//# e0=(u1-d0-(r2<<29))

//			e0 = (u1_p - r0) & (0x1FFFFFFF);
			r1 = u1_p - r0;
			e0 = r1 & (0x1FFFFFFF);
			//# print(((e0-r2)<<34))
			v1 = u1 - r2;
			if (e0 >= r2)
			{
				//# v0=(r-((e0-r2)<<34))
				//				v0 = (((r >> 32) - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
				//				v0 = ((R1 - (4 * (e0 - r2))) & 0xFFFFFFFF) << 32;
//				v0 = ((R1 - ((e0 - r2) << 2)) & 0xFFFFFFFF) << 32;
				r1 = e0 - r2;
				r1 <<= 2;
				r1 = R1 - r1;
				r1 &= (0xFFFFFFFF);
				v0 = r1;
				v0 <<= 32;
//				v1 = (u1 - r2 - 1);
				v1--;
			}
			else
			{
//				v0 = ((r2 - e0) << 34);
				r1 = (r2 - e0);
				v0 = r1;
				v0 <<= 34;

//				v1 = (u1 - r2);
			}
		}

		//########################
		//# d0>=u1

		//		if (r0 > (r1))
		//		if (!r3 && (r0 > (r1)))
		if ((!r3Stat) && (!r1Stat))
		{
			//# print(d0,u1)
			//			r2 = (u1 - d0) >> 29;
			//					r2 = 0;

			e0 = (r0 - r1);
//			asm("sub.u32 %0,%1,%2;":"=r"(e0):"r"(r0),"r"(r1));
			//# e0=(d0-u1-(r2<<29))
			//# e0 =r0-(r1+8*r2)
			//# r2=(e0>>29)-8*r3
			//# e0 = e0&(0xFFFFFFFF)
			//# e0=(u1-d0-(r2<<29))
			//			e0 = (u1 - d0) & (0x1FFFFFFF);
			//					if (e0 >= r2)
			{
				//						v0 = (e0 - r2) << 34;
				//						v1 = (r2 + u1);
				v0 = (e0);
				v0 <<= 34;
				v1 = (u1);
			}
			//			else
//			{
////				v0 = (R - (e1 << 34) + (e0 << 34));
//				v0 = (R1 - (e1 << 2) + (e0 <<2))<<32;
//				v1 = (e1 + d1);
//			}
		}
	}
	/******************************************************************/
//	v2 = 0;
//	s1_out=v1;
//	return;
	/******************************************************/
//	c=s;
//	c>>=64;
//	 c=(s >> 64);
//	u0 = (s % (1 << 64)) + c * RC;
//	u0 = (s);
//	u1=0;
//	if(y2)
//		u1+=x0;
//	if(x2)
//		u1+=y0;
////	u1 = x0 * y2 + x2 * y0;
//
//	s=0;
//	if(x2)
//		s+=y1;
//	if(y2)
//		s+=x1;
	u1 = 0;
	s = 0;

//	u2 = x0;
//	u3 = x1;
//	if (y2)
//	{
//
////		u1 += x0;
////		s += x1;
//		asm("add.u64 %0,%0,%1;":"+l"(u1):"l"(u2));
//		asm("add.u64 %0,%0,%1;":"+l"(s):"l"(u3));
//	}
//	u2 = y0;
//	u3 = y1;
//	if (x2)
//	{
////		u1 += y0;
////		s += y1;
//		asm("add.u64 %0,%0,%1;":"+l"(u1):"l"(u2));
//		asm("add.u64 %0,%0,%1;":"+l"(s):"l"(u3));
//	}

	r0 = 0;
	r1 = 0;
	if (y2)
	{

		//		u1 += x0;
		//		s += x1;
//			asm("add.u64 %0,%0,%1;":"+l"(u1):"l"(u2));
//			asm("add.u64 %0,%0,%1;":"+l"(s):"l"(u3));
		r0 += x0;
		r1 += x1;
	}
//		u2 = y0;
//		u3 = y1;
	if (x2)
	{
		//		u1 += y0;
		//		s += y1;
//			asm("add.u64 %0,%0,%1;":"+l"(u1):"l"(u2));
//			asm("add.u64 %0,%0,%1;":"+l"(s):"l"(u3));
		r0 += y0;
		r1 += y1;
	}

	u1 += r0;
	s += r1;

	u2 = 0;
	u3 = 0;
//	u1 = x0 * y2 + x2 * y0;
//	s = x2 * y1 + x1 * y2;

//	usfixn32 s0 = s;
////	s0 = s0 % R1;
//	usfixn32 s1 = s;
//	if(s>R1)
//	s1 = s>>31;
//
//	s0=s-(R<<s1);
//
//
////	u2=(s0<<32)-r0*s1;
//	u2 = (s0);
//	u2 <<= 32;
////	u2 -= R0 * s1;
//
////	u3 = x2 * y2 + s1;
//	u3 = s1+(x2&&y2);

//	usfixn32 s0 = s;
//	s0 = s0 % R1;
//	usfixn32 r0 = s;
	r0 = s;
	if (s > R1)
		r0 = s >> 31;
	u3 = r0 + (x2 && y2);

//	r0 = s;// - (R << r0);  // ?!?! to be fixed

	u2 = (r0);
//	u2=s;
	u2 <<= 32;
//	u2 -= R0 * s1;

//	u3 = x2 * y2 + s1;

	if (u0 >= R)
	{
		u0 -= R;
		u1++;
	}
//
//	v0 = 0;
//	v1 = 0;
//	v2 = 0;
//
////	v0, v1;
////	c = c2 * 2;
//	c = c2 << 1;
//	d1 = c >> 29;
////	d0 = c - (d1 << 29);
//	d0 = (d1 << 29);
//	d0 = c - (d0);
//
//	if (d0 >= d1)
//	{
//		v2 = d0 - d1;
//		e1 = v2 >> 29;
//		e0 = v2 - (e1 << 29);
//		e0 >= e1 ? (v0 = (e0 - e1) << 34) : (v0 = R - (e1 << 34) + (e0 << 34));
//		e0 >= e1 ? (v1 = e1 + d1) : (v1 = e1 + d1 - 1);
//		/*
//		 if(e0>=e1)
//		 {
//		 v16=(e0-e1)<<34;
//		 v17=e1+d1;
//		 }
//		 else
//		 {
//		 v17=e1+d1-1;
//		 v16=R-(e1<<34)+(e0<<34);
//		 }
//		 */
//	}
//	else
//	{
//		//d1>d0
//		v2 = d1 - d0;
//		e1 = v2 >> 29;
//		e0 = v2 - (e1 << 29);
//		e0 >= e1 ? (v0 = R - ((e0 - e1) << 34)) : (v0 = (e1 - e0) << 34);
//		e0 >= e1 ? (v1 = d1 - e1 - 1) : (v1 = d1 - e1);
//	}
//
////	e1 = (d0 - d1) >> 29;
//	e1 = (d0 - d1);
//	e1 = e1 >> 29;
//	e0 = (d0 - d1 - (e1 << 29));
//	v0 = (e0 - e1);
//	v0 <<= 34;
//	v1 = (e1 + d1);

//	usfixn32 cc;
//	asm("add.cc.u64 %0,%1,%2;":"=l"( s):"l"(u1),"l"(u2));
//	asm("addc.u32 %0,0,0;":"=r"(cc):);
////	v2 = 0;
////	s = u1 + u2;
//
//	if (cc > 1)
//	{
//		if (s >= R)
//			s = s - R;
//		else
//			s = s + RC;
//		u3++;
//	}

//	asm("add.cc.u64 %0,%1,%2;":"=l"( s):"l"(u1),"l"(u2));
//	asm("addc.u32 %0,0,0;":"=r"(r0):);

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.u32 %1,0,0;\n\t"
			"}"
			:"=l"(s),"=r"(r0):"l"(u1),"l"(u2)
	);

	//	v2 = 0;
	//	s = u1 + u2;

	if (r0 > 1)
	{
		if (s >= R)
			s = s - R;
		else
			s = s + RC;
		u3++;
	}

//	if (s < u1 || s < v2)
//	{
//		c = 1;
//		s = s + RC;
//	}
//	else if (s >= R)
//	{
//		c = 1;
//		s = s - R;
//	}
//	else
//	{
//		c = 0;
//	}
	u1 = s;
//  u3 += c;
//	c2 = c;

//	s = u1 + 2 * c2;
	c2 <<= 1;
//	s = u1 + c2;

//	if (s < u1 || s < (c2 * 2))

//	asm("add.cc.u64 %0,%1,%2;":"=l"(s):"l"(u1),"l"(c2));
//	asm("addc.u32 %0,0,0;":"=r"(cc):);
//
//	if (cc > 1)
//	{
//		if (s >= R)
//			s = s - R;
//		else
//			s = s + RC;
//		u3++;
//	}

//	asm("add.cc.u64 %0,%1,%2;":"=l"(s):"l"(u1),"l"(c2));
//	asm("addc.u32 %0,0,0;":"=r"(r0):);

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.u32 %1,0,0;\n\t"
			"}"
			:"=l"(s),"=r"(r0):"l"(u1),"l"(c2)
	);
	if (r0 > 1)
	{
		if (s >= R)
			s = s - R;
		else
			s = s + RC;
		u3++;
	}
//  if (s < u1 || s < (c2))
//    {
//      c = 1;
//      s = s + RC;
//    }
//  else if (s >= R)
//    {
//      c = 1;
//      s = s - R;
//    }
//  else
//    {
//      c = 0;
//    }
//  u3 += c;
	u1 = s;

	c = 0;
	s = v0 + c;

//	asm("sub.cc.u64 %0,%0,%1;":"+l"(u0):"l"(s));
//	asm("addc.u64 %0,0,0;":"=l"(c):);
	if (u0 < s)
	{
		c = 1;
		u0 = R - s + u0;
//		asm("sub.cc.u64 %0,%0,%1;":"+l"(u0):"l"(s));
//		u0+=R;
	}
	else
	{
		c = 0;

		u0 = u0 - s;
	}

//
//	v0 = 0xFFFFFFFFFFFFFFFF - v0;
//	v0++;
////	s = u0 + v0;
//	asm("{\n\t"
//			"add.cc.u64 %0,%0,%2;\n\t"
//			"addc.u64 %1,0,0;\n\t"
//			"}"
//			:"+l"(u0),"=l"(c)
//			:"l"(v0)
//	);
//
////	s = s & (0xFFFFFFFFFFFFFFFF);
////	if (s < v0 || s < u0)
////		c = 1;
////	else
////		c = 0;
//	c = !c;
//	if (c == 1)
//	{
////		s = s - RC;
//		u0 = u0 - RC;
//	}

//	u0 = s;

	s = v1 + c;
	if (u1 < s)
	{
		c = 1;
		u1 = R - s + u1;
	}
	else
	{
		c = 0;
		u1 = u1 - s;
	}

//	v1 += c;
//	v1 = 0xFFFFFFFFFFFFFFFF - v1;
//	v1++;
//	//	s = u0 + v0;
//	asm("{\n\t"
//			"add.cc.u64 %0,%0,%2;\n\t"
//			"addc.u64 %1,0,0;\n\t"
//			"}"
//			:"+l"(u1),"=l"(c)
//			:"l"(v1)
//	);

	//	s = s & (0xFFFFFFFFFFFFFFFF);
	//	if (s < v0 || s < u0)
	//		c = 1;
	//	else
	//		c = 0;
//	c = !c;
//	if (c == 1)
//	{
//		//		s = s - RC;
//		u1 = u1 - RC;
//	}

//	s = v2 + c;

//	v1 = c;
//	v1 = 0xFFFFFFFFFFFFFFFF - v1;
//	v1++;
//	//	s = u0 + v0;
//	asm("{\n\t"
//			"add.cc.u64 %0,%0,%2;\n\t"
//			"addc.u64 %1,0,0;\n\t"
//			"}"
//			:"+l"(u3),"=l"(c)
//			:"l"(v1)
//	);

	s = c;
	if (u3 < s)
	{
		c = 1;
		u3 = R - s + u3;
	}
	else
	{
		c = 0;
		u3 = u3 - s;
	}
////		c=!c;
//	if (c == 1)
//		u3 = u3 - RC;
	s0_out = u0;
//	u1 = u1;
	s1_out = u1;
	s2_out = u3;

}

/**********************************************/
//__device__ __inline__ void device_p3_mulLong_2_plain_2
__device__ __inline__ void
device_p3_mulLong_2_revised_1_working (usfixn64 const &x, usfixn64 const &y,
														 usfixn64 &s0_out, usfixn64 &s1_out,
														 usfixn64 &s2_out)
{
	usfixn64 u0, u1, u2, u3;
	usfixn64 r0 = R;
	usfixn64 r1 = R;
	usfixn64 rc0 = RC;
	usfixn64 rc1 = RC;
//	r0 = r0 % (1 << 32);

	r0 = r0 & (0xFFFFFFFF);
	r1 >>= (32);

	rc0 = rc0 & (0xFFFFFFFF);
	rc1 >>= (32);

	usfixn64 x0 = (x & 0xFFFFFFFF);
	usfixn64 x2 = (x >= R);
	usfixn64 x1 = (x - x0 - x2 * R);
	x1 = (x2 * R);
	x1 = x - x0 - x1;
	x1 >>= 32;

//	x1=(x-x0-r*x2)>>32

//	y0 = y % (1 << 32)
//	y2 = int(y >= r)
//	y1 = (y - y0 - r * y2) >> 32
	usfixn64 y0 = (y & 0xFFFFFFFF);
	usfixn64 y2 = (y >= R);
	usfixn64 y1 = (y - y0 - y2 * R);
	y1 = (y2 * R);
	y1 = y - y0 - y1;
	y1 >>= 32;

	usfixn64 s;
	usfixn64 c;
//	usfixn64 t=x1;
//	usfixn64 t1;
//	t<<=32;
////	s=(x0 + (x1 << 32));
//	t= x0 + t;
////	* (y0 + (y1 << 32));
//	t1=y1;
//	t1<<=32;
//	t1+=y0;
//	s=t1*t;
//	c=0;
//	if(s<t1||s<t)
//		c=1;

/////////////////////////////////////////////// problem is here, in mutliplication at the beginnging,
//	usfixn64 a0, a1, a2, a3;
//	c = 0;
//	a0 = (x0 * y0) & (0xFFFFFFFF);
//	a1 = (x0 * y1) & (0xFFFFFFFF);
//	a2 = (x1 * y0) & (0xFFFFFFFF);
//	a3 = (x1 * y1) & (0xFFFFFFFF);
//
//	a0 = (a0) & (0xFFFFFFFF);
//	a1 = (a1) & (0xFFFFFFFF);
//	a2 = (a2) & (0xFFFFFFFF);
//	a3 = (a3) & (0xFFFFFFFF);
//
//	a1 <<= 32;
//	s = a0 + (a1);
//	if (s < a0 || s < a1)
//		c++;
//	a0 = s;
//	a2 <<= 32;
//	s = a0 + a2;
//	if (s < a0 || s < a2)
//		c++;
//
//	a1 = (x0 * y1);
//	a2 = (x1 * y0);
//	a1 = (a1) & (0xFFFFFFFF);
//	a2 = (a2) & (0xFFFFFFFF);
//	a1 >>= 32;
//	a2 >>= 32;
//	c = c + a3 + a1 + a2;
	usfixn64 X, Y;
	X = x1;
	Y = y1;
	X <<= 32;
	Y <<= 32;
	X = X + x0;
	Y = Y + y0;

	device_p3_mult_64bit_with_carry (X, Y, s, c);

//	c=s;
//	c>>=64;
//	 c=(s >> 64);
//	u0 = (s % (1 << 64)) + c * RC;
	u0 = (s);

	u1 = x0 * y2 + x2 * y0;

	s = x2 * y1 + x1 * y2;
	usfixn64 s0 = s;
	s0 = s0 % r1;
	usfixn64 s1 = s;
	s1 = s1 / r1;

//	u2=(s0<<32)-r0*s1;
	u2 = (s0);
	u2 <<= 32;
	u2 -= r0 * s1;

	u3 = x2 * y2 + s1;

	if (u0 >= R)
	{
		u0 -= R;
		u1++;
	}

	usfixn64 v0, v1, v2;
	usfixn64 k0, k1;
	v0 = 0;
	v1 = 0;
	v2 = 0;
	usfixn64 c2 = c;

	usfixn64 d0, d1, e0, e1;
	usfixn64 r = R;

	usfixn64 aux;
//	v0, v1;
	c = c2 * 2;
	d1 = c >> 29;
//	d0 = c - (d1 << 29);
	d0 = (d1 << 29);
	d0 = c - (d0);

//	e1 = (d0 - d1) >> 29;
//	e1 = (d0 - d1);
////	aux=e1;
////	e1 = e1 >> 29;
//	e1 >>=29;
////	e0 = (d0 - d1 - (e1 << 29));
//	e0 = ((e1 << 29));
//	e0 = d0-d1-e0;
//	v0 = (e0 - e1);
//	v0 <<= 34;
//	v1 = (e1 + d1);
	usfixn64 d2;

	if (d0 >= d1)
	{
		d2 = d0 - d1;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v0 = (e0 - e1) << 34) : (v0 = R - (e1 << 34) + (e0 << 34));
		e0 >= e1 ? (v1 = e1 + d1) : (v1 = e1 + d1 - 1);
		/*
		 if(e0>=e1)
		 {
		 v16=(e0-e1)<<34;
		 v17=e1+d1;
		 }
		 else
		 {
		 v17=e1+d1-1;
		 v16=R-(e1<<34)+(e0<<34);
		 }
		 */
	}
	else
	{
		//d1>d0
		d2 = d1 - d0;
		e1 = d2 >> 29;
		e0 = d2 - (e1 << 29);
		e0 >= e1 ? (v0 = R - ((e0 - e1) << 34)) : (v0 = (e1 - e0) << 34);
		e0 >= e1 ? (v1 = d1 - e1 - 1) : (v1 = d1 - e1);
	}
	v2 = 0;

//	c = 0;
//	s = v0 + c;
//	if (0 < s)
//	{
//		v0 = r - v0;
//		c = 1;
//	}
//	else
//	{
//		c = 0;
//	}
//
//	s = v1 + c;
//	if (0 < s)
//	{
//		v1 = r - v1;
//		c = 1;
//	}
//	else
//	{
//		c = 0;
//	}
//	s = v2 + c;
//	if (0 < s)
//	{
//		v2 = r - v2;
//		c = 1;
//	}
//	else
//	{
//		c = 0;
//	}

//	s0_out = u0;
//		u1 = u1 + u2;
//		s1_out = u1;
//		s2_out = u3

	s = u1 + u2;
	if (s < u1 || s < v2)
	{
		c = 1;
		s = s + RC;
	}
	else if (s >= R)
	{
		c = 1;
		s = s - R;
	}
	else
	{
		c = 0;
	}
	u1 = s;
	u3 += c;
//	c2 = c;

	s = u1 + 2 * c2;
	if (s < u1 || s < (c2 * 2))
	{
		c = 1;
		s = s + RC;
	}
	else if (s >= R)
	{
		c = 1;
		s = s - R;
	}
	else
	{
		c = 0;
	}
	u1 = s;
	u3 += c;

	c = 0;
	s = v0 + c;
	if (u0 < s)
	{
		c = 1;
		u0 = r - s + u0;
	}
	else
	{
		c = 0;
		u0 = u0 - s;
	}

	s = v1 + c;
	if (u1 < s)
	{
		c = 1;
		u1 = r - s + u1;
	}
	else
	{
		c = 0;
		u1 = u1 - s;
	}

	s = v2 + c;
	if (u3 < s)
	{
		c = 1;
		u3 = r - s + u3;
	}
	else
	{
		c = 0;
		u3 = u3 - s;
	}

	s0_out = u0;
//	u1 = u1;
//	u1=v1;
	s1_out = u1;
//	s1_out = aux;
	s2_out = u3;

}

/**********************************************/
// Multiplication in Z/(R^8 + 1)Z (us = xs * ys)
//__device__ void bigMult(usfixn64 * xs, usfixn64 *ys, usfixn64 *us)
//__device__ void bigMult_revised(usfixn64 * xs,
//		usfixn64 * __restrict__ ys, usfixn64 *us)
//{
//	usfixn64 ts1[8];
//	usfixn64 ts2[8];
//	usfixn64 rs[3];
//	short c0, c1, c2, c3, c4, c5, c6, c7;
//	usfixn64 l0, l1, l2, l3, l4, l5, l6, l7, h0, h1, h2, h3, h4, h5, h6, h7;
//
//
//	usfixn64 r0,r1,r2;
//	//x0*y0
////	device_p3_mulLong_2_plain_2(xs[0], ys[0], rs[0], rs[1], rs[2]);
////	l0 = rs[0];
////	h0 = rs[1];
////	c0 = (short) rs[2];
//	device_p3_mulLong_2_plain_2(xs[0], ys[0], r0, r1, r2);
//	l0 = r0;
//	h0 = r1;
//	c0 = (short) r2;
//
//	//x0*y1+x1*y0
//	device_p3_mulLong_2_plain(xs[0], ys[1], rs[0], rs[1], rs[2]);
//	l1 = rs[0];
//	h1 = rs[1];
//	c1 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y2+x1*y1+x2*y0
//	device_p3_mulLong_2_plain(xs[0], ys[2], rs[0], rs[1], rs[2]);
//	l2 = rs[0];
//	h2 = rs[1];
//	c2 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y3+x1*y2+x2*y1+x3*y0
//	device_p3_mulLong_2_plain(xs[0], ys[3], rs[0], rs[1], rs[2]);
//	l3 = rs[0];
//	h3 = rs[1];
//	c3 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
//	device_p3_mulLong_2_plain(xs[0], ys[4], rs[0], rs[1], rs[2]);
//	l4 = rs[0];
//	h4 = rs[1];
//	c4 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
//	device_p3_mulLong_2_plain(xs[0], ys[5], rs[0], rs[1], rs[2]);
//	l5 = rs[0];
//	h5 = rs[1];
//	c5 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
//	device_p3_mulLong_2_plain(xs[0], ys[6], rs[0], rs[1], rs[2]);
//	l6 = rs[0];
//	h6 = rs[1];
//	c6 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
//
//	//x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
//	device_p3_mulLong_2_plain(xs[0], ys[7], rs[0], rs[1], rs[2]);
//	l7 = rs[0];
//	h7 = rs[1];
//	c7 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[1], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[2], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[0], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
//
//	// (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
//	//	ts1[0] = l0;
//	//	ts1[1] = h0;
//	//	ts1[2] = c0;
//	//	ts1[3] = c1;
//	//	ts1[4] = c2;
//	//	ts1[5] = c3;
//	//	ts1[6] = c4;
//	//	ts1[7] = c5;
//	//
//	//	ts2[0] = 0;
//	//	ts2[1] = l1;
//	//	ts2[2] = h1;
//	//	ts2[3] = h2;
//	//	ts2[4] = h3;
//	//	ts2[5] = h4;
//	//	ts2[6] = h5;
//	//	ts2[7] = h6;
//	//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	//	ts2[0] = 0;
//	//	ts2[1] = 0;
//	//	ts2[2] = l2;
//	//	ts2[3] = l3;
//	//	ts2[4] = l4;
//	//	ts2[5] = l5;
//	//	ts2[6] = l6;
//	//	ts2[7] = l7;
//	////	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//	//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//
//	//###############################################
//	ts1[0] = 0;
//	ts1[1] = 0;
//	ts1[2] = c0;
//	ts1[3] = c1;
//	ts1[4] = c2;
//	ts1[5] = c3;
//	ts1[6] = c4;
//	ts1[7] = c5;
//
//	ts2[0] = 0;
//	ts2[1] = h0;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
//
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	ts2[0] = l0;
//	ts2[1] = l1;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = l7;
//	//	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	//###############################################
//	ts2[0] = c6;
//	ts2[1] = c7;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	//	device_p3_bigSub(ts1, ts2, um);
//	device_p3_bigSub(ts1, ts2, ts1);
//	ts2[0] = h7;
//	ts2[1] = 0;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	//	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	//(x7*y7)r^6
//	device_p3_mulLong_2_plain(xs[7], ys[7], rs[0], rs[1], rs[2]);
//	l6 = rs[0];
//	h6 = rs[1];
//	c6 = (short) rs[2];
//
//	//(x6*y7+x7*y6)r^5
//	device_p3_mulLong_2_plain(xs[6], ys[7], rs[0], rs[1], rs[2]);
//	l5 = rs[0];
//	h5 = rs[1];
//	c5 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[7], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
//
//	//(x5*y7+x6*y6+x7*y5)r^4
//	device_p3_mulLong_2_plain(xs[5], ys[7], rs[0], rs[1], rs[2]);
//	l4 = rs[0];
//	h4 = rs[1];
//	c4 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[6], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
//
//	//(x4*y7+x5*y6+x6*y5+x7*y4)r^3
//	device_p3_mulLong_2_plain(xs[4], ys[7], rs[0], rs[1], rs[2]);
//	l3 = rs[0];
//	h3 = rs[1];
//	c3 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[5], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
//
//	//(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
//	device_p3_mulLong_2_plain(xs[3], ys[7], rs[0], rs[1], rs[2]);
//	l2 = rs[0];
//	h2 = rs[1];
//	c2 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[4], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
//
//	//(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
//	device_p3_mulLong_2_plain(xs[2], ys[7], rs[0], rs[1], rs[2]);
//	l1 = rs[0];
//	h1 = rs[1];
//	c1 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[3], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
//
//	//(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
//	device_p3_mulLong_2_plain(xs[1], ys[7], rs[0], rs[1], rs[2]);
//	l0 = rs[0];
//	h0 = rs[1];
//	c0 = (short) rs[2];
//	device_p3_mulLong_2_plain(xs[2], ys[6], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[3], ys[5], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[4], ys[4], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[5], ys[3], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[6], ys[2], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//	device_p3_mulLong_2_plain(xs[7], ys[1], rs[0], rs[1], rs[2]);
//	device_p3_smallAdd(&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
//
//	//	if(c6!=0)
//	//		printf("not zero \n");
//	////#####################################################
//	////(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
//	//	ts2[0] = l0;
//	//	ts2[1] = h0;
//	//	ts2[2] = c0;
//	//	ts2[3] = c1;
//	//	ts2[4] = c2;
//	//	ts2[5] = c3;
//	//	ts2[6] = c4;
//	//	ts2[7] = c5;
//	////	device_p3_bigSub(ts1, ts2);
//	//	device_p3_bigSub(ts1, ts2, ts1);
//	//
//	//	ts2[0] = 0;
//	//	ts2[1] = l1;
//	//	ts2[2] = h1;
//	//	ts2[3] = h2;
//	//	ts2[4] = h3;
//	//	ts2[5] = h4;
//	//	ts2[6] = h5;
//	//	ts2[7] = h6;
//	////	device_p3_bigSub(ts1, ts2);
//	//	device_p3_bigSub(ts1, ts2, ts1);
//	//
//	//	ts2[0] = 0;
//	//	ts2[1] = 0;
//	//	ts2[2] = l2;
//	//	ts2[3] = l3;
//	//	ts2[4] = l4;
//	//	ts2[5] = l5;
//	//	ts2[6] = l6;
//	//	ts2[7] = 0;
//	////	device_p3_bigSub(ts1, ts2);
//	//	device_p3_bigSub(ts1, ts2, ts1);
//	//
//	//	ts2[0] = c6;
//	//	ts2[1] = 0;
//	//	ts2[2] = 0;
//	//	ts2[3] = 0;
//	//	ts2[4] = 0;
//	//	ts2[5] = 0;
//	//	ts2[6] = 0;
//	//	ts2[7] = 0;
//	//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	//
//	//	us[0] = ts1[0];
//	//	us[1] = ts1[1];
//	//	us[2] = ts1[2];
//	//	us[3] = ts1[3];
//	//	us[4] = ts1[4];
//	//	us[5] = ts1[5];
//	//	us[6] = ts1[6];
//	//	us[7] = ts1[7];
//	//#####################################################
//	//(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = c0;
//	ts2[3] = c1;
//	ts2[4] = c2;
//	ts2[5] = c3;
//	ts2[6] = c4;
//	ts2[7] = c5;
//	//	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = h0;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
//	//	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = l0;
//	ts2[1] = l1;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = 0;
//	//	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	//c7=0, h7=0, l7=0 as there are no subtraction for coefficient of r^7
//	ts2[0] = c6;
//	ts2[1] = 0;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//
//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
//	//#####################################################
//}
/**************************************************************************/
/**************************************************************************/
//__device__ __inline__ void mulLong_register2(usfixn64 x0, usfixn64 y0,
//		usfixn64 &s0, usfixn64 &s1, usfixn64 &s2)
//{
//	short x1 = 0, y1 = 0;
////	usfixn64 l, h, c;
////		usfixn64 l;
////		usfixn64 h;	//replaced with s1
////		usfixn64 c; //
//
////	usfixn64 x0, y0;
////	usfixn64 v2, v5, v9, v10, v11, v14, v15, v16, v17;
//	usfixn64 v2, v5;
////	usfixn64 v10, v11;
////	usfixn64  v14, v15;
////	usfixn64 v16, v17;
////	usfixn64 q;
////	usfixn64 t;
//	usfixn64 t0, t1, t2;
//	usfixn64 a0, a1, b0, b1, c0, c1;
////	usfixn64 c1prime, d0, d1, d2, e0, e1;
////	usfixn64 d0, d1, e0, e1;
//
//	if (x0 <= SQRTR && y0 <= SQRTR)
//	{
////		s[0] = x * y;
////		s[1] = 0;
////		s[2] = 0;
//		s0 = x0 * y0;
//		s1 = 0;
//		s2 = 0;
//		return;
//	}
//
////	x1 = (x0 >= R ? 1 : 0);
////	x0 = (x1 > 0 ? x0 - R : x0);
////
//	if (x0 >= R)
//	{
//		x1 = 1;
//		x0 -= R;
//	}
//
////	y1 = (y >= R ? 1 : 0);
////	y0 = (y1 > 0 ? y - R : y);
//	if (y0 >= R)
//	{
//		y1 = 1;
//		y0 -= R;
//	}
//
//	if (x0 <= SQRTR && y0 <= SQRTR)
//	{
////		s[0] = x0 * y0;
////		s[1] = h;
////		s[2] = c;
//		s0 = x0 * y0;
////		s1 = h;
////		s2 = s2;
//
//		return;
//	}
//
////can use some if and else for here, as x1 a
//	v2 = x0 * y1;		//[0,v2,0];
//	v5 = x1 * y0;		//[0,v5,0];
//	s2 = x1 * y1;		//[0,0,1];
//
////	c = v9;
////	l = 0;
////	h = v5 + v2;
////	h = v5 + v2;
//	s1 = v5 + v2;
////	h < v5 || h < v2 ? (c = c + 1) : (c = c);
////	c > v9 ? (h = h + RC) : (h = h);
////	if (h < v5 || h < v2)
////	{
////		c++;
////		h += RC;
////	}
//	if (s1 < v5 || s1 < v2)
//	{
////			c++;
//		s2++;
//		s1 += RC;
//	}
//
////lhc
////x0*y0
////	a1 = x0 >> 32;
////	a0 = x0 - (a1 << 32);
////	b1 = y0 >> 32;
////	b0 = y0 - (b1 << 32);
//
//	c0 = 0;
//	a1 = x0 >> 32;
//	b1 = y0 >> 32;
//
////	b0 = y0 - ((y0 >> 32) << 32);
////	a0 = x0 - ((x0 >> 32) << 32);
//
//	a0 = x0 - (a1 << 32);
//	b0 = y0 - (b1 << 32);
//
//	c1 = a1 * b1;
//
////######################################################
//
////
////	t = (a0 * b1);
//////	q = t >> 32;
////	q = (a0 * b1) >> 32;
////	t = (t - (q << 32)) << 32;
////	c1 += q;
////	c0 += t;  //safe
////
////	t = a1 * b0;
//////	q = t >> 32;
////	q = (a1 * b0) >> 32;
////	t = (t - (q << 32)) << 32;
////	c1 += q;
////	q = c0 + t;               //here, is not related to r.
//////	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
////	if (q < c0 || q < t)
////		c1++;
////	c0 = q;
////
//////	t = a0 * b0;
//////	q = c0 + t;
////	t = a0 * b0;
////	q = c0 + (a0 * b0);
//////	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=x0*y0
////	if (q < c0 || q < t)
////		c1++;
////	c0 = q;
////
//////	c1prime = c1 << 1;
////	c1 <<= 1;
////
////	v5 = 0;
////	v2 = c0;
//////	c0 >= R ? (v11 = 1) : (v11 = 0);
//////	v11 > 0 ? (v10 = c0 - R) : (v10 = c0);
//////v12=0;
////	if (c0 >= R)
////	{
////		v5 = 1;
//////			v10=c0-R;
////		v2 -= R;
////	}
////
//////	s0=l; s1=h; b0=c;
////	s0 = 0;
//////	q = l + v2;  //[l,h,c] + [v10,v11,0]
////	q = s0 + v2;  //[l,h,c] + [v10,v11,0]
//////	q < l || q < v10 ? (v11 = v11 + 1) : (v11 = v11);
//////	q < l || q < v10 ? (l = q + RC) : (l = q);
////
////	if (q < s0 || q < v2)
////	{
////		v5++;
//////		l = q + RC;
////		s0 = q + RC;
////	}
////	else
////	{
//////		l = q;
////		s0 = q;
////	}
////
//////	if (l >= R)
//////	{
//////		l = l - R;
//////		v5++;
//////	}
////	if (s0 >= R)
////	{
////		s0 = s0 - R;
////		v5++;
////	}
////
////	q = s1 + v5;
//////	q < h || q < v11 ? (c = c + 1) : (c = c);
//////	q < h || q < v11 ? (h = q + RC) : (h = q);
////	if (q < s1 || q < v5)
////	{
//////		c++;
////		s2++;
////		s1 = q + RC;
////	}
////	else
////	{
////		s1 = q;
////	}
////
////	if (s1 >= R)
////	{
////		s1 = s1 - R;
//////		c++;
////		s2++;
////	}
////
////	q = s1 + v5;
////	//	q < h || q < v11 ? (c = c + 1) : (c = c);
////	//	q < h || q < v11 ? (h = q + RC) : (h = q);
////	if (q < s1 || q < v5)
////	{
//////			c++;
////		s2++;
////		s1 = q + RC;
////	}
////	else
////	{
////		s1 = q;
////	}
////
////	if (s1 >= R)
////	{
////		s1 = s1 - R;
//////			c++;
////		s2++;
////	}
//
////######################################################
//
//	t0 = (a0 * b1);
//	t1 = a1 * b0;
//	t2 = a0 * b0;
//
////	q0,q1,q2=a0,a1,b0
//
////	q = t >> 32;
////	q = (a0 * b1) >> 32;
//	a0 = (t0) >> 32;
//	a1 = t1 >> 32;
//	b0 = c0 + (t2);
//
//	t0 = (t0 - (a0 << 32)) << 32;
//	t1 = (t1 - (a1 << 32)) << 32;
//
//	c1 += a0;
//	c0 += t0;		//safe
//
////	q = t >> 32;
////	q = (a1 * b0) >> 32;
//
//	c1 += a1;
//	a1 = c0 + t1;		//here, is not related to r.
////	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
//	if (a1 < c0 || a1 < t1)
//		c1++;
//	c0 = a1;
//
////	t = a0 * b0;
////	q = c0 + t;
//
////	q = c0 + (a0 * b0);
//
////	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=x0*y0
//	if (b0 < c0 || b0 < t2)
//		c1++;
//	c0 = b0;
//
////	c1prime = c1 << 1;
//	c1 <<= 1;
//	v5 = 0;
//	v2 = c0;
////	c0 >= R ? (v11 = 1) : (v11 = 0);
////	v11 > 0 ? (v10 = c0 - R) : (v10 = c0);
////v12=0;
//	if (c0 >= R)
//	{
//		v5 = 1;
//		//			v10=c0-R;
//		v2 -= R;
//	}
////###########
////	s0=l; s1=h; b0=c;
//	s0 = 0;
////	q = l + v2;  //[l,h,c] + [v10,v11,0]
//	a0 = s0 + v2;		//[l,h,c] + [v10,v11,0]
////	q < l || q < v10 ? (v11 = v11 + 1) : (v11 = v11);
////	q < l || q < v10 ? (l = q + RC) : (l = q);
//
//	if (a0 < s0 || a0 < v2)
//	{
//		v5++;
//		//		l = q + RC;
//		s0 = a0 + RC;
//	}
//	else
//	{
//		//		l = q;
//		s0 = a0;
//	}
//
////	if (l >= R)
////	{
////		l = l - R;
////		v5++;
////	}
//	if (s0 >= R)
//	{
//		s0 = s0 - R;
//		v5++;
//	}
//
//	a1 = s1 + v5;
////	q < h || q < v11 ? (c = c + 1) : (c = c);
////	q < h || q < v11 ? (h = q + RC) : (h = q);
//	if (a1 < s1 || a1 < v5)
//	{
//		//		c++;
//		s2++;
//		s1 = a1 + RC;
//	}
//	else
//	{
//		s1 = a1;
//	}
//
//	if (s1 >= R)
//	{
//		s1 = s1 - R;
//		//		c++;
//		s2++;
//	}
//
//	a1 = s1 + v5;
////	q < h || q < v11 ? (c = c + 1) : (c = c);
////	q < h || q < v11 ? (h = q + RC) : (h = q);
//	if (a1 < s1 || a1 < v5)
//	{
//		//			c++;
//		s2++;
//		s1 = a1 + RC;
//	}
//	else
//	{
//		s1 = a1;
//	}
//
//	if (s1 >= R)
//	{
//		s1 = s1 - R;
//		//			c++;
//		s2++;
//	}
////#####################################
////v13=0;
////	c1 >= R ? (v15 = 1) : (v15 = 0);
////	v15 > 0 ? (v14 = c1 - R) : (v14 = c1); //v13=0;
//
////	v15=0;
//	v5 = 0;
//	v2 = c1;
//	if (c1 >= R)		//-> v15=1>0
//	{
//		v5 = 1;
//		v2 = v2 - R;
//	}
//
//	b0 = s1 + v2;  //[l,h,c]+[0,v14,v15]
////	q < h || q < v14 ? (c = c + v15 + 1) : (c = c + v15);
////	q < h || q < v14 ? (h = q + RC) : (h = q);
////	c = c + v5;
//	s2 += v5;
//	if (b0 < s1 || b0 < v2)
//	{
//		//		c++;
//		s2++;
//		s1 = b0 + RC;
//	}
//	else
//	{
//		s1 = b0;
//	}
//
//	if (s1 >= R)
//	{
//		s1 = s1 - R;
//		//		c++;
//		s2++;
//	}
////#################
////[l,h,c]
//
////	d1 = c1 >> 29;
//	a1 = c1 >> 29;
////	d0 = c1 - (d1 << 29);
////	d0 = c1 - ((c1 >> 29) << 29);
//	a0 = c1 - ((c1 >> 29) << 29);
////[l,h,c]-[v16,v17,0]
//	if (a0 >= a1)
//	{
//		//		d2 = a0 - a1;
//		b1 = (a0 - a1) >> 29;
//		//		b0 = d2 - (b1 << 29);
//		b0 = (a0 - a1) - (b1 << 29);
//
//		//		b0 >= b1 ? (v16 = (b0 - b1) << 34) : (v16 = R - (b1 << 34) + (b0 << 34));
//		//		b0 >= b1 ? (v17 = b1 + a1) : (v17 = b1 + a1 - 1);
//
//		v5 = b1 + a1;
//		if (b0 >= b1)
//		{
//			v2 = (b0 - b1) << 34;
//		}
//		else
//		{
//			v2 = R - (b1 << 34) + (b0 << 34);
//			v5--;
//		}
//		/*
//		 if(b0>=b1)
//		 {
//		 v16=(b0-b1)<<34;
//		 v17=b1+a1;
//		 }
//		 else
//		 {
//		 v17=b1+a1-1;
//		 v16=R-(b1<<34)+(b0<<34);
//		 }
//		 */
//	}
//	else
//	{
//		//a1>a0
//		//		d2 = a1 - a0;
//		//		b1 = d2 >> 29;
//		b1 = (a1 - a0) >> 29;
//		//		b0 = d2 - (b1 << 29);
//		b0 = (a1 - a0) - (b1 << 29);
//		//		b0 >= b1 ? (v16 = R - ((b0 - b1) << 34)) : (v16 = (b1 - b0) << 34);
//		//		b0 >= b1 ? (v17 = a1 - b1 - 1) : (v17 = a1 - b1);
//		v5 = a1 - b1;
//		if (b0 >= b1)
//		{
//			v2 = R - ((b0 - b1) << 34);
//			v5--;
//		}
//		else
//		{
//			v2 = (b1 - b0) << 34;
//		}
//
//		/*
//		 if(b0>=b1)
//		 {
//		 v16=R-((b0-b1)<<34);
//		 v17=a1-b1-1;
//		 }
//		 else
//		 {
//		 v16=(b1-b0)<<34;
//		 v17=a1-b1;
//		 }
//		 */
//	}
//
////q
//	a0 = 0;
//	if (s0 >= v2)
//	{
//		s0 = s0 - v2;
//	}
//	else
//	{
//		s0 = R - v2 + s0;
//		a0 = 1;
//	}
////t
//	if (s1 < a0 + v5)
//	{
////		c = c - 1;
////		c--;
//		s2--;
//		s1 = R - a0 - v5 + s1;
//	}
//	else
//	{
//		s1 = s1 - a0 - v5;
//	}
////	s[0] = l;
////	s[1] = h;
////	s[2] = c;
////	s0 = l;
////	s1 = h;
////	s2 = c;
//}
/**********************************************/

//x0 and y0 will be changed, can only trust values in s0, s1, s2
__device__ void
device_p3_mulLong_register2_ptr (usfixn64 * __restrict__ x0, usfixn64 * __restrict__ y0,
											 usfixn64 * __restrict__ s0, usfixn64 *__restrict__ s1,
											 usfixn64 *__restrict__ s2)
{
//	usfixn64 t0, t1, t2;
//	usfixn64 a0 = 0, a1 = 0, b0, b1, c0 = 0, c1 = 0;
//
//	if (*x0 <= SQRTR && *y0 <= SQRTR)
//	{
//		*s0 = *x0 * *y0;
//		*s1 = 0;
//		*s2 = 0;
//		return;
//	}
//
//	if (*x0 >= R)
//	{
//		c0 = 1;
//		*x0 -= R;
//	}
//
//	if (*y0 >= R)
//	{
//		c1 = 1;
//		*y0 -= R;
//	}
//
//	if (*x0 <= SQRTR && *y0 <= SQRTR)
//	{
//		*s0 = *x0 * *y0;
//		return;
//	}
//
////can use some if and else for here, as c0 a
////		a0 = *x0 * c1; //[0,v2,0];
////		a1 = c0 * *y0; //[0,v5,0];
////		*s2 = c0 * c1; //[0,0,1];
//
////		a0 = *x0 * c1; //[0,v2,0];
//
////	a1 = c0 * *y0; //[0,v5,0];
////		*s2 = c0 * c1; //[0,0,1];
//
//	*s1 = a1 + a0;
//	if (*s1 < a1 || *s1 < a0)
//	{
//		*s2++;
//		*s1 += RC;
//	}
//
////###################################################### sets a0, a1
//	c0 = 0;
//	a1 = *x0 >> 32;
//	b1 = *y0 >> 32;
//
//	a0 = *x0 - (a1 << 32);
//	b0 = *y0 - (b1 << 32);
//
//	c1 = a1 * b1;
//
////######################################################
//
//	t0 = (a0 * b1);
//	t1 = a1 * b0;
//	t2 = a0 * b0;
//
//	a0 = (t0) >> 32;
//	a1 = t1 >> 32;
//	b0 = c0 + (t2);
//
//	t0 = (t0 - (a0 << 32)) << 32;
//	t1 = (t1 - (a1 << 32)) << 32;
//
//	c1 += a0;
//	c0 += t0;		//safe
//
//	c1 += a1;
//	a1 = c0 + t1;		//here, is not related to r.
////	q < c0 || q < t ? (c1++) : (c1 = c1);  //c0=c0+t and carry, safe
////	if (a1 < c0 || a1 < t1)
////		c1++;
//
//	c1 += (a1 < c0 || a1 < t1);
//	c0 = a1;
//
////	q < c0 || q < t ? (c1++) : (c1 = c1);  //Now we finish [c0,c1]=*x0**y0
////	if (b0 < c0 || b0 < t2)
////		c1++;
//
//	c1 += (b0 < c0 || b0 < t2);
//	c0 = b0;
//
////	c1prime = c1 << 1;
//	c1 <<= 1;
////	*y0 = 0;
////	*x0 = c0;
////
////	if (c0 >= R)
////	{
////		*y0 = 1;
////		//			v10=c0-R;
////		*x0 -= R;
////	}
//
//	*y0 = (c0 >= R);
//	*x0 = c0;
//
//	if (c0 >= R)
//	{
////		*y0 = 1;
//		//			v10=c0-R;
//		*x0 -= R;
//	}
//
////########### //sets a0, a1,
//	*s0 = 0;
//	a0 = *s0 + *x0;		//[l,h,c] + [v10,v11,0]
//
////	if (a0 < *s0 || a0 < *x0)
////	{
////		*y0++;
////		*s0 = a0 + RC;
////	}
////	else
////	{
////		*s0 = a0;
////	}
//
////	if
////	{
//	*y0 += (a0 < *s0 || a0 < *x0);
//	*s0 = a0 + (a0 < *s0 || a0 < *x0) * RC;
////	}
////	else
////	{
////		*s0 = a0;
////	}
//
////
////	if (*s0 >= R)
////	{
////		*s0 = *s0 - R;
////		*y0++;
////	}
//
////			if
////			{
//	*y0 += (*s0 >= R);
//	*s0 = *s0 - (*s0 >= R) * R;
//
////			}
//
//	a1 = *s1 + *y0;
////	if (a1 < *s1 || a1 < *y0)
////	{
////		//		c++;
////		*s2++;
////		*s1 = a1 + RC;
////	}
////	else
////	{
////		*s1 = a1;
////	}
//
////	if (a1 < *s1 || a1 < *y0)
////		{
////		c++;
//	*s2 += (a1 < *s1 || a1 < *y0);
//	*s1 = a1 + (a1 < *s1 || a1 < *y0) * RC;
////		}
////		else
////		{
////			*s1 = a1;
////		}
//
////	if (*s1 >= R)
////	{
////		*s1 = *s1 - R;
////		*s2++;
////	}
//
////			if
////				{
//	*s2 += (*s1 >= R);
//	*s1 -= R * (*s1 >= R);
//
////				}
//
////	a1 = *s1 + *y0;
////	if (a1 < *s1 || a1 < *y0)
////	{
////		//			c++;
////		*s2++;
////		*s1 = a1 + RC;
////	}
////	else
////	{
////		*s1 = a1;
////	}
////
////	if (*s1 >= R)
////	{
////		*s1 = *s1 - R;
////		*s2++;
////	}
////#####################################
//	*y0 = 0;
//	*x0 = c1;
//
////	if (c1 >= R) //-> v15=1>0
////	{
////		*y0 = 1;
////		*x0 = *x0 - R;
////	}
//
////	if  //-> v15=1>0
//	{
//		*y0 = (c1 >= R);
//		*x0 -= (c1 >= R) * R;
//	}
//
//	b0 = *s1 + *x0;  //[l,h,c]+[0,v14,v15]
//	*s2 += *y0;
////	if (b0 < *s1 || b0 < *x0)
////	{
////		*s2++;
////		*s1 = b0 + RC;
////	}
////	else
////	{
////		*s1 = b0;
////	}
//
////	if (b0 < *s1 || b0 < *x0)
//	{
//		*s2 += (b0 < *s1 || b0 < *x0);
//		*s1 = b0 + (b0 < *s1 || b0 < *x0) * RC;
//	}
////		else
////		{
////			*s1 = b0;
////		}
//
////	if (*s1 >= R)
////	{
////		*s1 = *s1 - R;
////		*s2++;
////	}
//
////		if (*s1 >= R)
//	{
////				*s1 = *s1 - R;
//		*s2 += (*s1 >= R);
//		*s1 -= (*s1 >= R) * R;
//	}
////#################
////[l,h,c]
//	a1 = c1 >> 29;
//	a0 = c1 - ((c1 >> 29) << 29);
//	if (a0 >= a1)
//	{
//		b1 = (a0 - a1) >> 29;
//		b0 = (a0 - a1) - (b1 << 29);
//		*y0 = b1 + a1;
//		if (b0 >= b1)
//		{
//			*x0 = (b0 - b1) << 34;
//		}
//		else
//		{
//			*x0 = R - (b1 << 34) + (b0 << 34);
//			*y0--;
//		}
//	}
//	else
//	{
//		b1 = (a1 - a0) >> 29;
//		b0 = (a1 - a0) - (b1 << 29);
//		*y0 = a1 - b1;
//		if (b0 >= b1)
//		{
//			*x0 = R - ((b0 - b1) << 34);
//			*y0--;
//		}
//		else
//		{
//			*x0 = (b1 - b0) << 34;
//		}
//	}
//
//	a0 = 0;
////	if (*s0 >= *x0)
////	{
////		*s0 = *s0 - *x0;
////	}
////	else
////	{
////		*s0 = R - *x0 + *s0;
////		a0 = 1;
////	}
//	a0 = (*s0 < *x0);
//	*s0 = *s0 - *x0 + (*s0 < *x0) * R;
//
////	if (*s1 < a0 + *y0)
////	{
////		*s2--;
////		*s1 = R - a0 - *y0 + *s1;
////	}
////	else
////	{
////		*s1 = *s1 - a0 - *y0;
////	}
//	*s2 -= (*s1 < a0 + *y0);
//	*s1 = *s1 - a0 - *y0 + (*s1 < a0 + *y0) * R;

}
/**********************************************/
__device__ void
device_p3_mulLong_register4_4 (usfixn64* __restrict__ x0, usfixn64* __restrict__ x1,
										 usfixn64*__restrict__ y0, usfixn64* __restrict__ y1,
										 usfixn64 * __restrict__ s0, usfixn64 *__restrict__ s1,
										 usfixn64 *__restrict__ s2)
{

//	usfixn64 m0, m1, m2;
//	mulLong_register2_ptr(x0, y0, s0, s1, s2);
//	mulLong_register2_ptr(x1, y1, &m0, &m1, &m2);
//	smallAdd3_modified(s0, s1, s2, &m0, &m1, &m2);
//
////	mulLong_register2_ptr(x2, y2, &m0, &m1, &m2);
////	smallAdd3_modified(s0, s1, s2, &m0, &m1, &m2);
////
////	mulLong_register2_ptr(x3, y3, &m0, &m1, &m2);
////	smallAdd3_modified(s0, s1, s2, &m0, &m1, &m2);
}

/**********************************************/
//store in l0, h0, c0
__device__ void
device_p3_smallAdd (usfixn64 *l0, usfixn64 *h0, short *c0, usfixn64 *l1, usfixn64 *h1,
					usfixn64 *c1)
{
	short c = 0;
	usfixn64 s = 0;

	s = *l0 + *l1;
	s < *l0 || s < *l1 ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}

	*l0 = s;
	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
	s = *h0 + *h1;
	s < *h0 || s < *h1 ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*h0 = s;
	*c0 = *c0 + (short) *c1 + c;
}

/**********************************************/
//change name to smallAdd_ptx_v0
__device__ void
device_p3_smallAdd_plain (usfixn64 &l0, usfixn64 &h0, usfixn64 &c0, usfixn64 &l1,
								usfixn64 &h1, usfixn64 &c1)
{
//	usfixn64 c = 0;
//	short c = 0;
	usfixn16 c = 0;
	usfixn64 s = 0;

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.u16 %1,0x00,0x00;\n\t"
			"}"
			:"=l"(s),"=h"(c)
			:"l"(l0),"l"(l1));
//	s = l0 + l1;
//	s < l0 || s < l1 ? c = 1 : c = 0;

//	c > 0 ? s = s + RC : s = s;
	if (c > 0)
		s += RC;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}

	l0 = s;
	h1 = h1 + c;  //h1<r<2^64-1. This means no overflow
//	s = h0 + h1;
//	s < h0 || s < h1 ? c = 1 : c = 0;

	asm("{\n\t"
			"add.cc.u64 %0,%2,%3;\n\t"
			"addc.u16 %1,0x00,0x00;\n\t"
			"}"
			:"=l"(s),"=h"(c)
			:"l"(h0),"l"(h1));

//very risky!
//	asm("{\n\t"
//				"add.u64 %0,%2,%3;\n\t"
//				"addc.u16 %1,0,0;\n\t"
//				"}"
//				:"=l"(s),"=h"(c):"l"(h0),"l"(h1));

//	c > 0 ? s = s + RC : s = s;
	if (c > 0)
		s += RC;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	h0 = s;
	c0 = c0 + c1 + c;
}

/**********************************************/
__device__ void
device_p3_smallAdd_plain_working (usfixn64 &l0, usfixn64 &h0, usfixn64 &c0, usfixn64 &l1,
												usfixn64 &h1, usfixn64 &c1)
{
//	usfixn64 c = 0;
	short c = 0;
	usfixn64 s = 0;

	s = l0 + l1;
	s < l0 || s < l1 ? c = 1 : c = 0;

//	asm("{\n\t"
//				"add.u64 %0,%2,%3;\n\t"
//				"addc.u16 %1,0,0;\n\t"
//				"}"
//				:"=l"(s),"=h"(c):"l"(l0),"l"(l1));

	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}

	l0 = s;
	h1 = h1 + c;  //h1<r<2^64-1. This means no overflow
	s = h0 + h1;
	s < h0 || s < h1 ? c = 1 : c = 0;

//very risky!
//	asm("{\n\t"
//				"add.u64 %0,%2,%3;\n\t"
//				"addc.u16 %1,0,0;\n\t"
//				"}"
//				:"=l"(s),"=h"(c):"l"(h0),"l"(h1));

	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	h0 = s;
	c0 = c0 + c1 + c;
}
/**********************************************/
//(xm[0:2],c)=xm[0:2]+ym[0:2]+c
__device__ __inline__ void
device_p3_smallAdd2_plain (usfixn64 *xm, usfixn64*ym, short & c)
{
////	unsigned short c = 0;
//	short pos1;
//	usfixn64 num1, num2;
//
////	usfixn64 um[8];
//
//	num1 = 0;
//	num2 = R - 1;
//	short i;
//
//	for (i = 0; i < 2; i++)
//	{
//		num1 = xm[i] + ym[i] + c;
//		if (num1 < xm[i] || num1 < ym[i]) //there is overflow/truncation
//		{
//			xm[i] = num1 + RC;
//			c = 1;
//		}
//		else if (num1 >= R)
//		{
//			c = 1;
//			xm[i] = num1 - R;
//		}
//		else
//		{
//			xm[i] = num1;
//			c = 0;
//		}
//	}

//	usfixn64 c = 0;
	c = 0;
	usfixn64 s = 0;

	s = xm[0] + ym[0];
	s < xm[0] || s < ym[0] ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}

	xm[0] = s;
	ym[1] = ym[1] + c;  //h1<r<2^64-1. This means no overflow
	s = xm[1] + ym[1];
	s < xm[1] || s < ym[1] ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	xm[1] = s;
	xm[2] = xm[2] + ym[2] + c;
}

/**********************************************/
//store in l0, h0, c0 ->different from "device_p3_smallAdd" just in type of input arguments -> all are usfixn64
__device__ __inline__ void
smallAdd2 (usfixn64 * __restrict__ l0, usfixn64 * __restrict__ h0,
					 usfixn64 * __restrict__ c0, usfixn64 * __restrict__ l1,
					 usfixn64 * __restrict__ h1, usfixn64 * __restrict__ c1)
{
	short c = 0;
	usfixn64 s = 0;

	s = *l0 + *l1;
	s < *l0 || s < *l1 ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}

	*l0 = s;
	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
	s = *h0 + *h1;
	s < *h0 || s < *h1 ? c = 1 : c = 0;
	c > 0 ? s = s + RC : s = s;

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*h0 = s;
	*c0 = *c0 + *c1 + c;
}

/**********************************************/

__device__ void
smallAdd3_modified (usfixn64 * __restrict__ l0, usfixn64 * __restrict__ h0,
										usfixn64 * __restrict__ c0,
										const usfixn64 * __restrict__ l1,
										const usfixn64 * __restrict__ h1,
										const usfixn64 * __restrict__ c1)
{
	short c = 0;
	usfixn64 s = 0;

	s = *l0 + *l1;
	if (s < *l0 || s < *l1)
	{
		c = 1;
		s = s + RC;
	}

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*l0 = s;

	*c0 = *c0 + *c1;
//	s = *h0 + *h1;
//	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
	s = *h0 + *h1 + c;

//	c = 0;
	if (s < *h0 || s < (*h1 + c))
	{
//		c = 1;

//		s = s + RC;
		s = s + RC;
		*c0++;
	}
//	else
//	{
//		 c = 0;
//	}

	*h0 = s;
//	if (s >= R)
//	{
//		s = s - R;
////		c = 1;
//		*c0++;
//	}
	if (*h0 >= R)
	{
		*h0 = *h0 - R;
		//		c = 1;
		*c0++;
	}
//	*c0 = *c0 + *c1 + c;
//	*c0 +=c;
//	*h0 = s;
}

/**********************************************/

//__device__  void smallAdd3_modified(usfixn64* __restrict__ l0, usfixn64* __restrict__ h0, usfixn64* __restrict__ c0,
//		const usfixn64 * __restrict__ l1, const usfixn64 * __restrict__ h1, const usfixn64 * __restrict__ c1)
//{
//	short c = 0;
//	usfixn64 s = 0;
//
//	s = *l0 + *l1;
//	if (s < *l0 || s < *l1)
//	{
//		c = 1;
//		s = s + RC;
//	}
//
//	if (s >= R)
//	{
//		s = s - R;
//		c = 1;
//	}
//	*l0 = s;
//
//	*c0 = *c0 + *c1;
////	s = *h0 + *h1;
//	//	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
//	s = *h0 + *h1 + c;
//
////	c = 0;
//	if (s < *h0 || s < (*h1 + c))
//	{
////		c = 1;
//
////		s = s + RC;
//		s = s + RC;
//		*c0++;
//	}
////	else
////	{
////		 c = 0;
////	}
//
//	*h0 = s;
////	if (s >= R)
////	{
////		s = s - R;
//////		c = 1;
////		*c0++;
////	}
//	if (*h0 >= R)
//	{
//		*h0 = *h0 - R;
//		//		c = 1;
//		*c0++;
//	}
////	*c0 = *c0 + *c1 + c;
////	*c0 +=c;
////	*h0 = s;
//}
/**********************************************/

//
//__device__ __inline__ void smallAdd4_modified(usfixn64 * __restrict__ l0,
//		usfixn64 * __restrict__ h0, usfixn64 * __restrict__ c0, const usfixn64 l1,
//		const usfixn64 h1, const usfixn64 c1)
//{
//	short c = 0;
//	usfixn64 s = 0;
//
//	s = *l0 + l1;
//	if (s < *l0 || s < l1)
//	{
//		c = 1;
//		s = s + RC;
//	}
//
//	if (s >= R)
//	{
//		s = s - R;
//		c = 1;
//	}
//	*l0 = s;
//
//	*c0 = *c0 + c1;
////	s = *h0 + *h1;
////	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
//	s = *h0 + h1 + c;
//
////	c = 0;
//	if (s < *h0 || s < (h1 + c))
//	{
////		c = 1;
//
////		s = s + RC;
//		s = s + RC;
//		*c0++;
//	}
////	else
////	{
////		 c = 0;
////	}
//
//	*h0 = s;
////	if (s >= R)
////	{
////		s = s - R;
//////		c = 1;
////		*c0++;
////	}
//	if (*h0 >= R)
//	{
//		*h0 = *h0 - R;
//		//		c = 1;
//		*c0++;
//	}
////	*c0 = *c0 + *c1 + c;
////	*c0 +=c;
////	*h0 = s;
//}
/**********************************************/

//
////this code is supposed to do subtraction, but actually does addition, internals should be changed.
//__device__ __inline__ void smallSub2_modified(usfixn64 * __restrict__ l0,
//		usfixn64 * __restrict__ h0, usfixn64 * __restrict__ c0,
//		usfixn64 * __restrict__ l1, usfixn64 * __restrict__ h1,
//		usfixn64 * __restrict__ c1)
//{
//	short c = 0;
//	usfixn64 s = 0;
//
//	s = *l0 + *l1;
//	if (s < *l0 || s < *l1)
//	{
//		c = 1;
//		s = s + RC;
//	}
//
//	if (s >= R)
//	{
//		s = s - R;
//		c = 1;
//	}
//	*l0 = s;
//
//	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
//	s = *h0 + *h1;
//	c = 0;
//	if (s < *h0 || s < *h1)
//	{
//		c = 1;
//		s = s + RC;
//	}
////	else
////	{
////		 c = 0;
////	}
//
//	if (s >= R)
//	{
//		s = s - R;
//		c = 1;
//	}
//	*h0 = s;
//	*c0 = *c0 + *c1 + c;
//}
/**********************************************/

//store in l0, h0, c0 ->different from "device_p3_smallAdd" just in type of input arguments -> all are usfixn64
__device__ void
smallAdd4_v0 (usfixn64 * __restrict__ l0, usfixn64 * __restrict__ h0,
							usfixn64 * __restrict__ c0, usfixn64 * __restrict__ l1,
							usfixn64 * __restrict__ h1, usfixn64 * __restrict__ c1,
							usfixn64 * __restrict__ l2, usfixn64 * __restrict__ h2,
							usfixn64 * __restrict__ c2, usfixn64 * __restrict__ l3,
							usfixn64 * __restrict__ h3, usfixn64 * __restrict__ c3)
{
	short cSum0 = 0;
	short cSum2 = 0;
	usfixn64 s0 = 0;
	usfixn64 s2 = 0;

	s0 = *l0 + *l1;
	s2 = *l2 + *l3;
	if (s0 < *l0 || s0 < *l1)
	{
		cSum0 = 1;
		s0 = s0 + RC;
	}

	if (s0 >= R)
	{
		s0 = s0 - R;
		cSum0 = 1;
	}

	if (s2 < *l2 || s2 < *l3)
	{
		cSum2 = 1;
		s2 = s2 + RC;
	}

	if (s2 >= R)
	{
		s2 = s2 - R;
		cSum2 = 1;
	}
	*l0 = s0;
	*l2 = s2;

	*h1 = *h1 + cSum0;  //h1<r<2^64-1. This means no overflow
	*h3 = *h3 + cSum2;  //h3<r<2^64-1. This means no overflow

	s0 = *h0 + *h1;
	cSum0 = 0;

	if (s0 < *h0 || s0 < *h1)
	{
		cSum0 = 1;
		s0 = s0 + RC;
	}
	if (s0 >= R)
	{
		s0 = s0 - R;
		cSum0 = 1;
	}
	s2 = *h2 + *h3;
	cSum2 = 0;

	if (s2 < *h2 || s2 < *h3)
	{
		cSum2 = 1;
		s2 = s2 + RC;
	}

	if (s2 >= R)
	{
		s2 = s2 - R;
		cSum2 = 1;
	}

	*h0 = s0;
	*h2 = s2;
	*c0 = *c0 + *c1 + cSum0;
	*c2 = *c2 + *c3 + cSum2;

//############################
	s0 = *l0 + *l2;
	if (s0 < *l0 || s0 < *l2)
	{
		cSum0 = 1;
		s0 = s0 + RC;
	}

	if (s0 >= R)
	{
		s0 = s0 - R;
		cSum0 = 1;
	}
	*l0 = s0;

	*h2 = *h2 + cSum0;  //h2<r<2^64-1. This means no overflow
	s0 = *h0 + *h2;
	cSum0 = 0;
	if (s0 < *h0 || s0 < *h2)
	{
		cSum0 = 1;
		s0 = s0 + RC;
	}

	if (s0 >= R)
	{
		s0 = s0 - R;
		cSum0 = 1;
	}
	*h0 = s0;
	*c0 = *c0 + *c2 + cSum0;

}

/**********************************************/

__device__ __inline__ void
smallAdd3_v0 (usfixn64 * __restrict__ l0, usfixn64 * __restrict__ h0,
							usfixn64 * __restrict__ c0, usfixn64 * __restrict__ l1,
							usfixn64 * __restrict__ h1, usfixn64 * __restrict__ c1,
							usfixn64 * __restrict__ l2, usfixn64 * __restrict__ h2,
							usfixn64 * __restrict__ c2)
{
	short c = 0;
	usfixn64 s = 0;

	s = *l0 + *l1;
	if (s < *l0 || s < *l1)
	{
		c = 1;
		s = s + RC;
	}

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*l0 = s;

	*h1 = *h1 + c;  //h1<r<2^64-1. This means no overflow
	s = *h0 + *h1;
	c = 0;
	if (s < *h0 || s < *h1)
	{
		c = 1;
		s = s + RC;
	}
//	else
//	{
//		 c = 0;
//	}

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*h0 = s;
	*c0 = *c0 + *c1 + c;
//###########################
	s = *l0 + *l2;
	if (s < *l0 || s < *l2)
	{
		c = 1;
		s = s + RC;
	}

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*l0 = s;

	*h2 = *h2 + c;  //h2<r<2^64-1. This means no overflow
	s = *h0 + *h2;
	c = 0;
	if (s < *h0 || s < *h2)
	{
		c = 1;
		s = s + RC;
	}
//	else
//	{
//		 c = 0;
//	}

	if (s >= R)
	{
		s = s - R;
		c = 1;
	}
	*h0 = s;
	*c0 = *c0 + *c2 + c;

}

/**********************************************/
// Multiplication in Z/(R^8 + 1)Z (us = xs * ys)
//__device__ void bigMult(usfixn64 * xs, usfixn64 *ys, usfixn64 *us)
__device__ void
bigMult (usfixn64 * __restrict__ xs, usfixn64 * __restrict__ ys,
				 usfixn64 * __restrict__ us)
{
	usfixn64 ts1[8];
	usfixn64 ts2[8];
	usfixn64 rs[3];
	short c0, c1, c2, c3, c4, c5, c6, c7;
	usfixn64 l0, l1, l2, l3, l4, l5, l6, l7, h0, h1, h2, h3, h4, h5, h6, h7;

//x0*y0
	device_p3_mulLong (xs[0], ys[0], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];

//x0*y1+x1*y0
	device_p3_mulLong (xs[0], ys[1], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[0], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//x0*y2+x1*y1+x2*y0
	device_p3_mulLong (xs[0], ys[2], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[1], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[0], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//x0*y3+x1*y2+x2*y1+x3*y0
	device_p3_mulLong (xs[0], ys[3], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[2], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[1], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[0], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
	device_p3_mulLong (xs[0], ys[4], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[3], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[2], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[1], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[0], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
	device_p3_mulLong (xs[0], ys[5], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[4], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[3], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[2], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[1], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[0], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
	device_p3_mulLong (xs[0], ys[6], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[5], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[4], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[3], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[2], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[1], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[0], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);

//x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
	device_p3_mulLong (xs[0], ys[7], rs);
	l7 = rs[0];
	h7 = rs[1];
	c7 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[6], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[5], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[4], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[3], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[2], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[1], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[0], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);

// (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
//	ts1[0] = l0;
//	ts1[1] = h0;
//	ts1[2] = c0;
//	ts1[3] = c1;
//	ts1[4] = c2;
//	ts1[5] = c3;
//	ts1[6] = c4;
//	ts1[7] = c5;
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = l7;
////	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);

//###############################################
	ts1[0] = 0;
	ts1[1] = 0;
	ts1[2] = c0;
	ts1[3] = c1;
	ts1[4] = c2;
	ts1[5] = c3;
	ts1[6] = c4;
	ts1[7] = c5;

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;

	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = l7;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
//###############################################
	ts2[0] = c6;
	ts2[1] = c7;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2, um);
	device_p3_bigSub (ts1, ts2, ts1);
	ts2[0] = h7;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

//(x7*y7)r^6
	device_p3_mulLong (xs[7], ys[7], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];

//(x6*y7+x7*y6)r^5
	device_p3_mulLong (xs[6], ys[7], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[7], ys[6], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//(x5*y7+x6*y6+x7*y5)r^4
	device_p3_mulLong (xs[5], ys[7], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[6], ys[6], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[5], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//(x4*y7+x5*y6+x6*y5+x7*y4)r^3
	device_p3_mulLong (xs[4], ys[7], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[5], ys[6], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[5], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[4], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
	device_p3_mulLong (xs[3], ys[7], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[4], ys[6], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[5], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[4], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[3], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
	device_p3_mulLong (xs[2], ys[7], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[3], ys[6], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[5], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[4], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[3], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[2], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
	device_p3_mulLong (xs[1], ys[7], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];
	device_p3_mulLong (xs[2], ys[6], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[5], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[4], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[3], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[2], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[1], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);

//	if(c6!=0)
//		printf("not zero \n");
////#####################################################
////(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
//	ts2[0] = l0;
//	ts2[1] = h0;
//	ts2[2] = c0;
//	ts2[3] = c1;
//	ts2[4] = c2;
//	ts2[5] = c3;
//	ts2[6] = c4;
//	ts2[7] = c5;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = 0;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = c6;
//	ts2[1] = 0;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//
//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
//#####################################################
//(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
	ts2[0] = 0;
	ts2[1] = 0;
	ts2[2] = c0;
	ts2[3] = c1;
	ts2[4] = c2;
	ts2[5] = c3;
	ts2[6] = c4;
	ts2[7] = c5;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

//c7=0, h7=0, l7=0 as there are no subtraction for coefficient of r^7
	ts2[0] = c6;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);

	us[0] = ts1[0];
	us[1] = ts1[1];
	us[2] = ts1[2];
	us[3] = ts1[3];
	us[4] = ts1[4];
	us[5] = ts1[5];
	us[6] = ts1[6];
	us[7] = ts1[7];
//#####################################################
}

/**********************************************/
__device__ void
device_p3_bigMult_plain (usfixn64 * xs, usfixn64 * ys)
{
	usfixn64 ts1[8];
	usfixn64 ts2[8];
	usfixn64 rs[3];
	short c0, c1, c2, c3, c4, c5, c6, c7;
	usfixn64 l0, l1, l2, l3, l4, l5, l6, l7, h0, h1, h2, h3, h4, h5, h6, h7;

//x0*y0
	device_p3_mulLong (xs[0], ys[0], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];

//x0*y1+x1*y0
	device_p3_mulLong (xs[0], ys[1], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[0], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//x0*y2+x1*y1+x2*y0
	device_p3_mulLong (xs[0], ys[2], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[1], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[0], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//x0*y3+x1*y2+x2*y1+x3*y0
	device_p3_mulLong (xs[0], ys[3], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[2], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[1], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[0], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
	device_p3_mulLong (xs[0], ys[4], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[3], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[2], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[1], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[0], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
	device_p3_mulLong (xs[0], ys[5], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[4], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[3], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[2], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[1], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[0], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
	device_p3_mulLong (xs[0], ys[6], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[5], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[4], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[3], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[2], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[1], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[0], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);

//x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
	device_p3_mulLong (xs[0], ys[7], rs);
	l7 = rs[0];
	h7 = rs[1];
	c7 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[6], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[5], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[4], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[3], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[2], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[1], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[0], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);

// (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
//	ts1[0] = l0;
//	ts1[1] = h0;
//	ts1[2] = c0;
//	ts1[3] = c1;
//	ts1[4] = c2;
//	ts1[5] = c3;
//	ts1[6] = c4;
//	ts1[7] = c5;
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = l7;
////	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);

//###############################################
	ts1[0] = 0;
	ts1[1] = 0;
	ts1[2] = c0;
	ts1[3] = c1;
	ts1[4] = c2;
	ts1[5] = c3;
	ts1[6] = c4;
	ts1[7] = c5;

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;

	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = l7;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
//###############################################
	ts2[0] = c6;
	ts2[1] = c7;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2, um);
	device_p3_bigSub (ts1, ts2, ts1);
	ts2[0] = h7;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

//(x7*y7)r^6
	device_p3_mulLong (xs[7], ys[7], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];

//(x6*y7+x7*y6)r^5
	device_p3_mulLong (xs[6], ys[7], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[7], ys[6], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//(x5*y7+x6*y6+x7*y5)r^4
	device_p3_mulLong (xs[5], ys[7], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[6], ys[6], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[5], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//(x4*y7+x5*y6+x6*y5+x7*y4)r^3
	device_p3_mulLong (xs[4], ys[7], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[5], ys[6], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[5], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[4], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
	device_p3_mulLong (xs[3], ys[7], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[4], ys[6], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[5], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[4], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[3], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
	device_p3_mulLong (xs[2], ys[7], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[3], ys[6], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[5], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[4], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[3], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[2], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
	device_p3_mulLong (xs[1], ys[7], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];
	device_p3_mulLong (xs[2], ys[6], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[5], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[4], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[3], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[2], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[1], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);

//	if(c6!=0)
//		printf("not zero \n");
////#####################################################
////(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
//	ts2[0] = l0;
//	ts2[1] = h0;
//	ts2[2] = c0;
//	ts2[3] = c1;
//	ts2[4] = c2;
//	ts2[5] = c3;
//	ts2[6] = c4;
//	ts2[7] = c5;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = 0;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = c6;
//	ts2[1] = 0;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//
//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
//#####################################################
//(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
	ts2[0] = 0;
	ts2[1] = 0;
	ts2[2] = c0;
	ts2[3] = c1;
	ts2[4] = c2;
	ts2[5] = c3;
	ts2[6] = c4;
	ts2[7] = c5;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

//c7=0, h7=0, l7=0 as there are no subtraction for coefficient of r^7
	ts2[0] = c6;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);

	xs[0] = ts1[0];
	xs[1] = ts1[1];
	xs[2] = ts1[2];
	xs[3] = ts1[3];
	xs[4] = ts1[4];
	xs[5] = ts1[5];
	xs[6] = ts1[6];
	xs[7] = ts1[7];

//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
}

/**********************************************/
__device__ void
device_p3_bigMult_plain_r1 (usfixn64 * xs, usfixn64 * ys)
{
	usfixn64 ts1[8];
	usfixn64 ts2[8];

	usfixn64 lArray[8];
	usfixn64 lArraySub[8];

	usfixn64 rs[3];
	short c0, c1, c2, c3, c4, c5, c6, c7;
	usfixn64 l0, l1, l2, l3, l4, l5, l6, l7, h0, h1, h2, h3, h4, h5, h6, h7;

//x0*y0
	device_p3_mulLong (xs[0], ys[0], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];

//x0*y1+x1*y0
	device_p3_mulLong (xs[0], ys[1], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[0], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//x0*y2+x1*y1+x2*y0
	device_p3_mulLong (xs[0], ys[2], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[1], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[0], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//x0*y3+x1*y2+x2*y1+x3*y0
	device_p3_mulLong (xs[0], ys[3], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[2], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[1], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[0], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
	device_p3_mulLong (xs[0], ys[4], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[3], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[2], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[1], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[0], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
	device_p3_mulLong (xs[0], ys[5], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[4], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[3], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[2], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[1], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[0], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
	device_p3_mulLong (xs[0], ys[6], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[5], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[4], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[3], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[2], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[1], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[0], rs);
	device_p3_smallAdd (&l6, &h6, &c6, &rs[0], &rs[1], &rs[2]);

//x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
	device_p3_mulLong (xs[0], ys[7], rs);
	l7 = rs[0];
	h7 = rs[1];
	c7 = (short) rs[2];
	device_p3_mulLong (xs[1], ys[6], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[2], ys[5], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[4], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[3], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[2], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[1], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[0], rs);
	device_p3_smallAdd (&l7, &h7, &c7, &rs[0], &rs[1], &rs[2]);

// (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
//	ts1[0] = l0;
//	ts1[1] = h0;
//	ts1[2] = c0;
//	ts1[3] = c1;
//	ts1[4] = c2;
//	ts1[5] = c3;
//	ts1[6] = c4;
//	ts1[7] = c5;
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = l7;
////	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);

//###############################################
	ts1[0] = 0;
	ts1[1] = 0;
	ts1[2] = c0;
	ts1[3] = c1;
	ts1[4] = c2;
	ts1[5] = c3;
	ts1[6] = c4;
	ts1[7] = c5;

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;

	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = l7;

	memcpy (lArray, ts2, 8 * sizeof(usfixn64));

//	device_p3_bigPrimeAdd_correct(ts1, ts2, um);
//reset ts1 to zero
	memset (ts1, 0x00, 8 * sizeof(usfixn64));
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);
//###############################################
	ts2[0] = c6;
	ts2[1] = c7;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2, um);
//	device_p3_bigSub(ts1, ts2, ts1);
	ts2[0] = h7;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
//	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);

//(x7*y7)r^6
	device_p3_mulLong (xs[7], ys[7], rs);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];

//(x6*y7+x7*y6)r^5
	device_p3_mulLong (xs[6], ys[7], rs);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong (xs[7], ys[6], rs);
	device_p3_smallAdd (&l5, &h5, &c5, &rs[0], &rs[1], &rs[2]);

//(x5*y7+x6*y6+x7*y5)r^4
	device_p3_mulLong (xs[5], ys[7], rs);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong (xs[6], ys[6], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[5], rs);
	device_p3_smallAdd (&l4, &h4, &c4, &rs[0], &rs[1], &rs[2]);

//(x4*y7+x5*y6+x6*y5+x7*y4)r^3
	device_p3_mulLong (xs[4], ys[7], rs);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong (xs[5], ys[6], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[5], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[4], rs);
	device_p3_smallAdd (&l3, &h3, &c3, &rs[0], &rs[1], &rs[2]);

//(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
	device_p3_mulLong (xs[3], ys[7], rs);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong (xs[4], ys[6], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[5], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[4], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[3], rs);
	device_p3_smallAdd (&l2, &h2, &c2, &rs[0], &rs[1], &rs[2]);

//(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
	device_p3_mulLong (xs[2], ys[7], rs);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong (xs[3], ys[6], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[5], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[4], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[3], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[2], rs);
	device_p3_smallAdd (&l1, &h1, &c1, &rs[0], &rs[1], &rs[2]);

//(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
	device_p3_mulLong (xs[1], ys[7], rs);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];
	device_p3_mulLong (xs[2], ys[6], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[3], ys[5], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[4], ys[4], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[5], ys[3], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[6], ys[2], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);
	device_p3_mulLong (xs[7], ys[1], rs);
	device_p3_smallAdd (&l0, &h0, &c0, &rs[0], &rs[1], &rs[2]);

//	if(c6!=0)
//		printf("not zero \n");
////#####################################################
////(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
//	ts2[0] = l0;
//	ts2[1] = h0;
//	ts2[2] = c0;
//	ts2[3] = c1;
//	ts2[4] = c2;
//	ts2[5] = c3;
//	ts2[6] = c4;
//	ts2[7] = c5;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = l1;
//	ts2[2] = h1;
//	ts2[3] = h2;
//	ts2[4] = h3;
//	ts2[5] = h4;
//	ts2[6] = h5;
//	ts2[7] = h6;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = 0;
//	ts2[1] = 0;
//	ts2[2] = l2;
//	ts2[3] = l3;
//	ts2[4] = l4;
//	ts2[5] = l5;
//	ts2[6] = l6;
//	ts2[7] = 0;
////	device_p3_bigSub(ts1, ts2);
//	device_p3_bigSub(ts1, ts2, ts1);
//
//	ts2[0] = c6;
//	ts2[1] = 0;
//	ts2[2] = 0;
//	ts2[3] = 0;
//	ts2[4] = 0;
//	ts2[5] = 0;
//	ts2[6] = 0;
//	ts2[7] = 0;
//	device_p3_bigPrimeAdd_correct(ts1, ts2, ts1);
//
//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
//#####################################################
//(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
	ts2[0] = 0;
	ts2[1] = 0;
	ts2[2] = c0;
	ts2[3] = c1;
	ts2[4] = c2;
	ts2[5] = c3;
	ts2[6] = c4;
	ts2[7] = c5;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;
//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = 0;

	memset (lArray, 0x00, 8 * sizeof(usfixn64));
	device_p3_bigSub (lArray, ts2, lArray);
//	memcpy(lArraySub, ts2, 8 * sizeof(usfixn64));
	memcpy (lArray, ts2, 8 * sizeof(usfixn64));

	xs[0] = lArray[0];
	xs[1] = lArray[1];
	xs[2] = lArray[2];
	xs[3] = lArray[3];
	xs[4] = lArray[4];
	xs[5] = lArray[5];
	xs[6] = lArray[6];
	xs[7] = lArray[7];

	return;

//	device_p3_bigSub(ts1, ts2);
	device_p3_bigSub (ts1, ts2, ts1);

//c7=0, h7=0, l7=0 as there are no subtraction for coefficient of r^7
	ts2[0] = c6;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	device_p3_bigPrimeAdd_correct (ts1, ts2, ts1);

//	xs[0] = ts1[0];
//	xs[1] = ts1[1];
//	xs[2] = ts1[2];
//	xs[3] = ts1[3];
//	xs[4] = ts1[4];
//	xs[5] = ts1[5];
//	xs[6] = ts1[6];
//	xs[7] = ts1[7];

	device_p3_bigSub (lArray, lArraySub, lArray);

	xs[0] = lArray[0];
	xs[1] = lArray[1];
	xs[2] = lArray[2];
	xs[3] = lArray[3];
	xs[4] = lArray[4];
	xs[5] = lArray[5];
	xs[6] = lArray[6];
	xs[7] = lArray[7];

//	us[0] = ts1[0];
//	us[1] = ts1[1];
//	us[2] = ts1[2];
//	us[3] = ts1[3];
//	us[4] = ts1[4];
//	us[5] = ts1[5];
//	us[6] = ts1[6];
//	us[7] = ts1[7];
}

/**********************************************/
__device__ void
device_p3_bigMult_plain_revised_v0 (usfixn64 * xs, usfixn64 * ys)
{
	usfixn64 ts1[8];
	usfixn64 ts2[8];
	usfixn64 rs[3];
	usfixn64 c0, c1, c2, c3, c4, c5, c6, c7;
	usfixn64 l0, l1, l2, l3, l4, l5, l6, l7, h0, h1, h2, h3, h4, h5, h6, h7;

	//x0*y0
	//	device_p3_mulLong(xs[0], ys[0], rs);
	device_p3_mulLong_2_plain_2 (xs[0], ys[0], rs[0], rs[1], rs[2]);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];

	//x0*y1+x1*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[1], rs[0], rs[1], rs[2]);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);

	//x0*y2+x1*y1+x2*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[2], rs[0], rs[1], rs[2]);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);

	//x0*y3+x1*y2+x2*y1+x3*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[3], rs[0], rs[1], rs[2]);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);

	//x0*y4+x1*y3+x2*y2+x3*y1+x4*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[4], rs[0], rs[1], rs[2]);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);

	//x0*y5+x1*y4+x2*y3+x3*y2+x4*y1+x5*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[5], rs[0], rs[1], rs[2]);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);

	//x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[6], rs[0], rs[1], rs[2]);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l6, h6, c6, rs[0], rs[1], rs[2]);

	//x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0
	device_p3_mulLong_2_plain_2 (xs[0], ys[7], rs[0], rs[1], rs[2]);
	l7 = rs[0];
	h7 = rs[1];
	c7 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[1], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[2], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[0], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l7, h7, c7, rs[0], rs[1], rs[2]);

	// (c5+h6+l7)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1-c7)*r+(l0-c6-h7)
	//	ts1[0] = l0;
	//	ts1[1] = h0;
	//	ts1[2] = c0;
	//	ts1[3] = c1;
	//	ts1[4] = c2;
	//	ts1[5] = c3;
	//	ts1[6] = c4;
	//	ts1[7] = c5;
	//
	//	ts2[0] = 0;
	//	ts2[1] = l1;
	//	ts2[2] = h1;
	//	ts2[3] = h2;
	//	ts2[4] = h3;
	//	ts2[5] = h4;
	//	ts2[6] = h5;
	//	ts2[7] = h6;
	//	device_p3_bigPrimeAdd_plain_ptx_v0(ts1, ts2, ts1);
	//	ts2[0] = 0;
	//	ts2[1] = 0;
	//	ts2[2] = l2;
	//	ts2[3] = l3;
	//	ts2[4] = l4;
	//	ts2[5] = l5;
	//	ts2[6] = l6;
	//	ts2[7] = l7;
	////	device_p3_bigPrimeAdd_plain_ptx_v0(ts1, ts2, um);
	//	device_p3_bigPrimeAdd_plain_ptx_v0(ts1, ts2, ts1);

	//###############################################
	ts1[0] = 0;
	ts1[1] = 0;
	ts1[2] = c0;
	ts1[3] = c1;
	ts1[4] = c2;
	ts1[5] = c3;
	ts1[6] = c4;
	ts1[7] = c5;

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;

	device_p3_bigPrimeAdd_plain_ptx_v0 (ts1, ts2);
	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = l7;
	//	device_p3_bigPrimeAdd_plain_ptx_v0(ts1, ts2, um);
	device_p3_bigPrimeAdd_plain_ptx_v0 (ts1, ts2);
	//###############################################
	ts2[0] = c6;
	ts2[1] = c7;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2, um);
	device_p3_bigPrimeSub_plain_ptx_v0 (ts1, ts2);
	ts2[0] = h7;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	device_p3_bigPrimeSub_plain_ptx_v0 (ts1, ts2);

	//(x7*y7)r^6
	device_p3_mulLong_2_plain_2 (xs[7], ys[7], rs[0], rs[1], rs[2]);
	l6 = rs[0];
	h6 = rs[1];
	c6 = (short) rs[2];

	//(x6*y7+x7*y6)r^5
	device_p3_mulLong_2_plain_2 (xs[6], ys[7], rs[0], rs[1], rs[2]);
	l5 = rs[0];
	h5 = rs[1];
	c5 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[7], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l5, h5, c5, rs[0], rs[1], rs[2]);

	//(x5*y7+x6*y6+x7*y5)r^4
	device_p3_mulLong_2_plain_2 (xs[5], ys[7], rs[0], rs[1], rs[2]);
	l4 = rs[0];
	h4 = rs[1];
	c4 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[6], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l4, h4, c4, rs[0], rs[1], rs[2]);

	//(x4*y7+x5*y6+x6*y5+x7*y4)r^3
	device_p3_mulLong_2_plain_2 (xs[4], ys[7], rs[0], rs[1], rs[2]);
	l3 = rs[0];
	h3 = rs[1];
	c3 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[5], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l3, h3, c3, rs[0], rs[1], rs[2]);

	//(x3*y7+x4*y6+x5*y5+x6*y4+x7*y3)r^2
	device_p3_mulLong_2_plain_2 (xs[3], ys[7], rs[0], rs[1], rs[2]);
	l2 = rs[0];
	h2 = rs[1];
	c2 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[4], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l2, h2, c2, rs[0], rs[1], rs[2]);

	//(x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2)r
	device_p3_mulLong_2_plain_2 (xs[2], ys[7], rs[0], rs[1], rs[2]);
	l1 = rs[0];
	h1 = rs[1];
	c1 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[3], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l1, h1, c1, rs[0], rs[1], rs[2]);

	//(x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1)
	device_p3_mulLong_2_plain_2 (xs[1], ys[7], rs[0], rs[1], rs[2]);
	l0 = rs[0];
	h0 = rs[1];
	c0 = (short) rs[2];
	device_p3_mulLong_2_plain_2 (xs[2], ys[6], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[3], ys[5], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[4], ys[4], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[5], ys[3], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[6], ys[2], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);
	device_p3_mulLong_2_plain_2 (xs[7], ys[1], rs[0], rs[1], rs[2]);
	device_p3_smallAdd_plain (l0, h0, c0, rs[0], rs[1], rs[2]);

	//	if(c6!=0)
	//		printf("not zero \n");
	////#####################################################
	////(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
	//	ts2[0] = l0;
	//	ts2[1] = h0;
	//	ts2[2] = c0;
	//	ts2[3] = c1;
	//	ts2[4] = c2;
	//	ts2[5] = c3;
	//	ts2[6] = c4;
	//	ts2[7] = c5;
	////	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2, ts1);
	//
	//	ts2[0] = 0;
	//	ts2[1] = l1;
	//	ts2[2] = h1;
	//	ts2[3] = h2;
	//	ts2[4] = h3;
	//	ts2[5] = h4;
	//	ts2[6] = h5;
	//	ts2[7] = h6;
	////	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2, ts1);
	//
	//	ts2[0] = 0;
	//	ts2[1] = 0;
	//	ts2[2] = l2;
	//	ts2[3] = l3;
	//	ts2[4] = l4;
	//	ts2[5] = l5;
	//	ts2[6] = l6;
	//	ts2[7] = 0;
	////	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2, ts1);
	//
	//	ts2[0] = c6;
	//	ts2[1] = 0;
	//	ts2[2] = 0;
	//	ts2[3] = 0;
	//	ts2[4] = 0;
	//	ts2[5] = 0;
	//	ts2[6] = 0;
	//	ts2[7] = 0;
	//	device_p3_bigPrimeAdd_plain_ptx_v0(ts1, ts2, ts1);
	//
	//	us[0] = ts1[0];
	//	us[1] = ts1[1];
	//	us[2] = ts1[2];
	//	us[3] = ts1[3];
	//	us[4] = ts1[4];
	//	us[5] = ts1[5];
	//	us[6] = ts1[6];
	//	us[7] = ts1[7];
	//#####################################################
	//(c5+h6)*r^7+(c4+h5+l6)*r^6+(c3+h4+l5)*r^5+(c2+h3+l4)*r^4+(c1+h2+l3)*r^3+(c0+h1+l2)*r^2+(h0+l1)*r+(l0-c6)
	ts2[0] = 0;
	ts2[1] = 0;
	ts2[2] = c0;
	ts2[3] = c1;
	ts2[4] = c2;
	ts2[5] = c3;
	ts2[6] = c4;
	ts2[7] = c5;
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	device_p3_bigPrimeSub_plain_ptx_v0 (ts1, ts2);

	ts2[0] = 0;
	ts2[1] = h0;
	ts2[2] = h1;
	ts2[3] = h2;
	ts2[4] = h3;
	ts2[5] = h4;
	ts2[6] = h5;
	ts2[7] = h6;
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	device_p3_bigPrimeSub_plain_ptx_v0 (ts1, ts2);

	ts2[0] = l0;
	ts2[1] = l1;
	ts2[2] = l2;
	ts2[3] = l3;
	ts2[4] = l4;
	ts2[5] = l5;
	ts2[6] = l6;
	ts2[7] = 0;
	//	device_p3_bigPrimeSub_plain_ptx_v0(ts1, ts2);
	device_p3_bigPrimeSub_plain_ptx_v0 (ts1, ts2);

	//c7=0, h7=0, l7=0 as there are no subtraction for coefficient of r^7
	ts2[0] = c6;
	ts2[1] = 0;
	ts2[2] = 0;
	ts2[3] = 0;
	ts2[4] = 0;
	ts2[5] = 0;
	ts2[6] = 0;
	ts2[7] = 0;
	device_p3_bigPrimeAdd_plain_ptx_v0 (ts1, ts2);

	xs[0] = ts1[0];
	xs[1] = ts1[1];
	xs[2] = ts1[2];
	xs[3] = ts1[3];
	xs[4] = ts1[4];
	xs[5] = ts1[5];
	xs[6] = ts1[6];
	xs[7] = ts1[7];

	//	us[0] = ts1[0];
	//	us[1] = ts1[1];
	//	us[2] = ts1[2];
	//	us[3] = ts1[3];
	//	us[4] = ts1[4];
	//	us[5] = ts1[5];
	//	us[6] = ts1[6];
	//	us[7] = ts1[7];
}

/**********************************************/
__device__ void
device_p3_bigPrimeMultiplication_permutated_k (short k, usfixn64 * __restrict__ xs,
																						uConstArray8_align8 & yVector,
																						usfixn64* __restrict__ lVector,
																						usfixn64 * __restrict__ hVector,
																						usfixn64* __restrict__ cVector,
																						usfixn64 * __restrict__ lVectorSub,
																						usfixn64 * __restrict__ hVectorSub,
																						usfixn64 * __restrict__ cVectorSub,
																						usfixn64 & permutationStride)
{
	short step = 0;
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;

	usfixn64 offset = tid;
	usfixn64 m0, m1, m2;
	usfixn64 lhc[3];
	usfixn64 lhcSub[3];
	short c = 0;
	short j = 0;
//	short k = 0;
//	for (k = 0; k < 8; k++)
	{
		step = 8 - k;
		lhc[0] = 0;
		lhc[1] = 0;
		lhc[2] = 0;

		lhcSub[0] = 0;
		lhcSub[1] = 0;
		lhcSub[2] = 0;

		offset = tid + (permutationStride * (7 - k));
		lhc[0] = lVector[offset];
		lhc[1] = hVector[offset];
		lhc[2] = cVector[offset];
		lhcSub[0] = lVectorSub[offset];
		lhcSub[1] = hVectorSub[offset];
		lhcSub[2] = cVectorSub[offset];

		i = k;
		if (i > 0)
//			device_shiftRight_uvector8(yVector, offset);
		{
//			for (j = 0; j < i; j++)
			device_p3_shiftRight_uConstArray8 (yVector, offset);
		}
		offset = tid;
//	for (i = 0; i < step; i++)
//			if(i<step)
		{

			i = 0;
			//##########################
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[7], m0, m1, m2);

				//				device_p3_mulLong_2_plain(xStep, y.i[7 - i], m0, m1, m2);
				//				device_p3_smallAdd2_plain(lhc, tmp, c);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_smallAdd3_plain(lhc, tmp);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=1;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[6], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=2;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[5], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=3;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[4], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=4;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[3], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=5;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=6;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[1], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=7;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[0], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
		}

		offset = tid + (permutationStride * (7 - k));
		lVector[offset] = lhc[0];
		hVector[offset] = lhc[1];
		cVector[offset] = lhc[2];
		lVectorSub[offset] = lhcSub[0];
		hVectorSub[offset] = lhcSub[1];
		cVectorSub[offset] = lhcSub[2];
	}
}

/**********************************************/
__device__ void
device_p3_bigPrimeMultiplication_permutated (usfixn64 * __restrict__ xs,
																					uConstArray8_align8 & yVector,
																					usfixn64* __restrict__ lVector,
																					usfixn64 * __restrict__ hVector,
																					usfixn64* __restrict__ cVector,
																					usfixn64 * __restrict__ lVectorSub,
																					usfixn64 * __restrict__ hVectorSub,
																					usfixn64 * __restrict__ cVectorSub,
																					usfixn64 & permutationStride)
{
	short step = 0;
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;

	usfixn64 offset = tid;
	usfixn64 m0, m1, m2;
	usfixn64 lhc[3];
	usfixn64 lhcSub[3];
	short c = 0;
	short k = 0;
	for (k = 0; k < 8; k++)
	{
		step = 8 - k;
		lhc[0] = 0;
		lhc[1] = 0;
		lhc[2] = 0;

		lhcSub[0] = 0;
		lhcSub[1] = 0;
		lhcSub[2] = 0;

		i = k;
		if (i > 0)
//			device_shiftRight_uvector8(yVector, offset);
			device_p3_shiftRight_uConstArray8 (yVector, offset);
		offset = tid;
//	for (i = 0; i < step; i++)
//			if(i<step)
		{

			i = 0;
			//##########################
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[7], m0, m1, m2);

				//				device_p3_mulLong_2_plain(xStep, y.i[7 - i], m0, m1, m2);
				//				device_p3_smallAdd2_plain(lhc, tmp, c);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_smallAdd3_plain(lhc, tmp);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=1;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[6], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=2;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[5], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=3;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[4], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=4;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[3], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=5;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=6;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[1], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=7;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[0], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
		}

		offset = tid + (permutationStride * (7 - k));
		lVector[offset] = lhc[0];
		hVector[offset] = lhc[1];
		cVector[offset] = lhc[2];
		lVectorSub[offset] = lhcSub[0];
		hVectorSub[offset] = lhcSub[1];
		cVectorSub[offset] = lhcSub[2];
	}
}

/**********************************************/
__device__ void
device_p3_bigPrimeMultiplication_permutated_v1 (usfixn64 * __restrict__ xs,
																						 uConstArray8_align8 & yVector,
																						 usfixn64* __restrict__ lVector,
																						 usfixn64 * __restrict__ hVector,
																						 usfixn64* __restrict__ cVector,
																						 usfixn64 * __restrict__ lVectorSub,
																						 usfixn64 * __restrict__ hVectorSub,
																						 usfixn64 * __restrict__ cVectorSub,
																						 usfixn64 & permutationStride)
{
	short step = 0;
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;

	usfixn64 offset = tid;
	usfixn64 m0, m1, m2;
	usfixn64 lhc[3];
	usfixn64 lhcSub[3];
	short c = 0;
	short k = 0;
	for (k = 0; k < 8; k++)
	{
		step = 8 - k;
		lhc[0] = 0;
		lhc[1] = 0;
		lhc[2] = 0;

		lhcSub[0] = 0;
		lhcSub[1] = 0;
		lhcSub[2] = 0;

		i = k;
		if (i > 0)
//			device_shiftRight_uvector8(yVector, offset);
			device_p3_shiftRight_uConstArray8 (yVector, offset);
		offset = tid;

//	for (i = 0; i < step; i++)
//			if(i<step)
		{

			i = 0;
			//##########################
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[7], m0, m1, m2);

				//				device_p3_mulLong_2_plain(xStep, y.i[7 - i], m0, m1, m2);
				//				device_p3_smallAdd2_plain(lhc, tmp, c);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_smallAdd3_plain(lhc, tmp);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=1;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[6], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)				//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=2;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[5], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=3;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[4], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}

			//			if (i < step)		//i=4;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[3], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=5;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=6;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[1], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
			//			if (i < step)		//i=7;
			{
				device_p3_mulLong_2_plain_2 (xs[offset], yVector.i[0], m0, m1, m2);
				//				device_p3_smallAdd_plain(lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i < step)		//i=0;
					device_p3_smallAdd_plain (lhc[0], lhc[1], lhc[2], m0, m1, m2);
				if (i >= step)
					device_p3_smallAdd_plain (lhcSub[0], lhcSub[1], lhcSub[2], m0, m1, m2);
				offset += permutationStride;
				i++;
			}
		}

		offset = tid + (permutationStride * (7 - k));
		lVector[offset] = lhc[0];
		hVector[offset] = lhc[1];
		cVector[offset] = lhc[2];
		lVectorSub[offset] = lhcSub[0];
		hVectorSub[offset] = lhcSub[1];
		cVectorSub[offset] = lhcSub[2];
	}
}

/**********************************************/
__device__ void
device_p3_newMult15_join_r1 (usfixn64* __restrict__ xs,
													usfixn64* __restrict__ lVector,
													usfixn64 * __restrict__ hVector,
													usfixn64* __restrict__ cVector,
													usfixn64 * __restrict__ parameters)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
	memset (v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = lVector[tid + (i) * permutationStride];
		offset += permutationStride;
	}

	v1.i[0] = hVector[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0]=0
	v1.i[0] = 0;

	v1.i[0] = cVector[tid + 6 * permutationStride];
	v1.i[1] = cVector[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0:1]=0
	v1.i[0] = 0;
	v1.i[1] = 0;

//positive h's [1:7]
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = hVector[tid + (i - 1) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0:1], v1.i[2:7] wil be filled with cVector[0:5]
	v1.i[0] = 0;
	v1.i[1] = 0;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = cVector[tid + (i - 2) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

//	############################# writing back to g-memory
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = v0.i[i];
		//		xs[offset] = v0_sub.i[i];
		//		xs[offset] = lVectorSub[offset];
		//		xs[offset] = hVector[offset];
		//		xs[offset] = cVector[offset];
		offset += permutationStride;
	}
}

/**********************************************/
__device__ void
device_p3_newMult15_join_r2 (usfixn64* __restrict__ xs,
													usfixn64* __restrict__ lVectorSub,
													usfixn64 * __restrict__ hVectorSub,
													usfixn64* __restrict__ cVectorSub,
													usfixn64 * __restrict__ parameters)
{

	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################
//
	offset = tid;
//	//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = xs[offset];
		offset += permutationStride;
	}
//#################################################
	uConstArray8_align8 v0_sub;
//	uConstArray8_align8 v1_sub;
//	uConstArray8_align8 v2_sub;

//	memset(v0_sub.i, 0x00, 8 * sizeof(usfixn64));
	memset (v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2_sub.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0_sub.i[i] = lVectorSub[tid + (i) * permutationStride];
		offset += permutationStride;
	}
	v0_sub.i[7] = 0;

//	v1.i[0] = hVectorSub[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
//	device_p3_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = cVectorSub[tid + 6 * permutationStride];

//	v1.i[1] = cVectorSub[tid + 7 * permutationStride];

//	device_p3_bigPrimeAdd_plain(v1.i, v2_sub.i);
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
//positive h's [1:7]
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = hVectorSub[tid + (i - 1) * permutationStride];
		//		v0_sub.i[i] = hVectorSub[tid + (8-i) * permutationStride];
		//		v0_sub.i[i] = hVectorSub[tid + (6+i) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0_sub.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
	v1.i[1] = 0;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = cVectorSub[tid + (i - 2) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0_sub.i, v1.i);

//#################################################
//	device_p3_bigSub_plain(v0.i, v0_sub.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
//#################################################

//	############################# writing back to g-memory
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
		offset += permutationStride;
	}
}

/**********************************************/
__device__ void
device_p3_newMult15_join_r1_mulOmega (usfixn64* __restrict__ xs,
																	 usfixn64* __restrict__ lVector,
																	 usfixn64 * __restrict__ hVector,
																	 usfixn64* __restrict__ cVector,
																	 usfixn64 * __restrict__ parameters, short po)
{
	short rp, wp;
	if (po <= 0)
	{
		return;
	}
	rp = po >> 4; //w^16=r
	wp = po - (rp << 4);
	if (wp <= 0)
		return;

	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
	memset (v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = lVector[tid + (i) * permutationStride];
		offset += permutationStride;
	}

	v1.i[0] = hVector[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0]=0
	v1.i[0] = 0;

	v1.i[0] = cVector[tid + 6 * permutationStride];
	v1.i[1] = cVector[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0:1]=0
	v1.i[0] = 0;
	v1.i[1] = 0;

//positive h's [1:7]
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = hVector[tid + (i - 1) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//just clear v1.i[0:1], v1.i[2:7] wil be filled with cVector[0:5]
	v1.i[0] = 0;
	v1.i[1] = 0;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = cVector[tid + (i - 2) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

//	############################# writing back to g-memory
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = v0.i[i];
		//		xs[offset] = v0_sub.i[i];
		//		xs[offset] = lVectorSub[offset];
		//		xs[offset] = hVector[offset];
		//		xs[offset] = cVector[offset];
		offset += permutationStride;
	}
}

/**********************************************/
__device__ void
device_p3_newMult15_join_r2_mulOmega (usfixn64* __restrict__ xs,
																	 usfixn64* __restrict__ lVectorSub,
																	 usfixn64 * __restrict__ hVectorSub,
																	 usfixn64* __restrict__ cVectorSub,
																	 usfixn64 * __restrict__ parameters, short po)
{

	short rp, wp;
	if (po <= 0)
	{
		return;
	}
	rp = po >> 4; //w^16=r
	wp = po - (rp << 4);
	if (wp <= 0)
		return;

	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################
//
	offset = tid;
//	//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = xs[offset];
		offset += permutationStride;
	}
//#################################################
	uConstArray8_align8 v0_sub;
//	uConstArray8_align8 v1_sub;
//	uConstArray8_align8 v2_sub;

//	memset(v0_sub.i, 0x00, 8 * sizeof(usfixn64));
	memset (v1.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v2_sub.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0_sub.i[i] = lVectorSub[tid + (i) * permutationStride];
		offset += permutationStride;
	}
	v0_sub.i[7] = 0;

//	v1.i[0] = hVectorSub[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
//	device_p3_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = cVectorSub[tid + 6 * permutationStride];

//	v1.i[1] = cVectorSub[tid + 7 * permutationStride];

//	device_p3_bigPrimeAdd_plain(v1.i, v2_sub.i);
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
//positive h's [1:7]
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = hVectorSub[tid + (i - 1) * permutationStride];
		//		v0_sub.i[i] = hVectorSub[tid + (8-i) * permutationStride];
		//		v0_sub.i[i] = hVectorSub[tid + (6+i) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0_sub.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
	v1.i[1] = 0;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = cVectorSub[tid + (i - 2) * permutationStride];
		offset += permutationStride;
	}
//	device_p3_bigPrimeAdd_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeAdd_plain_inPlace (v0_sub.i, v1.i);

//#################################################
//	device_p3_bigSub_plain(v0.i, v0_sub.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
//#################################################

//	############################# writing back to g-memory
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
		offset += permutationStride;
	}
}

/**********************************************/

/**********************************************/
//# subtraction modulo r
//device_p3_subR:=proc(x,y,c0,r)
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
device_p3_subR (usfixn64 * __restrict__ xm, usfixn64 * __restrict__ ym)
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

/**********************************************/
//# addition modulo r
//device_p3_addR:=proc(x,y,c0,r)
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
device_p3_addR (usfixn64 *__restrict__ xm, usfixn64 *__restrict__ ym)
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
device_p3_add128 (usfixn64 x0, usfixn64 x1, usfixn64 y0, usfixn64 y1, usfixn64 c0,
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
__device__ void
device_p3_sub128 (usfixn64 x0, usfixn64 x1, usfixn64 y0, usfixn64 y1, usfixn64 & u0,
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
device_p3_sub192_v0 (usfixn64 &x0, usfixn64 &x1, usfixn64 &x2, usfixn64 &y0, usfixn64 &y1,
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
device_p3_sub192 (usfixn64 &x0, usfixn64 &x1, usfixn64 &x2, usfixn64 &y0, usfixn64 &y1,
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
__device__ __inline__ void
device_p3_divR (usfixn64 &xs0, usfixn64 &xs1, usfixn64 &q, usfixn64 &m)
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
				device_p3_sub128 (s0, s1, m0, 0, s0, s1);
				f2 = -1;
			}
			else
			{
//	      # s:=irem(s,(r-r0))-q0*r0;
//	      # s:=m0-q0*r0;
				device_p3_sub128 (m0, 0, s0, s1, s0, s1);
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
//device_p3_toLHC:=proc (s,r)
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
//  d:=device_p3_divR(t0,r);
//  t0q:=d[1];
//  t0r:=d[2];
//
//  # t1q:=iquo(t1,r);
//  # t1r:=irem(t1,r);
//  d:=device_p3_divR(t1,r);
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
//  p1,tmp:=device_p3_addR(p1,t1q,tmp,r);
//  p2:=p2+tmp;
//
//  tmp:=0;
//  p1,tmp:=device_p3_subR(p1,t0q,tmp,r);
//  p2:=p2+tmp;
//
//  tmp:=0;
//  # p1:=p1-iquo(t0,r)-t2;
//  p1,tmp:=device_p3_subR(p1,t2,tmp,r);
//  p2:=p2+tmp;
//  tmp:=0;
//  # print(p1,"after");
//  # p1:=2*s1-t0q-t2+t1q;
//  # print(p1-r,"before");
//  # p1:=p1-iquo(t0,r)-t2;
//  # p0:=p0-t0r+t1r;
//
//  tmp:=0;
//  p0,tmp:=device_p3_addR(p0,t1r,tmp,r);
//  # p1,tmp:=device_p3_addR(p1,tmp,0,r);
//  # p2:=p2+tmp;
//  p1:=p1+tmp;
//
//  tmp:=0;
//  p0,tmp:=device_p3_subR(p0,t0r,tmp,r);
//  # p1,tmp:=device_p3_addR(p1,tmp,0,r);
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
__device__ void
device_p3_toLHC_v0 (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
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
		device_p3_divR (t00, t01, t0q, t0r);

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

//	device_p3_divR (t0, t0q, t0r);
//	device_p3_divR (t1, t1q, t1r);
	device_p3_divR (t00, t01, t0q, t0r);
	device_p3_divR (t10, t11, t1q, t1r);

//	printf("tqs , %lu, %lu, %lu, %lu\n",t0q,t0r,t1q,t1r); //correct up to this point

//	tmp = 0;
//	device_p3_addR (p1, t1q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;

//	printf("ps before final %lu, %lu, %lu \n",p0,p1,p2);
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t1q;
	device_p3_bigPrimeAdd_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t0q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t0q;
	device_p3_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t2, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t2;
	device_p3_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

	a1[1] = p0;
	a1[2] = p1;
	a1[3] = p2;
	a2[1] = t1r;

	device_p3_bigPrimeAdd_plain (a1, a2);
	a2[1] = t0r;
	device_p3_bigSub_plain (a1, a2);
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
__device__ void
device_p3_toLHC_v1 (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
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

		device_p3_divR (t00, t01, t0q, t0r);

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

//	device_p3_divR (t0, t0q, t0r);
//	device_p3_divR (t1, t1q, t1r);
	device_p3_divR (t00, t01, t0q, t0r);
	device_p3_divR (t10, t11, t1q, t1r);

//	printf("tqs , %lu, %lu, %lu, %lu\n",t0q,t0r,t1q,t1r); //correct up to this point

//	tmp = 0;
//	device_p3_addR (p1, t1q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;

//	printf("ps before final %lu, %lu, %lu \n",p0,p1,p2);
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t1q;
	device_p3_bigPrimeAdd_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t0q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t0q;
	device_p3_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t2, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t2;
	device_p3_bigSub_plain (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

	a1[1] = p0;
	a1[2] = p1;
	a1[3] = p2;
	a2[1] = t1r;

	device_p3_bigPrimeAdd_plain (a1, a2);
	a2[1] = t0r;
	device_p3_bigSub_plain (a1, a2);
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
device_p3_toLHC (usfixn64&s0, usfixn64 & s1, usfixn64 &s2, usfixn64 &p0, usfixn64&p1,
			 usfixn64&p2)
{
//	usfixn64 p2, p1, p0;
	usfixn64 t0, t1;
	usfixn64 t2;
	usfixn64 t00, t01;
	usfixn64 t10, t11;
	usfixn64 l = 0, h = 0, c = 0, tmp = 0;
//	usfixn64 p0,p1,p2;
//	p2 = 4 * s2;
//	p1 = 2 * s1;
	p2 = s2 << 2;
//		p1 = s1<<1;
	p0 = s0;
	usfixn64 t0q = 0, t0r = 0, t1q = 0, t1r = 0;
//	usfixn64 a1[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 a2[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
	usfixn64 a1[3] =
		{ 0, 0, 0 };
	usfixn64 a2[3] =
		{ 0, 0, 0 };
	usfixn64 r1 = 0x8000000000000000;
//	if (s1 < r1)
	{
//		p1=2*s1;
		p1 = s1 << 1;
	}
//	else
	if (s1 >= r1)
	{

//		tmp = 2;
//		asm("{\n\t"
//				"mul.lo.u64 %0,%2,%3;\n\t"
//				"mul.hi.u64 %1,%2,%3;\n\t"
//				"}"
//				:"=l"(t00),"=l"(t01)
//				:"l"(s1),"l"(tmp)
//		);

		t00 = p1;
		t01 = s1;
//		t00 = s1<<1;
		t01 >>= 63;

		device_p3_divR (t00, t01, t0q, t0r);

		p2 = p2 + t0q;
		p1 = t0r;
		t0q = 0;
		t0r = 0;
		t00 = 0;
		t01 = 0;
	}
//	usfixn64 r0;
//	r0 = 1;
//	r0 <<= 34;
//	t0 = (2 * s1 * r0);
//	usfixn64 r02 = 2 * r0;
//	usfixn64 r02 = r0<<1;
	usfixn64 r02 = 0x800000000;

//	asm("{\n\t"
//			"mul.lo.u64 %0,%2,%3;\n\t"
//			"mul.hi.u64 %1,%2,%3;\n\t"
//			"}"
//			:"=l"(t00),"=l"(t01)
//			:"l"(s1),"l"(r02)
//	);
//	t00 = s1 & ((1 << 29) - 1);
	t00 = s1 & (0x1FFFFFFF);
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

//	device_p3_divR (t0, t0q, t0r);
//	device_p3_divR (t1, t1q, t1r);
	device_p3_divR (t00, t01, t0q, t0r);
	device_p3_divR (t10, t11, t1q, t1r);

//	printf("tqs , %lu, %lu, %lu, %lu\n",t0q,t0r,t1q,t1r); //correct up to this point

//	tmp = 0;
//	device_p3_addR (p1, t1q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;

//	printf("ps before final %lu, %lu, %lu \n",p0,p1,p2);
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t1q;
//	device_p3_bigPrimeAdd_plain (a1, a2);
	device_p3_addR (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t0q, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t0q;
//	device_p3_bigSub_plain (a1, a2);
	device_p3_subR (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

//	tmp = 0;
//	device_p3_subR (p1, t2, tmp, p1, tmp);		//suspicious with c
//	p2 = p2 + tmp;
	a1[0] = p1;
	a1[1] = p2;
	a2[0] = t2;
//	device_p3_bigSub_plain (a1, a2);
	device_p3_subR (a1, a2);
	p1 = a1[0];
	p2 = a1[1];

	a1[0] = p0;
	a1[1] = p1;
	a1[2] = p2;
	a2[0] = t1r;

//	device_p3_bigPrimeAdd_plain (a1, a2);
	device_p3_addR (a1, a2);
	a2[0] = t0r;
//	device_p3_bigSub_plain (a1, a2);
	device_p3_subR (a1, a2);
	p0 = a1[0];
	p1 = a1[1];
	p2 = a1[2];

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
__global__ void
kernel_p3_mult_revised_singleThread (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
{
//	usfixn64 tid = threadIdx.x+blockIdx.x*blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 tid = 0;

	usfixn64 u[3] =
		{ 0, 0, 0 }, uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 c = 0;
	short i = 0, j = 0;
	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
			cVector[COEFFICIENT_SIZE];
	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
			cVectorSub[COEFFICIENT_SIZE];
	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];

//	usfixn64 x[COEFFICIENT_SIZE]=
//	{4122069494891546112,
//	416448630046793856,
//	3768218334765421568,
//	4312036060060413440,
//	2597503149375243776,
//	6581614685471126528,
//	3507585816910166528,
//	7624981294075398144};
//
//	usfixn64 y[COEFFICIENT_SIZE]={
//	8244138994078059520,
//	832897264388555008,
//	7536436669530843136,
//	8624072120120826880,
//	5195006303045454848,
//	3939857316907610112,
//	7015171638115300352,
//	6026590529821184000};

	usfixn64 x[COEFFICIENT_SIZE], y[COEFFICIENT_SIZE];

//		for (i=0;i<8;i++)
//		{
//			x[i]=R-2;
//			y[i]=R-2;
//		}

//	n = 1;
	for (tid = 0; tid < n; tid++)
	{
		offset = tid;
		for (i = 0; i < 8; i++)
		{
			x[i] = xs[offset];
			y[i] = ys[offset];
			offset += permutationStride;
//			printf("x[i]=%lu, y[i]=%lu\n",x[i],y[i]);
		}
		for (i = 0; i < 8; i++)
		{
			lVectorSub[i] = 0;
			hVectorSub[i] = 0;
			cVectorSub[i] = 0;
			lVector[i] = 0;
			hVector[i] = 0;
			cVector[i] = 0;
		}

		for (i = 0; i < 8; i++)
		{
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
				device_p3_oneShiftRight (y);
			for (j = 8 - i; j < 8; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);
//				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
			}
			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
						 cVectorSub[7 - i]);
//			lVectorSub[8 - i] = uSub[0];
//			hVectorSub[8 - i] = uSub[1];
//			cVectorSub[8 - i] = uSub[2];

			for (j = 0; j < 8 - i; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);
//				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
				c = u[2];
				u[2] = 0;
				device_p3_add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
			}
//			results of device_p3_add128 doesn't match with maple implementation
//			printf("u0=%lu, u1=%lu, u2=%lu \n", u[0], u[1], u[2]);
			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
//			lVector[8 - i] = u[0];
//			hVector[8 - i] = u[1];
//			cVector[8 - i] = u[2];
		}

//		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//		       lVector[0],lVector[1],lVector[2],
//		       lVector[3],lVector[4],lVector[5],
//		       lVector[6],lVector[7]);
//
//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);
		device_p3_oneShiftRight (hVectorSub);
		device_p3_oneShiftRight (cVectorSub);
		device_p3_oneShiftRight (cVectorSub);
		for (i = 0; i < 8; i++)
		{
			v0[i] = 0;
			v1[i] = 0;
		}
		v0[0] = hVectorSub[0];
		hVectorSub[0] = 0;

		v1[0] = cVectorSub[0];
		v1[1] = cVectorSub[1];
		cVectorSub[0] = 0;
		cVectorSub[1] = 0;
		device_p3_bigPrimeAdd_plain (lVectorSub, hVectorSub);

		device_p3_bigPrimeAdd_plain (lVectorSub, cVectorSub);
		device_p3_bigSub_plain (lVectorSub, v0);
//		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//						       lVectorSub[0],lVectorSub[1],lVectorSub[2],
//						       lVectorSub[3],lVectorSub[4],lVectorSub[5],
//						       lVectorSub[6],lVectorSub[7]);
		device_p3_bigSub_plain (lVectorSub, v1);
//		printf("v10= %lu, v11=%lu \n", v1[0], v1[1]);
//		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//						       lVectorSub[0],lVectorSub[1],lVectorSub[2],
//						       lVectorSub[3],lVectorSub[4],lVectorSub[5],
//						       lVectorSub[6],lVectorSub[7]);

		device_p3_oneShiftRight (hVector);
		device_p3_oneShiftRight (cVector);
		device_p3_oneShiftRight (cVector);
		for (i = 0; i < 8; i++)
		{
			v0[i] = 0;
			v1[i] = 0;
		}
		v0[0] = hVector[0];
		hVector[0] = 0;

		v1[0] = cVector[0];
		v1[1] = cVector[1];
		cVector[0] = 0;
		cVector[1] = 0;
		device_p3_bigPrimeAdd_plain (lVector, hVector);
		device_p3_bigPrimeAdd_plain (lVector, cVector);
		device_p3_bigSub_plain (lVector, v0);
		device_p3_bigSub_plain (lVector, v1);

		device_p3_bigSub_plain (lVector, lVectorSub);
		offset = tid;
		for (i = 0; i < 8; i++)
		{
			xs[offset] = lVector[i];
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread_8lhc (usfixn64 * xs, usfixn64* ys,
																usfixn64* parameters)
{
//	usfixn64 tid = threadIdx.x+blockIdx.x*blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 tid = 0;

	usfixn64 u[3] =
		{ 0, 0, 0 }, uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 c = 0;
	short i = 0, j = 0;
	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
			cVector[COEFFICIENT_SIZE];
	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
			cVectorSub[COEFFICIENT_SIZE];
	usfixn64 signVector[COEFFICIENT_SIZE];
	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 complement[COEFFICIENT_SIZE] =
		{ 0, 0, 0, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE };

	usfixn64 l, h;
//	usfixn64 x[COEFFICIENT_SIZE]=
//	{4122069494891546112,
//	416448630046793856,
//	3768218334765421568,
//	4312036060060413440,
//	2597503149375243776,
//	6581614685471126528,
//	3507585816910166528,
//	7624981294075398144};
//
//	usfixn64 y[COEFFICIENT_SIZE]={
//	8244138994078059520,
//	832897264388555008,
//	7536436669530843136,
//	8624072120120826880,
//	5195006303045454848,
//	3939857316907610112,
//	7015171638115300352,
//	6026590529821184000};

	usfixn64 x[COEFFICIENT_SIZE], y[COEFFICIENT_SIZE];

//		for (i=0;i<8;i++)
//		{
//			x[i]=R-2;
//			y[i]=R-2;
//		}

//	n = 1024;
	for (tid = 0; tid < n; tid++)
	{
		offset = tid;
		for (i = 0; i < 8; i++)
		{
			x[i] = xs[offset];
			y[i] = ys[offset];
			offset += permutationStride;
//			printf("x[i]=%lu, y[i]=%lu\n",x[i],y[i]);
		}
		for (i = 0; i < 8; i++)
		{
			lVectorSub[i] = 0;
			hVectorSub[i] = 0;
			cVectorSub[i] = 0;
			lVector[i] = 0;
			hVector[i] = 0;
			cVector[i] = 0;
		}
		complement[0] = 0;
		complement[1] = 0;
		complement[2] = 0;

		for (i = 3; i < 8; i++)
			complement[i] = R_MINUS_ONE;
		for (i = 0; i < 8; i++)
		{
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
				device_p3_oneShiftRight (y);
			for (j = 8 - i; j < 8; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);
//				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
			}

//			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
//						 cVectorSub[7 - i]);
//			lVectorSub[8 - i] = uSub[0];
//			hVectorSub[8 - i] = uSub[1];
//			cVectorSub[8 - i] = uSub[2];

			for (j = 0; j < 8 - i; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);
//				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
				c = u[2];
				u[2] = 0;
				device_p3_add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
			}
//			results of device_p3_add128 doesn't match with maple implementation
//			printf("u0=%lu, u1=%lu, u2=%lu \n", u[0], u[1], u[2]);
//			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
			usfixn64 sign = 0; //if negative =1 else =0
			device_p3_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
//			printf ("device_p3_sub192=  %lu, %lu, %lu , sign=%lu \n", u[0], u[1], u[2], sign);
			l = 0;
			h = 0;
			c = 0;
			device_p3_toLHC (u[0], u[1], u[2], l, h, c);

			for (j = 0; j < COEFFICIENT_SIZE; j++)
			{
				v0[j] = 0;
				v1[j] = 0;
			}
			v0[0] = l;
			v0[1] = h;
			v0[2] = c;
			if (sign == 1)
			{
//				device_p3_bigSubZero_3 (l, h, c);
				device_p3_bigSub_plain (v1, v0);
				l = v1[0];
				h = v1[1];
				c = v1[2];
			}
//			printf ("l=%lu, h=%lu, c=%lu \n", l, h, c);
//			printf ("===========================\n");

			lVector[7 - i] = l;
			hVector[7 - i] = h;
			cVector[7 - i] = c;
			signVector[7 - i] = sign;
//			printf("end of loop i=%d \n",i);
		}

//		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//		       lVector[0],lVector[1],lVector[2],
//		       lVector[3],lVector[4],lVector[5],
//		       lVector[6],lVector[7]);
//
//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (i > 0)
				device_p3_cyclicShift_plain (complement, 1);
			if (signVector[i] == 1)
			{
				device_p3_bigPrimeAdd_plain (lVector, complement);
			}
		}

//		device_p3_oneShiftRight (hVectorSub);
//		device_p3_oneShiftRight (cVectorSub);
//		device_p3_oneShiftRight (cVectorSub);
//		for (i = 0; i < 8; i++)
//		{
//			v0[i] = 0;
//			v1[i] = 0;
//		}
//		v0[0] = hVectorSub[0];
//		hVectorSub[0] = 0;
//
//		v1[0] = cVectorSub[0];
//		v1[1] = cVectorSub[1];
//		cVectorSub[0] = 0;
//		cVectorSub[1] = 0;
//		device_p3_bigPrimeAdd_plain (lVectorSub, hVectorSub);
//
//		device_p3_bigPrimeAdd_plain (lVectorSub, cVectorSub);
//		device_p3_bigSub_plain (lVectorSub, v0);
////		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
////						       lVectorSub[0],lVectorSub[1],lVectorSub[2],
////						       lVectorSub[3],lVectorSub[4],lVectorSub[5],
////						       lVectorSub[6],lVectorSub[7]);
//		device_p3_bigSub_plain (lVectorSub, v1);
////		printf("v10= %lu, v11=%lu \n", v1[0], v1[1]);
////		printf("l \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
////						       lVectorSub[0],lVectorSub[1],lVectorSub[2],
////						       lVectorSub[3],lVectorSub[4],lVectorSub[5],
////						       lVectorSub[6],lVectorSub[7]);

		device_p3_oneShiftRight (hVector);
		device_p3_oneShiftRight (cVector);
		device_p3_oneShiftRight (cVector);
		for (i = 0; i < 8; i++)
		{
			v0[i] = 0;
			v1[i] = 0;
		}
		v0[0] = hVector[0];
		hVector[0] = 0;

		v1[0] = cVector[0];
		v1[1] = cVector[1];
		cVector[0] = 0;
		cVector[1] = 0;
		device_p3_bigPrimeAdd_plain (lVector, hVector);
		device_p3_bigPrimeAdd_plain (lVector, cVector);
		device_p3_bigSub_plain (lVector, v0);
		device_p3_bigSub_plain (lVector, v1);

//		device_p3_bigSub_plain (lVector, lVectorSub);
		offset = tid;
		for (i = 0; i < 8; i++)
		{
			xs[offset] = lVector[i];
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_v0 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 c = 0;
	short i = 0, j = 0;
	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
			cVector[COEFFICIENT_SIZE];
	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
			cVectorSub[COEFFICIENT_SIZE];
	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 x[COEFFICIENT_SIZE], y[COEFFICIENT_SIZE];

//	for (tid = 0; tid < n; tid++)
	{
		offset = tid;
		for (i = 0; i < 8; i++)
		{
			x[i] = xs[offset];
			y[i] = ys[offset];
			offset += permutationStride;
		}
		for (i = 0; i < 8; i++)
		{
			lVectorSub[i] = 0;
			hVectorSub[i] = 0;
			cVectorSub[i] = 0;
			lVector[i] = 0;
			hVector[i] = 0;
			cVector[i] = 0;
		}

		for (i = 0; i < 8; i++)
		{
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
				device_p3_oneShiftRight (y);
			for (j = 8 - i; j < 8; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);

				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
			}
			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
						 cVectorSub[7 - i]);
			for (j = 0; j < 8 - i; j++)
			{
//				m=xs[j]*y[8-j];
				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[7-j])
				);
				c = u[2];
				u[2] = 0;
				device_p3_add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
// 				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
			}
			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);

		}

//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);
		device_p3_oneShiftRight (hVectorSub);
		device_p3_oneShiftRight (cVectorSub);
		device_p3_oneShiftRight (cVectorSub);
		for (i = 0; i < 8; i++)
		{
			v0[i] = 0;
			v1[i] = 0;
		}
		v0[0] = hVectorSub[0];
		hVectorSub[0] = 0;

		v1[0] = cVectorSub[0];
		v1[1] = cVectorSub[1];
		cVectorSub[0] = 0;
		cVectorSub[1] = 0;
		device_p3_bigPrimeAdd_plain (lVectorSub, hVectorSub);
		device_p3_bigPrimeAdd_plain (lVectorSub, cVectorSub);
		device_p3_bigSub_plain (lVectorSub, v0);
		device_p3_bigSub_plain (lVectorSub, v1);

		device_p3_oneShiftRight (hVector);
		device_p3_oneShiftRight (cVector);
		device_p3_oneShiftRight (cVector);
		for (i = 0; i < 8; i++)
		{
			v0[i] = 0;
			v1[i] = 0;
		}
		v0[0] = hVector[0];
		hVector[0] = 0;

		v1[0] = cVector[0];
		v1[1] = cVector[1];
		cVector[0] = 0;
		cVector[1] = 0;
		device_p3_bigPrimeAdd_plain (lVector, hVector);
		device_p3_bigPrimeAdd_plain (lVector, cVector);
		device_p3_bigSub_plain (lVector, v0);
		device_p3_bigSub_plain (lVector, v1);

		device_p3_bigSub_plain (lVector, lVectorSub);

		offset = tid;
		for (i = 0; i < 8; i++)
		{
			xs[offset] = lVector[i];
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_step1 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
										usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
										usfixn64* lVectorSub, usfixn64* hVectorSub,
										usfixn64* cVectorSub)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 c = 0;
	short i = 0, j = 0;
//	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
//			cVector[COEFFICIENT_SIZE];
//	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
//			cVectorSub[COEFFICIENT_SIZE];
//	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 y[COEFFICIENT_SIZE];

	usfixn64 lAdd, hAdd, cAdd;
	usfixn64 lSub, hSub, cSub;
//	for (tid = 0; tid < n; tid++)
	{
		offset = tid;
		for (i = 0; i < 8; i++)
		{
//			x[i] = xs[offset];
			y[7 - i] = ys[offset];
			offset += permutationStride;
		}
//		for (i = 0; i < 8; i++)
//		{
//			lVectorSub[i] = 0;
//			hVectorSub[i] = 0;
//			cVectorSub[i] = 0;
//			lVector[i] = 0;
//			hVector[i] = 0;
//			cVector[i] = 0;
//		}

		offset = tid + (permutationStride << 3);
		usfixn64 xOffset;
		for (i = 0; i < 8; i++)
		{
			offset -= permutationStride;
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
//				device_p3_oneShiftRight (y);
				device_p3_oneShiftLeft (y);

			xOffset = tid + 7 * permutationStride;
//			for (j = 8 - i; j < 8; j++)
//			for (j = 7; j >= 8 - i; j--)
			for (j = 7; j >= 0; j--)
			{
//				m=xs[j]*y[8-j];
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[7-j])
//				);

				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(xs[xOffset]),"l"(y[j])
				);
				xOffset -= permutationStride;

				if (j > 7 - i)
					device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
				else
					device_p3_add128 (u[0], u[1], m0, m1, u[2], u[0], u[1], u[2]);

			}

//			//			for (j = 0; j < 8 - i; j++)
//			for (j = 7 - i; j >= 0; j--)
//			{
//				//				m=xs[j]*y[8-j];
////				asm("{\n\t"
////						"mul.lo.u64 %0,%2,%3;\n\t"
////						"mul.hi.u64 %1,%2,%3;\n\t"
////						"}"
////						:"=l"(m0),"=l"(m1)
////						:"l"(x[j]),"l"(y[7-j])
////				);
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);
//				xOffset -= permutationStride;
//				c = u[2];
//				u[2] = 0;
//				device_p3_add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
//				// 				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
//			}

//			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
//						 cVectorSub[7 - i]);

//			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
//			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lSub, hSub, cSub);
//			device_p3_toLHC (u[0], u[1], u[2], lAdd, hAdd, cAdd);
//			offset = tid + (7 - i) * permutationStride;
			lVector[offset] = u[0];
			hVector[offset] = u[1];
			cVector[offset] = u[2];
			lVectorSub[offset] = uSub[0];
			hVectorSub[offset] = uSub[1];
			cVectorSub[offset] = uSub[2];
//			offset=offset-permutationStride;
		}

//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);

	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_step2 (usfixn64 * xs, usfixn64* ys, usfixn64* parameters,
										usfixn64* lVector, usfixn64* hVector, usfixn64* cVector,
										usfixn64* lVectorSub, usfixn64* hVectorSub,
										usfixn64* cVectorSub)
{

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];

	usfixn64 s0, s1, s2;
	usfixn64 l, h, c;

	offset = tid;
	short i = 0;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		s0 = lVector[offset];
		s1 = hVector[offset];
		s2 = cVector[offset];
		device_p3_toLHC (s0, s1, s2, l, h, c);
		lVector[offset] = l;
		hVector[offset] = h;
		cVector[offset] = c;
		offset += permutationStride;
	}

	offset = tid;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{

		s0 = lVectorSub[offset];
		s1 = hVectorSub[offset];
		s2 = cVectorSub[offset];
		device_p3_toLHC (s0, s1, s2, l, h, c);
		lVectorSub[offset] = l;
		hVectorSub[offset] = h;
		cVectorSub[offset] = c;
		offset += permutationStride;
	}

}

/**********************************************/
__global__ void
kernel_p3_mult_revised_step3 (usfixn64* __restrict__ xs, usfixn64* __restrict__ lVector,
										usfixn64 * __restrict__ hVector,
										usfixn64* __restrict__ cVector,
										usfixn64* __restrict__ lVectorSub,
										usfixn64 * __restrict__ hVectorSub,
										usfixn64* __restrict__ cVectorSub,
										usfixn64 * __restrict__ parameters)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	for (i = 0; i < 8; i++)
	{
//		v0.i[i]=0;
		v1.i[i] = 0;
	}

//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = lVector[offset];
		offset += permutationStride;
	}

	offset = tid + (permutationStride << 3) - permutationStride;

//	v1.i[0] = hVector[tid + 7 * permutationStride];
	v1.i[0] = hVector[offset];
//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}
//	v1.i[1] = cVector[tid + 7 * permutationStride];
//	v1.i[0] = cVector[tid + 6 * permutationStride];
	v1.i[1] = cVector[offset];
	offset -= permutationStride;
	v1.i[0] = cVector[offset];

//	device_p3_bigSub_plain(v0.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}

	v1.i[0] = 0;
//positive h's [1:7]
	offset = tid;
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
//		v1.i[i] = hVector[tid + (i - 1) * permutationStride];
		v1.i[i] = hVector[offset];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v0.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//			{
//				v1.i[i]=0;
//			}
	v1.i[0] = 0;
	v1.i[1] = 0;
	offset = tid;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
//		v1.i[i] = cVector[tid + (i - 2) * permutationStride];
		v1.i[i] = cVector[offset];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v0.i, v1.i);

//#################################################
	uConstArray8_align8 v0_sub;
//	uConstArray8_align8 v1_sub;
//	uConstArray8_align8 v2_sub;

//	memset(v0_sub.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	for (i = 0; i < 8; i++)
	{
//				v0_sub.i[i]=0;
		v1.i[i] = 0;
	}
//	memset(v2_sub.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

	offset = tid;
//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
//		v0_sub.i[i] = lVectorSub[tid + (i) * permutationStride];
		v0_sub.i[i] = lVectorSub[offset];
		offset += permutationStride;
	}
	v0_sub.i[7] = 0;

//	v1.i[0] = hVectorSub[tid + 7 * permutationStride];
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
//	device_p3_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = cVectorSub[tid + 6 * permutationStride];
//	v1.i[1] = cVectorSub[tid + 7 * permutationStride];

//	device_p3_bigPrimeAdd_plain(v1.i, v2_sub.i);
//	device_p3_bigSub_plain(v0_sub.i, v1.i);
	device_p3_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);

	offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
//positive h's [1:7]
	offset = tid;
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
//		v1.i[i] = hVectorSub[tid + (i - 1) * permutationStride];
		v1.i[i] = hVectorSub[offset];
//		v0_sub.i[i] = hVectorSub[tid + (8-i) * permutationStride];
//		v0_sub.i[i] = hVectorSub[tid + (6+i) * permutationStride];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);

	offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
	v1.i[0] = 0;
	v1.i[1] = 0;
	offset = tid;
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
//		v1.i[i] = cVectorSub[tid + (i - 2) * permutationStride];
		v1.i[i] = cVectorSub[offset];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);

//#################################################
//	device_p3_bigSub_plain(v0.i, v0_sub.i);
	device_p3_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
//#################################################

//	############################# writing back to g-memory
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
		offset += permutationStride;
	}
}

/************************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step1 (const usfixn64 * __restrict__ xs,
												 const usfixn64* __restrict__ ys, usfixn64* parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 usfixn32* __restrict__ signVector)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];

	usfixn64 u[3] =
		{ 0, 0, 0 };
	usfixn64 uSub[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	usfixn64 s = 0, c = 0;
	short i = 0, j = 0;
	usfixn64 sign = 0;
//	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
//			cVector[COEFFICIENT_SIZE];
//	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
//			cVectorSub[COEFFICIENT_SIZE];
//	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];
	usfixn64 y[COEFFICIENT_SIZE];
	usfixn64 x[COEFFICIENT_SIZE];

//	usfixn64 lAdd, hAdd, cAdd;
//	usfixn64 lSub, hSub, cSub;
//	usfixn64 m0Array[COEFFICIENT_SIZE], m1Array[COEFFICIENT_SIZE];

	usfixn64 h0, h1, h2;
	usfixn64 xOffset, xOffset_init = tid + 7 * permutationStride;
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			x[i] = xs[offset];
			y[7 - i] = ys[offset];
			offset += permutationStride;
		}
//		for (i = 0; i < 8; i++)
//		{
//			lVectorSub[i] = 0;
//			hVectorSub[i] = 0;
//			cVectorSub[i] = 0;
//			lVector[i] = 0;
//			hVector[i] = 0;
//			cVector[i] = 0;
//		}

		offset = tid + (permutationStride << 3);

		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			offset -= permutationStride;
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			uSub[0] = 0;
			uSub[1] = 0;
			uSub[2] = 0;
			if (i > 0)
//				device_p3_oneShiftRight (y);
				device_p3_oneShiftLeft (y);

			xOffset = xOffset_init;
//			for (j = 8 - i; j < 8; j++)
//			for (j = 7; j >= 8 - i; j--)
			for (j = 7; j >= 0; j--)
			{
//				m=xs[j]*y[8-j];
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(x[j]),"l"(y[7-j])
//				);

//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);

				asm("{\n\t"
						"mul.lo.u64 %0,%2,%3;\n\t"
						"mul.hi.u64 %1,%2,%3;\n\t"
						"}"
						:"=l"(m0),"=l"(m1)
						:"l"(x[j]),"l"(y[j])
				);
//
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0Array[j]),"=l"(m1Array[j])
//						:"l"(xs[xOffset]),"l"(y[j])
//				);

				xOffset -= permutationStride;

//				if (j > 7 - i)
//					device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1], uSub[2]);
//				else
//					device_p3_add128 (u[0], u[1], m0, m1, u[2], u[0], u[1], u[2]);

//			}
//
//			for (j = 7; j >= 0; j--)
//			{
				h0 = u[0];
				h1 = u[1];
				h2 = u[2];
				if (j > 7 - i)
				{
					h0 = uSub[0];
					h1 = uSub[1];
					h2 = uSub[2];
				}

				/*****************************/
//				device_p3_add128 (h0, h1, m0, m1, h2, h0, h1, h2);
//				device_p3_add128 (h0, h1, m0Array[j], m1Array[j], h2, h0, h1, h2);
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.cc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(s),"=l"(c):"l"(h0),"l"(m0));

				h0 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%1,%2;\n\t"
						"addc.cc.u64 %1,0,0;\n\t"
						"}"
						:"=l"(s),"+l"(c):"l"(m1));
				usfixn64 u1 = s;
				asm("{\n\t"
						"add.cc.u64 %0,%2,%3;\n\t"
						"addc.cc.u64 %1,%1,0;\n\t"
						"}"
						:"=l"(s),"+l"(c):"l"(h1),"l"(u1));
				h1 = s;
				h2 += c;

				/*****************************/
				if (j > 7 - i)
				{
					uSub[0] = h0;
					uSub[1] = h1;
					uSub[2] = h2;
				}
				else
				{
					u[0] = h0;
					u[1] = h1;
					u[2] = h2;
				}
			}

//			//			for (j = 0; j < 8 - i; j++)
//			for (j = 7 - i; j >= 0; j--)
//			{
//				//				m=xs[j]*y[8-j];
////				asm("{\n\t"
////						"mul.lo.u64 %0,%2,%3;\n\t"
////						"mul.hi.u64 %1,%2,%3;\n\t"
////						"}"
////						:"=l"(m0),"=l"(m1)
////						:"l"(x[j]),"l"(y[7-j])
////				);
//
//				asm("{\n\t"
//						"mul.lo.u64 %0,%2,%3;\n\t"
//						"mul.hi.u64 %1,%2,%3;\n\t"
//						"}"
//						:"=l"(m0),"=l"(m1)
//						:"l"(xs[xOffset]),"l"(y[j])
//				);
//				xOffset -= permutationStride;
//				c = u[2];
//				u[2] = 0;
//				device_p3_add128 (u[0], u[1], m0, m1, c, u[0], u[1], u[2]);
//				// 				printf("m0=%lu, m1=%lu, j=%d \n", m0, m1, j);
//			}

//			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i], hVectorSub[7 - i],
//						 cVectorSub[7 - i]);

//			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i], cVector[7 - i]);
//			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lSub, hSub, cSub);
//			device_p3_toLHC (u[0], u[1], u[2], lAdd, hAdd, cAdd);
//			offset = tid + (7 - i) * permutationStride;
			sign = 0;
			device_p3_sub192 (u[0], u[1], u[2], uSub[0], uSub[1], uSub[2], sign);
			lVector[offset] = u[0];
			hVector[offset] = u[1];
			cVector[offset] = u[2];
			signVector[offset] = sign;
//			hVectorSub[offset] = uSub[1];
//			cVectorSub[offset] = uSub[2];
//			offset=offset-permutationStride;
		}

//		printf("h \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       hVector[0],hVector[1],hVector[2],
//				       hVector[3],hVector[4],hVector[5],
//				       hVector[6],hVector[7]);
//
//		printf("c \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n %lu \n",
//				       cVector[0],cVector[1],cVector[2],
//				       cVector[3],cVector[4],cVector[5],
//				       cVector[6],cVector[7]);

	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step2 (const usfixn64* __restrict__ parameters,
												 usfixn64* __restrict__ lVector,
												 usfixn64* __restrict__ hVector,
												 usfixn64* __restrict__ cVector,
												 const usfixn32* __restrict__ signVector)
{

	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];

//	usfixn64 s0, s1, s2;
//	usfixn64 sign;
	usfixn64 l, h, c;
//	usfixn64 v0[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v1[COEFFICIENT_SIZE] =
//		{ 0, 0, 0, 0, 0, 0, 0, 0 };
//	usfixn64 v0[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};
//		usfixn64 v1[COEFFICIENT_SIZE] =
//			{ 0, 0, 0};
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
		offset = tid;
		short i = 0;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
//		memset (v0, 0x00, 3 * sizeof(usfixn64));
//		memset (v1, 0x00, 3* sizeof(usfixn64));
//		s0 = lVector[offset];
//		s1 = hVector[offset];
//		s2 = cVector[offset];
//		sign = lVectorSub[offset];
			l = 0;
			h = 0;
			c = 0;
//		device_p3_toLHC (s0, s1, s2, l, h, c);
			device_p3_toLHC (lVector[offset], hVector[offset], cVector[offset], l, h, c);

//		v0[0] = l;
//		v0[1] = h;
//		v0[2] = c;
//		if (sign == 1)
			if (signVector[offset] == 1)
			{
				device_p3_bigSubZero_3 (l, h, c);

//			device_p3_bigSub_plain (v1, v0);
//			l = v1[0];
//			h = v1[1];
//			c = v1[2];
			}
			lVector[offset] = l;
			hVector[offset] = h;
			cVector[offset] = c;
			offset += permutationStride;
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_mult_revised_8lhc_step3 (usfixn64* __restrict__ xs,
												 const usfixn64* __restrict__ lVector,
												 const usfixn64 * __restrict__ hVector,
												 const usfixn64* __restrict__ cVector,
												 const usfixn32* __restrict__ signVector,
												 const usfixn64 * __restrict__ parameters)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
//	short op = short(parameters[15]);
	usfixn64 offset = tid;
//	short c = 0;
	i = 0;
	usfixn64 n = parameters[7];

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
	usfixn64 complement[COEFFICIENT_SIZE] =
		{ 0, 0, 0, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE, R_MINUS_ONE };

	short j = 0;
	for (tid; tid < n; tid += gridDim.x * blockDim.x)
	{
//	uConstArray8_align8 v2;

//	memset(v0.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
		for (i = 0; i < 8; i++)
		{
//		v0.i[i]=0;
			v1.i[i] = 0;
		}

//	memset(v2.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

		offset = tid;
//all l values [0:7]
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			v0.i[i] = lVector[offset];
			offset += permutationStride;
		}

		complement[0] = 0;
		complement[1] = 0;
		complement[2] = 0;
		for (j = 3; j < COEFFICIENT_SIZE; j++)
			complement[j] = R_MINUS_ONE;
		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			if (i > 0)
//			device_3_cyclicShift_plain (complement, 1);
				device_p3_cyclicShift_permutated_12 (complement);
			if (signVector[offset] == 1)
			{
				device_p3_bigPrimeAdd_plain_inPlace (v0.i, complement);
//			device_p3_bigPrimeAdd_plain_ptx_v0(v0.i,complement);
			}
			offset += permutationStride;
		}

		offset = tid + (permutationStride << 3) - permutationStride;

//	v1.i[0] = hVector[tid + 7 * permutationStride];
		v1.i[0] = hVector[offset];
//	device_p3_bigSub_plain(v0.i, v1.i);
		device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}
//	v1.i[1] = cVector[tid + 7 * permutationStride];
//	v1.i[0] = cVector[tid + 6 * permutationStride];
		v1.i[1] = cVector[offset];
		offset -= permutationStride;
		v1.i[0] = cVector[offset];

//	device_p3_bigSub_plain(v0.i, v1.i);
		device_p3_bigPrimeSub_plain_inPlace (v0.i, v1.i);

		offset = tid;
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//		{
//			v1.i[i]=0;
//		}

		v1.i[0] = 0;
//positive h's [1:7]
		offset = tid;
		for (i = 1; i < COEFFICIENT_SIZE; i++)
		{
//		v1.i[i] = hVector[tid + (i - 1) * permutationStride];
			v1.i[i] = hVector[offset];
			offset += permutationStride;
		}
//		device_p3_bigPrimeAdd_plain (v0.i, v1.i);
		device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

		offset = tid;
//positive c's [2:7]
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for(i=0;i<8;i++)
//			{
//				v1.i[i]=0;
//			}
		v1.i[0] = 0;
		v1.i[1] = 0;
		offset = tid;
		for (i = 2; i < COEFFICIENT_SIZE; i++)
		{
//		v1.i[i] = cVector[tid + (i - 2) * permutationStride];
			v1.i[i] = cVector[offset];
			offset += permutationStride;
		}
//		device_p3_bigPrimeAdd_plain (v0.i, v1.i);
		device_p3_bigPrimeAdd_plain_inPlace (v0.i, v1.i);

//#################################################
//	uConstArray8_align8 v0_sub;
//	uConstArray8_align8 v1_sub;
//	uConstArray8_align8 v2_sub;

//	memset(v0_sub.i, 0x00, 8 * sizeof(usfixn64));
//	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	for (i = 0; i < 8; i++)
//	{
////				v0_sub.i[i]=0;
//		v1.i[i] = 0;
//	}
//	memset(v2_sub.i, 0x00, 8 * sizeof(usfixn64));
//#################################################

//	offset = tid;
////all l values [0:7]
//	for (i = 0; i < COEFFICIENT_SIZE; i++)
//	{
////		v0_sub.i[i] = lVectorSub[tid + (i) * permutationStride];
//		v0_sub.i[i] = lVectorSub[offset];
//		offset += permutationStride;
//	}
//	v0_sub.i[7] = 0;
//
////	v1.i[0] = hVectorSub[tid + 7 * permutationStride];
////	device_p3_bigSub_plain(v0_sub.i, v1.i);
////	device_p3_bigPrimeSub_plain_inPlace(v0_sub.i, v1.i);
//
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = cVectorSub[tid + 6 * permutationStride];
////	v1.i[1] = cVectorSub[tid + 7 * permutationStride];
//
////	device_p3_bigPrimeAdd_plain(v1.i, v2_sub.i);
////	device_p3_bigSub_plain(v0_sub.i, v1.i);
//	device_p3_bigPrimeSub_plain_inPlace (v0_sub.i, v1.i);
//
//	offset = tid;
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = 0;
////positive h's [1:7]
//	offset = tid;
//	for (i = 1; i < COEFFICIENT_SIZE; i++)
//	{
////		v1.i[i] = hVectorSub[tid + (i - 1) * permutationStride];
//		v1.i[i] = hVectorSub[offset];
////		v0_sub.i[i] = hVectorSub[tid + (8-i) * permutationStride];
////		v0_sub.i[i] = hVectorSub[tid + (6+i) * permutationStride];
//		offset += permutationStride;
//	}
//	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);
//
//	offset = tid;
////positive c's [2:7]
////	memset(v1.i, 0x00, 8 * sizeof(usfixn64));
//	v1.i[0] = 0;
//	v1.i[1] = 0;
//	offset = tid;
//	for (i = 2; i < COEFFICIENT_SIZE; i++)
//	{
////		v1.i[i] = cVectorSub[tid + (i - 2) * permutationStride];
//		v1.i[i] = cVectorSub[offset];
//		offset += permutationStride;
//	}
//	device_p3_bigPrimeAdd_plain (v0_sub.i, v1.i);
//
////#################################################
////	device_p3_bigSub_plain(v0.i, v0_sub.i);
//	device_p3_bigPrimeSub_plain_inPlace (v0.i, v0_sub.i);
////#################################################

//	############################# writing back to g-memory
		offset = tid;
		for (i = 0; i < COEFFICIENT_SIZE; i++)
		{
			xs[offset] = v0.i[i];
//		xs[offset] = v0_sub.i[i];
//		xs[offset] = lVectorSub[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
			offset += permutationStride;
		}
	}
}

/************************************************/
__global__ void
kernel_p3_lhc (usfixn64 * xs, usfixn64* ys, usfixn64* parameters)
{
//	usfixn64 tid = threadIdx.x+blockIdx.x*blockDim.x;
	usfixn64 n = parameters[7];
	usfixn64 tid = 0;
	usfixn64 x[COEFFICIENT_SIZE], y[COEFFICIENT_SIZE];
	usfixn64 u[3] =
		{ 0, 0, 0 }, t[3] =
		{ 0, 0, 0 };
	usfixn64 offset = 0;
	usfixn64 permutationStride = parameters[5];
	usfixn64 m0 = 0, m1 = 0;
	short i = 0, j = 0;
	usfixn64 lVector[COEFFICIENT_SIZE], hVector[COEFFICIENT_SIZE],
			cVector[COEFFICIENT_SIZE];
	usfixn64 lVectorSub[COEFFICIENT_SIZE], hVectorSub[COEFFICIENT_SIZE],
			cVectorSub[COEFFICIENT_SIZE];
	usfixn64 v0[COEFFICIENT_SIZE], v1[COEFFICIENT_SIZE];

//	u[0]=R-2;
//	u[1]=R-2;
//	u[2]=R-2;

//	u[2] = 2;
//	u[1] = 137438953583;
//	u[0] = 18446743523953737760;
//	 should give  (l,h,c)= 32, 9223372054034644960, 7

//	u[0]=18446743523953737760; u[1]=137438953583; u[2]=2 ;
	u[0] = 18446743592673214492, u[1] = 13835058175541248097, u[2] = 1;
//	u[0]=18446743661392691224, u[1]=9223372139933990995, u[2]=1 ;
//	u[0]=18446743730112167956, u[1]=4611686104326733893, u[2]=1 ;
//	u[0]=18446743798831644688, u[1]=68719476791, u[2]=1 ;
//	u[0]=18446743867551121420, u[1]=13835058106821771305, u[2]=0 ;
//	u[0]=18446743936270598152, u[1]=9223372071214514203, u[2]=0 ;
//	u[0]=18446744004990074884, u[1]=4611686035607257101, u[2]=0 ;

//	device_p3_sub192=
	u[0] = 1210760581564858368, u[1] = 2053551088978343176, u[2] = 0;
//	l=7252122177670479873, h=5116269883728032574, c=9223372054034644991
//	===========================
//	device_p3_sub192=  4566277235430883328, 4967893369498670099, 0 , sign=1
//	l=2028246543331524609, h=8510957387578794724, c=9223372054034644990
//	===========================
//	device_p3_sub192=  4483291542710910976, 8784722215583078071, 18446744073709551615 , sign=0
//	l=441273785808322560, h=8346072447485085807, c=9223372019674906621
//	===========================
//	device_p3_sub192=  6977625007207743488, 6233633731900117096, 0 , sign=1
//	l=5546970682761412609, h=5979476667491151013, c=9223372054034644990
//
	device_p3_toLHC (u[0], u[1], u[2], t[0], t[1], t[2]);

//	printf("u0=%lu, u1=%lu, u2=%lu \n ", u[0], u[1], u[2]);
//	printf ("t0=%lu, t1=%lu, t2=%lu \n ", t[0], t[1], t[2]);
//	device_p3_divR(u[0],u[1],t[1],t[0]);
//		offset = tid;
//		for (i = 0; i < 8; i++)
//	xs[0] = t[0];
//	xs[permutationStride] = t[1];
//	xs[2 * permutationStride] = t[2];
//			offset += permutationStride;
}

/**********************************************/
__global__ void
kernel_p3_bigMult_plain (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
{
	//
	short operation = parameters[0];
	short iterations = parameters[1];
	short paddingMethod = parameters[2];
	short dynamicMemSize = parameters[3];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	short shiftNo = 3; //parameters ?
	short i;
//	, pos, pos1, t;
	unsigned short c = 0;
//	usfixn64 num1, num2;
//	usfixn64 *xd, *yd, *ud, *xm, *ym, *um;
	usfixn64 *xd, *yd, *xm, *ym;
//	usfixn64 *xd, *xm;

//	num1 = 0;
//	num2 = R - 1;

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * 8);
	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * 8);
//	ud = (usfixn64*) ((char*) us + tid * sizeof(usfixn64) * 8);

//	device_3_cyclicShift_plain(xd,shiftNo);
//	bigMult(xs, ys, us);
	device_p3_bigMult_plain (xd, yd);
//	xd[0]=tid;
}

/**********************************************/
__global__ void
kernel_p3_bigMult_plain_2 (usfixn64 * xs, usfixn64* ys, usfixn64 *parameters)
//__global__ void plain(const usfixn64 * __restrict__ xs, const usfixn64 * __restrict__  ys, usfixn64 *us, const short * __restrict__ parameters)
{
	//
	short operation = parameters[0];
	short iterations = parameters[1];
	short paddingMethod = parameters[2];
	short dynamicMemSize = parameters[3];

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	short shiftNo = 3; //parameters ?
	short i;
//	, pos, pos1, t;
	unsigned short c = 0;
//	usfixn64 num1, num2;
//	usfixn64 *xd, *yd, *ud, *xm, *ym, *um;
	usfixn64 *xd, *yd, *xm, *ym;
//	usfixn64 *xd, *xm;

//	num1 = 0;
//	num2 = R - 1;

	xd = (usfixn64*) ((char*) xs + tid * sizeof(usfixn64) * 8);
	yd = (usfixn64*) ((char*) ys + tid * sizeof(usfixn64) * 8);
//	ud = (usfixn64*) ((char*) us + tid * sizeof(usfixn64) * 8);

//	device_3_cyclicShift_plain(xd,shiftNo);
//	bigMult(xs, ys, us);
	device_p3_bigMult_plain_r1 (xd, yd);
//	xd[0]=tid;
}

/**********************************************/
__global__ void
kernel_p3_multiplication_permutated (usfixn64 *xs, usfixn64 *parameters)
{

//	implement one of newMult functions here

//	short operation = parameters[0];
//	usfixn64 permutationStride = parameters[5];
//	short shuffle = parameters[6];
//	short padding = 0;//parameters[?]
//	short nShift= 3;//parameters [?]
//	usfixn64 idx;
//	usfixn64 tid = (threadIdx.x + blockIdx.x * blockDim.x);
//
//	//idx = (tid / permutationBlockSize) * 8 * permutationBlockSize + (tid % permutationBlockSize);
//	//following indexing is slightly faster than above indexing
//
//	idx = tid;
//
//	if(padding==0)
//		device_cyclicShift_permutated_2(&xs[idx], nShift, permutationStride);

}


/**********************************************/
//xs = xs*ys
//most efficient big multiplication
//not verified
//uses both shared memory, lmem, register, global mem
//try to fix negate and bigAddPlain to have static indexing
__global__ void
kernel_p3_newMult14 (usfixn64* __restrict__ xs, const usfixn64* __restrict__ ys,
					 usfixn64 * __restrict__ parameters, usfixn64* lVector,
					 usfixn64 *hVector, usfixn64* cVector)
{
//	usfixn64 n = parameters[0];
	usfixn64 permutationStride = parameters[5];
	short op = short (parameters[15]);
	short step = 0;
//	step=7;
	//read x0, x1
	//do one column
	//pass y to upper thread
	//8 threads for each coefficient

	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	uConstArray8_align8 y;
	usfixn64 offset = tid;
	usfixn64 m0, m1, m2;
//	unsigned short c = 0;
	usfixn64 lhc[8];
	usfixn64 tmp[8];
	usfixn64 lhcSub[8];
	usfixn64 lArraySub;

//	unsigned short cSub = 0;
	short c = 0;
	offset = tid;
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		y.i[i] = ys[offset];
		offset += permutationStride;
	}
//	shIdx = threadIdx.x;
	offset = tid;
	short k = 0;

	if (op == 0)	//addition part
	{
		for (k = 0; k < COEFFICIENT_SIZE; k++)
		{
			step = 8 - k;
			offset = tid;
//		l = 0, h = 0, c = 0;
			lhc[0] = 0;
			lhc[1] = 0;
			lhc[2] = 0;
//			memset(lhc,0x00,8*sizeof(usfixn64));
//			c = 0;
			//shift right y for #step times
			if (k > 0)
			{
				device_p3_oneShiftRight (y.i);
			}

			offset = tid;
			c = 0;
			for (i = 0; i < step; i++)
			{

//				m0 = 0;
//				m1 = 0;
//				m2 = 0;
				//##########################
				device_p3_mulLong_2_plain (xs[offset], y.i[7 - i], m0, m1, m2);

				tmp[0] = m0;
				tmp[1] = m1;
				tmp[2] = m2;
				device_p3_smallAdd2_plain (lhc, tmp, c);
//				device_smallAdd3_plain(lhc, tmp);
				offset += permutationStride;
			}
//			bigPrimeAdd_check(lhc, c, 8);
//			memset(tmp,0x00,8*sizeof(usfixn64));
//			device_p3_bigPrimeAdd_plain(lhc, tmp);
			step = k;
			offset = tid;
//			offset = tid + (permutationStride * k);
			offset = tid + (permutationStride * (7 - k));
//			lVector[offset] = lhc[0];
//			hVector[offset] = lhc[1];
//			cVector[offset] = lhc[2];

			lVector[offset] = lhc[0];
			hVector[offset] = lhc[1];
			cVector[offset] = lhc[2];
		}
//		return;
	}
	if (op == 1)
	{
		for (k = 0; k < COEFFICIENT_SIZE; k++)
		{
			step = 8 - k;
			offset = tid;
			//		l = 0, h = 0, c = 0;
			lhcSub[0] = 0;
			lhcSub[1] = 0;
			lhcSub[2] = 0;
//			cSub = 0;
			//shift right y for #step times //should reset value of y.i
			//will  be a problem in repetition
			if (k > 0)
			{
				device_p3_oneShiftRight (y.i);
			}

//			offset = tid + (k) * permutationStride;
			offset = tid + (7 - k) * permutationStride;
			lhc[0] = lVector[offset];
			lhc[1] = hVector[offset];
			lhc[2] = cVector[offset];
//			c = cVector[tid + k * permutationStride];
			offset = tid;
			c = 0;
			for (i = step; i < 8; i++)
			{

				m0 = 0;
				m1 = 0;
				m2 = 0;
				//##########################
				device_p3_mulLong_2_plain (xs[offset], y.i[7 - i], m0, m1, m2);

				tmp[0] = m0;
				tmp[1] = m1;
				tmp[2] = m2;
				device_p3_smallAdd2_plain (lhcSub, tmp, c);
//				device_smallAdd3_plain(lhcSub, tmp);
				offset += permutationStride;
			}
//			bigPrimeAdd_check(lhcSub, c, 3);
//			c-=cSub;
//			device_smallSub2_plain(lhc, lhcSub, c);
			c = 0;
//			device_smallSub3_plain(lhc, lhcSub,c);
//			bigPrimeSub_check(lhc, c, 3);
//			memset(&lhc[3], 0x00, 5 * sizeof(usfixn64));
//			memset(&lhcSub[3], 0x00, 5 * sizeof(usfixn64));
			memset (lhc, 0x00, 8 * sizeof(usfixn64));
//			device_p3_bigSub_plain(lhc, lhcSub);
			for (short x = 3; x < 8; x++)
			{
				lhc[x] = 0;
			}
//			if(k==0)
//			{
//				lhcSub[0]=0;
//				lhcSub[1]=0;
//				lhcSub[2]=0;
//			}
//			device_p3_bigSub(lhc, lhcSub,lhc);
//			device_p3_bigPrimeAdd_correct(lhc, lhcSub,lhc);
//			device_smallSub3_plain(lhc, lhcSub, c);
//			device_p3_smallAdd2_plain(lhc, lhcSub, c);

//			if(c==1)
//			bigPrimeSub_check(lhc, c, 8);
//			{
////				memset(&lhc[8], 0x00, 8 * sizeof(usfixn64));
//				lhc[0]=c;
//				lhc[1]=c;
//				lhc[2]=c;
//			}

//			offset = tid + (k * permutationStride);
			offset = tid + ((7 - k) * permutationStride);
			lVector[offset] = lhc[0];
			hVector[offset] = lhc[1];
			cVector[offset] = lhc[2];
		}
	}
}

/**********************************************/
__global__ void
kernel_p3_newMult14_incStep (usfixn64 * __restrict__ parameters)
{
	parameters[15]++;
}

/**********************************************/
__global__ void
kernel_p3_newMult14_join2 (usfixn64* __restrict__ xs, usfixn64 * __restrict__ parameters,
								 usfixn64* __restrict__ lVector,
								 usfixn64 * __restrict__ hVector,
								 usfixn64* __restrict__ cVector)
{
	usfixn64 permutationStride = parameters[5];
	short i = 0;
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	short op = short (parameters[15]);
	usfixn64 offset = tid;
	short c = 0;
	i = 0;

	uConstArray8_align8 v0;
	uConstArray8_align8 v1;
	uConstArray8_align8 v2;

	memset (v0.i, 0x00, 8 * sizeof(usfixn64));
	memset (v1.i, 0x00, 8 * sizeof(usfixn64));
	memset (v2.i, 0x00, 8 * sizeof(usfixn64));
	//#################################################

	offset = tid;
	//all l values [0:7]
	for (i = 0; i < COEFFICIENT_SIZE; i++)
	{
		v1.i[i] = lVector[tid + (i) * permutationStride];
//		v1.i[i] = lVector[tid + (7-i) * permutationStride];
		offset += permutationStride;
	}

//	c=1;
////	if(cVector[tid]==1)
//		bigPrimeSub_check(v1.i, c);
	memset (v0.i, 0x00, 8 * sizeof(usfixn64));
	v0.i[0] = hVector[tid + 7 * permutationStride];
//	v0.i[0] = hVector[tid + 0 * permutationStride];
//	device_p3_bigSub_plain(v1.i, v0.i);

	memset (v2.i, 0x00, 8 * sizeof(usfixn64));
//	v0.i[0] = cVector[tid + 6 * permutationStride];
//	v0.i[1] = cVector[tid + 7 * permutationStride];
	v2.i[0] = cVector[tid + 6 * permutationStride];
	v2.i[1] = cVector[tid + 7 * permutationStride];

//	v0.i[0] = cVector[tid + 1 * permutationStride];
//		v0.i[1] = cVector[tid + 0* permutationStride];
//	device_p3_bigSub_plain(v1.i, v0.i);
	device_p3_bigPrimeAdd_plain (v2.i, v0.i);

	offset = tid;
	memset (v0.i, 0x00, 8 * sizeof(usfixn64));
	//positive h's [1:7]
	for (i = 1; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = hVector[tid + (i - 1) * permutationStride];
//		v0.i[i] = hVector[tid + (8-i) * permutationStride];
//		v0.i[i] = hVector[tid + (6+i) * permutationStride];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v1.i, v0.i);

	offset = tid;
	//positive c's [2:7]
	memset (v0.i, 0x00, 8 * sizeof(usfixn64));
	for (i = 2; i < COEFFICIENT_SIZE; i++)
	{
		v0.i[i] = cVector[tid + (i - 2) * permutationStride];
//		v0.i[i] = cVector[tid + (9-i) * permutationStride];
		offset += permutationStride;
	}
	device_p3_bigPrimeAdd_plain (v1.i, v0.i);

//	device_p3_bigSub_plain(v1.i, v2.i);
	//#################################################
	offset = tid;
//	############################# writing back to g-memory
//#pragma unroll 8
	if (op == 0)
	{
		for (i = 0; i < 8; i++)
		{
			//		xs[offset] = l.i[j];

			lVector[offset] = v1.i[i];
			hVector[offset] = 0;
			cVector[offset] = 0;
//			xs[offset]=lVector[offset];
//				xs[offset] = v1.i[i];
			//		xs[offset] = cVector[offset];
			offset += permutationStride;
		}
	}
	offset = tid;

	if (op == 1)
		for (i = 0; i < 8; i++)
		{
//		xs[offset] = l.i[j];
//			xs[offset] = v1.i[7 - i];
			xs[offset] = v1.i[i];
			xs[offset] = lVector[offset];
//		xs[offset] = hVector[offset];
//		xs[offset] = cVector[offset];
			offset += permutationStride;
		}
}

/**********************************************/
__global__ void
kernel_p3_mulLong_revised (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;

	usfixn64 offset = 0;
	short i = 0;
	usfixn64 l[8], h[8], c[8];
	usfixn64 permutationStride = parameters[5];
	usfixn64 x[8], y[8];

	offset = tid;
	for (i = 0; i < 8; i++)
	{
		x[i] = xs[offset];
		y[i] = ys[offset];
		offset += permutationStride;
	}
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		device_p3_mulLong_2_plain_2 (x[i], y[i], l[i], h[i], c[i]);
//		device_p3_mulLong_2_plain_2(xs[offset], ys[offset], l[i], h[i], c[i]);
//		offset += permutationStride;
	}

	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = l[i];
		xs[offset] = h[i];
//			xs[offset]=c[i];
		offset += permutationStride;
	}
}

/**********************************************/
__global__ void
kernel_p3_mulLong_plain (usfixn64 * xs, usfixn64 *ys, usfixn64* parameters)
{
	usfixn64 tid = threadIdx.x + blockIdx.x * blockDim.x;
	usfixn64 offset = 0;
	short i = 0;
	usfixn64 lhc[8][3];
	usfixn64 permutationStride = parameters[5];

	usfixn64 xt, yt;
	offset = tid;
	for (i = 0; i < 8; i++)
	{
		device_p3_mulLong (xs[offset], ys[offset], lhc[i]);
//		xt=xs[offset];
//		yt=ys[offset];
//		xt=xt%R;
//		yt=yt%R;
//		lhc[i][0]=xt*yt;
		offset += permutationStride;
	}

	offset = tid;
	for (i = 0; i < 8; i++)
	{
		xs[offset] = lhc[i][0];
		xs[offset] = lhc[i][1];
//			xs[offset]=lhc[i][2];
		offset += permutationStride;
	}
}
/**********************************************/
#endif
