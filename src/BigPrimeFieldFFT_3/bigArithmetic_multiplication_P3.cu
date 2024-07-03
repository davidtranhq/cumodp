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


#ifndef BIG_ARITHMETIC_MULTIPLICATION_CU_
#define BIG_ARITHMETIC_MULTIPLICATION_CU_

//#include "BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "../../include/BigPrimeFieldFFT_3/bigPrimeField_P3.h"
#include "bigArithmetic_multiplication_P3_device.h"

/**********************************************/
__global__ void
kernel_p3_mult_revised_singleThread (usfixn64 * xs, usfixn64* ys,
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
				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1],
													uSub[2]);
			}
			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i],
											 hVectorSub[7 - i], cVectorSub[7 - i]);
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
			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i],
											 cVector[7 - i]);
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
				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1],
													uSub[2]);
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

				device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1],
													uSub[2]);
			}
			device_p3_toLHC (uSub[0], uSub[1], uSub[2], lVectorSub[7 - i],
											 hVectorSub[7 - i], cVectorSub[7 - i]);
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
			device_p3_toLHC (u[0], u[1], u[2], lVector[7 - i], hVector[7 - i],
											 cVector[7 - i]);

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
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub)
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
					device_p3_add128 (uSub[0], uSub[1], m0, m1, uSub[2], uSub[0], uSub[1],
														uSub[2]);
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
															usfixn64* lVector, usfixn64* hVector,
															usfixn64* cVector, usfixn64* lVectorSub,
															usfixn64* hVectorSub, usfixn64* cVectorSub)
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
kernel_p3_mult_revised_step3 (usfixn64* __restrict__ xs,
															usfixn64* __restrict__ lVector,
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
																	 const usfixn64* __restrict__ ys,
																	 usfixn64* parameters,
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
			device_p3_toLHC (lVector[offset], hVector[offset], cVector[offset], l, h,
											 c);

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
kernel_p3_newMult14_join2 (usfixn64* __restrict__ xs,
													 usfixn64 * __restrict__ parameters,
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
