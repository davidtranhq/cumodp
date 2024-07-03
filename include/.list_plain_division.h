#ifndef _LIST_PLAIN_DIVISION_H_
#define _LIST_PLAIN_DIVISION_H_
/*******************************************************************
* This code is intended for for computing remainder i.e f rem m    *
* The degree of m is T*W-1 that is 511 and the degree of f is 1020
********************************************************************/

#include<iostream>
#include <ctime>
#include<cmath>

#include "inlines.h"
#include "printing.h"
#include "types.h"

/*************************************************************
*                 T*W <= 2^9                                 *
* For bigger polynomial, make T = 32 then W will be W = 16   *
*************************************************************/
const sfixn T = 128;
/***********************************************************
* one thread is responsible for computing "W" coefficients *
***********************************************************/
const sfixn W = 4; 


#define BASE_1 31
__global__ void list_divCUDA(sfixn *M, sfixn *F, sfixn start, sfixn length, sfixn threadsPerDiv, sfixn DivPerBlock, sfixn polyNum, sfixn P);

void list_div(sfixn *M, sfixn *F, sfixn n, sfixn m, sfixn  start_offset, sfixn length_poly, sfixn poly_on_layer, sfixn p );

#endif
