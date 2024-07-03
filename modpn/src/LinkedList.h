/* This file is part of the MODPN library

    MODPN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MODPN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Copyright, Marc Moreno Maza <moreno@csd.uwo.ca>
*/


/* Authors: Xin Li <xli96@csd.uwo.ca>, Marc Moreno Maza <moreno@csd.uwo.ca> */
/* Copyright (c) 2009 by Marc Moreno Maza.  All rights reserved             */
#ifndef __LinkedList_h
#define __LinkedList_h

#include "Types.h"
#include "MPMMTS.h"

RegularPair *EX_RegularPair_Init(preFFTRep * poly, TriSet *ts);

void EX_RegularPair_Print(void *element);

void EX_RegularPair_Free(void *element);

RegularListPair *EX_RegularListPair_Init(sfixn no, preFFTRep **polyList, TriSet *ts);

void EX_RegularListPair_Print(void *element);

void EX_RegularListPair_Free(void *element);

void EX_Poly_Print(void *element);

void EX_Poly_Free(void *element);

TaskPair2 *EX_TaskPair2_Init(int i, int d, TriSet *ts);

void EX_TaskPair2_Print(void *element);

void EX_TaskPair2_Free(void *element);

TaskPair *EX_TaskPair_Init(int index, TriSet *ts);

void EX_TaskPair_Print(void *element);

void EX_TaskPair_Free(void *element);

LinearNode *EX_LinearNode_Init(void *element);

void EX_LinearNode_Print(LinearNode *node, void (*printELement)(void *));


void EX_LinearNode_Free(LinearNode *node, void (*freeElement)(void *));

LinkedQueue *EX_LinkedQueue_Init();

int EX_LinkedQueue_IsEmpty(LinkedQueue *queue);

void EX_LinkedQueue_Print(LinkedQueue *queue, void (*printELement)(void *));

void FILE_LinkedQueue_Print(LinkedQueue *queue, FILE *file);

void EX_LinkedQueue_Free(LinkedQueue *queue, void (*freeElement)(void *));

void EX_LinkedQueue_Enqeue(LinkedQueue *queue, void *element);

void *EX_LinkedQueue_Deqeue(LinkedQueue *queue);

LinkedQueue *EX_LinkedQueue_Copy(LinkedQueue *queue, void *(*copyElement)(void *) );

int EX_LinkedQueue_Size(LinkedQueue *queue);

LinkedQueue * EX_LinkedQueue_Concat_1(LinkedQueue *queue1, LinkedQueue *queue2);

void **LinkedQueue2Array(LinkedQueue *queue, void *(*copyElement)(void *));

void *EX_CopyRegularPair(void *pair);

void *EX_CopyRegularListPair(void *pair);

void *EX_CopyPoly(void *element);

void EX_RegularListPair_List_Free(void *element);

void EX_RegularPair_List_Free(void *element);

#endif
