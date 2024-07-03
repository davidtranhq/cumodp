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
#ifndef __HashTable_h
#define __HashTable_h 

#include "Types.h"
#include "UniHensel.h"

#define  HashModuli 101

typedef struct node_struct NODE; 
typedef struct hashtable_struct HASHTABLE;

struct node_struct {
  operand opnd;
  struct node_struct *next;	
};


// array index is the key value.
struct hashtable_struct
{
        NODE ** table;
	unsigned int length;
};


NODE *list_create(operand opnd);

NODE *list_insert_after(NODE *node, operand opnd);

NODE *list_insert_beginning(NODE *list, operand opnd);

int list_remove(NODE *list, NODE *node);

int compareOperand(operand opnd1, operand opnd2);

NODE *list_find(NODE *node, operand opnd);

void list_free(NODE *list);

int hashFunction(operand opnd);


// yes, return the existing node.
// no, return a NULL;
operand searchNodeinHashTable(operand opnd, HASHTABLE *hash_table);



HASHTABLE *newHashTable(unsigned int length);


void freeHashTable(HASHTABLE *hash_table);


operand  searchNodeinHashTable(operand opnd, HASHTABLE *hash_table);

#endif
