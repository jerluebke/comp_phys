#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3

/* TODO: print function for nodes */

void insert( const Node *, const key_t *, const Value * );
void delete( Node * );
Node *search( Node *, key_t );
DArray_Node *search_children( const Node *, const Node *, DArray_Node *);
