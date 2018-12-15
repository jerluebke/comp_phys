#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3

void insert( const Node *, const key_t *, Value * );
void delete( Node * );
Node *search( Node *head, key_t key );
DArray_Node *search_children( Node *head, const Node *ref, DArray_Node *res );
