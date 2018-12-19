#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3


#define FIND_NEIGHBOURS_RETURN_NODE 1


Node *build_tree( const Item * );
void cleanup( Node * );

void insert( const Node *, const Item * );
Node *search( key_t , Node *, lvl_t );

#if FIND_NEIGHBOURS_RETURN_NODE
void find_neighbours( key_t , Node *, DArray_Node * );
#else
void find_neighbours( key_t , Node *, DArray_Value * );
#endif

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
