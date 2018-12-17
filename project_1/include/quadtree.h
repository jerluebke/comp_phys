#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3


Node *build_tree( const Item * );
void find_neighbours( key_t , Node *, DArray_Value * );

void insert( const Node *, const Item * );
void delete( Node * );
Node *search( Node *, key_t , lvl_t );
void print_node( Node * );

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
