#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3


Node *build_tree( const Item * );
void cleanup( Node * );

void insert( const Node *, const Item * );
Node *search( key_t , Node *, lvl_t );

void find_neighbours( key_t , Node *, DArray_Item * );

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
