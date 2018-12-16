#pragma once

#include "types.h"

/* dimensions */
#define DIM 2
/* number of children */
#define NOC 4
/* significant bits */
#define MASK 0x3

/* TODO: print function for nodes */

Node *build_tree( const Item * );
DArray_Value *find_neighbours( key_t , DArray_Value * );

void insert( const Node *, const Item * );
void delete( Node * );
Node *search( Node *, key_t );
void search_children( const Node *, const Node *, DArray_Node * );
void print_node( Node * );
