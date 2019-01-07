#pragma once

#include "types.h"
#include "morton.h"
#include "quadtree.h"


void search_naive( const unsigned int *, size_t, double );
void search_fast( const unsigned int *, size_t, double );
void search_fastfast( const unsigned int *, size_t, double );

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
