#pragma once

#include "types.h"
#include "morton.h"
#include "quadtree.h"

#define SQUARE(x) (x)*(x)

typedef struct fargs_s {
    const unsigned int *data;
    size_t size;
    double r_sq;
} fargs_t;

int search_naive( const fargs_t * );
int search_fast( const fargs_t * );
int search_fastfast( const fargs_t * );

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
