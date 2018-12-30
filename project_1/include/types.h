#pragma once

#include <stdint.h>
#include "dynamic_array.h"
/* #define key_t __quadtree_key_t */

/* for 256x256 Resolution */
typedef uint16_t key_t;
typedef uint8_t lvl_t;
static const uint8_t maxlvl = 8;

typedef struct Value Value;
typedef struct Item Item;
typedef struct Node Node;

DARRAY_TYPEDEF(const Item *, Item)


/* structure implementations */

struct Value {
    lvl_t x, y;
};

struct Item {
    key_t key;
    const Value *val;
    size_t idx;     /* original index */
    uint8_t last;
};

struct Node {
    key_t key;
    const Item *i;  /* content */
    lvl_t lvl;
    Node **c;       /* children */
    uint8_t allocated;
};

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
