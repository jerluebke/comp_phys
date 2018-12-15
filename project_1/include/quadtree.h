#pragma once

#include <stdint.h>
#include "dynamic_array.h"

/* number of children */
#define NOC 4

/* for 256x256 Resolution */
typedef uint16_t key_t;
typedef uint8_t lvl_t;
static const uint8_t maxlvl = 8;

typedef struct Value Value;
typedef struct Node Node;

DARRAY(Value *, Value)
DARRAY(Node *, Node)

void insert(Node *, key_t *, Value *);
void delete(Node *);
// Res *search(Node *, key_t, Res *);
