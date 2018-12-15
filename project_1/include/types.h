#pragma once

#include <stdint.h>
#include "dynamic_array.h"

/* for 256x256 Resolution */
typedef uint16_t key_t;
typedef uint8_t lvl_t;
static const uint8_t maxlvl = 8;

typedef struct Value Value;
typedef struct Node Node;

DARRAY(Value *, Value)
DARRAY(Node *, Node)

