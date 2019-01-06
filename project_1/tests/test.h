#pragma once

#include <stdio.h>
#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"
#include "../include/morton.h"
#include "../include/quadtree.h"


typedef struct {
    Item *i; key_t *k; DArray_Item *da; size_t s;
    insert_fptr_t ifunc;
} DataStruct;

typedef struct {
    Value *vals; size_t size; key_t *exp;
} TreeInput;

typedef struct {
    key_t key; key_t *val;
} KeyValueInput;


/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
