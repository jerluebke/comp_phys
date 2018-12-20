#include <stdio.h>
#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"
#include "itoa.h"
#include "../include/morton.h"
#include "../include/quadtree.h"


typedef struct {
    Item *i; key_t *k; DArray_Node *da; size_t s;
} DataStruct;

typedef struct {
    Value *vals; size_t size; key_t *exp;
} TreeInput;

typedef struct {
    key_t key; key_t *val;
} KeyValueInput;


void print_node(Node *, FILE *);
void print_tree(Node *, char *);
void print_leafs(Node *, Item *, char *);
void print_neighbours(Node *, key_t, FILE *);
void print_all_neighbours(Node *, Item *, char *);
void print_morton(Item *, char *);

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
