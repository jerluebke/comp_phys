#pragma once
#include <graphviz/gvc.h>
#include "types.h"

typedef struct QuadtreeEnv {
    size_t idx;
    Value *vals;
    Item *items;
    Node *head;
    GVC_t *gvc;
} QuadtreeEnv;

unsigned int qtenv_get_key(QuadtreeEnv *, unsigned int);
unsigned int qtenv_insert(QuadtreeEnv *, double *);
int qtenv_is_last(QuadtreeEnv *);
QuadtreeEnv *qtenv_setup(const unsigned int *, size_t, unsigned int *);
void qtenv_free(QuadtreeEnv *);

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */