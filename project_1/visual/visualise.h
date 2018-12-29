#include "types.h"

typedef struct QuadtreeEnv {
    size_t idx;
    Value *vals;
    Item *items;
    Node *head;
} QuadtreeEnv;

lvl_t qtenv_insert(QuadtreeEnv *, lvl_t *);
QuadtreeEnv *qtenv_setup(const lvl_t *, size_t);
void qtenv_free(QuadtreeEnv *);

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
