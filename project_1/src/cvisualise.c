#include <math.h>
#include "cvisualise.h"
#include "morton.h"
#include "quadtree.h"

static inline key_t bap( key_t k, lvl_t j, lvl_t l )
{
    return ( k >> DIM * (l - j) ) & MASK;
}

/* TODO: draw graph for each step */
lvl_t qtenv_insert(QuadtreeEnv *this, double *res)
{
    lvl_t i, nl, rl;    /* index, new levels, relevant level */
    Node *n, *m;        /* temporary nodes */
    Item *item;
    Value v;            /* temporary coordinate */
    double r;           /* resolution */
    key_t key;
    n       = this->head;
    item    = &(this->items[this->idx]);
    key     = item->key;
    ++this->idx;

    nl  = insert(n, item);
    m   = search(key, n, maxlvl);
    rl  = m->lvl - nl;

    while ( n->lvl < rl-1 )
        n = n->c[bap(m->key, n->lvl+1, m->lvl)];
    for ( i = 0; n->lvl < m->lvl; i+=3) {
        n = n->c[bap(m->key, n->lvl+1, m->lvl)];
        v = coords2(n->key);
        r = pow(2, n->lvl-2);
        res[i]      = ((double)v.x) / r;
        res[i+1]    = ((double)v.y) / r;
        res[i+2]    = 1.0f / r;
    }

    return nl;
}

int qtenv_is_last(QuadtreeEnv *this)
{
    return this->items[this->idx].last;
}

QuadtreeEnv *qtenv_setup(const unsigned int *in, size_t size, unsigned int *si)
{
    size_t i, j;
    Value *vals;
    Item *items;
    Node *head;
    QuadtreeEnv *this;

    vals = xmalloc(sizeof(Value) * size);
    for ( i = 0, j = 0; i < 2*size; i+=2, j++ ) {
        vals[j].x = in[i];
        vals[j].y = in[i+1];
    }
    items = xmalloc(sizeof(Item)*size);
    items = build_morton(vals, items, size);
    for ( i = 0; i < size; ++i ) si[i] = items[i].idx;
    head = xmalloc(sizeof(Node));
    head->allocated = 1;
    head->key = 0;
    head->lvl = 0;
    head->i = NULL;
    head->c = xmalloc(sizeof(Node *)*NOC);
    for ( i = 0; i < NOC; ++i ) head->c[i] = NULL;
    this = xmalloc(sizeof(QuadtreeEnv));
    this->idx   = 0;
    this->vals  = vals;
    this->items = items;
    this->head  = head;

    return this;
}

void qtenv_free(QuadtreeEnv *this)
{
    free(this->vals);
    free(this->items);
    cleanup(this->head);
    free(this);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
