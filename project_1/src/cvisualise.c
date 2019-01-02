#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "itoa.h"
#include "cvisualise.h"
#include "morton.h"
#include "quadtree.h"

#define KEYSIZE     16u
#define NAMESIZE    10u
#define LABELSIZE   10u

#define FILENAMETEMPLATE    "data/graphs/out-%.3lu.png"
#define FILENAMESIZE        30u

#define LABELTEMPLATE       "current node: 0x%.4X = 0b%s"
#define GRAPHLABELSIZE      64u


static inline key_t bap( key_t k, lvl_t j, lvl_t l )
{
    return ( k >> DIM * (l - j) ) & MASK;
}

#ifndef _WIN32
static void build_graph(Agraph_t *g, Node *head, Agnode_t *prev, char *buf)
{
    size_t i;
    Agnode_t *n;

    snprintf(buf, NAMESIZE, "%u-%x", head->lvl, head->key);
    n = agnode(g, buf, 1);
    agsafeset(n, "label", "", "");
    agsafeset(n, "style", "filled", "");

    if ( prev )
        agedge(g, prev, n, 0, 1);

    if ( !head->c )
        agsafeset(n, "fillcolor", "green", "");

    else {
        agsafeset(n, "fillcolor", "grey", "");
        for ( i = 0; i < NOC; ++i )
            if ( head->c[i] )
                build_graph(g, head->c[i], n, buf);
    }
}
#endif

unsigned int qtenv_get_key(QuadtreeEnv *this, unsigned int idx)
{
    return this->items[idx].key;
}

unsigned int qtenv_insert(QuadtreeEnv *this, double *res)
{
    lvl_t i, nl, rl, sb;/* index, new levels, relevant level, significant bytes */
    Node *n, *m;        /* temporary nodes */
    Item *item;
    Value v;            /* temporary coordinate */
    double r;           /* resolution */
    char buf1[GRAPHLABELSIZE], buf2[LABELSIZE];
#ifndef _WIN32
    Agraph_t *g;
    Agnode_t *gn;
#endif

    key_t key;
    n       = this->head;
    item    = &(this->items[this->idx]);
    key     = item->key;
    ++this->idx;

    nl  = insert(n, item);
    m   = search(key, n, maxlvl);
    rl  = m->lvl - nl;
#ifndef _WIN32
    g   = agopen("g", Agundirected, NULL);
    build_graph(g, n, NULL, buf1);
#endif

    while ( n->lvl < rl-1 ) {
        sb = bap(m->key, n->lvl+1, m->lvl);
        n = n->c[sb];
#ifndef _WIN32
        snprintf(buf1, NAMESIZE, "%u-%x", n->lvl, n->key);
        gn = agnode(g, buf1, 0);
        snprintf(buf1, LABELSIZE, "0b%.*s%s", sb > 1 ? 0 : 1, "0",
                                              my_itoa(sb, buf2, 2));
        agsafeset(gn, "label", buf1, "");
#endif
    }
    for ( i = 0; n->lvl < m->lvl; i+=3 ) {
        sb = bap(m->key, n->lvl+1, m->lvl);
        n = n->c[sb];
#ifndef _WIN32
        snprintf(buf1, NAMESIZE, "%u-%x", n->lvl, n->key);
        gn = agnode(g, buf1, 0);
        snprintf(buf1, LABELSIZE, "0b%.*s%s", sb > 1 ? 0 : 1, "0",
                                              my_itoa(sb, buf2, 2));
        agsafeset(gn, "label", buf1, "");
        agsafeset(gn, "labelfontsizel", "64", "");
        if ( i )
            agsafeset(gn, "fillcolor", "red", "");
#endif
        v = coords2(n->key);
        r = pow(2, n->lvl-2);
        res[i]      = ((double)v.x) / r;
        res[i+1]    = ((double)v.y) / r;
        res[i+2]    = 1.0f / r;
    }
#ifndef _WIN32
    agsafeset(gn, "style", "bold,filled", "");
    agsafeset(gn, "fillcolor", "red", "");
    /* snprintf(buf1, GRAPHLABELSIZE, LABELTEMPLATE, m->key, my_itoa(m->key, buf2, 2)); */
    /* agsafeset(g, "label", buf1, ""); */
    /* agsafeset(g, "labelloc", "b", ""); */

    FILE *tf = fopen("data/graphs/key.file", "a+");
    fprintf(tf, "0b%s\n", my_itoa(m->key, buf2, 2));
    fclose(tf);

    char fname[FILENAMESIZE];
    snprintf(fname, FILENAMESIZE, FILENAMETEMPLATE, this->idx);
    FILE *fp = fopen(fname, "w+");
    gvLayout(this->gvc, g, "dot");
    gvRender(this->gvc, g, "png", fp);
    gvFreeLayout(this->gvc, g);
    fclose(fp);

    agclose(g);
#endif

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
#ifndef _WIN32
    this->gvc   = gvContext();
#endif

    return this;
}

void qtenv_free(QuadtreeEnv *this)
{
    free(this->vals);
    free(this->items);
    cleanup(this->head);
#ifndef _WIN32
    gvFreeContext(this->gvc);
#endif
    free(this);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
