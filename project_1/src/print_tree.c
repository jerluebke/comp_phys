#define _CRT_SECURE_NO_WARNINGS 1
#include <math.h>
#include <stdio.h>
#include <graphviz/gvc.h>
#include "morton.h"
#include "quadtree.h"


#define NAMESZ  10u

#define any_values_1_size 16u
static Value any_values_1[] = {
    { 0x1C, 0x68 }, { 0x7F, 0x7F }, { 0x20, 0xA0 }, { 0x68, 0x86 },
    { 0x0, 0x0 }, { 0x1C, 0x5E }, { 0x2A, 0x54 }, { 0x24, 0x5A },
    { 0x54, 0x8A }, { 0x55, 0x8A }, { 0x55, 0x8B }, { 0xFF, 0xFF },
    { 0x60, 0x20 }, { 0x80, 0x7F }, { 0x3F, 0xC0 }, { 0x40, 0xFF }
};


void build_graph(Agraph_t *g, Node *head, Agnode_t *prev, char *buf)
{
    size_t i;
    Agnode_t *n/*, *m*/;

    snprintf(buf, NAMESZ, "%u-%x", head->lvl, head->key);
    n = agnode(g, buf, 1);

    if ( prev )
        agedge(g, prev, n, 0, 1);

    if ( !head->c )
        agsafeset(n, "color", "blue", "");

    else {
        agsafeset(n, "color", "grey", "");
        for ( i = 0; i < NOC; ++i ) {
            if ( head->c[i] )
                build_graph(g, head->c[i], n, buf);
        }
    }
}


int render_tree(Node *head, FILE *fp)
{
    char buf[NAMESZ];
    Agraph_t *g;
    GVC_t *gvc;

    gvc = gvContext();
    g = agopen("g", Agundirected, NULL);

    build_graph(g, head, NULL, &buf[0]);

    gvLayout(gvc, g, "dot");
    gvRender(gvc, g, "svg", fp);
    /* gvRender(gvc, g, "dot", stdout); */
    gvFreeLayout(gvc, g);

    agclose(g);

    return gvFreeContext(gvc);
}


int main()
{
    char buf[120], filename[128];
    size_t size = any_values_1_size;
    Value *vals = any_values_1;
    Item *items = xmalloc( sizeof(Item) * size );
    items = build_morton( vals, items, size );
    Node *head = build_tree( items );

    printf("Enter name for output file: ");
    scanf("%120s", buf);
    snprintf(filename, 128, "data/%s.svg", buf);

    FILE *fp = fopen(filename, "w+");
    if (render_tree(head, fp))
        return EXIT_FAILURE;

    fclose(fp);
    cleanup( head );
    free( items );

    return EXIT_SUCCESS;
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */