#define _CRT_SECURE_NO_WARNINGS 1
/* #include <math.h> */
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
    Agnode_t *n, *m;

    snprintf(buf, NAMESZ, "%u-%x", head->lvl, head->key);
    n = agnode(g, buf, 1);

    if ( prev )
        agedge(g, prev, n, 0, 1);

    if ( !head->c )
        agsafeset(n, "color", "blue", "");

    else {
        agsafeset(n, "color", "gray", "");
        for ( i = 0; i < NOC; ++i ) {
            if ( head->c[i] )
                build_graph(g, head->c[i], n, buf);
            /* else {
             *     snprintf(buf, NAMESZ, "%u-%x",
             *              head->lvl+1, (key_t)((head->key << 2) | i));
             *     m = agnode(g, buf, 1);
             *     agedge(g, n, m, NULL, 1);
             *     agsafeset(m, "color", "gray", "");
             * } */
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
    gvRender(gvc, g, "dot", stdout);
    gvFreeLayout(gvc, g);

    agclose(g);

    return gvFreeContext(gvc);
}


#if 0
void print_node(Node *head, FILE *fp)
{
    int i = 0;
    Coords2d_8bit c = coords2(head->key);

    char ws[maxlvl];
    while ( i < head->lvl ) ws[i++] = ' ';
    ws[i] = '\0';

    fprintf(fp, "%s[(%u, %u), %lf, '%c', [\n",
            ws, c.x, c.y, 1./pow(2, head->lvl), head->c ? 'l' : 'b');
    if ( !head->c )
        fprintf(fp, "%s 'none'", ws);
    else
        for ( i = 0; i < 4; ++i )
            if ( head->c[i] )
                print_node(head->c[i], fp);
    fprintf(fp, "\n%s]],", ws);
}
#endif


int main()
{
    size_t size = any_values_1_size;
    Value *vals = any_values_1;
    Item *items = xmalloc( sizeof(Item) * size );
    items = build_morton( vals, items, size );
    Node *head = build_tree( items );

#if 0
    FILE *fp = fopen("data/res.py", "w+");
    fputs("tree = ", fp);

    print_node(head, fp);
    fputs("\n", fp);

    fclose(fp);
#endif

    FILE *fp = fopen("data/out.svg", "w+");
    if (render_tree(head, fp))
        return EXIT_FAILURE;
    fclose(fp);

    cleanup( head );
    free( items );

    return EXIT_SUCCESS;
}
