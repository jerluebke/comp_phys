#include <stdio.h>
#include "itoa.h"
#include "../include/types.h"
#include "../include/quadtree.h"

void print_node(Node *head, FILE *fp)
{
    size_t i;
    char ws[maxlvl], buf[2*maxlvl];

    for ( i = 0; i < head->lvl; ++i) ws[i] = ' ';
    ws[++i] = '\0';

    fprintf(fp, "%s[[%u, %u], %u, %u, {\n", ws, head->val_arr->p[0]->x,
            head->val_arr->p[0]->y, head->lvl, head->allocated);
    for ( i = 0; i < 4; ++i ) {
        fprintf(fp, "%s0b%s :\n", ws, __itoa(head->c[i]->key, buf, 2));
        print_node(head->c[i], fp);
    }
    fprintf(fp, "%s}]\n", ws);
}

void print_tree(Node *head, char *filename)
{
    char buf[2];
    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "{0b%s :\n", __itoa(head->key, buf, 2));
    print_node(head, fp);
    fprintf(fp, "}\n");
    fclose(fp);
}

void print_leafs(Node *head, Item *items, char *filename)
{
    char buf[2*maxlvl];
    Node *node;
    FILE *fp = fopen(filename, "w+");

    fputc('{', fp);
    while ( items != NULL ) {
        node = search(items->key, head, maxlvl);
        fprintf(fp, "0b%s : [[%u, %u], %u]\n", __itoa(node->key, buf, 2),
                node->val_arr->p[0]->x, node->val_arr->p[0]->y, node->lvl);
        ++items;
    }
    fputc('}', fp);

    fclose(fp);
}

void print_neighbours(Node *head, key_t key, FILE *fp)
{
    size_t i;
    char buf[2*maxlvl];
    DArray_Node ns = {NULL, 0, 0};

    DArray_Node_init(&ns, 8);
    find_neighbours(key, head, &ns);

    fprintf(fp, "ref: [0b%s\t%u]\n", __itoa(ns.p[0]->key, buf, 2), head->lvl);
    for ( i = 1; i < ns._used; ++i )
        fprintf(fp, "\t[0b%s,\t%u]\n", __itoa(ns.p[i]->key, buf, 2), ns.p[i]->lvl);
    fputs("\n\n", fp);
}

void print_all_neighbours(Node *head, Item *items, char *filename)
{
    FILE *fp = fopen(filename, "w+");

    while ( items != NULL ) {
        print_neighbours(head, items->key, fp);
        ++items;
    }

    fclose(fp);
}

void print_morton(Item *items, char *filename)
{
    char buf[2*maxlvl];
    FILE *fp = fopen(filename, "w+");

    while ( items != NULL ) {
        fprintf(fp, "[0b%s,\t[%u, %u]]\n", __itoa(items->key, buf, 2),
                items->val->x, items->val->y);
        ++items;
    }
}

/* vim: set ff=dos tw=79 sw=4 ts=4 et ic ai : */
