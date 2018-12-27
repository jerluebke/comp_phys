#include <math.h>
#include <stdio.h>
#include "../tests/test_data.h"


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


int main()
{
    size_t size = __any_values_1_size;
    Value *vals = __any_values_1;
    Item *items = xmalloc( sizeof(Item) * size );
    items = build_morton( vals, items, size );
    Node *head = build_tree( items );

    FILE *fp = fopen("data/res.py", "w+");
    fputs("tree = ", fp);

    print_node(head, fp);
    fputs("\n", fp);

    fclose(fp);

    cleanup( head );
    free( items );
}
