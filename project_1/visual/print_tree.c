#include <math.h>
#include <stdio.h>
#include "../tests/test_data.h"

extern inline lvl_t decode(key_t);
extern inline Value coords(key_t);

static const key_t B[] = {0x5555, 0x3333, 0x0F0F};
static const key_t S[] = {1, 2, 4};

inline lvl_t decode(key_t k)
{
    k &= B[0];
    k = (k ^ (k >> S[0])) & B[1];
    k = (k ^ (k >> S[1])) & B[2];
    k = (k ^ (k >> S[2]));
    return k;
}

inline Value coords(key_t k)
{
    Value val = { .x = decode(k), .y = decode(k >> 1) };
    return val;
}

void print_node(Node *head, FILE *fp)
{
    int i = 0;
    Value c = coords(head->key);

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
