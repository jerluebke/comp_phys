#include "test_data.h"
#include "csearch.h"


int main()
{

    const double r_sq = 8.0f;
    const size_t size = 16;
    static const unsigned int in[] = {
        0x1C, 0x68, 0x7F, 0x7F , 0x20, 0xA0 , 0x68, 0x86 ,
        0x0, 0x0 , 0x1C, 0x5E , 0x2A, 0x54 , 0x24, 0x5A ,
        0x54, 0x8A , 0x55, 0x8A , 0x55, 0x8B , 0xFF, 0xFF ,
        0x60, 0x20 , 0x80, 0x7F , 0x3F, 0xC0 , 0x40, 0xFF
    };

#if 0
    fputs("call with `valgrind bin/mem_test.out`\n\n", stderr);

    size_t size = __any_values_1_size;
    Value *vals = __any_values_1;

    fputs("building morton keys...\n\n", stderr);

    Item *items = xmalloc( sizeof(Item) * size );
    items = build_morton( vals, items, size );

    fputs("done.\nbuilding quadtree...\n\n", stderr);

    Node *head = build_tree( items, insert_fast );

    fputs("done.\ncleaning up...\n\n", stderr);

    cleanup( head );
    free( items );

    /**/

    double r_sq = 8.0f;
    unsigned int *in = xmalloc(sizeof(lvl_t) * 2*size);
    size_t i, j;
    for ( i = 0, j = 0; i < size; i++, j+=2 ) {
        in[j] = vals[i].x;
        in[j+1] = vals[i].y;
    }
#endif

    fputs("running neighbour search for every value...\n\n", stderr);
    fputs("naive search...\n", stderr);
    search_naive(in, size, r_sq);
    fputs("quadtree search with naive insert...\n", stderr);
    search_fast(in, size, r_sq);
    fputs("quadtree search with fast insert...\n", stderr);
    search_fastfast(in, size, r_sq);
    fputs("\ndone searching.\n", stderr);

    fputs("done.\n", stderr);

    return EXIT_SUCCESS;
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
