#include "../tests/test_data.h"

int main()
{
    fputs("call with `valgrind bin/mem_test.out`\n\n", stderr);

    size_t size = __any_values_1_size;
    Value *vals = __any_values_1;

    fputs("building morton keys...\n\n", stderr);

    Item *items = xmalloc( sizeof(Item) * size );
    items = build_morton( vals, items, size );

    fputs("done.\nbuilding quadtree...\n\n", stderr);

    Node *head = build_tree( items );

    fputs("done.\ncleaning up...\n\n", stderr);

    cleanup( head );
    free( items );

    fputs("done.\n", stderr);

    return EXIT_SUCCESS;
}
