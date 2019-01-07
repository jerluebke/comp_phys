#include "types.h"
#include "morton.h"
#include "quadtree.h"

#define SQUARE(x) (x)*(x)
#define SD(v, w, a) SQUARE((v)->a - (w)->a)         /* squared difference */
#define METRIC(v, w) (SD(v, w, x) + SD(v, w, y))

#define SEARCH_FUNC(TYPE, CHECK_QUERY, VALUE_ACCESS)                        \
int search_##TYPE( Value *query, DArray_##TYPE *vals, double r_sq,          \
                   DArray_##TYPE *res )                                     \
{                                                                           \
    int num;                                                                \
    TYPE##Iterator *it, *end;                                               \
    it = DArray_##TYPE##_start(vals);                                       \
    end = DArray_##TYPE##_end(vals);                                        \
    for ( num = 0; it != end; ++it ) {                                      \
        CHECK_QUERY                                                         \
        if ( METRIC(query, (*it)VALUE_ACCESS) < r_sq ) {                    \
            DArray_##TYPE##_append(res, *it);                               \
            ++num;                                                          \
        }                                                                   \
    }                                                                       \
    return num;                                                             \
}

#define CHECK_QUERY if (*it==query) continue;
#define VALUE_ACCESS ->val
SEARCH_FUNC(Value, CHECK_QUERY,)
SEARCH_FUNC(Item,, VALUE_ACCESS)


void search_naive( const unsigned int *in, size_t size, double r_sq )
{
    size_t i, j;
    Value *vals;
    DArray_Value val_arr, res;

    vals = xmalloc(sizeof(Value) * size);
    for ( i = 0, j = 0; i < 2*size; i+=2, j++ ) {
        vals[j].x = in[i];
        vals[j].y = in[i+1];
    }

    val_arr = (DArray_Value) { (const Value **)&vals, size, size };
    DArray_Value_init(&res, 8);

    for ( i = 0; i < size; ++i ) {
        res._used = 0;
        search_Value( &vals[i], &val_arr, r_sq, &res );
    }

    DArray_Value_free(&res);
    free(vals);
}


void search_tree( const unsigned int *in, size_t size, double r_sq )
{
    size_t i, j;
    Value *vals;
    Item *items;
    Node *head;
    DArray_Item tmp, res;

    vals = xmalloc(sizeof(Value) * size);
    for ( i = 0, j = 0; i < 2*size; i+=2, j++ ) {
        vals[j].x = in[i];
        vals[j].y = in[i+1];
    }
    items = xmalloc(sizeof(Item)*size);
    items = build_morton(vals, items, size);
    head = build_tree(items, insert_simple);

    DArray_Item_init(&tmp, 8);
    DArray_Item_init(&res, 8);

    for ( i = 0; i < size; ++i ) {
        res._used = 0;
        find_neighbours( items[i].key, head, &tmp );
        search_Item( (Value *)items[i].val, &tmp, r_sq, &res );
    }

    DArray_Item_free(&tmp);
    DArray_Item_free(&res);
    free(items);
    free(vals);
}


/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
