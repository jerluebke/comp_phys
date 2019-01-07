#include "types.h"
#include "morton.h"
#include "quadtree.h"

#define SQUARE(x) (x)*(x)
#define SD(v, w, a) SQUARE((v)->a - (w)->a)         /* squared difference */
#define METRIC(v, w) (SD(v, w, x) + SD(v, w, y))

#define SEARCH_FUNC(TYPE, CHECK_QUERY, VALUE_ACCESS)                        \
static int search_##TYPE( Value *query, DArray_##TYPE *vals, double r_sq,   \
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


#define SEARCH_SETUP(NAME, TYPE, DECL, INIT, PREP, PARAM, FREE)             \
void search_##NAME( const unsigned int *in, size_t size, double r_sq )      \
{                                                                           \
    size_t i, j;                                                            \
    Value *vals;                                                            \
    DECL                                                                    \
    DArray_##TYPE tmp, res;                                                 \
    vals = xmalloc(sizeof(Value) * size);                                   \
    for ( i = 0, j = 0; i < 2*size; i+=2, j++ ) {                           \
        vals[j].x = in[i];                                                  \
        vals[j].y = in[i+1];                                                \
    }                                                                       \
    INIT                                                                    \
    DArray_##TYPE##_init(&res, 8);                                          \
    for ( i = 0; i < size; ++i ) {                                          \
        res._used = 0;                                                      \
        PREP                                                                \
        search_##TYPE(PARAM, &tmp, r_sq, &res);                             \
    }                                                                       \
    FREE                                                                    \
    DArray_##TYPE##_free(&res);                                             \
    free(vals);                                                             \
}


#define FAST_DECL   \
    Item *items;    \
    Node *head;
#define FAST_INIT(INSERT)                       \
    items = xmalloc(sizeof(Item)*size);         \
    items = build_morton(vals, items, size);    \
    head = build_tree(items, INSERT);           \
    DArray_Item_init(&tmp, 8);
#define FAST_PREP find_neighbours( items[i].key, head, &tmp );
#define FAST_PARAM (Value *)items[i].val
#define FAST_FREE           \
    DArray_Item_free(&tmp); \
    free(items);

SEARCH_SETUP(fast, Item, FAST_DECL, FAST_INIT(insert_simple), FAST_PREP,
        FAST_PARAM, FAST_FREE)
SEARCH_SETUP(fastfast, Item, FAST_DECL, FAST_INIT(insert_fast), FAST_PREP,
        FAST_PARAM, FAST_FREE)


#define SIMPLE_INIT \
    tmp = (DArray_Value) { (const Value **)&vals, size, size };
#define SIMPLE_PARAM &vals[i]

SEARCH_SETUP(naive, Value, , SIMPLE_INIT, , SIMPLE_PARAM, )



#if 0
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


void search_fast( const unsigned int *in, size_t size, double r_sq )
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
#endif


/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
