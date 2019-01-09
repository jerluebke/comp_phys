/* This file contains functions which wrap all necessary steps of setting up the
 * requiered data structures and perform the neighbour search for each given item
 *
 * Used for performance and memory testing.
 */
#include "search.h"


void print_num( const Value *val, int n )
{
    printf("(%u, %u)\tfound: %d\n", val->x, val->y, n);
}

#define DO_PRINT 0
#if DO_PRINT
#define PRINT_NUM(v, n) print_num(v, n);
#define PRINT_INFO(s) printf("\n\n"#s"\nquery\tfound\n");
#else
#define PRINT_NUM(v, n) n;
#define PRINT_INFO(s)
#endif


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

#define EMPTY
#define CHECK_QUERY if (*it==query) continue;
#define VALUE_ACCESS ->val
SEARCH_FUNC(Value, CHECK_QUERY, EMPTY)
SEARCH_FUNC(Item, EMPTY, VALUE_ACCESS)


#define SEARCH_SETUP(NAME, TYPE, DECL, INIT, PREP, PARAM, FREE)             \
int search_##NAME( const fargs_t *fargs )                                   \
{                                                                           \
    const unsigned int *in = fargs->data;                                   \
    size_t size = fargs->size;                                              \
    double r_sq = fargs->r_sq;                                              \
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
    PRINT_INFO(NAME)                                                        \
    for ( i = 0; i < size; ++i ) {                                          \
        res._used = 0;                                                      \
        PREP                                                                \
        PRINT_NUM(PARAMS, search_##TYPE(PARAM, &tmp, r_sq, &res))           \
    }                                                                       \
    FREE                                                                    \
    DArray_##TYPE##_free(&res);                                             \
    free(vals);                                                             \
    return 0;                                                               \
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
    cleanup(head);          \
    free(items);

SEARCH_SETUP(fast, Item, FAST_DECL, FAST_INIT(insert_simple), FAST_PREP,
        FAST_PARAM, FAST_FREE)
SEARCH_SETUP(fastfast, Item, FAST_DECL, FAST_INIT(insert_fast), FAST_PREP,
        FAST_PARAM, FAST_FREE)


#define SIMPLE_INIT                 \
    DArray_Value_init(&tmp, size);  \
    for ( i = 0; i < size; ++i )    \
        tmp.p[i] = &vals[i];        \
    tmp._used = size;
#define SIMPLE_PARAM &vals[i]
#define SIMPLE_FREE DArray_Value_free(&tmp);

SEARCH_SETUP(naive, Value, EMPTY, SIMPLE_INIT, EMPTY, SIMPLE_PARAM, SIMPLE_FREE)

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
