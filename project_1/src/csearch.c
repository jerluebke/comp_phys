#if defined(unix) || defined(__unix__) || defined (__unix)
#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#if !_POSIX_TIMERS
#error This file requieres the POSIX timer implementation, which is not \
    supported on your machine!
#endif
#else
#error This file was written for a unix OS with POSIX support! 
#endif

#include <time.h>
#include <limits.h>
#include <math.h>
#include "types.h"
#include "morton.h"
#include "quadtree.h"

#include "sample_data.h"


typedef struct fargs_s {
    const unsigned int *data;
    size_t size;
    double r_sq;
} fargs_t;


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

#define EMPTY
#define CHECK_QUERY if (*it==query) continue;
#define VALUE_ACCESS ->val
SEARCH_FUNC(Value, CHECK_QUERY, EMPTY)
SEARCH_FUNC(Item, EMPTY, VALUE_ACCESS)


#define SEARCH_SETUP(NAME, TYPE, DECL, INIT, PREP, PARAM, FREE)             \
void search_##NAME( const fargs_t *fargs )                                  \
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

/**/

void timeit(void(*func)(const fargs_t *), const fargs_t *fargs,
            unsigned int iter, char *name)
{
    unsigned int i;
    double avg = 0, std = 0;
    long min = LONG_MAX, max = LONG_MIN;
    long *res = xmalloc(sizeof(long) * iter);
    struct timespec tp_a, tp_b, tp_c;

    clock_gettime(CLOCK_REALTIME, &tp_a);

    for ( i = 0; i < iter; ++i ) {
        clock_gettime(CLOCK_REALTIME, &tp_b);
        (*func)(fargs);
        clock_gettime(CLOCK_REALTIME, &tp_c);

        res[i] = tp_c.tv_nsec - tp_b.tv_nsec;
        avg += (double) res[i];
        if ( res[i] < min )
            min = res[i];
        if ( res[i] > max)
            max = res[i];
    }

    clock_gettime(CLOCK_REALTIME, &tp_c);

    avg /= iter;
    for ( i = 0; i < iter; ++i )
        std += (double)(SQUARE(res[i] - avg));
    std = sqrt(std / ((double) iter-1));

    printf("\nTIMEIT %s - %u runs\n"        \
           "average     :   %lf +- %lf ns\n"\
           "    max/min :   %lu / %lu\n"    \
           "total time  :   %lu ns\n\n",
           name, iter, avg, std, max, min, tp_c.tv_nsec - tp_a.tv_nsec);
}


#if 1
int main()
{
    unsigned int iter = 10u;
    const fargs_t fargs_small = { .data=input_data_256, .size=size_256, .r_sq=16.0f };
    timeit(search_naive, &fargs_small, iter, "naive - 256");
    timeit(search_fast, &fargs_small, iter, "fast - 256");
    timeit(search_fastfast, &fargs_small, iter, "fastfast - 256");

    const fargs_t fargs = { .data=input_data_1265, .size=size_1265, .r_sq=16.0f };
    timeit(search_naive, &fargs, iter, "naive - 1265");
    timeit(search_fastfast, &fargs, iter, "fastfast - 1265");
    return 0;
}
#endif

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
