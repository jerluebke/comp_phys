#include "types.h"

#define SQUARE(x) (x)*(x)
#define SD(v, w, a) SQUARE((v)->(a) - (w)->(a))     /* squared difference */
#define METRIC(v, w) (SD(v, w, x) + SD(v, w, y))

#define SEARCH_FUNC(TYPE, CHECK_QUERY, IS_VALUE)                            \
int search_##TYPE( Value *query, DArray_##TYPE *vals, double r_sq,          \
                   DArray_##TYPE *res )                                     \
{                                                                           \
    int num;                                                                \
    TYPE##Iterator *it, *end;                                               \
    it = DArray_##TYPE##_start(vals);                                       \
    end = DArray_##TYPE##_end(vals);                                        \
    for ( num = 0; it != end; ++it ) {                                      \
        CHECK_QUERY                                                         \
        if ( METRIC(query, (*it)IS_VALUE) < r_sq ) {                        \
            DArray_##TYPE##_append(res, *it);                               \
            ++num;                                                          \
        }                                                                   \
    }                                                                       \
    return num;                                                             \
}

SEARCH_FUNC(Value, if(*it==query)continue;, ->val)
SEARCH_FUNC(Item,,)
