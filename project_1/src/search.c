/* #include "types.h" */

#define SQUARE(x) (x)*(x)
#define SD(v, w, a) SQUARE((v)->(a) - (w)->(a))     /* squared difference */
#define METRIC(v, w) (SD(v, w, x) + SD(v, w, y))

#define SEARCH_FUNC(TYPE, CHECK_QUERY, IS_VALUE) \
int search_##TYPE( Value *query, DArray_##TYPE *vals, double r, DArray_##TYPE *res ) \
{ \
    int num; \
    TYPE##Iterator *it, *end; \
    it = DArray_##TYPE##_start(vals); \
    end = DArray_##TYPE##_end(vals); \
    for ( num = 0; it != end; ++it ) { \
        CHECK_QUERY \
        if ( METRIC(query, (*it)IS_VALUE) < r ) { \
            DArray_##TYPE##_append(res, *it); \
            ++num; \
        } \
    } \
    return num; \
}

/* #define CHECK_QUERY if ( *it == query ) continue;
#define IS_VALUE ->val */
SEARCH_FUNC(Value, if(*it==query)continue;, ->val)

/* #define CHECK_QUERY
#define IS_VALUE */
SEARCH_FUNC(Item,,)


#if 0
int search_candidates( Value *query, DArray_Item *candidates, double r, DArray_Item *res )
{
    int num;
    ItemIterator *it, *end;

    it = DArray_Item_start(candidates);
    end = DArray_Item_end(candidates);
    for ( num = 0; it != end; ++it ) {
        if ( METRIC(query, (*it)->val) < r ) {
            DArray_Item_append(res, *it);
            ++num;
        }
    }

    return num;
}
#endif
