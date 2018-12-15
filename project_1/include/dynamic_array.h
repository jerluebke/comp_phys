#pragma once

#include "xmalloc.h"

#define FACTOR 2ul

#ifdef __cplusplus
extern "C" {
#endif


#define DARRAY(TYPE, NAME)  \
DARRAY_STRUCT(TYPE, NAME)   \
DARRAY_GET(TYPE, NAME)      \
DARRAY_INIT(TYPE, NAME)     \
DARRAY_INSERT(TYPE, NAME)   \
DARRAY_FREE(NAME)


/* struct DArray
 * wrapper for regular array with dynamic insertion
 * for reference see: https://stackoverflow.com/a/3536261/9133910
 *
 * Members
 * =======
 * _data, TYPE *   :   array of data
 * _used, size_t   :   space in use
 * _size, size_t   :   space available
 *  */
#define DARRAY_STRUCT(TYPE, NAME)                                           \
typedef struct {                                                            \
    TYPE *_data;                                                            \
    size_t _used, _size;                                                    \
} DArray_##NAME;                                                            \


/* DArray_get
 * subscript wrapper to access some DArray's data
 * NO BOUNDCHECKING!
 *
 * Params
 * ======
 * da, DArray *    :   DArray whichs _data to access
 * i, size_t       :   index to access
 *  */
#define DARRAY_GET(TYPE, NAME)                                              \
TYPE DArray_##NAME##_get (DArray_##NAME *da, size_t i)                      \
{                                                                           \
    return da->_data[i];                                                    \
}                                                                           \


/* DArray_init
 * allocate _data to given size
 *
 * Params
 * ======
 * da, DArray *    :   DArray to initialize
 * initial, size_t :   its initial size
 *  */
#define DARRAY_INIT(TYPE, NAME)                                             \
void DArray_##NAME##_init(DArray_##NAME *da, size_t initial)                \
{                                                                           \
    da->_data = xmalloc(sizeof(TYPE) * initial);                            \
    da->_used = 0;                                                          \
    da->_size = initial;                                                    \
}


/* DArray_insert
 * insert value of TYPE TYPE at the end of given DArray
 * realloc its _data if necessary by factor FACTOR
 *
 * Params
 * ======
 * da, DArray *    :   DArray in which to insert
 * elem, TYPE      :   value to insert
 *  */
#define DARRAY_INSERT(TYPE, NAME)                                           \
void DArray_##NAME##_insert(DArray_##NAME *da, TYPE elem)                   \
{                                                                           \
    if ( da->_used == da->_size )  {                                        \
        da->_size *= FACTOR;                                                \
        da->_data = xrealloc(da->_data, sizeof(TYPE) * da->_size);          \
    }                                                                       \
    da->_data[da->_used++] = elem;  /* set elem and increment da->_used */  \
}


/* DArray_free
 * free and set NULL _data of given DArray
 * set _used and _size to 0
 *
 * Params
 * ======
 * da, DArray *    :   DArray to free
 *  */
#define DARRAY_FREE(NAME)                                                   \
void DArray_##NAME##_free(DArray_##NAME *da)                                \
{                                                                           \
    free(da->_data);                                                        \
    da->_data = NULL;                                                       \
    da->_used = da->_size = 0;                                              \
}


#ifdef __cpluplus
}
#endif
