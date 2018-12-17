/* for the splitting and interleaving, see:
 *     http://www-graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
 *     https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
 *
 * for calculating left, right, top, bottom and for more information regarding
 * morton keys, see:
 *     https://en.wikipedia.org/wiki/Z-order_curve#Coordinate_values
 *  */
#include "../include/morton.h"

/* lookup tables for `split2` */
static const uint16_t B[] = {0x5555, 0x3333, 0x0F0F};
static const uint16_t S[] = {1, 2, 4};


/* split2
 * split up binary representation of given x such that between each two bits
 * a 0 is inserted (right bound), e.g.: 0b000111 -> 0b010101
 *
 * Params
 * ======
 * x, uint16_t     :   number to split up, must be less than 256
 *
 * Returns
 * =======
 * uint16_t
 *  */
static inline uint16_t split2( uint16_t x )
{
    x = (x | (x << S[2])) & B[2];
    x = (x | (x << S[1])) & B[1];
    x = (x | (x << S[0])) & B[0];
    return x;
}


/* interleave8
 * interleave two 8-bit integers by splitting up their binary representation
 * and shuffeling their bits, i.e. y7 x7 ... y1 x1 y0 x0
 *
 * Params
 * ======
 * x, y (uint8_t)  :   numbers to interleave
 *
 * Returns
 * =======
 * uint16_t
 *  */
static inline uint16_t interleave8( uint8_t x, uint8_t y )
{
    return split2(x) | (split2(y) << 1);
}


/* cmp_keys
 * compare two pointers to Item structs by comparing their interger-keys
 *
 * passed to `qsort`
 *
 * Params
 * ======
 * a, b (const void *) :   casted to `const Item *` and compared
 *
 * Returns
 * =======
 *  1, if a > b
 * -1, if a < b
 *  0, if a == b
 *  */
static inline int cmp_keys( const void *a, const void *b )
{
    const Item *arg1, *arg2;
    arg1 = (const Item *)a;
    arg2 = (const Item *)b;

    return (arg1->key > arg2->key) - (arg1->key < arg2->key);
}


/* build_morton
 * calculates and sortes morton keys for given data
 *
 * Params
 * ======
 * vals, Value *   :   values to hash with `size` elements
 * keys, Item *     :   array to write result with `size+1` elements
 * size, size_t    :   size of `vals`, `keys` needs to have size+1 (for
 *                      terminating character, for look-ahead during tree
 *                      building)
 *
 * Returns
 * =======
 * sorted array of keys
 *  */
Item *build_morton( const Value *vals, Item *items, size_t size )
{
    size_t i;
    Item *end;

    /* calculate keys from coordinates and set reference to corresp. value */
    for ( i = 0; i < size; ++i ) {
        items[i].key = interleave8(vals[i].x, vals[i].y);
        items[i].val = &vals[i];
    }

    qsort(items, size, sizeof(Item), cmp_keys);

    /* set last element to NULL (terminating character) */
    end = items + size;
    end = NULL;

    return items;
}


/*
 * convenience functions to retreive adjacent keys
 *  */

inline uint16_t left( uint16_t key )
{
    return ((key & 0x5555) - 1 & 0x5555) | (key & 0xAAAA);
}

inline uint16_t right( uint16_t key )
{
    return ((key | 0xAAAA) + 1 & 0x5555) | (key & 0xAAAA);
}

inline uint16_t top( uint16_t key )
{
    return ((key & 0xAAAA) - 1 & 0xAAAA) | (key & 0x5555);
}

inline uint16_t bot( uint16_t key )
{
    return ((key | 0x5555) + 1 & 0xAAAA) | (key & 0x5555);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
