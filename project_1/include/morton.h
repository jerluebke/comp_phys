#pragma once

#include "types.h"

/* TODO: allow 32- and 64-bit integers as keys */

Item *build_morton( const Value *, Item *, size_t );


/********************************************************************/
/* Implementation of extern inline functions                        */
/********************************************************************/

extern const uint16_t B[];
extern const uint16_t S[];

typedef struct {
    uint8_t x, y;
} Coords2d_8bit;


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
 *
 */
inline uint16_t split2( uint16_t x )
{
    x = (x | (x << S[2])) & B[2];
    x = (x | (x << S[1])) & B[1];
    x = (x | (x << S[0])) & B[0];
    return x;
}

/* interleave
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
 *
 */
inline uint16_t interleave( uint8_t x, uint8_t y )
{
    return split2(x) | (split2(y) << 1);
}

/* decode
 * from a given integer x extract every second bit (starting at index 0)
 * i.e. x coordinate from morton key; to retreive the y coordinate, first
 *     right-shift by one (k >> 1)
 *
 * Params
 * ======
 * k, uint16_t     :   morton key to decode
 *
 * Returns
 * =======
 * uint8_t, x component of given key
 *
 */
inline uint8_t decode(uint16_t k)
{
    k &= B[0];
    k = (k ^ (k >> S[0])) & B[1];
    k = (k ^ (k >> S[1])) & B[2];
    k = (k ^ (k >> S[2]));
    return k;
}

/* coords2
 * from a given key extract x- and y-component and return as struct Coords2d_8bit:
 *     typedef struct {
 *         uint8_t x, y
 *     } Coords2d_8bit;
 *
 * Params
 * ======
 * k, uint16_t     :   morton key to decode
 *
 * Returns
 * =======
 * Coords2d_8bit, struct holding decoded x and y component of given key
 *
 */
inline Coords2d_8bit coords2(uint16_t k)
{
    Coords2d_8bit c = { .x = decode(k), .y = decode(k >> 1) };
    return c;
}


inline uint16_t left( uint16_t key )
{
    return (((key & 0x5555) - 1) & 0x5555) | (key & 0xAAAA);
}

inline uint16_t right( uint16_t key )
{
    return (((key | 0xAAAA) + 1) & 0x5555) | (key & 0xAAAA);
}

inline uint16_t top( uint16_t key )
{
    return (((key & 0xAAAA) - 1) & 0xAAAA) | (key & 0x5555);
}

inline uint16_t bot( uint16_t key )
{
    return (((key | 0x5555) + 1) & 0xAAAA) | (key & 0x5555);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
