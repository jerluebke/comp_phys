#pragma once

#include "types.h"

/* TODO: allow 32- and 64-bit integers as keys */

Item *build_morton( const Value *, Item *, size_t );


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
