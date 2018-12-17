#pragma once

#include "types.h"

/* TODO: allow 32- and 64-bit integers as keys */

Item *build_morton( const Value *, Item *, size_t );
inline key_t left( key_t );
inline key_t right( key_t );
inline key_t top( key_t );
inline key_t bot( key_t );

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
