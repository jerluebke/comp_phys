#pragma once

#include "types.h"

/* TODO: allow 32- and 64-bit integers as keys */

Key *build_morton( const Value *, Key *, size_t );
inline key_t left( key_t );
inline key_t right( key_t );
inline key_t top( key_t );
inline key_t bot( key_t );
