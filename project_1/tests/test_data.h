#include "test.h"


/**********/
/* MORTON */
/**********/

Value morton_build_values[] = {
    { 1, 1}
};

key_t morton_build_expected[] = {
    0x2
};

size_t morton_build_size[] = {
    10
};

static MunitParameterEnum morton_build_params[] = {
    { "values", (char **)&morton_build_values },
    { "expected", (char **)&morton_build_expected},
    { "size", (char **)&morton_build_size },
    { NULL, NULL },
};

// TODO

/* vim: set ff=dos tw=79 sw=4 ts=4 et ic ai : */
