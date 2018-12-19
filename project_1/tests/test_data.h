#include "test.h"

/* TODO: fill params, add terminating NULL */


/**********/
/* MORTON */
/**********/

static Value __morton_build_values_raw[] = {
    { 1, 1 }
};

static key_t __morton_build_expected_raw[] = {
    0x2
};

static size_t __morton_build_size_raw[] = {
    10
};

/*********************************************************************/

static char *morton_build_values[] = {
    (char *)__morton_build_values_raw,
    NULL
};

static char *morton_build_expected[] = {
    (char *)__morton_build_expected_raw,
    NULL
};

static char *morton_build_size[] = {
    (char *)__morton_build_size_raw,
    NULL
};

static MunitParameterEnum morton_build_params[] = {
    { "values",     morton_build_values },
    { "expected",   morton_build_expected},
    { "size",       morton_build_size },
    { NULL, NULL },
};


/*********************************************************************/
/*********************************************************************/


static key_t __morton_left_in_raw[] = {
    0x11
};

static key_t __morton_left_exp_raw[] = {

};


static key_t __morton_right_in_raw[] = {
    0x11
};

static key_t __morton_right_exp_raw[] = {

};


static key_t __morton_top_in_raw[] = {
    0x11
};

static key_t __morton_top_exp_raw[] = {

};


static key_t __morton_bot_in_raw[] = {
    0x11
};

static key_t __morton_bot_exp_raw[] = {

};

/*********************************************************************/

#define MORTON_DIRECTION_PARAMS(DIR) \
static char *morton_##DIR##_in[] = { \
    (char *)__morton_##DIR##_in_raw, \
    NULL \
}; \
\
static char *morton_##DIR##_exp[] = { \
    (char *)__morton_##DIR##_exp_raw, \
    NULL \
}; \
\
static MunitParameterEnum morton_##DIR##_params[] = { \
    { "in", morton_##DIR##_in }, \
    { "exp", morton_##DIR##_exp } \
};

MORTON_DIRECTION_PARAMS(left)
MORTON_DIRECTION_PARAMS(right)
MORTON_DIRECTION_PARAMS(top)
MORTON_DIRECTION_PARAMS(bot)


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/


/************/
/* QUADTREE */
/************/

static MunitParameterEnum quadtree_build_params[] = {
// values, size, expected
};

static MunitParameterEnum quadtree_neighbours_params[] = {
// values, size, expected, reference
};


/* vim: set ff=dos tw=79 sw=4 ts=4 et ic ai : */
