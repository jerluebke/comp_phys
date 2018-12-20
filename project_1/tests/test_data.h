#include "test.h"

/* TODO:
 *     fill test values
 *     consider order in memory */

static const size_t __any_values_1_size = 10;
static Value __any_values_1[] = {

};


#define BUILD_INPUT(WHICH, NO) \
static TreeInput __##WHICH##_build_input_##NO[] = { \
    { __any_values_##NO, __any_values_##NO##_size, __##WHICH##_build_expected_##NO } \
};

/*********************************************************************/
/*********************************************************************/

/**********/
/* MORTON */
/**********/

static key_t __morton_build_expected_1[] = {

};


BUILD_INPUT(morton, 1)

static char *morton_build_input[] = {
    (char *)__morton_build_input_1,
    // (char *)__morton_build_input_2,
    // etc.
    NULL
};

static MunitParameterEnum morton_build_params[] = {
    { "input", morton_build_input },
    { NULL, NULL },
};


/*********************************************************************/
/*********************************************************************/

static const key_t __morton_any_key_1 = 0x1;

static key_t __morton_left_val_1[] = {

};

static key_t __morton_right_val_1[] = {

};

static key_t __morton_top_val_1[] = {

};

static key_t __morton_bot_val_1[] = {

};


#define MORTON_DIRECTION_ELEM(DIR, NO) \
static KeyValueInput __morton_##DIR##_input_##NO[] = { \
    { __morton_any_key_##NO, __morton_##DIR##_val_##NO } \
};

#define MORTON_DIRECTION_PARAMS(DIR) \
MORTON_DIRECTION_ELEM(DIR, 1) \
/* MORTON_DIRECTION_ELEM(DIR, 2), etc... */ \
static char *morton_##DIR##_kv[] = { \
    (char *)__morton_##DIR##_input_1, \
    /* (char *)__morton_##DIR##_input_2, */ \
    NULL \
}; \
\
static MunitParameterEnum morton_##DIR##_params[] = { \
    { "input", morton_##DIR##_kv }, \
    { NULL, NULL } \
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

static key_t __quadtree_build_expected_1[] = {

};


BUILD_INPUT(quadtree, 1)

static char *quadtree_build_input[] = {
    (char *)__quadtree_build_input_1,
    /* etc. */
    NULL
};

static MunitParameterEnum quadtree_build_params[] = {
    { "setup", quadtree_build_input },
    {  NULL, NULL }
};


/*********************************************************************/
/*********************************************************************/

static const key_t __quadtree_neighbours_ref_1 = 0x1;

static key_t __quadtree_neighbours_exp_1[] = {

};


#define QUADTREE_NEIGHBOURS_INPUT(NO) \
static KeyValueInput __quadtree_neighbours_input_##NO[] = { \
    { __quadtree_neighbours_ref_##NO, __quadtree_neighbours_exp_##NO} \
};

QUADTREE_NEIGHBOURS_INPUT(1)

static char *quadtree_neighbours_input[] = {
    (char *)__quadtree_neighbours_input_1,
    /* etc. */
    NULL
};

static MunitParameterEnum quadtree_neighbours_params[] = {
    { "setup", quadtree_build_input },
    { "input", quadtree_neighbours_input },
    { NULL, NULL }
};


/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
