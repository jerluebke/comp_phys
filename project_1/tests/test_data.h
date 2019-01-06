#pragma once
#include "test.h"


#define __any_values_1_size 16u
static Value __any_values_1[] = {
    /* 1 */     { 0x1C, 0x68 }, { 0x7F, 0x7F }, { 0x20, 0xA0 }, { 0x68, 0x86 },
    /* 0 */     { 0x0, 0x0 }, { 0x1C, 0x5E }, { 0x2A, 0x54 }, { 0x24, 0x5A },
    /* 2 */     { 0x54, 0x8A }, { 0x55, 0x8A }, { 0x55, 0x8B }, { 0xFF, 0xFF },
    { 0x60, 0x20 }, { 0x80, 0x7F }, { 0x3F, 0xC0 }, { 0x40, 0xFF }
};


#define __any_values_2_size 3u
static Value __any_values_2[] = {
    { 0x40, 0xC0 }, { 0xA0, 0x20 }, { 0xE0, 0x60 }
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
    0x0, 0x1C00, 0x23F8, 0x2664,
    0x2698, 0x29D0, 0x3FFF, 0x6AAA,
    0x8C00, 0x9198, 0x9199, 0x919B,
    0x9468, 0xA555, 0xBAAA, 0xFFFF
};

static key_t __morton_build_expected_2[] = {
    0x4C00, 0x7C00, 0xB000
};


BUILD_INPUT(morton, 1)
BUILD_INPUT(morton, 2)

static char *morton_build_input[] = {
    (char *)__morton_build_input_1,
    (char *)__morton_build_input_2,
    /* etc. */
    NULL
};

static MunitParameterEnum morton_build_params[] = {
    { "input", morton_build_input },
    { NULL, NULL },
};


/*********************************************************************/
/*********************************************************************/

#define __morton_any_key_1  0x3ull
static key_t __morton_left_val_1        = 0x2;
static key_t __morton_right_val_1       = 0x6;
static key_t __morton_top_val_1         = 0x1;
static key_t __morton_bot_val_1         = 0x9;

#define __morton_any_key_2  0x0ull
static key_t __morton_left_val_2        = 0x5555;
static key_t __morton_right_val_2       = 0x1;
static key_t __morton_top_val_2         = 0xAAAA;
static key_t __morton_bot_val_2         = 0x2;

#define __morton_any_key_3  0x3FFFull
static key_t __morton_left_val_3        = 0x3FFE;
static key_t __morton_right_val_3       = 0x6AAA;
static key_t __morton_top_val_3         = 0x3FFD;
static key_t __morton_bot_val_3         = 0x9555;

#define __morton_any_key_4  0x4000ull
static key_t __morton_left_val_4        = 0x1555;
static key_t __morton_right_val_4       = 0x4001;
static key_t __morton_top_val_4         = 0xEAAA;
static key_t __morton_bot_val_4         = 0x4002;


#define MORTON_DIRECTION_ELEM(DIR, NO) \
static KeyValueInput __morton_##DIR##_input_##NO[] = { \
    { __morton_any_key_##NO, &__morton_##DIR##_val_##NO } \
};

#define MORTON_DIRECTION_PARAMS(DIR) \
MORTON_DIRECTION_ELEM(DIR, 1) \
MORTON_DIRECTION_ELEM(DIR, 2) \
MORTON_DIRECTION_ELEM(DIR, 3) \
MORTON_DIRECTION_ELEM(DIR, 4) \
static char *morton_##DIR##_kv[] = { \
    (char *)__morton_##DIR##_input_1, \
    (char *)__morton_##DIR##_input_2, \
    (char *)__morton_##DIR##_input_3, \
    (char *)__morton_##DIR##_input_4, \
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
    0x0, 0x1, 0x8, 0x99,
    0x9A, 0xA, 0x3, 0x1,
    0x8, 0x9198, 0x9199, 0x919B,
    0x25, 0xA, 0xB, 0x3
};


static key_t __quadtree_build_expected_2[] = {
    0x4, 0x7, 0x2
};


BUILD_INPUT(quadtree, 1)
BUILD_INPUT(quadtree, 2)

static char *quadtree_build_input_1[] = {
    (char *)__quadtree_build_input_1,
    /* etc. */
    NULL
};

static char *insert_func_ids[] = {
    "f",
    NULL
};


static MunitParameterEnum quadtree_build_params_1[] = {
    { "setup", quadtree_build_input_1 },
    { "ifunc", insert_func_ids },
    {  NULL, NULL }
};


/* seperate set of params to avoid conflict with neighbour tests (since all
 * permutations of one given set of params are tested) */
static char *quadtree_build_input_2[] = {
    (char *)__quadtree_build_input_2,
    /* etc. */
    NULL
};

static MunitParameterEnum quadtree_build_params_2[] = {
    { "setup", quadtree_build_input_2 },
    { "ifunc", insert_func_ids },
    {  NULL, NULL }
};


/*********************************************************************/
/*********************************************************************/

/* first element is number of expected key */

#define __quadtree_neighbours_ref_1     0x0ull
static key_t __quadtree_neighbours_exp_1[] = {
    3, 0x1, 0x8, 0x3
};

#define __quadtree_neighbours_ref_2     0x23F8ull
static key_t __quadtree_neighbours_exp_2[] = {
    3, 0x9A, 0x0, 0xA
};

#define __quadtree_neighbours_ref_3     0x2664ull
static key_t __quadtree_neighbours_exp_3[] = {
    1, 0x9A
};

#define __quadtree_neighbours_ref_4     0x3FFFull
static key_t __quadtree_neighbours_exp_4[] = {
    6, 0x1, 0x1, 0x25, 0x0, 0x8, 0x3
};

#define __quadtree_neighbours_ref_5     0x9199ull
static key_t __quadtree_neighbours_exp_5[] = {
    2, 0x9198, 0x919B
};

#define __quadtree_neighbours_ref_6     0xFFFFull
static key_t __quadtree_neighbours_exp_6[] = {
    4, 0x25, 0xB, 0x1, 0x3
};

#define __quadtree_neighbours_ref_7     0x9468ull
static key_t __quadtree_neighbours_exp_7[] = {
    3, 0x3, 0x3, 0x1
};

#define __quadtree_neighbours_ref_8     0xA555ull
static key_t __quadtree_neighbours_exp_8[] = {
    2, 0xB, 0x8
};

#define __quadtree_neighbours_ref_9     0x6AAAull
static key_t __quadtree_neighbours_exp_9[] = {
    4, 0x1, 0x3, 0x3, 0x25
};


#define __quadtree_neighbours_ref_10    0x4C00ull
static key_t __quadtree_neighbours_exp_10[] = {
    1, 0x7
};

#define __quadtree_neighbours_ref_11    0xB000ull
static key_t __quadtree_neighbours_exp_11[] = {
    0
};


#define QUADTREE_NEIGHBOURS_INPUT(NO) \
static KeyValueInput __quadtree_neighbours_input_##NO[] = { \
    { __quadtree_neighbours_ref_##NO, __quadtree_neighbours_exp_##NO } \
};

QUADTREE_NEIGHBOURS_INPUT(1)
QUADTREE_NEIGHBOURS_INPUT(2)
QUADTREE_NEIGHBOURS_INPUT(3)
QUADTREE_NEIGHBOURS_INPUT(4)
QUADTREE_NEIGHBOURS_INPUT(5)
QUADTREE_NEIGHBOURS_INPUT(6)
QUADTREE_NEIGHBOURS_INPUT(7)
QUADTREE_NEIGHBOURS_INPUT(8)
QUADTREE_NEIGHBOURS_INPUT(9)

QUADTREE_NEIGHBOURS_INPUT(10)
QUADTREE_NEIGHBOURS_INPUT(11)


/* find_neighbours: general */
static char *quadtree_neighbours_input_1[] = {
    (char *)__quadtree_neighbours_input_1,
    (char *)__quadtree_neighbours_input_2,
    (char *)__quadtree_neighbours_input_3,
    (char *)__quadtree_neighbours_input_4,
    (char *)__quadtree_neighbours_input_5,
    (char *)__quadtree_neighbours_input_6,
    (char *)__quadtree_neighbours_input_7,
    (char *)__quadtree_neighbours_input_8,
    (char *)__quadtree_neighbours_input_9,
    /* etc. */
    NULL
};

static MunitParameterEnum quadtree_neighbours_params_1[] = {
    { "setup", quadtree_build_input_1 },
    { "ifunc", insert_func_ids },
    { "input", quadtree_neighbours_input_1 },
    { NULL, NULL }
};


/* find_neighbours: node without neighbours and multiple items */
static char *quadtree_neighbours_input_2[] = {
    (char *)__quadtree_neighbours_input_10,
    (char *)__quadtree_neighbours_input_11,
    NULL
};

static MunitParameterEnum quadtree_neighbours_params_2[] = {
    { "setup", quadtree_build_input_2 },
    { "ifunc", insert_func_ids },
    { "input", quadtree_neighbours_input_2 },
    { NULL, NULL }
};

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
