#include "test_data.h"
#include "test.h"

typedef struct {
    Item *i; key_t *k; DArray_Node *da; size_t s;
} data_struct;


/*********************************************************************/


/********************/
/*      MORTON      */
/********************/

static MunitResult
test_morton_build(const MunitParameter params[], void *data)
{
   (void ) data;

    Value *vals;
    key_t *exp;
    size_t size, i;

    vals = (Value *)munit_parameters_get(params, "values");
    exp = (key_t *)munit_parameters_get(params, "expected");
    size = *(size_t *)munit_parameters_get(params, "size");

    Item items[size+1], *res;
    key_t keys[size];

    res = build_morton(vals, items, size);
    for ( i = 0; i < size; ++i)
        keys[i] = res[i].key;

    assert_memory_equal(sizeof(key_t)*size, (void *)keys, (void *)exp);

    return MUNIT_OK;
}


#define TEST_MORTON_DIRECTION(DIR) \
static MunitResult \
test_morton_##DIR(const MunitParameter params[], void *data) \
{ \
    (void) data; \
    key_t in, out, exp; \
    in = *(key_t *)munit_parameters_get(params, "in"); \
    exp = *(key_t *)munit_parameters_get(params, "exp"); \
    out = DIR(in); \
    assert_ullong(out, ==, exp); \
    return MUNIT_OK; \
}

TEST_MORTON_DIRECTION(left)
TEST_MORTON_DIRECTION(right)
TEST_MORTON_DIRECTION(top)
TEST_MORTON_DIRECTION(bot)


/*********************************************************************/


/********************/
/*      QUADTREE    */
/********************/

static void *
quadtree_setup(const MunitParameter params[], void *data)
{
    (void) data;

    Value *vals;
    size_t size;
    Item *items;
    key_t *res;
    DArray_Node *neighbours;
    data_struct *data_ptr;

    vals = (Value *)munit_parameters_get(params, "values");
    size = *(size_t *)munit_parameters_get(params, "size");

    items = xmalloc(sizeof(Item) * (size+1));
    items = build_morton(vals, items, size);

    res = xmalloc(sizeof(key_t));

    neighbours = xmalloc(sizeof(DArray_Node));
    DArray_Node_init(neighbours, 8);

    data_ptr = xmalloc(sizeof(data_struct));
    data_ptr->i = items;
    data_ptr->k = res;
    data_ptr->da = neighbours;
    data_ptr->s = size;

    return (void *)data_ptr;
}

static void
quadtree_teardown(void *data)
{
    DArray_Node *ns;
    data_struct *dp;

    dp = (data_struct *)data;
    ns = dp->da;

    free(dp->i);
    free(dp->k);
    DArray_Node_free(ns);
    free(ns);
    free(dp);
}

static MunitResult
test_quadtree_build(const MunitParameter params[], void *data)
{
    size_t size, i = 0;
    key_t *exp;
    Item *items;
    Node *head, *tmp;
    data_struct *dp;

    exp = (key_t *)munit_parameters_get(params, "expected");

    dp = (data_struct *)data;
    size = dp->s;
    items = dp->i;
    head = build_tree(items);

    key_t leaf_keys[size];
    while ( items != NULL ) {
        tmp = search(items->key, head, maxlvl);
        leaf_keys[i++] = tmp->key;
    }

    cleanup(head);

    assert_memory_equal(sizeof(key_t)*size, (void *)leaf_keys, (void *)exp);

    return MUNIT_OK;
}


static MunitResult
test_quadtree_neighbours(const MunitParameter params[], void *data)
{
    size_t number, i;
    key_t refk, *exp, *res;
    Item *items;
    Node *head, *refn;
    DArray_Node *neighbours;
    data_struct *dp;

    refk = *(key_t *)munit_parameters_get(params, "reference");
    exp = (key_t *)munit_parameters_get(params, "expected");

    dp = (data_struct *)data;
    number = dp->s;     /* number of neighbours */
    neighbours = dp->da;
    res = dp->k;
    items = dp->i;
    head = build_tree(items);
    refn = search(refk, head, maxlvl);

    find_neighbours(refn->key, head, neighbours);
    res = xrealloc(res, sizeof(key_t) * neighbours->_used);
    for (i = 0; i < neighbours->_used; ++i)
        res[i] = DArray_Node_get(neighbours, i)->key;

    cleanup(head);

    assert_memory_equal(sizeof(key_t)*neighbours->_used, (void *)res, (void *)exp);

    return MUNIT_OK;
}


/*********************************************************************/


/****************************/
/*      MAIN and SUITES     */
/****************************/

#define TEST_MORTON_DIRECTION_CONFIG(DIR) \
    { "/test_morton_##DIR", test_morton_##DIR, NULL, NULL, \
        MUNIT_TEST_OPTION_NONE, morton_##DIR##_params }

MunitTest tests[] = {
    { "/test_morton_build", test_morton_build, NULL, NULL,
        MUNIT_TEST_OPTION_NONE, morton_build_params},
    TEST_MORTON_DIRECTION_CONFIG(left),
    TEST_MORTON_DIRECTION_CONFIG(right),
    TEST_MORTON_DIRECTION_CONFIG(top),
    TEST_MORTON_DIRECTION_CONFIG(bot),

    { "/test_quadtree_build", test_quadtree_build, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_build_params },
    { "/test_quadtree_neighbours", test_quadtree_neighbours, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_neighbours_params },

    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};


static const MunitSuite test_suite = {
    "/quadtree and morton tests", tests, NULL, 10, MUNIT_SUITE_OPTION_NONE
};


int main(int argc, char *argv[MUNIT_ARRAY_PARAM(argc+1)])
{
    return munit_suite_main(&test_suite, NULL, argc, argv);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
