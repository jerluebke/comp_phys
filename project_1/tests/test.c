#include "test_data.h"

#define MAX(a, b) (a > b ? a : b)


/********************/
/*      MORTON      */
/********************/

static MunitResult
test_morton_build(const MunitParameter params[], void *data)
{
    (void ) data;

    size_t i;

    TreeInput *in   = (TreeInput *)munit_parameters_get(params, "input");
    Value *vals     = in->vals;
    size_t size     = in->size;
    key_t *exp      = in->exp;

    Item items[size+1], *res;
    key_t keys[size];

    res = build_morton( vals, items, size );
    fprintf(stderr, "res\texpected\n\n");
    for ( i = 0; i < size; ++i ) {
        keys[i] = res[i].key;
        fprintf(stderr, "0x%X\t0x%X\n", keys[i], exp[i]);
    }

    assert_memory_equal( sizeof(key_t)*size,
                         (void *)keys, (void *)exp );

    return MUNIT_OK;
}


#define TEST_MORTON_DIRECTION(DIR)                                          \
static MunitResult                                                          \
test_morton_##DIR(const MunitParameter params[], void *data)                \
{                                                                           \
    (void) data;                                                            \
    key_t given, out, exp;                                                  \
    KeyValueInput                                                           \
        *inp    = (KeyValueInput *)munit_parameters_get(params, "input");   \
    given       = inp->key;                                                 \
    exp         = *(inp->val);                                              \
    out         = DIR(given);                                               \
    assert_ullong(out, ==, exp);                                            \
    return MUNIT_OK;                                                        \
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

    munit_log(MUNIT_LOG_DEBUG, "setting up quadtree...");

    TreeInput *in   = (TreeInput *)munit_parameters_get(params, "setup");
    Value *vals     = in->vals;
    size_t size     = in->size;

    Item *items     = xmalloc( sizeof(Item) * size );
    items           = build_morton( vals, items, size );

    key_t *res      = xmalloc( sizeof(key_t) );

    DArray_Item
        *neighbours = xmalloc( sizeof(DArray_Item) );
    DArray_Item_init(neighbours, 8);

    DataStruct
        *data_ptr   = xmalloc( sizeof(DataStruct) );
    data_ptr->i     = items;
    data_ptr->k     = res;
    data_ptr->da    = neighbours;
    data_ptr->s     = size;

    munit_log(MUNIT_LOG_DEBUG, "done setting up quadtree.");

    return (void *)data_ptr;
}

static void
quadtree_teardown(void *data)
{
    DataStruct *dp  = (DataStruct *)data;
    DArray_Item *ns = dp->da;

    free(dp->i);
    free(dp->k);
    DArray_Item_free(ns);
    free(ns);
    free(dp);
}

static MunitResult
test_quadtree_build(const MunitParameter params[], void *data)
{
    TreeInput *in   = (TreeInput *)munit_parameters_get(params, "setup");
    key_t *exp      = in->exp;

    insert_fptr_t ifunc;
    const char *ifunc_str = munit_parameters_get(params, "ifunc");
    switch ( ifunc_str[0] ) {
        case 's': ifunc = insert_simple; break;
        default : ifunc = insert_fast;  /* f */
    }

    DataStruct *dp  = (DataStruct *)data;
    size_t size     = dp->s;
    Item *items     = dp->i;
    munit_log(MUNIT_LOG_DEBUG, "building tree...");
    Node *head      = build_tree(items, ifunc);
    munit_log(MUNIT_LOG_DEBUG, "done building tree.");

    key_t leaf_keys[size-1];
    Node *tmp;
    size_t i = 0;
    while ( 1 ) {
        munit_logf(MUNIT_LOG_DEBUG, "searching for key 0x%X", items->key);
        tmp = search( items->key, head, maxlvl );
        munit_logf(MUNIT_LOG_DEBUG, "found! actual key: 0x%X", tmp->key);
        leaf_keys[i++] = tmp->key;
        if ( items->last )
            break;
        ++items;
    }

    munit_log(MUNIT_LOG_DEBUG, "deleting tree...");
    cleanup(head);
    munit_log(MUNIT_LOG_DEBUG, "done deleting tree.");

    assert_memory_equal( sizeof(key_t)*size,
                         (void *)leaf_keys, (void *)exp );

    return MUNIT_OK;
}


static MunitResult
test_quadtree_neighbours(const MunitParameter params[], void *data)
{
    size_t i;
    Node *tmp;

    KeyValueInput *in   = (KeyValueInput *)munit_parameters_get(params, "input");
    key_t refk          = in->key;
    key_t num           = *in->val++;
    key_t *exp          = in->val;

    DataStruct *dp      = (DataStruct *)data;
    DArray_Item
        *neighbours     = dp->da;
    key_t *res          = dp->k;
    Item *items         = dp->i;
    Node *head          = build_tree(items, insert_fast);
    Node *refn          = search(refk, head, maxlvl);
    munit_logf(MUNIT_LOG_DEBUG,
               "\ngiven key: 0x%X\nfound key: 0x%X\n", refk, refn->key);

    find_neighbours( refk, head, neighbours );

    res = xrealloc( res, sizeof(key_t) * (neighbours->_used+1) );
    for ( i = 0; i < neighbours->_used; ++i ) {
        tmp = search(neighbours->p[i]->key, head, maxlvl);
        res[i] = tmp->key;
    }

    fprintf(stderr, "res\texpected\n\n");
    for ( i = 0; i < MAX(neighbours->_used, num); ++i )
        fprintf(stderr, "0x%X\t0x%X\n", res[i], exp[i]);

    cleanup(head);

    assert_ullong(neighbours->_used, ==, num);
    if ( num == 0)
        return MUNIT_OK;

    assert_memory_equal( sizeof(key_t)*neighbours->_used,
                         (void *)res, (void *)exp );

    return MUNIT_OK;
}


/*********************************************************************/


/****************************/
/*      MAIN and SUITE      */
/****************************/

#define TEST_MORTON_DIRECTION_CONFIG(DIR) \
    { "/test_morton_"#DIR, test_morton_##DIR, NULL, NULL, \
        MUNIT_TEST_OPTION_NONE, morton_##DIR##_params }

MunitTest tests[] = {
    { "/test_morton_build", test_morton_build, NULL, NULL,
        MUNIT_TEST_OPTION_NONE, morton_build_params},
    TEST_MORTON_DIRECTION_CONFIG(left),
    TEST_MORTON_DIRECTION_CONFIG(right),
    TEST_MORTON_DIRECTION_CONFIG(top),
    TEST_MORTON_DIRECTION_CONFIG(bot),

    { "/test_quadtree_build", test_quadtree_build, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_build_params_1 },
    { "/test_quadtree_neighbours", test_quadtree_neighbours, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_neighbours_params_1 },

    { "/test_quadtree_build", test_quadtree_build, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_build_params_2 },
    { "/test_quadtree_neighbours", test_quadtree_neighbours, quadtree_setup,
        quadtree_teardown, MUNIT_TEST_OPTION_NONE, quadtree_neighbours_params_2 },

    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};


static const MunitSuite test_suite = {
    "/tests", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE
};


int main(int argc, char *argv[])
{
    return munit_suite_main(&test_suite, NULL, argc, argv);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
