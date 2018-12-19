#include "test_data.h"
#include "test.h"


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

    vals = (Value *)munit_parameters_get(params, "values");
    size = *(size_t *)munit_parameters_get(params, "size");

    items = xmalloc(sizeof(Item) * (size+1));
    items = build_morton(vals, items, size);

    return (void *)items;
}

static void
quadtree_teardown(void *data)
{
    Item *items;
    items = (Item *)data;
    free(items);
}

static MunitResult
test_quadtree_build(const MunitParameter params[], void *data)
{
    size_t size, i = 0;
    key_t *exp;
    Item *items;
    Node *head, *tmp;

    size = *(size_t *)munit_parameters_get(params, "size");
    exp = (key_t *)munit_parameters_get(params, "expected");

    items = (Item *)data;
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

    DArray_Node neighbours = { NULL, 0, 0 };
    DArray_Node_init(&neighbours, 8);

    number = *(size_t *)munit_parameters_get(params, "number");
    refk = *(key_t *)munit_parameters_get(params, "reference");
    exp = (key_t *)munit_parameters_get(params, "expected");

    items = (Item *)data;
    head = build_tree(items);
    refn = search(refk, head, maxlvl);

    find_neighbours(refn->key, head, &neighbours);
    res = xmalloc(sizeof(key_t) * neighbours._used);
    for (i = 0; i < neighbours._used; ++i)
        res[i] = DArray_Node_get(&neighbours, i)->key;

    cleanup(head);

    assert_memory_equal(sizeof(key_t)*neighbours._used, (void *)res, (void *)exp);

    DArray_Node_free(&neighbours);
    free(res);

    return MUNIT_OK;
}


/*********************************************************************/


int main(int argc, char *argv[MUNIT_ARRAY_PARAM(argc+1)])
{
    return munit_suite_main(&test_suite, (void *)"user_data", argc, argv);
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
