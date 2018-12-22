#include <assert.h>
#include "../include/quadtree.h"
#include "../include/morton.h"


/* lookup table for msb */
static const lvl_t log_table[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};


/* lookup tables for searching neighbours
 *
 *     index   |   direction
 *     --------+---------------
 *         0   |   left
 *         1   |   right
 *         2   |   top
 *         3   |   bottom
 *         4   |   top-left
 *         5   |   top-right
 *         6   |   bottom-left
 *         7   |   bottom-right
 *
 * bnds: boundaries
 * suffixes: in which direction to go when searching children of neighbours
 *
 * 0xDEAD is terminator
 *
 */

/* check bounds with `(key & bnds[i][0]) == bnds[i][1]` */
static const key_t bnds[8][2] = {
    { 0x5555, 0x0 },    { 0x5555, 0x5555 },
    { 0xAAAA, 0x0 },    { 0xAAAA, 0xAAAA },
    { 0xFFFF, 0x0 },    { 0xFFFF, 0x5555 },
    { 0xFFFF, 0xAAAA }, { 0xFFFF, 0xFFFF}
};

static const key_t suffixes[8][3] = {
    { 0x1, 0x3, 0xDEAD },   { 0x0, 0x2, 0xDEAD },
    { 0x2, 0x3, 0xDEAD },   { 0x0, 0x1, 0xDEAD },
    { 0x3, 0xDEAD },        { 0x2, 0xDEAD },
    { 0x1, 0xDEAD },        { 0x0, 0xDEAD }
};



/* msb - most significant bit
 *
 * lookup table method, see:
 *      http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
 *
 * Params
 * ======
 * k, key_t    :   Value of which to get the msb
 *
 * Returns
 * =======
 * msb of k, i.e. log_2( floor(k) )
 *
 */
static inline lvl_t msb( key_t k )
{
    key_t t;    /* temporaries */
    return (t = k >> 8) ? 8 + log_table[t] : log_table[k];
}


/* bap - bits at position
 *
 * Params
 * ======
 * k, key_t    :   key of which to read the bits
 * j, lvl_t    :   position in key at which to read the bits
 * l, lvl_t    :   length of key / DIM, i.e. maxlevel (might differ when
 *                 considering incomplete branches)
 *
 * Returns
 * =======
 * bits x_j, y_j at position j in key
 *
 */
static inline key_t bap( key_t k, lvl_t j, lvl_t l )
{
    return ( k >> DIM * (l - j) ) & MASK;
}


/* build_branch
 * build whole branch until the node specified by key is reached, starting
 * at head with a depth of nl
 *
 * Params
 * ======
 * head, Node *    :   node at which to start building the branch
 * nl, lvl_t       :   number of new levels
 * key, key_t      :   key of final node of new branch
 * val, Value *    :   content of final node
 *
 */
static Node *build_branch( lvl_t cl, lvl_t nl, key_t key, const Value *val )
{
    lvl_t i, j;
    Node *nn;           /* new nodes */
    DArray_Value *va;   /* value array for new node */

    nn = xmalloc(sizeof(Node) * nl);

    /* set children references between new Nodes */
    /* head = nn; */
    for ( i = 0; i < nl-1; ++i ) {
        /* allocate and init children */
        nn[i].c = xmalloc(sizeof(Node *) * NOC);
        for ( j = 0; j < NOC; ++j ) nn[i].c[j] = NULL;

        /* which direction is the next Node? -> look at x_j, y_j of key,
         * where j = head->lvl + i + 1;
         * head->lvl + i: level of current Node (i is 0-indexed)
         * +1 gives next level */
        nn[i].c[bap(key, cl+i+1, maxlvl)] = &nn[i+1];
    }
    nn[nl-1].c = NULL;

    /* init remaining attributes */
    for ( i = 0; i < nl; ++i ) {
        nn[i].val_arr   = NULL;
        nn[i].lvl       = cl + i;   /* TODO: works correctly for larger branches? */
        nn[i].key       = key >> DIM * (maxlvl - nn[i].lvl);
        nn[i].allocated = 0;
    }

    /* set value to last of the new Nodes */
    va = xmalloc(sizeof(DArray_Value));
    DArray_Value_init(va, 1);
    DArray_Value_append(va, val);
    nn[nl-1].val_arr = va;

    /* mark nn[0] as beginning of allocated Nodes for cleanup */
    nn[0].allocated = 1;

    return nn;
}


/* scr - search children recursivly
 *
 * Params
 * ======
 * head, Node*         :   node at which to start searching
 * suffix, key_t *     :   relevant search directions, terminated by 0xDEAD
 * res, DArray_Value * :   Array in which to write result
 *
 */
#if FIND_NEIGHBOURS_RETURN_NODE
static void scr( const Node *head, const key_t *suffix, DArray_Node *res )
#else
static void scr( const Node *head, const key_t *suffix, DArray_Value *res )
#endif
{
    const key_t *suffix_orig = suffix;

    /* if ( head->val_arr ) { */
    /*     assert( head->c == NULL ); */
    if ( !head->c ) {
#if FIND_NEIGHBOURS_RETURN_NODE
        DArray_Node_append(res, head);
#else
        DArray_Value_extend(res, head->val_arr);
#endif
    }
    else
        while ( *suffix != 0xDEAD ) {
            if ( head->c[*suffix] )
                scr( head->c[*suffix], suffix_orig, res );
            ++suffix;
        }
}


/* build_tree
 * convenience function: creates root node and inserts each element from given
 * items-array
 *
 * Params
 * ======
 * items, Item*    :   array for items to insert in tree
 *
 * Returns
 * =======
 * Node pointer to root node of newly created tree (don't forget to free after
 *     usage!)
 *
 */
Node *build_tree( const Item *items )
{
    int i;
    Node *head;

    head            = xmalloc(sizeof(Node));
    head->key       = 0;
    head->lvl       = 0;
    head->val_arr   = NULL;
    head->allocated = 1;
    head->c         = xmalloc(sizeof(Node *) * NOC);
    for ( i = 0; i < NOC; ++i ) head->c[i] = NULL;

    while ( !items->last ) {
        insert( head, items );
        ++items;
    }

    return head;
}


/* cleanup
 *
 * Params
 * ======
 * head, Node *    :   Node, whichs subtree to delete;
 *                     is set to NULL afterwards
 *
 */
void cleanup( Node *head )
{
    uint8_t i;

    if ( head->c ) {
        for ( i = 0; i < NOC; ++i )
            if ( head->c[i] )
                cleanup(head->c[i]);
        free(head->c);
    }
    if ( head->val_arr ) {
        DArray_Value_free(head->val_arr);
        free(head->val_arr);
    }
    if ( head->allocated ) {
        free(head);
    }
    head = NULL;
}


/* insert
 *
 * NOTICE: for fast tree building, the key following the currently considered
 * one is examined; the key past the end of the array should therefor be NULL
 * to indicate the end and prevent memory corruption
 *
 * Params
 * ======
 * head, Node *    :   starting Node of tree, in which to insert key
 * keys, key_t *   :   array of keys of which the first one is to be inserted,
 *                     its last element should be NULL
 *
 */
void insert( const Node *head, const Item *items )
{
    lvl_t nl;           /* number of new levels */
    key_t sb, lcl;      /* significant bits x_i, y_i; lowest common level */
    sb = bap(items->key, head->lvl+1, maxlvl);

    /* reached lowest level, Node already exists and is occupied:
     *     repeating keys, save in same Node */
    if ( head->lvl == maxlvl ) {
        DArray_Value_append(head->val_arr, items->val);
    }

    /* Node exists, insert in its child */
    else if ( head->c[sb] != NULL ) {
        insert(head->c[sb], items);
    }

    /* Node does not exist, create whole branch until lowest requiered level */
    else {  /* head->c[sb] == NULL */
        /* LOOK-AHEAD */
        if ( (items+1)->last ) {   /* current item is last */
            lcl = 0;
        } else {
            /* `keys[0] XOR keys[1]` sets to one the bits which differ between current and next key */
            /* lowest common level of current and next key */
            lcl = maxlvl - msb(items[0].key ^ items[1].key) / 2;
        }
        /* number of new levels */
        nl = (lcl - head->lvl > 0) ? lcl - head->lvl : 1;
        head->c[sb] = build_branch( head->lvl+1, nl, items->key, items->val );
    }
}


/* search - traverse tree until node without children or with given key is found
 *
 * Params
 * ======
 * head, Node *    :   node at which to start searching
 * key, key_t      :   key to look for
 * lvl, lvl_t      :   level of searched node, i.e. 2*(length of key)
 *
 * Returns
 * =======
 * Node pointer with the desired key or its deepest existing anchestor
 *
 */
Node *search( key_t key, Node *head, lvl_t lvl )
{
    key_t sb;
    while ( head->c
            && head->c[(sb = bap(key, head->lvl+1, lvl))]
            && lvl != head->lvl ) {
        head = head->c[sb];
    }
    return head;
}


/* find_neighbours
 * to given key, compute keys of potential neighbours.
 * search quadtree for those candidates to see if they exists, otherwise take
 *     their deepest existing ancestor.
 * if the candidates itself has children, search for their end nodes facing into
 *     the direction of the current node (i.e. its neighbours).
 *
 * Params
 * ======
 * key, key_t          :   key of node, whichs neighbours to search for
 * head, Node *        :   head of quadtree to search in
 * res, DArray_Value * :   Value array to write results into
 *                         Its first value is the Node of the given key (i.e. the
 *                             reference node), its neighbours start at index 1
 *
 * NOTICE: You get the Values of neighbouring nodes. Depending on the actual shape
 *     of the quadtree, the distances between the reference and its neighbours might
 *     vary significantly (i.e. values in opposite corners of large nodes).
 *     Perhaps you want to filter out those values which are too far away
 *
 */
#if FIND_NEIGHBOURS_RETURN_NODE
void find_neighbours( key_t key, Node *head, DArray_Node *res )
#else
void find_neighbours( key_t key, Node *head, DArray_Value *res )
#endif
{
    Node *c, *tmp;      /* current, temporary */
    int i;

    /* find current node given by key, actual existing key is c->key */
    c = search( key, head, maxlvl );

    /* candidate keys */
    key_t cand_keys[8] = {
        left(c->key),       right(c->key),
        top(c->key),        bot(c->key),
        top(left(c->key)),  top(right(c->key)),
        bot(left(c->key)),  bot(right(c->key))
    };

    /* overwrite res */
    res->_used = 0;

    /* iterate over directions */
    for ( i = 0; i < 8; ++i ) {
        /* if node is on boundary: skip */
        if ( (c->key & bnds[i][0]) == bnds[i][1] )
            continue;

        /* find neighbour candidate node */
        tmp = search( cand_keys[i], head, c->lvl );
        /* if ( tmp->key == c->key ) */
        /*     continue; */

        /* if it has further children: search them (findings will be written in
         *     res);
         * else: write tmp into res */
        if ( tmp->c )
            scr( tmp, suffixes[i], res );
        /* else if ( tmp->val_arr ) */
        else
#if FIND_NEIGHBOURS_RETURN_NODE
            DArray_Node_append(res, tmp);
#else
            DArray_Value_extend(res, tmp->val_arr);
#endif
    }
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
