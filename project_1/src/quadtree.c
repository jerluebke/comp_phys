#include "../include/quadtree.h"
#include "../include/morton.h"


/* lookup table for msb */
static const lvl_t log_table[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};


struct Value {
    /* TODO */
};

struct Node {
    key_t key;
    DArray_Value *val_arr;
    lvl_t lvl;
    Node **c;           /* children */
    DArray_Node *n_arr; /* neighbours */
    uint8_t allocated;
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
 *  */
static inline lvl_t msb( key_t k )
{
    key_t t;    /* temporaries */
    return (t = k >> 8) ? 8 + log_table[t] : log_table[k];
}


/* bap - bits at position
 *
 * Params
 * ======
 * key, key_t  :   key of which to read the bits
 * j, lvl_t    :   position in key at which to read the bits
 *
 * Returns
 * =======
 * bits x_j, y_j at position j in key
 *  */
static inline key_t bap( key_t key, lvl_t j )
{
    return ( key >> DIM * (maxlvl - j) ) & MASK;
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
 *  */
static void build_branch( const Node *head, lvl_t nl, key_t key, Value *val )
{
    lvl_t i, j;
    Node *nn;           /* new nodes */
    DArray_Value *va;   /* value array for new node */

    nn = xmalloc(sizeof(Node) * nl);

    /* set children references between new Nodes */
    head = nn;
    for ( i = 0; i < nl-1; ++i ) {
        /* allocate and init children */
        nn[i].c = xmalloc(sizeof(Node *) * NOC);
        for ( j = 0; j < NOC; ++j ) nn[i].c[j] = NULL;

        /* which direction is the next Node? -> look at x_j, y_j of key,
         * where j = head->lvl + i + 1;
         * head->lvl + i: level of current Node (i is 0-indexed)
         * +1 gives next level */
        nn[i].c[bap(key, head->lvl+i+1)] = &nn[i+1];
    }
    nn[nl-1].c = NULL;

    /* init remaining attributes */
    for ( i = 0; i < nl; ++i ) {
        nn[i].val_arr   = NULL;
        nn[i].n_arr     = NULL;
        nn[i].lvl       = head->lvl + i + 1;
        nn[i].key       = key >> DIM * (maxlvl - nn[i].lvl);
        nn[i].allocated = 0;
    }

    /* set value to last of the new Nodes */
    va = xmalloc(sizeof(DArray_Value));
    DArray_Value_init(va, 1);
    DArray_Value_insert(va, val);
    nn[nl-1].val_arr = va;

    /* mark nn[0] as beginning of allocated Nodes for cleanup */
    nn[0].allocated = 1;
}


/* insert
 *
 * Params
 * ======
 * head, Node *    :   starting Node of tree, in which to insert key
 * keys, key_t *   :   array of keys of which the first one is to be inserted,
 *                     its last element should be NULL
 *  */
void insert( const Node *head, const key_t *keys, Value *val )
{
    lvl_t nl;           /* number of new levels */
    key_t sb, lcl;      /* significant bits x_i, y_i; lowest common level */
    sb = bap(*keys, head->lvl+1);

    /* reached lowest level, Node already exists and is occupied:
     *     repeating keys, save in same Node */
    if ( head->lvl == maxlvl ) {
        DArray_Value_insert(head->val_arr, val);
    }

    /* Node exists, insert in its child */
    else if ( head->c[sb] != NULL ) {
        insert(head->c[sb], keys, val);
    }

    /* Node does not exist, create whole branch until lowest requiered level */
    else {  /* head->c[sb] == NULL */
        if ( (keys+1) == NULL ) {   /* current key is last */
            lcl = 0;
        } else {
            /* `keys[0] XOR keys[1]` sets to one the bits which differ between current and next key */
            /* lowest common level of current and next key */
            lcl = maxlvl - msb(keys[0] ^ keys[1]) / 2;
        }
        /* number of new levels */
        nl = (lcl - head->lvl > 0) ? lcl - head->lvl : 1;
        build_branch( head->c[sb], nl, *keys, val );
    }
}


/* delete
 *
 * Params
 * ======
 * head, Node *    :   Node, whichs subtree to delete;
 *                     is set to NULL afterwards
 *  */
void delete( Node *head )
{
    uint8_t i;

    if ( head->c ) {
        for ( i = 0; i < NOC; ++i )
            if ( head->c[i] )
                delete(head->c[i]);
        free(head->c);
    }
    if ( head->n_arr ) {
        DArray_Node_free(head->n_arr);
        free(head->n_arr);
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


/* search - traverse tree until node without children or with given key is found
 *
 * Params
 * ======
 * head, Node *    :   node at which to start searching
 * key, key_t      :   key to look for
 *
 * Returns
 * =======
 * Node * with the desired key or its deepest existing anchestor
 *  */
Node *search( Node *head, key_t key )
{
    if ( key == head->key || !head->c )
        return head;
    else
        return search( head->c[bap(key, head->lvl+1)], key );
}


/* scr - search children recursivly
 *
 * Params
 * ======
 * head, Node*         :   node at which to start searching
 * suffixes, key_t *   :   relevant search directions, terminated by 0xffff
 * res, DArray_Node *  :   Array in which to write result
 *  */
static void scr( Node *head, const key_t *suffixes, DArray_Node *res )
{
    if ( !head->c )
        DArray_Node_insert(res, head);

    while ( *suffixes != 0xffff ) {
        scr( head->c[*suffixes], suffixes, res );
        ++suffixes;
    }
}


/* search_children
 * search for children of head, which are facing ref (i.e. its possible neighbours)
 * calculates search directions and calls `scr`
 *
 * Params
 * ======
 * head, Node *        :   node at which to start searching
 * ref, Node *         :   node into whichs direction to go
 * res, DArray_Node *  :   Array in which to write result
 *
 * Returns
 * =======
 * pointer to res, where the result was written
 *  */
DArray_Node *search_children( Node *head, const Node *ref, DArray_Node *res )
{
    key_t suffixes[3];
    suffixes[1] = suffixes[2] = 0xffff;

    if ( head->key == left(ref->key) ) {
        suffixes[0] = 0x1;
        suffixes[1] = 0x3;
    } else if ( head->key == right(ref->key) ) {
        suffixes[0] = 0x0;
        suffixes[1] = 0x2;
    } else if ( head->key == top(ref->key) ) {
        suffixes[0] = 0x2;
        suffixes[1] = 0x3;
    } else if ( head->key == bot(ref->key) ) {
        suffixes[0] = 0x0;
        suffixes[1] = 0x1;
    } else if ( head->key == top(left(ref->key)) ) {
        suffixes[0] = 0x3;
    } else if ( head->key == top(right(ref->key)) ) {
        suffixes[0] = 0x2;
    } else if ( head->key == bot(left(ref->key)) ) {
        suffixes[0] = 0x1;
    } else if ( head->key == bot(right(ref->key)) ) {
        suffixes[0] = 0x0;
    }

    scr( head, suffixes, res );

    return res;
}

