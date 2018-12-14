#include <stdlib.h>
#include <stdint.h>

/* for 256x256 resolution */
typedef uint16_t key_t;
typedef uint8_t level_t;
const uint8_t maxlevel = 8;

typedef struct node node;
struct node {
    /* TODO: add the node's data */
    key_t address;
    level_t level;
    node **c;       /* children */
    node **n;       /* neighbours */
    uint8_t allocated;
};

typedef struct {
    size_t size;        /* number of found nodes */
    struct node *n;     /* array of found nodes */
} finding;

/* void insert(node *head, key_t key); */
/* void delete(node *head); */
/* finding *search(node *head, key_t key);  */

void insert_all(node *head, key_t *keys, size_t size)
{
    key_t *current, *end;
    current = keys;
    end = keys + size;

    while (current != end) {

    }
}

void insert(node *head, key_t *key)
{
    level_t i, j, shift;
    key_t sb;   /* significant bits */
    sb = *key >> (2*(maxlevel - head->level - 1)) & 0b11;

    if (head->c[sb] == NULL) {
        /* TODO: make node */
        i = head->level+1;
        shift = 2 * (maxlevel - i);
        while ( (key[0] >> shift & 0b11) == (key[1] >> shift & 0b11)
                && i != maxlevel )
            ++i, shift-=2;
        i -= head->level;

        head->c[sb] = malloc(sizeof(node) * i);
        head->c[sb][0].allocated = 1;
        for (j = 0; j < i; ++j)
            head->c[j]->c[/*TODO*/0] = head->c[j+1];
    }
    else {
        insert(head->c[sb], key);
    }
}
