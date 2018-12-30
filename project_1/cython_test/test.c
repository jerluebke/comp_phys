#include <stdio.h>
#include "cvisualise.h"
#include "data.h"

#if 0
compile with:
    gcc -Wall -Wextra -std=c11 test.c cvisualise.c ../src/*.c -I. -I../include/
                                                           -lm -g -o test.out
#endif

int main()
{
    lvl_t nl, res[16];
    QuadtreeEnv *qtenv;
    
    qtenv = qtenv_setup(data, 64);
    for (;;) {
        if ( qtenv_is_last(qtenv) )
            break;
        nl = qtenv_insert(qtenv, &res[0]);
        printf("new levels: %u\n", nl);
    }

    qtenv_free(qtenv);

    printf("done\n");
    return 0;
}
