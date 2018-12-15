/* This header provides wrapper around malloc and realloc including
 * convenience makros.
 * It is checked, whether malloc/realloc returns NULL, meaning failure to
 * perform the requested operation; the file and line in which the incident
 * occured are then printed and the program is aborted
 *
 * Taken from "Programming Projects in C" by Rouben Rostamian
 * ISBN 978-1-611973-49-5
 *  */
#pragma once

#include <stdlib.h>
#include <stdio.h>


void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
    void *res;
    if ( (res = malloc(nbytes)) == NULL ) {
        fprintf(stderr, "%s:line %d: malloc() of %zu bytes failed ...\n",
                file, line, nbytes);
        exit(EXIT_FAILURE);
    } else {
        return res;
    }
}


void *realloc_or_exit(void *ptr, size_t nbytes, const char *file, int line)
{
    void *tp;
    if ( (tp = realloc(ptr, nbytes)) == NULL ) {
        fprintf(stderr, "%s:line %d: realloc() of %zu bytes failed ...\n",
                file, line, nbytes);
        exit(EXIT_FAILURE);
    } else {
        return tp;
    }
}


// convenience makros
#define xmalloc(nbytes) malloc_or_exit((nbytes), __FILE__, __LINE__)
#define xrealloc(ptr, nbytes)   realloc_or_exit((ptr), (nbytes),    \
                                                __FILE__, __LINE__)

