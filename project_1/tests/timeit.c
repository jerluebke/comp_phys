/* This file contains a naive function to measure execution time of a given
 * function pointer.
 * args are passed via a fargs struct.
 * the used data is hardcoded in "sample_data.h" (for reproducibility and
 * sanity-checks, i.e. obeys the 16-bit-key-boundaries).
 * the results in nano seconds are printed to stdout.
 *
 * TODO: integer overflow might occure, perhaps print also seconds
 *
 * REQUIRES POSIX
 *
 */
#if defined(unix) || defined(__unix__) || defined (__unix)
#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#if !_POSIX_TIMERS
#error This file requieres the POSIX timer implementation, which is not \
    supported on your machine!
#endif
#else
#error This file was written for a unix OS with POSIX support!
#endif

#include <time.h>
#include <limits.h>
#include <math.h>
#include "search.h"
#include "sample_data.h"


void timeit(int(*func)(const fargs_t *), const fargs_t *fargs,
            unsigned int iter, char *name)
{
    unsigned int i;
    double avg = 0, std = 0;
    long min = LONG_MAX, max = LONG_MIN;
    long *res = xmalloc(sizeof(long) * iter);
    struct timespec tp_a, tp_b, tp_c;

    clock_gettime(CLOCK_REALTIME, &tp_a);

    for ( i = 0; i < iter; ++i ) {
        clock_gettime(CLOCK_REALTIME, &tp_b);
        (*func)(fargs);
        clock_gettime(CLOCK_REALTIME, &tp_c);

        res[i] = tp_c.tv_nsec - tp_b.tv_nsec;
        avg += (double) res[i];
        if ( res[i] < min )
            min = res[i];
        if ( res[i] > max)
            max = res[i];
    }

    clock_gettime(CLOCK_REALTIME, &tp_c);

    avg /= iter;
    for ( i = 0; i < iter; ++i )
        std += (double)(SQUARE(res[i] - avg));
    std = sqrt(std / ((double) iter-1));

    printf("\nTIMEIT %s - %u runs\n"        \
           "average     :   %lf +- %lf ns\n"\
           "    max/min :   %lu / %lu\n"    \
           "total time  :   %lu ns\n\n",
           name, iter, avg, std, max, min, tp_c.tv_nsec - tp_a.tv_nsec);

    free(res);
}


int main()
{
    unsigned int iter = 100u;
    const fargs_t fargs_small = { .data=input_data_256, .size=size_256, .r_sq=16.0f };
    timeit(search_naive, &fargs_small, iter, "naive - 256");
    timeit(search_fast, &fargs_small, iter, "fast - 256");
    timeit(search_fastfast, &fargs_small, iter, "fastfast - 256");

    const fargs_t fargs = { .data=input_data_1265, .size=size_1265, .r_sq=16.0f };
    timeit(search_naive, &fargs, iter, "naive - 1265");
    timeit(search_fast, &fargs, iter, "fast - 1265");
    timeit(search_fastfast, &fargs, iter, "fastfast - 1265");
    return 0;
}

/* vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : */
