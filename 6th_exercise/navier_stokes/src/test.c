#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "../include/navier_stokes.h"


void print_real_array(double *arr, size_t x, size_t y)
{
    size_t i, j;

    for (i = 0; i < y; ++i)
        for(j = 0; j < x; ++j)
            printf("%s%s%f%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]\n" : ","));
    puts("\n");
}


void print_complex_array(fftw_complex *arr, size_t x, size_t y)
{
    size_t i, j;

    printf("real = \n");
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            printf("%s%s%f%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x][REAL],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]\n" : ","));
    printf("imag = \n");
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            printf("%s%s%f%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x][IMAG],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]\n" : ","));
    puts("\n");
    puts("\n");
}


double inital_func(double x, double y)
{
    return exp(-4*(RSQUARE(x-1) + RSQUARE(y))) - exp(-4*(RSQUARE(x+1) + RSQUARE(y)));
}

static inline
void rfftshift(double *arr, size_t n0, size_t n1)
{
    size_t i, j, k;

    for (i = 0, k = 0; i < n1; ++i)
        for (j = 0; j < n0; ++j, ++k) {
            /* center zero-mode in the middle and flip along x-axis */
            arr[j+i*n0] *= (i+j) % 2 == 0 ? -1 : 1;
            /* revert x-axis flip */
            arr[k] *= k % 2 == 0 ? -1 : 1;
        }
}

void normalize(double *arr, size_t size, size_t n)
{
    double *end;
    end = arr+size;

    while (arr != end)
        *arr++ /= n;
}


int main()
{
    size_t i, j, Nx, Ny, Nkx, Nky, Ntot, Nktot;
    double xmin, xmax, ymin, ymax;
    double *x, *y, *z, *z_bu;
    fftw_complex *z_hat, *z_hat_bu;

    Nx = 16, Ny = 16;
    Nkx = Nx / 2 + 1, Nky = Ny;
    Ntot = Nx * Ny;
    Nktot = Nkx * Nky;
    xmin = -M_PI, xmax = M_PI, ymin = -M_PI, ymax = M_PI;

    x = malloc(sizeof(*x) * Nx);
    y = malloc(sizeof(*y) * Ny);
    z = malloc(sizeof(*z) * Ntot);

    rlinspace(xmin, xmax, Nx, x);
    rlinspace(ymin, ymax, Ny, y);
    for (i = 0; i < Ny; ++i)
        for (j = 0; j < Nx; ++j)
            z[j+i*Nx] = -inital_func(x[j], y[i]);

    /* [> to measure execution time in ns: <]
     * struct timespec tp0, tp_a, tp_b;
     * tp0.tv_sec = 0, tp0.tv_nsec = 0;
     * clock_settime(0, &tp0);
     * clock_gettime(0, &tp_a);
     * [> do something... <]
     * clock_gettime(0, &tp_b);
     * printf("time:\t%ld ns\n\n", tp_b.tv_nsec-tp_a.tv_nsec); */


    Params p = {Nx, Ny, .05, .0};
    Workspace *pde = init(p, z);

    printf("t = 0\n");
    print_real_array(pde->o, pde->Nx, pde->Ny);
    time_step(1, pde);
    printf("t = 0.05\n");
    print_real_array(pde->o, pde->Nx, pde->Ny);
    /* printf("o_hat = \n"); */
    /* print_complex_array(pde->ohat, pde->Nky, pde->Nkx); */

    cleanup(pde);

    free(x);
    free(y);
    free(z);

    return 0;
}
