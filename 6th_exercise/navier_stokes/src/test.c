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
            printf("%s%s%e%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x][REAL],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]\n" : ","));
    printf("imag = \n");
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            printf("%s%s%e%s%s",
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
    return exp(-4*(SQUARE(x-1) + SQUARE(y))) - exp(-4*(SQUARE(x+1) + SQUARE(y)));
}

static inline
void rfftshift(double *arr, size_t n0, size_t n1)
{
    size_t i, j;

    for (i = 0; i < n1; ++i)
        for (j = 0; j < n0; ++j)
            arr[j+i*n0] *= i % 2 == 0 ? 1. : -1.;
}


/* [> to measure execution time in ns: <]
 * struct timespec tp0, tp_a, tp_b;
 * tp0.tv_sec = 0, tp0.tv_nsec = 0;
 * clock_settime(0, &tp0);
 * clock_gettime(0, &tp_a);
 * [> do something... <]
 * clock_gettime(0, &tp_b);
 * printf("time:\t%ld ns\n\n", tp_b.tv_nsec-tp_a.tv_nsec); */
int main()
{
    size_t i, j, Nx, Ny, Nkx, Nky, Ntot, Nktot;
    double xmin, xmax, ymin, ymax;
    double *x, *y, *z;
    fftw_complex *z_hat, *res;

    Nx = 16, Ny = 16;
    Nkx = Nx / 2 + 1, Nky = Ny;
    Ntot = Nx * Ny;
    Nktot = Nkx * Nky;

    xmin = -M_PI, xmax = M_PI;
    ymin = -M_PI, ymax = M_PI;

    x = malloc(sizeof(*x) * Nx);
    y = malloc(sizeof(*y) * Ny);
    z = fftw_alloc_real(Ntot);
    z_hat = fftw_alloc_complex(Nktot);

    fftw_plan rfft = fftw_plan_dft_r2c_2d(Nx, Ny, z, z_hat, FFTW_MEASURE);


    x = linspace(xmin, xmax, Nx, x);
    y = linspace(ymin, ymax, Ny, y);

    for (i = 0; i < Ny; ++i)
        for (j = 0; j < Nx; ++j)
            z[j+i*Nx] = -inital_func(x[j], y[i]);


    rfftshift(z, Nx, Ny);
    fftw_execute(rfft);
    rfftshift(z, Nx, Ny);  

    /* print_real_array(z, Nx, Ny); */
    /* print_complex_array(z_hat, Nkx, Nky); */

    Params p = {Nx, Ny, .05, .0};
    Workspace *pde = init(p, z);

    res = rhs(z_hat, pde);
    
    /* print_complex_array(res, Nkx, Nky); */

    /* printf("t = 0\n"); */
    /* print_real_array(pde->o, pde->Nx, pde->Ny); */
    /* time_step(1, pde); */
    /* printf("t = 0.05\n"); */
    /* print_real_array(pde->o, pde->Nx, pde->Ny); */
    /* printf("o_hat = \n"); */
    /* print_complex_array(pde->ohat, pde->Nky, pde->Nkx); */

    free(x);
    free(y);
    /* free(z); */
    fftw_free(z);
    /* fftw_free(z_hat); */
    fftw_destroy_plan(rfft);
    /* fftw_cleanup(); */

    cleanup(pde);

    return 0;
}
