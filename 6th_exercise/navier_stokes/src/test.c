#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
void rfftshift(double *arr, size_t n0, size_t n1, char sign)
{
    size_t i, j, k;
    double dsign;
    dsign = sign > 0 ? 1. : -1.;

    /* for (i = 0; i < n0; ++i) */
    /*     for (j = 0; j < n1; ++j) */
    /*         arr[j+i*n1] *= dsign * ((i + j) % 2 == 0 ? -1 : 1); */
    /* for (i = 0; i < n0*n1; ++i)
     *     arr[i] *= i % 2 == 0 ? -1 : 1; */
    for (i = 0, k = 0; i < n1; ++i)
        for (j = 0; j < n0; ++j, ++k) {
            arr[j+i*n0] *= (i+j) % 2 == 0 ? -1 : 1;
            arr[k] *= k % 2 == 0 ? -1 : 1;
        }
}

int main()
{
    size_t i, j, Nx, Ny, Nkx, Nky, Ntot, Nktot;
    double xmin, xmax, ymin, ymax;
    double *x, *y, *z;
    fftw_complex *z_hat;

    Nx = 16, Ny = 16;
    Nkx = Nx / 2 + 1, Nky = Ny;
    Ntot = Nx * Ny;
    Nktot = Nkx * Nky;
    xmin = -M_PI, xmax = M_PI, ymin = -M_PI, ymax = M_PI;

    x = malloc(sizeof(*x) * Nx);
    y = malloc(sizeof(*y) * Ny);
    z = fftw_alloc_real(Ntot);
    z_hat = fftw_alloc_complex(Nktot);

    fftw_plan rfft2 = fftw_plan_dft_r2c_2d(Nx, Ny, z, z_hat, FFTW_MEASURE);

    rlinspace(xmin, xmax, Nx, x);
    rlinspace(ymin, ymax, Ny, y);
    for (i = 0; i < Ny; ++i)
        for (j = 0; j < Nx; ++j)
            z[j+i*Nx] = -inital_func(x[j], y[i]);

    /* print_real_array(z, Nx, Ny); */

    rfftshift(z, Nx, Ny, 1);
    fftw_execute(rfft2);
    /* rfftshift(z, Nx, Ny, 1); */

    /* print_real_array(z, Nx, Ny); */

    print_complex_array(z_hat, Nkx, Nky);

    /* Params p = {Nx, Ny, .05, .0}; */
    /* Workspace *pde = init(p, z); */

    /* printf("t = 0\n"); */
    /* print_real_array(pde->o, pde->Nx, pde->Ny); */
    /* time_step(1, pde); */
    /* printf("t = 0.05\n"); */
    /* print_real_array(pde->o, pde->Nx, pde->Ny); */
    /* printf("o_hat = \n"); */
    /* print_complex_array(pde->ohat, pde->Nky, pde->Nkx); */

    /* cleanup(pde); */

    free(x);
    free(y);
    fftw_free(z);
    fftw_free(z_hat);
    fftw_destroy_plan(rfft2);
    fftw_cleanup();

    return 0;
}
