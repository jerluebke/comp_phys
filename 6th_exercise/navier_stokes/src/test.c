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


int main()
{
    size_t i, j, Nx, Ny, Ntot;
    double xmin, xmax, ymin, ymax;
    double *x, *y, *z;

    Nx = 16, Ny = 16;
    Ntot = Nx * Ny;
    xmin = -M_PI, xmax = M_PI, ymin = -M_PI, ymax = M_PI;

    x = malloc(sizeof(*x) * Nx);
    y = malloc(sizeof(*y) * Ny);
    z = malloc(sizeof(*z) * Ntot);

    rlinspace(xmin, xmax, Nx, x);
    rlinspace(ymin, ymax, Ny, y);
    for (i = 0; i < Ny; ++i)
        for (j = 0; j < Nx; ++j)
            z[j+i*Nx] = inital_func(x[j], y[i]);

    Params p = {Nx, Ny, .05, .0};
    Workspace *pde = init(p, z);

    printf("t = 0\n");
    print_real_array(pde->o, pde->Nx, pde->Ny);
    time_step(1, pde);
    printf("t = 0.05\n");
    print_real_array(pde->o, pde->Nx, pde->Ny);
    printf("o_hat = \n");
    print_complex_array(pde->ohat, pde->Nky, pde->Nkx);

    cleanup(pde);

    free(x);
    free(y);
    free(z);

    return 0;
}
