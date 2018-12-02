#pragma once

#include "fftw3.h"

#define SQUARE(x) (x)*(x)
#define REAL 0
#define IMAG 1


typedef struct Params {
    const size_t Nx, Ny;
    const double dt, nu;
} Params;

typedef struct PDE {
    size_t Nx, Ny,          /* sizes of real and complex arrays */
           Nkx, Nky,
           Ntot, ktot;
    double dt;              /* time step */
    double nu;              /* viscosity parameter */
    double *o,              /* omega */
           *kx, *ky,        /* wave numbers, k-space */
           *ksq,            /* k spared, i.e. |kx|^2 + |ky|^2 */
           *utmp,           /* helper arrays for `double *rhs` */
           *prop_full,      /* helper array for propagator */
           *prop_pos_half,
           *prop_neg_half;
    fftw_complex *ohat,     /* fourier transform of omega */
                 *_ohat,    /* backup of ohat, since reverse dft doesn't preserve input */
                 *uhat,     /* fourier transform of utmp */
                 *res;      /* helper array for `double *rhs` */
    unsigned char *mask;    /* mask array for anti-aliasing */

    fftw_plan ohat_to_o,
              uhat_to_u,
              o_to_ohat,
              u_to_uhat;

} PDE;

PDE *init(const Params, const double *);
void cleanup(PDE *);
double *time_step(unsigned int, PDE *);
double *linspace(double, double, size_t, double *);
