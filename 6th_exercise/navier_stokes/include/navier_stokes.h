#pragma once

#include <fftw3.h>

#define RSQUARE(x) (x)*(x)
#define CSQUARE(x) (x)[0]*(x)[0] + (x)[1]*(x)[1]
#define REAL 0
#define IMAG 1


typedef struct Params {
    size_t Nx, Ny;
    double dt, nu;
} Params;

typedef struct Workspace {
    size_t Nx, Ny,          /* sizes of real and complex arrays */
           Nkx, Nky,
           Ntot, ktot;
    double dt;              /* time step */
    double nu;              /* viscosity parameter */
    double *o,              /* omega */
           *ksq,            /* k spared, i.e. |kx|^2 + |ky|^2 */
           *utmp,           /* helper arrays for `double *rhs` */
           *prop_full,      /* helper array for propagator */
           *prop_pos_half,
           *prop_neg_half;
    fftw_complex *kx, *ky,   /* wave numbers, k-space */
                 *ohat,      /* fourier transform of omega */
                 *_ohat,
                 *uhat,      /* fourier transform of utmp */
                 *res;       /* helper array for `double *rhs` */
    unsigned char *mask;    /* mask array for anti-aliasing */

    fftw_plan ohat_to_o,
              uhat_to_u,
              o_to_ohat,
              u_to_uhat;

} Workspace;

Workspace *init(Params, double *);
void cleanup(Workspace *);
double *time_step(unsigned short, Workspace *);

fftw_complex *clinspace(fftw_complex, fftw_complex,
                        size_t, fftw_complex *);
double *rlinspace(double, double, size_t, double *);
