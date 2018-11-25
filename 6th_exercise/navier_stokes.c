#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


#define CSQUARE(x) creal(x)*creal(x) + cimag(x)*cimag(x)


typedef struct Params {
    size_t Nx, Ny;
    double dt, nu;
} Params;

typedef struct Workspace {
    size_t Nx, Ny, Nkx, Nky, ktot;  /* sizes of real and complex arrays */
    double dt;              /* time step */
    double nu;              /* viscosity parameter */
    double *o,              /* omega */
           *ksq,            /* k spared, i.e. |kx|^2 + |ky|^2 */
           *utmp,           /* helper arrays for `double *rhs` */
           *prop_full,      /* helper array for propagator */
           *prop_pos_half,
           *prop_neg_half;
    double _Complex *kx, *ky,   /* wave numbers, k-space */
                    *ohat,      /* fourier transform of omega */
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
double _Complex *rhs(double _Complex *, Workspace *);
void scheme(Workspace *);

void fftshift(double _Complex *);
void ifftshift(double _Complex *);
double _Complex *clinspace(double _Complex, double _Complex,
                           size_t, double _Complex *);


/* set up workspace
 * allocate memory
 * init fftw plans
 *
 * call `cleanup` on the pointer returned by this function to free all
 *     allocated memory properly */
Workspace *init(Params params, double *iv)
{
    size_t i, j;
    Workspace *ws;

    /* allocate Workspace struct */
    ws = malloc(sizeof(*ws));

    /* init parameter */
    ws->Nx   = params.Nx;
    ws->Ny   = params.Ny;
    ws->Nkx  = params.Nx / 2 + 1;
    ws->Nky  = params.Ny;
    ws->ktot = ws->Nkx * ws->Nky;
    ws->dt   = params.dt;
    ws->nu   = params.nu;

    /* k-space: allocate and init as linspace */
    ws->kx = fftw_malloc_complex(ws->Nkx);
    ws->ky = fftw_malloc_complex(ws->Nky);
    ws->kx = clinspace(0., (double _Complex) (ws->Nx/2),
                       ws->Nkx, ws->kx);
    ws->ky = clinspace((double _Complex) (-ws->Nx/2),
                       (double _Complex) (ws->Nx/2-1),
                       ws->Nky, ws->ky);

    /* compute k^2 */
    ws->ksq = fftw_malloc_real(ws->ktot);
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j)
            /* ksq[x][y] = ksq[x + y*Nx] */
            ws->ksq[j+i*ws->Nkx] = CSQUARE(ws->kx[j]) + CSQUARE(ws->ky[i]);

    /* allocate omega and u and there transforms */
    ws->o    = fftw_malloc_real(ws->Nx * ws->Ny);
    ws->ohat = fftw_malloc_complex(ws->ktot);
    ws->utmp = fftw_malloc_real(ws->Nx * ws->Ny);
    ws->uhat = fftw_malloc_complex(ws->ktot);

    /* set up fftw plan */
    ws->o_to_ohat = fftw_plan_dft_r2c_2d(ws->Nx, ws->Ny,
                                         ws->o, ws->ohat,
                                         FFTW_MEASURE);
    ws->ohat_to_o = fftw_plan_dft_c2r_2d(ws->Nkx, ws->Nky,
                                         ws->ohat, ws->o,
                                         FFTW_MEASURE);
    ws->u_to_uhat = fftw_plan_dft_r2c_2d(ws->Nx, ws->Ny,
                                         ws->utmp, ws->uhat,
                                         FFTW_MEASURE);
    ws->uhat_to_u = fftw_plan_dft_c2r_2d(ws->Nkx, ws->Nky,
                                         ws->uhat, ws->utmp,
                                         FFTW_MEASURE);

    /* set inital value and compute FT */
    for (i = 0; i < ws->Nx * ws->Ny; ++i)
        ws->o[i] = iv[i];
    fftw_execute(ws->o_to_ohat);
    fftshift(ws->ohat);

    /* allocate and init mask for anti-aliasing */
    ws->mask = malloc(sizeof(unsigned char) * ws->ktot);
    double threshold = ws->Nx * ws->Ny / 9.;
    for (i = 0; i < ws->ktot; ++i)
        ws->mask[i] = ws->ksq[i] < threshold ? 1 : 0;

    /* allocate and init propagator lookup tables */
    ws->prop_full        = malloc(sizeof(double) * ws->ktot);
    ws->prop_pos_half    = malloc(sizeof(double) * ws->ktot);
    ws->prop_neg_half    = malloc(sizeof(double) * ws->ktot);
    for (i = 0; i < ws->ktot; ++i) {
        ws->prop_full[i]     = exp(-ws->nu * ws->ksq[i] * ws->dt);
        ws->prop_pos_half[i] = exp(-.5 *ws->nu * ws->ksq[i] * ws->dt);
        ws->prop_neg_half[i] = exp(.5 * ws->nu * ws->ksq[i] * ws->dt);
    }

    /* allocate helper array for `double *rhs` (is initilized there) */
    ws->res  = fftw_malloc_complex(ws->ktot);

    /* done */
    return ws;
}


/* void  */
void cleanup(Workspace *ws)
{
    fftw_free(ws->kx);
    fftw_free(ws->ky);
    fftw_free(ws->ksq);
    fftw_free(ws->o);
    fftw_free(ws->ohat);
    fftw_free(ws->utmp);
    fftw_free(ws->uhat);
    fftw_free(ws->res);
    free(ws->mask);
    free(ws->prop_full);
    free(ws->prop_pos_half);
    free(ws->prop_neg_half);
    free(ws);
}


/* calculate the next frame */
double *time_step(unsigned short steps, Workspace *ws)
{
    unsigned short i;

    for (i = 0; i < steps; ++i)
        scheme(ws);

    ifftshift(ws->ohat);
    fftw_execute(ws->ohat_to_o);

    return ws->o;
}


/* right hand side of equation */
double _Complex *rhs(double _Complex *ohat, Workspace *ws)
{
    size_t i;

    /* anti aliasing */
    for (i = 0; i < ws->ktot; ++i)
        ohat[i] = ws->mask[i] ? ohat[i] : 0.;

    /* iFT of ohat, yielding o */
    ifftshift(ws->ohat);
    fftw_execute(ws->ohat_to_o);


    /********************/
    /* x component of u */
    /********************/
    /* uhat_x = I * ky * ohat / k^2 */
    for (i = 0; i < ws->ktot; ++i)
        ws->uhat[i] = I * ws->ky[i] * ws->ohat[i] / ws->ksq[i];

    /* compute iFT of uhat, yielding utmp */
    ifftshift(ws->uhat);
    fftw_execute(ws->uhat_to_u);

    /* u_x * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    fftw_execute(ws->u_to_uhat);
    fftshift(ws->uhat);

    /* write into result */
    for (i = 0; i < ws->ktot; ++i)
        ws->res[i] = ws->ohat - I * ws->kx[i] * ws->uhat[i] * ws->dt;


    /********************/
    /* y component of u */
    /********************/
    /* uhat_y = I * kx * ohat / k^2 */
    for (i = 0; i < ws->ktot; ++i)
        ws->uhat[i] = I * ws->kx[i] * ws->ohat[i] / ws->ksq[i];

    /* compute iFT of uhat, yielding utmp */
    ifftshift(ws->uhat);
    fftw_execute(ws->uhat_to_u);

    /* u_y * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    fftw_execute(ws->u_to_uhat);
    fftshift(ws->uhat);

    /* write into result */
    for (i = 0; i < ws->ktot; ++i)
        ws->res[i] -= I * ws->ky[i] * ws->uhat[i] * ws->dt;


    return ws->res;
}


/* shu-osher scheme: 3-step runge-kutta */
void scheme(Workspace *ws)
{
    size_t i;
    double _Complex *otmp;

    /* step one */
    otmp = rhs(ws->ohat, ws);
    for (i = 0; i < ws->ktot; ++i)
        otmp[i] *= ws->prop_full[i];

    /* step two */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i)
        otmp[i] *= .25 * ws->prop_neg_half[i];

    /* step three */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i)
        ws->ohat[i] = 1./3. * (2. * otmp[i] * ws->prop_pos_half[i] \
                               + ws->ohat[i] * ws->prop_pos_half[i]);
}


double _Complex *clinspace(
        double _Complex start, double _Complex end,
        size_t np, double _Complex *dst)
{
    double _Complex z, dz, *startptr, *endptr;

    z = start;
    dz = (end - start) / ((double _Complex) np);
    startptr = dst;
    endptr = dst + np;

    while (dst != endptr)
        *dst++ = z, z += dz;

    return startptr;
}
