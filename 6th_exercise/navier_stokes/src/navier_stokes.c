#include <stdlib.h>
#include <math.h>
#include "../include/navier_stokes.h"


fftw_complex *rhs(fftw_complex *, Workspace *);
void scheme(Workspace *);

static inline void normalize(double *, size_t, size_t);
static inline void rfftshift(double *, size_t, size_t, char);
static inline void cfftshift(fftw_complex *, size_t, size_t, char);


/* set up workspace
 * allocate memory
 * init fftw plans
 *
 * call `cleanup` on the pointer returned by this function to free all
 *     allocated memory properly */
Workspace *init(Params params, double *iv)
{
    size_t i, j;
    fftw_complex kxmin, kxmax, kymin, kymax;
    Workspace *ws;

    /* allocate Workspace struct */
    ws = malloc(sizeof(*ws));

    /* init parameter */
    ws->Nx   = params.Nx;
    ws->Ny   = params.Ny;
    ws->Nkx  = params.Nx / 2 + 1;
    ws->Nky  = params.Ny;
    ws->Ntot = ws->Nx * ws->Ny;
    ws->ktot = ws->Nkx * ws->Nky;
    ws->dt   = params.dt;
    ws->nu   = params.nu;

    /* k-space: allocate and init as linspace */
    ws->kx = fftw_alloc_complex(ws->Nkx);
    ws->ky = fftw_alloc_complex(ws->Nky);
	kxmin[REAL] = 0.; kxmin[IMAG] = 0.;
	kxmax[REAL] = 0.; kxmax[IMAG] = (double) (ws->Nx/2);
	kymin[REAL] = 0.; kymin[IMAG] = -(double) (ws->Ny/2);
	kymax[REAL] = 0.; kymax[IMAG] = (double) (ws->Ny/2-1);
    ws->kx = clinspace(kxmin, kxmax, ws->Nkx, ws->kx);
    ws->ky = clinspace(kymin, kymax, ws->Nky, ws->ky);

    /* compute k^2 */
    ws->ksq = fftw_alloc_real(ws->ktot);
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j)
            /* ksq[x][y] = ksq[x + y*Nx] */
            ws->ksq[j+i*ws->Nkx] = CSQUARE(ws->kx[j]) + CSQUARE(ws->ky[i]);

    /* allocate omega and u and there transforms */
    ws->o    = fftw_alloc_real(ws->Ntot);
    ws->ohat = fftw_alloc_complex(ws->ktot);
    ws->utmp = fftw_alloc_real(ws->Ntot);
    ws->uhat = fftw_alloc_complex(ws->ktot);

    /* set up fftw plan */
    ws->o_to_ohat = fftw_plan_dft_r2c_2d(ws->Nx, ws->Ny,
                                         ws->o, ws->ohat,
                                         FFTW_MEASURE);
    ws->ohat_to_o = fftw_plan_dft_c2r_2d(ws->Nx, ws->Ny,
                                         ws->ohat, ws->o,
                                         FFTW_MEASURE);
    ws->u_to_uhat = fftw_plan_dft_r2c_2d(ws->Nx, ws->Ny,
                                         ws->utmp, ws->uhat,
                                         FFTW_MEASURE);
    ws->uhat_to_u = fftw_plan_dft_c2r_2d(ws->Nx, ws->Ny,
                                         ws->uhat, ws->utmp,
                                         FFTW_MEASURE);

    /* set inital value and compute FT */
    for (i = 0; i < ws->Ntot; ++i)
        ws->o[i] = iv[i];
    rfftshift(ws->o, ws->Nx, ws->Ny, 1);
    fftw_execute(ws->o_to_ohat);
    rfftshift(ws->o, ws->Nx, ws->Ny, -1);

    /* allocate and init mask for anti-aliasing */
    ws->mask = malloc(sizeof(unsigned char) * ws->ktot);
    double threshold = ws->Ntot / 9.;
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
    ws->res  = fftw_alloc_complex(ws->ktot);

    /* done */
    return ws;
}


/* free all allocated memory in a workspace struct
 * the pointer is afterwards NULL */
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
    fftw_destroy_plan(ws->o_to_ohat);
    fftw_destroy_plan(ws->ohat_to_o);
    fftw_destroy_plan(ws->u_to_uhat);
    fftw_destroy_plan(ws->uhat_to_u);
    fftw_cleanup();
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

    cfftshift(ws->ohat, ws->Nky, ws->Nkx, -1);
    fftw_execute(ws->ohat_to_o);
    normalize(ws->o, ws->Ntot, ws->ktot);
    cfftshift(ws->ohat, ws->Nky, ws->Nkx, 1);

    return ws->o;
}


/* right hand side of equation */
fftw_complex *rhs(fftw_complex *ohat, Workspace *ws)
{
    size_t i, j, idx;

    /* anti aliasing */
    for (i = 0; i < ws->ktot; ++i)
        if (! ws->mask[i]) {
            ohat[i][REAL] = 0.;
            ohat[i][IMAG] = 0.;
        }

    /* iFT of ohat, yielding o */
    cfftshift(ws->ohat, ws->Nky, ws->Nkx, -1);
    fftw_execute(ws->ohat_to_o);
    normalize(ws->o, ws->Ntot, ws->ktot);
    cfftshift(ws->ohat, ws->Nky, ws->Nkx, 1);


    /********************/
    /* x component of u */
    /********************/
    /* uhat_x = I * ky * ohat / k^2 */
    for (i = 0; i < ws->Nky; ++i)
        /* for (j = 0; j < ws->Nkx; ++j)
         *     ws->uhat[j+i*ws->Nkx] = \
         *         I * ws->ky[i] * ws->ohat[j+i*ws->Nkx] / ws->ksq[j+i*ws->Nkx]; */
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->uhat[idx][REAL] = \
                ws->ky[i][IMAG] * ws->ohat[idx][IMAG] / ws->ksq[idx];
            ws->uhat[idx][IMAG] = \
                ws->ky[i][REAL] * ws->ohat[idx][REAL] / ws->ksq[idx];
        }

    /* compute iFT of uhat, yielding utmp */
    cfftshift(ws->uhat, ws->Nky, ws->Nkx, -1);
    fftw_execute(ws->uhat_to_u);
    normalize(ws->utmp, ws->Ntot, ws->ktot);

    /* u_x * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    rfftshift(ws->utmp, ws->Nx, ws->Ny, 1);
    fftw_execute(ws->u_to_uhat);

    /* write into result */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->res[idx][REAL] = ws->ohat[idx][REAL] \
                - ws->kx[j][IMAG] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] = ws->ohat[idx][IMAG] \
                - ws->kx[j][REAL] * ws->uhat[idx][REAL] * ws->dt;
        }


    /********************/
    /* y component of u */
    /********************/
    /* uhat_y = I * kx * ohat / k^2 */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->uhat[idx][REAL] = \
                ws->kx[j][IMAG] * ws->ohat[idx][IMAG] / ws->ksq[idx];
            ws->uhat[idx][IMAG] = \
                ws->kx[j][REAL] * ws->ohat[idx][REAL] / ws->ksq[idx];
        }

    /* compute iFT of uhat, yielding utmp */
    cfftshift(ws->uhat, ws->Nky, ws->Nkx, -1);
    fftw_execute(ws->uhat_to_u);
    normalize(ws->utmp, ws->Ntot, ws->ktot);

    /* u_y * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    rfftshift(ws->utmp, ws->Nx, ws->Ny, 1);
    fftw_execute(ws->u_to_uhat);

    /* write into result */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->res[idx][REAL] -= \
                ws->ky[i][IMAG] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] -= \
                ws->ky[i][REAL] * ws->uhat[idx][REAL] * ws->dt;
        }


    return ws->res;
}


/* shu-osher scheme: 3-step runge-kutta */
void scheme(Workspace *ws)
{
    size_t i;
    fftw_complex *otmp;

    /* step one */
    otmp = rhs(ws->ohat, ws);
    for (i = 0; i < ws->ktot; ++i) {
        otmp[i][REAL] *= ws->prop_full[i];
        otmp[i][IMAG] *= ws->prop_full[i];
    }

    /* step two */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i) {
        otmp[i][REAL] *= .25 * ws->prop_neg_half[i];
        otmp[i][IMAG] *= .25 * ws->prop_neg_half[i];
    }

    /* step three */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i) {
        ws->ohat[i][REAL] = 1./3. * (2. * otmp[i][REAL] * ws->prop_pos_half[i] \
                + ws->ohat[i][REAL] * ws->prop_pos_half[i]);
        ws->ohat[i][IMAG] = 1./3. * (2. * otmp[i][IMAG] * ws->prop_pos_half[i] \
                + ws->ohat[i][IMAG] * ws->prop_pos_half[i]);
    }
}


/* a = a / n, n: size of transform */
static inline
void normalize(double *arr, size_t size, size_t n)
{
    double *end;
    end = arr+size;

    while (arr != end)
        *arr++ /= n;
}


static inline
void rfftshift(double *arr, size_t n0, size_t n1, char sign)
{
    size_t i, j;
    double dsign;
    dsign = sign > 0 ? 1. : -1.;

    for (i = 0; i < n0; ++i)
        for (j = 0; j < n1; ++j)
            arr[j+i*n1] *= dsign * ((i + j) % 2 == 0 ? 1 : -1);
}

static inline
void cfftshift(fftw_complex *arr, size_t n0, size_t n1, char sign)
{
    size_t i, j;
    double dsign;
    dsign = sign > 0 ? 1. : -1.;

    for (i = 0; i < n0; ++i) {
        for (j = 0; j < n1; ++j) {
            arr[j+i*n1][REAL] *= dsign * ((i + j) % 2 == 0 ? 1 : -1);
            arr[j+i*n1][IMAG] *= dsign * ((i + j) % 2 == 0 ? 1 : -1);
        }
    }
}


inline fftw_complex *clinspace(fftw_complex start, fftw_complex end,
                        size_t npoints, fftw_complex *dst)
{
    size_t i;
    fftw_complex z, dz;

    z[REAL] = start[REAL];
    z[IMAG] = start[IMAG];
    dz[REAL] = (end[REAL] - start[REAL]) / ((double) npoints);
    dz[IMAG] = (end[IMAG] - start[IMAG]) / ((double) npoints);

    for (i = 0; i < npoints; ++i) {
        dst[i][REAL] = z[REAL];
        dst[i][IMAG] = z[IMAG];
        z[REAL] += dz[REAL];
        z[IMAG] += dz[IMAG];
    }

    return dst;
}


inline double *rlinspace(double start, double end, size_t np, double *dst)
{
    double x, dx, *startptr, *endptr;

    x = start;
    dx = (end - start) / ((double) np);
    startptr = dst;
    endptr = dst + np;

    while (dst != endptr)
        *dst++ = x, x += dx;

    return startptr;
}
