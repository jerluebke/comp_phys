#include <stdlib.h>
#include <string.h>     /* memcpy */
#include <math.h>
#include "../include/navier_stokes.h"

/* TODO
 * check declarations in header for necessity
 * check function attributes (static, inline)
 * check parameter attributes (const) */


/* fftw_complex *rhs(fftw_complex *, Workspace *); */
void scheme(Workspace *);

void rfftshift(double *, size_t, size_t);
void rfft2(fftw_plan *, double *, Workspace *);
void irfft2(fftw_plan *, fftw_complex *, fftw_complex *, double *, Workspace *);


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
    ws->ksq[ws->Nkx*ws->Nky/2] = 1.;

    /* allocate omega and u and there transforms
     * include backup memory for inverse transforms */
    ws->o    = fftw_alloc_real(ws->Ntot);
    ws->ohat = fftw_alloc_complex(ws->ktot);
    ws->_ohat= fftw_alloc_complex(ws->ktot);
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
    rfft2(&ws->o_to_ohat, ws->o, ws);

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
    fftw_free(ws->_ohat);
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

    irfft2(&ws->ohat_to_o, ws->ohat, ws->_ohat, ws->o, ws);

    return ws->o;
}


/*********************************************************************/
/* TODO: remove later */
static
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


static
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


/*********************************************************************/


/* right hand side of equation */
fftw_complex *rhs(fftw_complex *ohat, Workspace *ws)
{
    size_t i, j, idx;
    fftw_complex *tmp;

    /* anti aliasing */
    for (i = 0; i < ws->ktot; ++i)
        if (! ws->mask[i])
            ohat[i][IMAG] = 0.;


    /* TODO: do switching outside of irfft */
    /* iFT of ohat, yielding o */
    irfft2(&ws->ohat_to_o, ws->ohat, ws->_ohat, ws->o, ws);
    tmp = ws->ohat;
    ws->ohat = ws->_ohat;
    ws->_ohat = tmp;


    /* TODO:
     * make linspace include max value */


    /********************/
    /* x component of u */
    /********************/
    /* uhat_x = I * ky * ohat / k^2 */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            /* ws->uhat[idx][REAL] = \ */
            /*     - ws->ky[i][IMAG] * ws->ohat[idx][IMAG] / ws->ksq[idx]; */
            ws->uhat[idx][IMAG] = \
                ws->ky[i][IMAG] * ws->ohat[idx][REAL] / ws->ksq[idx];
            /* FIXME */
            /* ws->uhat[idx][IMAG] = \
             *     ws->ky[i][IMAG] * ws->ohat[i][IMAG] / ws->ksq[idx]; */
        }

    /**************************/
    /* TODO: work in progress */
    /**************************/



    /* compute iFT of uhat, yielding utmp */
    irfft2(&ws->uhat_to_u, ws->uhat, ws->uhat, ws->utmp, ws);

    /* u_x * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    rfft2(&ws->u_to_uhat, ws->utmp, ws);

    /* write into result */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->res[idx][REAL] = ws->ohat[idx][REAL] \
                + ws->kx[j][IMAG] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] = ws->ohat[idx][IMAG] \
                - ws->kx[j][IMAG] * ws->uhat[idx][REAL] * ws->dt;
            /* FIXME: + or - ? */
            /* ws->res[idx][IMAG] = ws->ohat[idx][IMAG] \ */
            /*     + ws->kx[j][IMAG] * ws->uhat[i][IMAG] * ws->dt; */
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
                - ws->kx[j][IMAG] * ws->ohat[idx][REAL] / ws->ksq[idx];
            /* ws->uhat[idx][IMAG] = \ */
            /*     ws->kx[j][IMAG] * ws->ohat[idx][IMAG] / ws->ksq[idx]; */
        }

    /* compute iFT of uhat, yielding utmp */
    irfft2(&ws->uhat_to_u, ws->uhat, ws->uhat, ws->utmp, ws);

    /* u_y * o */
    for (i = 0; i < ws->ktot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* compute FT of u * o, yielding uhat */
    rfft2(&ws->u_to_uhat, ws->utmp, ws);

    /* write into result */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->res[idx][REAL] += \
                ws->ky[i][IMAG] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] -= \
                ws->ky[i][IMAG] * ws->uhat[idx][REAL] * ws->dt;
            /* ws->res[idx][IMAG] -= \ */
            /*     ws->ky[i][IMAG] * ws->uhat[idx][IMAG] * ws->dt; */
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
    for (i = 0; i < ws->ktot; ++i)
        otmp[i][IMAG] *= ws->prop_full[i];

    /* step two */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i)
        otmp[i][IMAG] = .25 * otmp[i][IMAG] * ws->prop_neg_half[i] \
                         + 3. * ws->ohat[i][IMAG] * ws->prop_pos_half[i];

    /* step three */
    otmp = rhs(otmp, ws);
    for (i = 0; i < ws->ktot; ++i)
        ws->ohat[i][IMAG] = 1./3. * (2. * otmp[i][IMAG] * ws->prop_pos_half[i] \
                            + ws->ohat[i][IMAG] * ws->prop_pos_half[i]);
}


/* perform real DFT
 *
 * Params
 * ======
 * plan    :   pointer to fftw_plan to execute
 * orig    :   array which will be transformed
 * ws      :   Workspace pointer holing auxillary values
 *  */
void rfft2(fftw_plan *plan, double *orig, Workspace *ws)
{
    /* prepare input for result with origin shifted to center */
    rfftshift(orig, ws->Nx, ws->Ny);
    /* do dft */
    fftw_execute(*plan);
    /* reverse changes */
    rfftshift(orig, ws->Nx, ws->Ny);
}


/* perform inverse DFT and normalize result
 *
 * Params
 * ======
 * plan    :   pointer to fftw_plan to execute
 * io      :   inout-struct containing pointer to memory for array to
 *             transform and to store its backup
 * res     :   double array, in which the result of the iDFT will be written
 * ws      :   Workspace pointer, which holds auxillary values
 *  */
void irfft2(fftw_plan *plan, fftw_complex *orig, fftw_complex *backup,
            double *res, Workspace *ws)
{
    double *end;
    fftw_complex *tmp;

    /* write backup of input, since c2r dft doesn't preserve input ... */
    if (orig != backup)
        memcpy((void *)backup, (void *)orig, sizeof(fftw_complex) * ws->ktot);

    /* do dft and center result */
    fftw_execute(*plan);
    rfftshift(res, ws->Nx, ws->Ny);

    /* normalize */
    end = res + ws->Ntot;
    while (res != end)
        *res++ /= ws->Ntot;

    /* set backup as original */
    tmp = orig;
    orig = backup;
    backup = tmp;
}


/* prepare input such that its DFT will have the origin frequence in the center
 *
 * Params
 * ======
 * arr     :   double array, input
 * nx, ny  :   int, dimension sizes
 *  */
void rfftshift(double *arr, size_t nx, size_t ny)
{
    size_t i, j;

    for (i = 0; i < ny; ++i)
        for (j = 0; j < nx; ++j)
            arr[j+i*nx] *= i % 2 == 0 ? 1. : -1.;
}


/* TODO: remove */
inline
fftw_complex *clinspace(fftw_complex start, fftw_complex end,
                        size_t npoints, fftw_complex *dst)
{
    size_t i;
    fftw_complex z, dz;

    z[REAL] = start[REAL];
    z[IMAG] = start[IMAG];
    dz[REAL] = (end[REAL] - start[REAL]) / ((double) npoints);
    dz[IMAG] = (end[IMAG] - start[IMAG]+1.) / ((double) npoints);

    for (i = 0; i < npoints; ++i) {
        dst[i][REAL] = z[REAL];
        dst[i][IMAG] = z[IMAG];
        z[REAL] += dz[REAL];
        z[IMAG] += dz[IMAG];
    }

    return dst;
}


/* TODO: include end */
inline
double *rlinspace(double start, double end, size_t np, double *dst)
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
