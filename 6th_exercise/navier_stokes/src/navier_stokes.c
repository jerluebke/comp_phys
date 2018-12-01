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
void irfft2(fftw_plan *, double *, Workspace *);


/* set up workspace
 * allocate memory
 * init fftw plans
 *
 * call `cleanup` on the pointer returned by this function to free all
 *     allocated memory properly */
Workspace *init(Params params, double *iv)
{
    size_t i, j;
    double kxmin, kxmax, kymin, kymax;
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
    ws->kx = fftw_alloc_real(ws->Nkx);
    ws->ky = fftw_alloc_real(ws->Nky);
	kxmin = 0.;
	kxmax = (double) (ws->Nx/2);
	kymin = -(double) (ws->Ny/2);
	kymax = (double) (ws->Ny/2-1);
    ws->kx = linspace(kxmin, kxmax, ws->Nkx, ws->kx);
    ws->ky = linspace(kymin, kymax, ws->Nky, ws->ky);

    /* compute k^2 */
    ws->ksq = fftw_alloc_real(ws->ktot);
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j)
            /* ksq[x][y] = ksq[x + y*Nx] */
            ws->ksq[j+i*ws->Nkx] = SQUARE(ws->kx[j]) + SQUARE(ws->ky[i]);
    /* FIX: this value is zero and would yield NaNs in division */
    ws->ksq[ws->Nkx*ws->Nky/2] = 1.;

    /* allocate omega and u and their transforms
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
                                         ws->_ohat, ws->o,
                                         FFTW_MEASURE);
    ws->u_to_uhat = fftw_plan_dft_r2c_2d(ws->Nx, ws->Ny,
                                         ws->utmp, ws->uhat,
                                         FFTW_MEASURE);
    ws->uhat_to_u = fftw_plan_dft_c2r_2d(ws->Nx, ws->Ny,
                                         ws->uhat, ws->utmp,
                                         FFTW_MEASURE);

    /* set inital value and compute FT */
    memcpy((void *)ws->o, (void *)iv, sizeof(double) * ws->Ntot);
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

    /* allocate helper array for `double *rhs` (is initialized there) */
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
    /* fftw_free(ws->res); */
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
double *time_step(int steps, Workspace *ws)
{
    int i;

    for (i = 0; i < steps; ++i)
        scheme(ws);

    memcpy((void *)ws->_ohat, (void *)ws->ohat,
            sizeof(fftw_complex) * ws->ktot);
    irfft2(&ws->ohat_to_o, ws->o, ws);

    return ws->o;
}


/*********************************************************************/
/* TODO: remove later */
static
void print_complex_array(fftw_complex *arr, size_t x, size_t y,
        char *name)
{
    size_t i, j;

    printf("%s_real = np.array(", name);
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            printf("%s%s%e%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x][REAL],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]" : ",")); 
    puts(")\n\n");
    printf("%s_imag = np.array(", name);
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            printf("%s%s%e%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x][IMAG],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]" : ","));
    puts(")\n\n");
}


static
void print_real_array(double *arr, size_t x, size_t y,
        char *name)
{
    size_t i, j;

    printf("%s = np.array(", name);
    for (i = 0; i < y; ++i)
        for(j = 0; j < x; ++j)
            printf("%s%s%e%s%s",
                    (i == 0 && j == 0 ? "[\n" : ""),
                    (j == 0 ? "[" : ""),
                    arr[j+i*x],
                    (j == x-1 ? "]\n" : ""),
                    (i == y-1 && j == x-1 ? "]" : ","));
    puts(")\n\n");
}


/*********************************************************************/


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

    /* print_complex_array(ohat, ws->Nkx, ws->Nky); */

    /* iFT of ohat, yielding o */
    memcpy((void *)ws->_ohat, (void *)ohat,
            sizeof(fftw_complex) * ws->ktot);
    irfft2(&ws->ohat_to_o, ws->o, ws);


    /* print_real_array(ws->kx, ws->Nkx, 1, "kx"); */
    /* print_real_array(ws->ky, 1, ws->Nky, "ky"); */
    /* print_real_array(ws->ksq, ws->Nkx, ws->Nky, "ksq"); */
    /* puts("# input ohat, after anti aliasing"); */
    /* print_complex_array(ohat, ws->Nkx, ws->Nky, */
    /*         "ohat_c"); */
    /* puts("# original ohat"); */
    /* print_complex_array(ws->ohat, ws->Nkx, ws->Nky, */
    /*         "pde_ohat_c"); */
    /* puts("# irfft(ohat)"); */
    /* print_real_array(ws->o, ws->Nx, ws->Ny, */
    /*         "o_c"); */


    /* printf("# *ohat_c = %p, *pde_ohat_c = %p\n", ohat, ws->ohat); */


    /********************/
    /* x component of u */
    /********************/
    /* uhat_x = I * ky * ohat / k^2 */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            /* multiplying by I switches real and imag parts
             * result of rfft has real part, but quite small (~1e-8)
             * but it is significant!
             * we only care about the imag part of uhat here */
            ws->uhat[idx][REAL] = \
                - ws->ky[i] * ws->ohat[idx][IMAG] / ws->ksq[idx];
            ws->uhat[idx][IMAG] = \
                ws->ky[i] * ws->ohat[idx][REAL] / ws->ksq[idx];
        }


    /* CORRECT */
    /* puts("# uhat = I * ky * ohat / k^2"); */
    /* print_complex_array(ws->uhat, ws->Nkx, ws->Nky, "uhat_x"); */


    /* compute iFT of uhat, yielding utmp */
    irfft2(&ws->uhat_to_u, ws->utmp, ws);

    /* puts("# irfft(uhat)"); */
    /* print_real_array(ws->utmp, ws->Nx, ws->Ny, "ux"); */
    /* print_real_array(ws->o, ws->Nx, ws->Ny, "o"); */

    /* u_x * o */
    for (i = 0; i < ws->Ntot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* puts("# u * o"); */
    /* print_real_array(ws->utmp, ws->Nx, ws->Ny, "uxo_c"); */

    /* compute FT of u * o, yielding uhat */
    rfft2(&ws->u_to_uhat, ws->utmp, ws);

    /* puts("# rfft(u*o)"); */
    /* print_complex_array(ws->uhat, ws->Nkx, ws->Nky, */
    /*         "uo_hat_c"); */
    /* print_real_array(ws->kx, 1, ws->Nkx, "kx"); */
    /* print_complex_array(ws->uhat, ws->Nkx, ws->Nky, "uhat_x"); */

    /* write into result */
    /* res = ohat - I * kx * uhat * dt */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            /* multiplying by I switches real and imag in uhat
             * again, only the imag part is relevant */
            ws->res[idx][REAL] = ohat[idx][REAL] \
                + ws->kx[j] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] = ohat[idx][IMAG] \
                - ws->kx[j] * ws->uhat[idx][REAL] * ws->dt;
        }


    /* TODO: there are a few values which devate too far */
    /* puts("# res = ohat - I * kx * uhat * dt"); */
    /* print_complex_array(ws->res, ws->Nkx, ws->Nky, */
    /*         "res_c"); */


    /********************/
    /* y component of u */
    /********************/
    /* uhat_y = -I * kx * ohat / k^2 */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->uhat[idx][REAL] = \
                + ws->kx[j] * ws->ohat[idx][IMAG] / ws->ksq[idx];
            ws->uhat[idx][IMAG] = \
                - ws->kx[j] * ws->ohat[idx][REAL] / ws->ksq[idx];
        }

    /* print_complex_array(ws->uhat, ws->Nkx, ws->Nky, "uhat_y"); */

    /* compute iFT of uhat, yielding utmp */
    irfft2(&ws->uhat_to_u, ws->utmp, ws);

    /* print_real_array(ws->utmp, ws->Nx, ws->Ny, "uy"); */
    /* print_real_array(ws->o, ws->Nx, ws->Ny, "o"); */

    /* u_y * o */
    for (i = 0; i < ws->Ntot; ++i)
        ws->utmp[i] *= ws->o[i];

    /* print_real_array(ws->utmp, ws->Nkx, ws->Nky, "uyo_c"); */

    /* compute FT of u * o, yielding uhat */
    rfft2(&ws->u_to_uhat, ws->utmp, ws);

    /* print_real_array(ws->ky, 1, ws->Nky, "ky"); */
    /* print_complex_array(ws->uhat, ws->Nkx, ws->Nky, "uhat_y"); */

    /* write into result */
    /* res -= I * ky * dt */
    for (i = 0; i < ws->Nky; ++i)
        for (j = 0; j < ws->Nkx; ++j) {
            idx = j + i * ws->Nkx;
            ws->res[idx][REAL] += \
                ws->ky[i] * ws->uhat[idx][IMAG] * ws->dt;
            ws->res[idx][IMAG] -= \
                ws->ky[i] * ws->uhat[idx][REAL] * ws->dt;
        }


    return ws->res;
}


void scheme(Workspace *ws)
{
    /* TODO: MEMORY LECK */
    ws->ohat = rhs(ws->ohat, ws);
/*     size_t i;
 *     fftw_complex *otmp;
 * 
 *     otmp = rhs(ws->ohat, ws);
 *     for (i = 0; i < ws->ktot; ++i)
 *         ws->ohat[i][IMAG] = otmp[i][IMAG]; */
}

/* shu-osher scheme: 3-step runge-kutta */
/* void scheme(Workspace *ws)
 * {
 *     size_t i;
 *     fftw_complex *otmp;
 * 
 *     [> print_complex_array(ws->ohat, ws->Nkx, ws->Nky, "_ohat"); <]
 * 
 *     [> step one <]
 *     otmp = rhs(ws->ohat, ws);
 * 
 *     [> print_complex_array(ws->ohat, ws->Nkx, ws->Nky, "ohat"); <]
 *     [> print_complex_array(otmp, ws->Nkx, ws->Nky, "_o1"); <]
 * 
 *     for (i = 0; i < ws->ktot; ++i)
 *         otmp[i][IMAG] *= ws->prop_full[i];
 * 
 *     [> print_complex_array(otmp, ws->Nkx, ws->Nky, "o1"); <]
 * 
 *     [> step two <]
 *     otmp = rhs(otmp, ws);
 * 
 *     [> print_complex_array(otmp, ws->Nkx, ws->Nky, "_o2"); <]
 *     [> print_complex_array(ws->ohat, ws->Nkx, ws->Nky, "_ohat"); <]
 * 
 *     for (i = 0; i < ws->ktot; ++i)
 *         otmp[i][IMAG] = .25 * (otmp[i][IMAG] * ws->prop_neg_half[i] \
 *                         + 3. * ws->ohat[i][IMAG] * ws->prop_pos_half[i]);
 * 
 *     [> print_complex_array(otmp, ws->Nkx, ws->Nky, "o2_"); <]
 * 
 * 
 *     [> step three <]
 *     otmp = rhs(otmp, ws);
 * 
 *     [> print_complex_array(otmp, ws->Nkx, ws->Nky, "_o3"); <]
 * 
 *     for (i = 0; i < ws->ktot; ++i)
 *         ws->ohat[i][IMAG] = 1./3. * (2. * otmp[i][IMAG] * ws->prop_pos_half[i] \
 *                             + ws->ohat[i][IMAG] * ws->prop_pos_half[i]);
 * 
 *     [> print_complex_array(ws->ohat, ws->Nkx, ws->Nky, "ohat_"); <]
 * } */


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
 * plan     :   pointer to fftw_plan to execute
 * res      :   double array, in which the result of the iDFT will be written
 * ws       :   Workspace pointer, which holds auxillary values
 *  */
void irfft2(fftw_plan *plan, double *res, Workspace *ws)
{
    double *end;

    /* do dft and center result */
    fftw_execute(*plan);
    rfftshift(res, ws->Nx, ws->Ny);

    /* normalize */
    end = res + ws->Ntot;
    while (res != end)
        *res++ /= ws->Ntot;
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

    /* multiply every second (odd) row with -1 */
    for (i = 0; i < ny; ++i)
        for (j = 0; j < nx; ++j)
            arr[j+i*nx] *= i % 2 == 0 ? 1. : -1.;
}


inline double *linspace(double start, double end, size_t np, double *dst)
{
    double x, dx, *startptr, *endptr;

    x = start;
    dx = (end - start + 1) / ((double) np);
    startptr = dst;
    endptr = dst + np;

    while (dst != endptr)
        *dst++ = x, x += dx;

    return startptr;
}
