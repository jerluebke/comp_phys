#include <stdlib.h>
#include <string.h>     /* memcpy */
#include <math.h>
#include "../include/navier_stokes.h"


/* fftw_complex *rhs(fftw_complex *, PDE *); */
static void scheme(PDE *);
static fftw_complex *rhs(fftw_complex *, PDE *);
static inline void rfftshift(double *, size_t, size_t);
static void rfft2(const fftw_plan *, double *, PDE *);
static void irfft2(const fftw_plan *, double *, PDE *);


/* set up workspace
 * allocate memory
 * init fftw plans
 *
 * call `cleanup` on the pointer returned by this function to free all
 *     allocated memory properly */
PDE *init(const Params params, const double *iv)
{
    size_t i, j;
    double kxmin, kxmax, kymin, kymax;
    PDE *this;

    /* allocate PDE struct */
    this = malloc(sizeof(*this));

    /* init parameter */
    this->Nx   = params.Nx;
    this->Ny   = params.Ny;
    this->Nkx  = params.Nx / 2 + 1;
    this->Nky  = params.Ny;
    this->Ntot = this->Nx * this->Ny;
    this->ktot = this->Nkx * this->Nky;
    this->dt   = params.dt;
    this->nu   = params.nu;

    /* k-space: allocate and init as linspace */
    this->kx = fftw_alloc_real(this->Nkx);
    this->ky = fftw_alloc_real(this->Nky);
	kxmin = 0.;
	kxmax = (double) (this->Nx/2);
	kymin = -(double) (this->Ny/2);
	kymax = (double) (this->Ny/2-1);
    this->kx = linspace(kxmin, kxmax, this->Nkx, this->kx);
    this->ky = linspace(kymin, kymax, this->Nky, this->ky);

    /* compute k^2 */
    this->ksq = fftw_alloc_real(this->ktot);
    for (i = 0; i < this->Nky; ++i)
        for (j = 0; j < this->Nkx; ++j)
            /* ksq[x][y] = ksq[x + y*Nx] */
            this->ksq[j+i*this->Nkx] = SQUARE(this->kx[j]) + SQUARE(this->ky[i]);
    /* FIX: this value is zero and would yield NaNs in division */
    this->ksq[this->Nkx*this->Nky/2] = 1.;

    /* allocate omega and u and their transforms
     * include backup memory for inverse transforms */
    this->o    = fftw_alloc_real(this->Ntot);
    this->ohat = fftw_alloc_complex(this->ktot);
    this->_ohat= fftw_alloc_complex(this->ktot);
    this->utmp = fftw_alloc_real(this->Ntot);
    this->uhat = fftw_alloc_complex(this->ktot);

    /* set up fftw plan */
    this->o_to_ohat = fftw_plan_dft_r2c_2d(this->Nx, this->Ny,
                                         this->o, this->ohat,
                                         FFTW_MEASURE);
    this->ohat_to_o = fftw_plan_dft_c2r_2d(this->Nx, this->Ny,
                                         this->_ohat, this->o,
                                         FFTW_MEASURE);
    this->u_to_uhat = fftw_plan_dft_r2c_2d(this->Nx, this->Ny,
                                         this->utmp, this->uhat,
                                         FFTW_MEASURE);
    this->uhat_to_u = fftw_plan_dft_c2r_2d(this->Nx, this->Ny,
                                         this->uhat, this->utmp,
                                         FFTW_MEASURE);

    /* set inital value and compute FT */
    memcpy((void *)this->o, (void *)iv, sizeof(double) * this->Ntot);
    rfft2(&this->o_to_ohat, this->o, this);

    /* allocate and init mask for anti-aliasing */
    this->mask = malloc(sizeof(unsigned char) * this->ktot);
    double threshold = this->Ntot / 9.;
    for (i = 0; i < this->ktot; ++i)
        this->mask[i] = this->ksq[i] < threshold ? 1 : 0;

    /* allocate and init propagator lookup tables */
    this->prop_full        = malloc(sizeof(double) * this->ktot);
    this->prop_pos_half    = malloc(sizeof(double) * this->ktot);
    this->prop_neg_half    = malloc(sizeof(double) * this->ktot);
    for (i = 0; i < this->ktot; ++i) {
        this->prop_full[i]     = exp(-this->nu * this->ksq[i] * this->dt);
        this->prop_pos_half[i] = exp(-.5 *this->nu * this->ksq[i] * this->dt);
        this->prop_neg_half[i] = exp(.5 * this->nu * this->ksq[i] * this->dt);
    }

    /* allocate helper array for `double *rhs` (is initialized there) */
    this->res  = fftw_alloc_complex(this->ktot);

    /* done */
    return this;
}


/* free all allocated memory in a workspace struct
 * the pointer is afterwards NULL */
void cleanup(PDE *this)
{
    fftw_free(this->kx);
    fftw_free(this->ky);
    fftw_free(this->ksq);
    fftw_free(this->o);
    fftw_free(this->ohat);
    fftw_free(this->_ohat);
    fftw_free(this->utmp);
    fftw_free(this->uhat);
    fftw_free(this->res);
    fftw_destroy_plan(this->o_to_ohat);
    fftw_destroy_plan(this->ohat_to_o);
    fftw_destroy_plan(this->u_to_uhat);
    fftw_destroy_plan(this->uhat_to_u);
    fftw_cleanup();
    free(this->mask);
    free(this->prop_full);
    free(this->prop_pos_half);
    free(this->prop_neg_half);
    free(this);
}


/* calculate the next frame */
double *time_step(unsigned int steps, PDE *this)
{
    unsigned int i;

    for (i = 0; i < steps; ++i)
        scheme(this);

    memcpy((void *)this->_ohat, (void *)this->ohat,
            sizeof(fftw_complex) * this->ktot);
    irfft2(&this->ohat_to_o, this->o, this);

    return this->o;
}


/* right hand side of equation */
static fftw_complex *rhs(fftw_complex *ohat, PDE *this)
{
    size_t i, j, idx;

    /* anti aliasing */
    for (i = 0; i < this->ktot; ++i)
        if (! this->mask[i]) {
            ohat[i][REAL] = 0.;
            ohat[i][IMAG] = 0.;
        }

    /* iFT of ohat, yielding o */
    memcpy((void *)this->_ohat, (void *)ohat,
            sizeof(fftw_complex) * this->ktot);
    irfft2(&this->ohat_to_o, this->o, this);


    /********************/
    /* x component of u */
    /********************/
    /* uhat_x = I * ky * ohat / k^2 */
    for (i = 0; i < this->Nky; ++i)
        for (j = 0; j < this->Nkx; ++j) {
            idx = j + i * this->Nkx;
            /* multiplying by I switches real and imag parts
             * result of rfft has real part, but quite small (~1e-8)
             * but it is significant!
             * Both parts (real and imag) are needed to calculate time
             * development! */
            this->uhat[idx][REAL] = \
                - this->ky[i] * this->ohat[idx][IMAG] / this->ksq[idx];
            this->uhat[idx][IMAG] = \
                this->ky[i] * this->ohat[idx][REAL] / this->ksq[idx];
        }

    /* compute iFT of uhat, yielding utmp */
    irfft2(&this->uhat_to_u, this->utmp, this);

    /* u_x * o */
    for (i = 0; i < this->Ntot; ++i)
        this->utmp[i] *= this->o[i];

    /* compute FT of u * o, yielding uhat */
    rfft2(&this->u_to_uhat, this->utmp, this);

    /* write into result */
    /* res = ohat - I * kx * uhat * dt */
    for (i = 0; i < this->Nky; ++i)
        for (j = 0; j < this->Nkx; ++j) {
            idx = j + i * this->Nkx;
            /* multiplying by I switches real and imag in uhat */
            this->res[idx][REAL] = ohat[idx][REAL] \
                + this->kx[j] * this->uhat[idx][IMAG] * this->dt;
            this->res[idx][IMAG] = ohat[idx][IMAG] \
                - this->kx[j] * this->uhat[idx][REAL] * this->dt;
        }


    /********************/
    /* y component of u */
    /********************/
    /* uhat_y = -I * kx * ohat / k^2 */
    for (i = 0; i < this->Nky; ++i)
        for (j = 0; j < this->Nkx; ++j) {
            idx = j + i * this->Nkx;
            this->uhat[idx][REAL] = \
                + this->kx[j] * this->ohat[idx][IMAG] / this->ksq[idx];
            this->uhat[idx][IMAG] = \
                - this->kx[j] * this->ohat[idx][REAL] / this->ksq[idx];
        }

    /* compute iFT of uhat, yielding utmp */
    irfft2(&this->uhat_to_u, this->utmp, this);

    /* u_y * o */
    for (i = 0; i < this->Ntot; ++i)
        this->utmp[i] *= this->o[i];

    /* compute FT of u * o, yielding uhat */
    rfft2(&this->u_to_uhat, this->utmp, this);

    /* write into result */
    /* res -= I * ky * dt */
    for (i = 0; i < this->Nky; ++i)
        for (j = 0; j < this->Nkx; ++j) {
            idx = j + i * this->Nkx;
            this->res[idx][REAL] += \
                this->ky[i] * this->uhat[idx][IMAG] * this->dt;
            this->res[idx][IMAG] -= \
                this->ky[i] * this->uhat[idx][REAL] * this->dt;
        }


    return this->res;
}


/* shu-osher scheme: 3-step runge-kutta */
static void scheme(PDE *this)
{
    size_t i;
    fftw_complex *otmp;

    /* step one */
    otmp = rhs(this->ohat, this);
    for (i = 0; i < this->ktot; ++i) {
        otmp[i][REAL] *= this->prop_full[i];
        otmp[i][IMAG] *= this->prop_full[i];
    }

    /* step two */
    otmp = rhs(otmp, this);
    for (i = 0; i < this->ktot; ++i) {
        otmp[i][REAL] = .25 * (otmp[i][REAL] * this->prop_neg_half[i] \
                        + 3. * this->ohat[i][REAL] * this->prop_pos_half[i]);
        otmp[i][IMAG] = .25 * (otmp[i][IMAG] * this->prop_neg_half[i] \
                        + 3. * this->ohat[i][IMAG] * this->prop_pos_half[i]);
    }

    /* step three */
    otmp = rhs(otmp, this);
    for (i = 0; i < this->ktot; ++i) {
        this->ohat[i][REAL] = 1./3. * (2. * otmp[i][REAL] * this->prop_pos_half[i] \
                            + this->ohat[i][REAL] * this->prop_pos_half[i]);
        this->ohat[i][IMAG] = 1./3. * (2. * otmp[i][IMAG] * this->prop_pos_half[i] \
                            + this->ohat[i][IMAG] * this->prop_pos_half[i]);
    }
}  


/* perform real DFT
 *
 * Params
 * ======
 * plan    :   pointer to fftw_plan to execute
 * orig    :   array which will be transformed
 * this      :   PDE pointer holing auxillary values
 *  */
static void rfft2(const fftw_plan *plan, double *orig, PDE *this)
{
    /* prepare input for result with origin shifted to center */
    rfftshift(orig, this->Nx, this->Ny);
    /* do dft */
    fftw_execute(*plan);
    /* reverse changes */
    rfftshift(orig, this->Nx, this->Ny);
}


/* perform inverse DFT and normalize result
 *
 * Params
 * ======
 * plan     :   pointer to fftw_plan to execute
 * res      :   double array, in which the result of the iDFT will be written
 * this       :   PDE pointer, which holds auxillary values
 *  */
static void irfft2(const fftw_plan *plan, double *res, PDE *this)
{
    double *end;

    /* do dft and center result */
    fftw_execute(*plan);
    rfftshift(res, this->Nx, this->Ny);

    /* normalize */
    end = res + this->Ntot;
    while (res != end)
        *res++ /= this->Ntot;
}


/* prepare input such that its DFT will have the origin frequence in the center
 *
 * Params
 * ======
 * arr     :   double array, input
 * nx, ny  :   int, dimension sizes
 *  */
static inline void rfftshift(double *arr, size_t nx, size_t ny)
{
    size_t i, j;

    /* multiply every second (odd) row with -1 */
    for (i = 0; i < ny; ++i)
        for (j = 0; j < nx; ++j)
            arr[j+i*nx] *= i % 2 == 0 ? 1. : -1.;
}


double *linspace(double start, double end, size_t np, double *dst)
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
