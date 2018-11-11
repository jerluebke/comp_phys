# -*- coding: utf-8 -*-

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import cg
sys.path.append('..\\3rd_exercise')
import multigrid as mg


def rho_func(x):
    return np.sin(x) * np.exp(-x**2)

N = 512
gridkwds = dict(rho_func=rho_func,
                N=N, xmin=-5, xmax=5,
                levels=int(math.log(N, 2)+1e-2)-1)

rho = rho_func(np.linspace(-5, 5, N))
L = sparse.diags([1, -2, 1], [-1, 0, 1], (N, N), format='csc')


def compare(imax_cg, imax_mg):
    x_cg = np.arange(imax_cg)
    x_mg = np.arange(0, imax_cg, imax_cg//imax_mg)
    err_cg          = cg.cg(L, rho, imax=imax_cg)[2]
    err_jacobi      = mg.err(solver='jacobi', imax=imax_mg, **gridkwds)
    err_omegajac    = mg.err(solver='omega_jacobi', imax=imax_mg, **gridkwds)
    err_gaussseidel = mg.err(solver='gauss_seidel', imax=imax_mg, **gridkwds)
    err_redblack    = mg.err(solver='red_black', imax=imax_mg, **gridkwds)

    ax = plt.gca()
    ax.semilogy(x_cg, err_cg, label='cg')
    ax.semilogy(x_mg, err_jacobi, label='mg - jacobi')
    ax.semilogy(x_mg, err_omegajac, label='mg - omega-jacobi')
    ax.semilogy(x_mg, err_gaussseidel, label='mg - gauss-seidel')
    ax.semilogy(x_mg, err_redblack, label='mg - red-black')
    ax.title('comparision - Multigrid and Conjugate Gradient')
    ax.legend()


compare(500, 30)
