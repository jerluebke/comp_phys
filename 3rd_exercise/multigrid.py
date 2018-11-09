# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse


class Grid:
    def __init__(self, n, xmin, xmax):
        self.x  = np.linspace(xmin, xmax, n)
        self.dx = self.x[1] - self.x[0]
        self.f  = np.zeros(n)
        self.rho= np.zeros(n)
        self.d  = np.zeros(n)
        self.L  = sparse.diags([1, -2, 1], [-1, 0, 1], (n-2, n-2), \
                               format='csc') / self.dx**2
        self.n = n

    @property
    def defect(self):
        n = self.n
        self.d[1:n-1] = self.rho[1:n-1] - self.L @ self.f[1:n-1]
        #  self.d[1:n-1] = self.rho[1:n-1] \
        #          - (self.f[2:n] - 2. * self.f[1:n-1] + self.f[0:n-2]) \
        #          / self.dx**2
        return self.d


def init(rho_func, N, xmin, xmax, levels):
    grids = [Grid(N // 2**i, xmin, xmax) for i in range(levels)]
    grids[0].rho[1:N-2] = rho_func(grids[0].x[1:N-2])
    return grids

def smooth(g, **kwds):
    n = g.n
    # jacobi
    g.f[1:n-1] = .5 * (g.f[0:n-2] + g.f[2:n] - g.dx**2 * g.rho[1:n-1])
    # omega-jacobi
    #  omega = kwds.get('omega', .5)
    #  g.f[1:n-1] = .5 * omega * (g.f[0:n-2] + g.f[2:n] - g.dx**2 * g.rho[1:n-1])\
    #          + (1. - omega) * g.f[1:n-1]
    # gauss-seidel
    #  g.f[1:n-1] = np.array([.5 * (g.f[i-1] + g.f[i+1] - g.dx**2 * g.rho[i])
    #                         for i in range(1, n-1)])
    # red-black gs
    #  g.f[1:n-1:2] = .5 * (g.f[0:n-2:2] + g.f[2:n:2] - g.dx**2 * g.rho[1:n-1:2])
    #  g.f[2:n-2:2] = .5 * (g.f[1:n-3:2] + g.f[3:n-1:2] - g.dx**2 * g.rho[2:n-2:2])
    # von Neumann boundary condition
    g.f[0] = g.f[-1] = 0

def restrict(arr):
    # injection
    #  return arr[2:arr.size-2:2]
    # full-weighting
    nf = arr.size
    nc = nf // 2
    res = np.zeros(nc)
    res[1:nc-1] = .5 * arr[2:nf-2:2] + .25 * (arr[3:nf-1:2] + arr[1:nf-3:2])
    return res

def prolong(arr):
    nc = arr.size
    nf = 2 * nc
    res = np.zeros(nf)
    res[2:nf-2:2] = arr[1:nc-1]
    res[1:nf-1:2] = (arr[0:nc-1] + arr[1:nc]) * .5
    return res

def solve_one_v(grids, level=0):
    fine = grids[level]
    smooth(fine)
    if level < len(grids)-1:
        coarse = grids[level+1]
        coarse.rho = restrict(fine.defect)
        solve_one_v(grids, level+1)
        fine.f += prolong(coarse.f)
    smooth(fine)


def solve(imax=100, tol=1e-7, ax=None, **gridkwds):
    if not ax:
        ax = plt.gca()
    grids = init(**gridkwds)
    err = np.zeros(imax)
    i = 0
    while i < imax:
        solve_one_v(grids)
        err[i] = np.max(np.abs(grids[0].defect))
        if err[i] < tol:
            print('convergence after %d iterations' % i)
            break
        i += 1
    ax.semilogy(err)
    return grids

N = 256
gridkwds = dict(rho_func=lambda x: np.sin(x)*np.exp(-x**2),
                N=N, xmin=-5, xmax=5, levels=int(math.log(N, 2)+1e-2)-1)
g = solve(**gridkwds)
