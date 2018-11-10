# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse


# runtime parameters
n = 256     # finest resolution
xmin = -5.0  # left border of domain
xmax =  5.0  # right border of domain
iterations = 30 # number of repititions of multigrid algorithm

# initialization for rho
init_rho = lambda x: np.sin(x) * np.exp(-x**2)

class Grid:
    def __init__(self, n, xmin, xmax):
        '''grid-data for one level; n + 1 gridpoints, range from xmin to xmax'''
        self.n = n + 1
        self.f = np.zeros(n + 1)
        self.rho = np.zeros(n + 1)
        self.defect = np.zeros(n + 1)
        self.dx = (xmax - xmin) / n
        self.x = np.linspace(xmin, xmax, n + 1) # g.x = xmin + (0:n)' * g.dx;

    def __repr__(self):
        return "Grid(" + str(self.n + 1) + ", " \
                + str(self.x[0]) + ", " + str(self.x[-1])+ ")"

def defect(grid):
    '''update defect array of grid'''
    n = grid.n
    grid.defect[1:n-1] = grid.rho[1:n-1] - \
                (grid.f[2:n] - 2.*grid.f[1:n-1] + grid.f[0:n-2]) / grid.dx**2

def prolong(coarse, fine):
    '''coarse to fine, fine is modified'''
    nc = coarse.n
    nf = fine.n
    fine.f[2:nf-2:2] += coarse.f[1:nc-1]
    fine.f[1:nf-1:2] += 0.5 * (coarse.f[0:nc-1] + coarse.f[1:nc])

def restrict(fine, coarse):
    '''fine to coarse, coarse is modified'''
    nc = coarse.n
    nf = fine.n
    coarse.rho[0:nc] = 0.
    # injektion
    #coarse.rho[1:nc-1] = fine.defect[2:nf-2:2]
    # full weighting
    coarse.rho[1:nc-1] = 0.5*fine.defect[2:nf-2:2]\
                        + 0.25*(fine.defect[3:nf-1:2] + fine.defect[1:nf-3:2])
    coarse.f[0:nc] = 0.

def smooth(grid):
    n = grid.n
    # Jacobi:
    # grid.f[1:n-1] = 0.5 * (grid.f[0:n-2] + grid.f[2:n] - grid.dx**2 * grid.rho[1:n-1])
    # omega-Jacobi:
    # omega = 0.5
    # grid.f[1:n-1] = 0.5 * omega * ( grid.f[0:n-2] + grid.f[2:n] - grid.dx**2 * grid.rho[1:n-1] ) \
    #              + (1.-omega) * grid.f[1:n-1];
    # Gauss-Seidel:
    # for i in range(1, n-1):
    #     grid.f[i] = 0.5 * (grid.f[i-1] + grid.f[i+1] - grid.dx**2 * grid.rho[i])
    # Red-Black Gauss-Seidel:
    grid.f[1:n-1:2] = 0.5 * (grid.f[0:n-2:2] + grid.f[2:n:2] - grid.dx**2 * grid.rho[1:n-1:2])
    grid.f[2:n-2:2] = 0.5 * (grid.f[1:n-3:2] + grid.f[3:n-1:2] - grid.dx**2 * grid.rho[2:n-2:2])

def solve_one(grids, level):
    '''one iteration on the specified level, calls itself recursively'''
    smooth(grids[level])
    if level < len(grids) - 1:
        defect(grids[level])
        restrict(grids[level], grids[level+1])
        solve_one(grids, level+1)
        prolong(grids[level+1], grids[level])
    smooth(grids[level])

# initialize grids
eps = 1.e-2 # for circumventing discrete floating point effects
grids = [Grid(n//2**i, xmin, xmax) for i in range(int(math.log(n, 2) + eps) - 1)]
grids[0].rho[1:grids[0].n-1] = init_rho(grids[0].x[1:grids[0].n-1])


def solve(imax=100, tol=1e-7, ax=None, **gridkwds):
    if not ax:
        ax = plt.gca()
    #  grids = init(**gridkwds)
    err = []
    i = 0
    while i < imax:
        solve_one(grids, 0)
        err.append(np.max(grids[0].defect))
        i += 1
    ax.semilogy(np.arange(i), err)
    return grids

N = 256
gridkwds = dict(rho_func=lambda x: np.sin(x)*np.exp(-x**2),
                N=N, xmin=-5, xmax=5, levels=int(math.log(N, 2)+1e-2)-1)
g = solve(**gridkwds)
