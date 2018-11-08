# -*- coding: utf-8 -*-

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
        self.L  = sparse.diags([1, -2, 1], [-1, 0, 1], (n, n), \
                               format='csc') / self.dx**2

def init(rho_func, N, xmin, xmax, levels):
    grids = []
    for i in range(levels):
        n = N // 2**i
        grids.append(Grid(n, xmin, xmax))
    grids[0].rho = rho_func(grids[0].x)
    return grids

def smooth(g):
    # jacobi
    g.f = (np.roll(g.f, 1) + np.roll(g.f, -1) - g.rho*g.dx**2) * .5
    g.f[0] = g.f[-1] = 0

def restrict(arr):
    # injection
    return arr[1::2]
    # full-weighting
    # wrong result...
    #  res = np.zeros(arr.size//2)
    #  res[:res.size-1] += (arr[:arr.size-2:2] + arr[2::2]) / 4.
    #  res += arr[1::2] / 2.
    #  return res

def prolong(arr):
    res = np.zeros(2*arr.size)
    res[2::2] = arr[1:]
    res[1:res.size-1:2] = (arr[:arr.size-1] + arr[1:]) * .5
    return res

def solve_one_v(grids, level=0):
    g = grids[level]
    smooth(g)
    if level < len(grids)-1:
        g.d = g.rho - g.L @ g.f
        grids[level+1].rho = restrict(g.d)
        solve_one_v(grids, level+1)
        g.f += prolong(grids[level+1].f)
    smooth(g)


def solve(imax=100, tol=1e-7, ax=None, **gridkwds):
    if not ax:
        ax = plt.gca()
    grids = init(**gridkwds)
    err = np.zeros(imax)
    i = 0
    while i < imax:
        solve_one_v(grids)
        err[i] = np.max(np.abs(
            grids[0].rho - grids[0].L @ grids[0].f)[1:-1])
        if err[i] < tol:
            print('convergence after %d iterations' % i)
            break
        i += 1
    ax.semilogy(err)
    return grids

N = 256
gridkwds = dict(rho_func=lambda x: np.sin(x)*np.exp(-x**2),
                N=N, xmin=-5, xmax=5, levels=int(np.log(N)))
g = solve(**gridkwds)
