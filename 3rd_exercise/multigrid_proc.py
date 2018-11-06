# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from copy import copy


class Grid:
    def __init__(self, xmin, xmax, n):
        if not n:
            pass
        else:
            self.f = np.zeros(n)    # RHS: Laplacian(f)
            self.rho = np.zeros(n)  # LHS: rho
            self.d = np.zeros(n)    # defect
            self.x = np.linspace(xmin, xmax, n)
            self.dx = self.x[1]-self.x[0]   # grid resolution

    def __copy__(self):
        copied = type(self)(0, 0, 0)
        for k, v in self.__dict__.items():
            if hasattr(v, 'copy'):
                copied.__dict__[k] = v.copy()
            else:
                copied.__dict__[k] = v
        return copied


def init(rho, xmin, xmax, N, levels):
    grids = []
    for i in range(levels):
        n = N // 2**i
        grids.append(Grid(xmin, xmax, n))
    grids[0].rho = rho if not callable(rho) \
                    else rho(grids[0].x)
    return grids


def Lap(g):
    # 1D
    return (np.roll(g.f, 1) - 2*g.f + np.roll(g.f, -1)) / g.dx**2

def Lap_arr(f, dx):
    return (np.roll(f, 1) - 2*f + np.roll(f, -1)) / dx**2

def smooth(g, steps):
    for i in range(steps):
        g.f = SOLVER(g.f, g.rho, g.dx)
        g.f[0] = g.f[-1] = 0

def _jacobi(f, rho, dx):
    return (np.roll(f, 1) + np.roll(f, -1) - rho*dx**2)/2.

def _omega_jacobi(f, rho, dx, omega):
    return (1-omega)*f+omega*(np.roll(f, 1)+np.roll(f, -1)-dx**2*rho)/2.

def _gauss_seidel(f, rho, dx):
    res = np.zeros(f.size)
    for i in range(1, f.size-1):
        res[i] = (f[i+1] + res[i-1] - rho[i]*dx**2)/2.
    return res

def _injection(f):
    return np.array([f[2*i] for i in range(f.size//2)])

def _weighting(f):
    return np.array([(f[2*i-1] + f[2*i+1])/4. + f[2*i]/2.
                     for i in range(1, f.size//2-1)])

def _interpolation(f):
    res = np.zeros(f.size*2)
    for i in range(f.size-1):
        res[2*i] = f[i]
        res[2*i+1] = (f[i]+f[i+1])/2.
    return res
    #  return np.reshape(np.array([(f[i], (f[i]+f[i+1])/2.)
    #                              for i in range(f.size-1)]), (f.size*2,))


RESTRICT = _injection
#  RESTRICT = _weighting
PROLONG = _interpolation
SOLVER = _jacobi


def solve_one_v(grids, level, steps):
    g = grids[level]
    smooth(g, steps[level])
    iters = steps[level]
    g.d = g.rho - Lap(g)
    if level < len(grids)-1 and steps[level+1]:
        grids[level+1].rho = RESTRICT(g.d)
        iters += solve_one_v(grids, level+1, steps)
        g.f += PROLONG(grids[level+1].f)
    smooth(g, steps[level])
    return iters + steps[level]


def rho_1(x):
    return np.sin(x)*np.exp(-x**2)

def rho_2(x, lo=-10, hi=11, L=1):
    return np.sum([np.sin(2*np.pi*x/L*2**i) for i in range(lo, hi)],
                  axis=0)

def rho_3(x, d=1):
    return rho_2(x)*np.exp(-x**2*np.pi*d)


#  %timeit solve(init(rho, -5, 5, 64, 4), [100, 30, 0], [1, 1, 1])
#  %timeit solve(init(rho, -5, 5, 64, 3), [50, 20, 0], [1, 1, 1])
#  %timeit solve(init(rho, -5, 5, 64, 3), [1, 0, 0], [1, 1, 1])

def test_mg(N, rho_func,
            solver=(_jacobi, 'jacobi'),
            steps=[2, 50],
            solve_func=solve_one_v,
            tol=1e-7, maxsteps=1000,
            ax=None, label='v'):
    levels = len(steps)

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.set(xlabel='iterations', ylabel='error mean',
           title='Multigrid - datapoints = %d, solver = %s' % (N, solver[1]))

    g = init(rho_func, 0, 1, N, levels)
    res_grid = None

    iters = np.zeros(maxsteps)
    err = np.zeros(maxsteps)
    conv = False
    for i in range(1, maxsteps):
        iters[i] = iters[i-1] + solve_one_v(g, 0, steps)
        err[i] = np.mean((g[0].rho - Lap(g[0]))[1:-1])
        if not conv and abs(err[i]) < tol:
            res_grid = copy(g[0])
            ax.plot([iters[i]], [err[i]], 'ro', zorder=2)
            conv = True
            print('convergence after %d iterations' % iters[i])

    ax.plot(iters, err, label='%s - %s' % (label, str(steps)), zorder=1)
    ax.plot([0, iters[-1]], [tol, tol], 'k:')
    ax.legend()

    return res_grid


def test_solvers(N, rho_func, solvers, tol=1e-7, maxiter=10000):
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='iterations', ylabel='error mean',
                         title='Comparison of solvers - datapoints=%d' % N)

    x = np.linspace(0, 1, N)
    dx = x[1] - x[0]
    rho = rho_func(x)

    for solver, label in solvers:
        f = np.zeros(N)
        err = np.zeros(maxiter)
        conv = False
        for i in range(maxiter):
            f = solver(f, rho, dx)
            f[0] = f[-1] = 0
            err[i] = np.mean((rho-Lap_arr(f, dx))[1:-1])
            if conv and err[i] < tol:
                ax.plot([i], [err[i]], 'ro', zorder=2)
                conv = True
                print('%s reached tol after %d iterations' % (label, i))
        ax.loglog(err, label=label, zorder=1)

    ax.plot([0, maxiter], [tol, tol], 'k:')
    ax.legend()

solvers = [(_jacobi, 'jacobi'),
           (lambda f, r, d: _omega_jacobi(f, r, d, 2./3.), 'jacobi, omega=2/3'),
           (_gauss_seidel, 'gauss-seidel')]

#  test_solvers(100, rho_1, solvers, maxiter=20000)
