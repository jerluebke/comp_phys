# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from typing import List

class Grid:
    def __init__(self, xmin, xmax, n):
        self.f = np.zeros(n)    # RHS: Laplacian(f)
        self.rho = np.zeros(n)  # LHS: rho
        self.d = np.zeros(n)    # defect
        self.x = np.linspace(xmin, xmax, n)
        self.dx = self.x[1]-self.x[0]   # grid resolution

GridList = List[Grid]
IntList = List[int]


def init(rho, xmin, xmax, N, levels):
    grids = []
    for i in range(levels):
        n = N // 2**i
        grids.append(Grid(xmin, xmax, n))
    grids[0].rho = rho if not callable(rho) \
                    else rho(np.linspace(xmin, xmax, N))
    return grids


def Lap(g):
    # 1D
    return (np.roll(g.f, 1) - 2*g.f + np.roll(g.f, -1)) / g.dx**2

def smooth(g, steps):
    for i in range(steps):
        g.f = solver(g.f, g.rho, g.dx)
        g.f[0] = g.f[-1] = 0

def _jacobi(f, rho, dx):
    return (np.roll(f, 1) + np.roll(f, -1) - rho*dx**2)/2.

def _gaus_seidel(f, rho, dx):
    res = np.zeros(f.size)
    for i in range(1, f.size):
        res[i] = (f[i+1] + res[i-1] - rho*dx**2)/2.

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


restrict = _injection
prolong = _interpolation
solver = _jacobi


def solve_one(grids, level, steps, repeat):
    g = grids[level]
    smooth(g, steps[level])
    g.d = g.rho - Lap(g)
    for _ in range(repeat[level]):
        if level < len(grids)-1 and steps[level+1]:
            grids[level+1].rho = restrict(g.d)
            solve_one(grids, level+1, steps, repeat)
            g.f += prolong(grids[level+1].f)
        smooth(g, steps[level])


def solve(grids, steps, repeat, tol=1e-7, maxiter=1000):
    for i in range(maxiter):
        solve_one(grids, 0, steps, repeat)
        d = grids[0].rho - Lap(grids[0])
        if np.mean(np.abs(d[1:-1])) < tol:
            #  print('convergence after %d iterations' % i)
            return i
    return -1


def rho(x):
    return np.sin(x)*np.exp(-x**2)


%timeit solve(init(rho, -5, 5, 64, 4), [100, 30, 0], [1, 1, 1])
%timeit solve(init(rho, -5, 5, 64, 3), [50, 20, 0], [1, 1, 1])
%timeit solve(init(rho, -5, 5, 64, 3), [1, 0, 0], [1, 1, 1])
