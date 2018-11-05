# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from collections import named_tuple
from typing import List

Grid = named_tuple('Grid', [
    'n',    # size
    'f',    # RHS: Laplace(f)
    'rho',  # LHS: rho
    'd',    # defect
    'dx',   # grid resolution
    'x'     # data
])

GridArray = List[Grid]
IntList = List[int]


def Lap(f : np.ndarray, dx : float) -> np.ndarray:
    # 1D
    return (np.roll(f, 1) - 2*f + np.roll(f, -1)) / dx**2

def defect(g : Grid) -> None:
    g.d = g.rho - Lap(g.f, g.dx)

def smooth(g : Grid, steps : int) -> None:
    pass

def _jacobi(f : np.array,
            rho : np.array,
            dx : float) -> np.array():
    return (np.roll(f, 1) + np.roll(f, -1) - dx**2*rho)/2.

def _injection(f : np.ndarray) -> np.ndarray:
    res = np.zeros(f.size//2)
    for i in range(res.size):
        res[i] = f[2*i]
    return res
    #  return np.array([f[2*i] for i in range(f.size//2)])

def _weighting(f : np.ndarray) -> np.ndarray:
    res = np.zeros(f.size//2)
    for i in range(1, res.size-1):
        res[i] = (f[2*i-1] + f[2*i+1])/4. + f[2*i]/2.
    return res

def _interpolation(f : np.ndarray) -> np.ndarray:
    res = np.zeros(f.size*2)
    for i in range(f.size-1):
        res[2*i] = f[i]
        res[2*i+1] = (f[i] + f[i+1])/2.
    return res


restrict = _injection
prolong = _interpolation


def solve_one(grids : GridArray,
              level : int,
              steps : IntList) -> None:
    g = grids[level]
    smooth(g, steps[level])
    defect(g)
    if level < len(grids)-1 and steps[level]:
        grids[level+1].rho = restrict(g.defect)
        solve_one(grids, level+1, steps)
        g.f += prolong(grids[level+1].f)
    smooth(g, steps[level])
