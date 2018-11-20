# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg


#  Schrödingers equation with crank-nicolson

def crank_nicolson(f, H, dt, steps, *fargs, **fkwds):
    I = sparse.identity(f.size, format='csc')
    A = I - 1j * dt / 2. * H
    B = I + 1j * dt / 2. * H
    B_LU = linalg.splu(B)
    for i in range(steps):
        f = B_LU.solve(A.dot(f))
        yield f

a = -50
b = 50
N = 1000

s = 4.
x0 = 0.     # inital position
k0 = 6.     # inital momentum
x = np.linspace(a, b, N)
dx_sq = (x[1] - x[0])**2
f0 = np.exp(1j*k0*x-(x-x0)**2/(2*s**2)) / (s*math.sqrt(2*np.pi))

L = sparse.diags([1., -2., 1.], [-1, 0, 1], shape=(N, N),
                 format='csc') / dx_sq
