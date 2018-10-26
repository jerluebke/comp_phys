# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt

def crank_nicolson(f, D, dt, steps):
    res = np.zeros((f.size, steps))
    res[:,0] = f
    I = sparse.identity(f.size, format='csc')
    A = (I + dt * D)
    B = linalg.inv(I - dt * D)
    for i in range(1, steps):
        res[:,i] = (B @ A) @ res[:,i-1]
    return res


def test(s=1., N=100, steps=200):
    k = 1
    x = np.linspace(0, 10, N)
    dx_sq = (x[1] - x[0])**2
    L = sparse.diags([1., -2., 1.], [-1, 0, 1], shape=(N, N), format='csc') / dx_sq
    get_f0 = lambda s: np.exp(-x**2/(2*s**2)) / (s*np.sqrt(2*np.pi))
    f0 = get_f0(s)

    fig, ax = plt.subplots()
    ax.set(xlabel='time', ylabel='x', title='heat diffusion in 1D over time')
    im = ax.imshow(crank_nicolson(f0, L*k/2., dx_sq/(2.*k), steps))
    return im
