# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt

def crank_nicolson(f, D, dt, steps):
    res = np.zeros((f.size, steps))
    res[:,0] = f
    # sparse.linalg.inv is more efficient with csc format
    I = sparse.identity(f.size, format='csc')
    A = (I + dt * D)
    B = linalg.inv(I - dt * D)
    for i in range(1, steps):
        res[:,i] = (B @ A) @ res[:,i-1]
    return res

def fwd_euler(f, D, dt, steps):
    res = np.zeros((f.size, steps))
    res[:,0] = f
    I = sparse.identity(f.size, format='csc')
    for i in range(1, steps):
        res[:,i] = (I + dt* D) @ res[:,i-1]
    return res

def leapfrog(f0, f1, D, dt, steps):
    res = np.zeros((f0.size, steps))
    res[:,0] = f0
    res[:,1] = f1
    for i in range(2, steps):
        res[:,i] = res[:,i-2] + 2 * dt * D @ res[:,i-1]
    return res


def test(s=1., N=100, steps=200):
    k = 1
    x = np.linspace(-6, 6, N)
    dx_sq = (x[1] - x[0])**2
    dt = dx_sq / (2. * k)   # cfl condition, nyquist wavelength
    L = sparse.diags([1., -2., 1.], [-1, 0, 1], shape=(N, N), format='csc') / dx_sq
    #  get_f0 = lambda s: np.exp(-x**2/(2*s**2)) / (s*np.sqrt(2*np.pi))
    #  f0 = get_f0(s)
    f0 = np.piecewise(x, [np.abs(x) < 1.], [1.])
    f1 = fwd_euler(f0, L*k, dt, 2)[:,1]

    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.suptitle('heat diffusion in 1D over time')
    ax[0].set(ylabel='x', title='crank-nicolson')
    ax[1].set(xlabel='time', ylabel='x', title='forward euler')
    im0 = ax[0].imshow(crank_nicolson(f0, L*k/2., dt, steps))
    # weird results...
    #  im1 = ax[1].imshow(leapfrog(f0, f1, k*L, dt/1024., steps))
    im1 = ax[1].imshow(fwd_euler(f0, k*L, dt/2, steps))
    return im0, im1
