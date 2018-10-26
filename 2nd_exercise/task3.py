# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt

def rk4(F, x, t, dt, steps, *fargs, **fkwds):
    """
    classical runge-kutta (4 steps)

    fn = f + a1*k1 + a2*k2 + a3*k3 + a4*k4
    k1 = dt*F(f, t)
    k2 = dt*F(f + nu21*k1, t + nu21*dt)
    k3 = dt*F(f + nu31*k1 + nu32*k2, t + nu31*dt + nu32*dt)
    k4 = dt*F(f + nu41*k1 + nu42*k2 + nu43*k3, t + dt*(nu41 + nu42 + nu43))
    ...
    kn = dt*F(f + sum(nu_ni*ki, i=1,n-1), t + dt*sum(nu_ni, i=1,n-1))
    """
    res = np.zeros((x.size, steps))
    i = 0
    while i < steps:
        res[:,i] = x
        k1 = dt * F(x, t, *fargs, **fkwds)
        k2 = dt * F(x + k1/2., t + dt/2., *fargs, **fkwds)
        k3 = dt * F(x + k2/2., t + dt/2., *fargs, **fkwds)
        k4 = dt * F(x + k3, t + dt, *fargs, **fkwds)
        x = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6.
        t += dt
        i += 1
    return res

#  def F_rk(x, t, k):
#      return k * L @ x

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
