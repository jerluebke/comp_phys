# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse


def cg(A, b, x, imax=1000, tol=1e-7):
    i = 0
    err_arr = np.zeros(imax)
    r = b - A @ x
    p = r
    dn = np.dot(r, r)
    d0 = dn
    while i < imax:
        if dn < tol**2*d0:
            print('convergence after %d iterations' % i)
            break
        q = A @ p
        alpha = dn / np.dot(p, q)
        x = x + alpha*p
        r = r - alpha*q
        d_tmp = dn
        dn = np.dot(r, r)
        beta = dn / d_tmp
        p = r + beta*p
        err_arr[i] = dn / d0
        i += 1
    return x, i, err_arr


N = 1000
y = np.zeros(N)
x = np.linspace(-5, 5, N)
b = np.sin(x) * np.exp(-x**2)
L = sparse.diags([1, -2, 1], [-1, 0, 1], (N, N), format='csc')

res, i, err = cg(L, b, y, N+1)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x, b, label='b')
ax1.plot(x, L@res, label='result')
ax1.legend()

ax2.semilogy(np.arange(i), err[:i], label='error')
ax2.legend()
