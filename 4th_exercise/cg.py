# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse


def cg(A, b, x=None, imax=1000, tol=1e-7):
    i = 0
    err_arr = np.zeros(imax)
    if x is None:
        x = np.zeros_like(b)
    r = b - A @ x
    p = r
    dn = np.dot(r, r)
    d0 = dn
    while i < imax:
        #  if dn < tol**2*d0:
        #      print('convergence after %d iterations' % i)
        #      break
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
