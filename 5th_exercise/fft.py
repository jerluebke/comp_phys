# -*- coding: utf-8 -*-

import numpy as np

def dft(x):
    """
    direct discrete FT
    """
    N = x.size
    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)

def recursive_ft(x):
    """
    recursive discrete FT
    calls dft on lowest level
    """
    N = x.size
    if N <= 2:
        return dft(x)
    else:
        X_even = recursive_ft(x[::2])
        X_odd = recursive_ft(x[1::2])
        W = np.exp(-2j * np.pi * np.arange(N//2) / N)
        return np.concatenate([X_even + W * X_odd,
                               X_even - W * X_odd])
