# -*- coding: utf-8 -*-

import math
import numpy as np

def direct_ft(x):
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

    x.size needs to be power of 2
    """
    N = x.size
    if N <= 2:
        return direct_ft(x)
    else:
        X_even = recursive_ft(x[::2])
        X_odd = recursive_ft(x[1::2])
        W = np.exp(-2j * np.pi * np.arange(N//2) / N)
        # X_i = E_i + W_i * O_i
        # X_i+N/2 = E_i - W_i * O_i
        return np.concatenate([X_even + W * X_odd,
                               X_even - W * X_odd])

def bitrev_ft(x):
    return _bitrev_ft(x.astype(complex))

def _bitrev_ft(x):
    """
    discret fast fourier transform

    x needs to be complex (e.g.: x=x.astype(complex))
       and x.size a power of 2

    iterative and in-place implementation of the recursive approach:
        the elements are reordered and grouped into "recursion layers",
        mulitplied with the corresponding factors and summed over
    """
    def reverse(x, s):
        """
        bit reversal of input x with wordlength s
        """
        y = 0
        for i in range(s):
            y = (y << 1) | (x & 1)
            x >>= 1
        return y

    n           = x.size
    layers      = int(math.log(n, 2))
    x           = x[[reverse(i, layers) for i in range(n)]]     # bit reversed
    exptable    = np.exp([-2j * np.pi * i / n for i in range(n // 2)])

    # outer loop
    # goes through log(n, 2) "recursion" layers
    # on the l-th layer (l=log(s,2)) the fft is divided into (n-1) / s
    #  segments, where s is the size of each segment on the given layer
    s = 2
    while s <= n:
        half_s      = s // 2
        tablestep   = n // s
        # middle loop
        # goes through the fft's segments at layer log(s, 2)
        for i in range(0, n, s):
            # inner loop
            # goes through the segment's s elements and performs the
            #  "butterfly" operation
            k = 0
            for j in range(i, i+half_s):
                tmp             = x[j + half_s] * exptable[k]
                # x_right = x_left - x_right * exp(-i*pi*element / (2**layer))
                x[j + half_s]   = x[j] - tmp
                # x_left = x_left + x_right * exp(-i*pi*element / (2**layer))
                x[j]    += tmp
                k       += tablestep
        s *= 2
    return x
