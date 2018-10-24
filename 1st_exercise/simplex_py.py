# -*- coding: utf-8 -*-

import numpy as np

A = 1.
G = 2.
Q = .5
S = .5


def step(x, f, tol):
    if np.all(np.std(x, axis=0) < tol):
        raise StopIteration

    z = f(x)
    ix = np.argsort(z)
    z = z[ix]
    x = x[ix]

    m = np.sum(x[:2], axis=1)/2.
    r = m+A*(m-x[2])
    fr = f(r)
    if z[0] <= fr <= z[1]:
        x[2] = r
        return x
    if fr < z[0]:
        e = m+G*(r-m)
        if f(e) < fr:
            x[2] = e
        else:
            x[2] = r
        return x

    c = m+Q*(x[2]-m)
    if f(c) < z[2]:
        x[2] = c
        return x

    x[1:3] = x[0]+S*(x[1:3]-x[0])
    return x


def find(x, f, tol, maxiter):
    steps = 0
    while 1:
        steps += 1
        if np.all(np.std(x, axis=0) < tol):
            return 1, x
        elif steps > maxiter:
            return 0, x

        z = f(x)
        ix = np.argsort(z)
        z = z[ix]
        x = x[ix]

        m = np.sum(x[:2], axis=1)/2.
        r = m+A*(m-x[2])
        fr = f(r)
        if z[0] <= fr <= z[1]:
            x[2] = r
            continue
        if fr < z[0]:
            e = m+G*(r-m)
            if f(e) < fr:
                x[2] = e
            else:
                x[2] = r
            continue

        c = m+Q*(x[2]-m)
        if f(c) < z[2]:
            x[2] = c
            continue

        x[1:3] = x[0]+S*(x[1:3]-x[0])
        continue
