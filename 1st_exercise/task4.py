# -*- coding: utf-8 -*-


def bisek(f, a, b, tol=1e-19):
    if abs(f(a)) <  tol or abs(f(b)) < tol:
        return a, b
    c = (a + b) / 2
    if (f(a) > 0  and f(c) > 0) or (f(a) < 0 and f(c) < 0):
        a = c
    else:
        b = c
    #  print(c)
    return bisek(f, a, b, tol)
