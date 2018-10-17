# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


#################
# ALGORITHMS    #
#################

def bisek(f, a, b, tol=1e-12, n=0):
    """
    RECURSIVE

    Parameters
    ==========
        f   :   callable, function to find roots of
        a, b:   intervall, in which to search
        tol :   convergence tolerance
        n   :   number steps

    Returns
    =======
        a   :   left border of intervall (in case of convergence, a and b are
                close together)
        n   :   steps until result

    """
    #  if abs(f(a)) <  tol or abs(f(b)) < tol:
    if abs(a-b) < tol or abs(min(f(a), f(b))) < tol:
        return a, n
    c = (a + b) / 2
    if np.sign(f(a)) == np.sign(f(c)):
        return bisek(f, c, b, tol, n+1)
    else:
        return bisek(f, a, c, tol, n+1)


def sekant(f, xn, xm, tol=1e-12, n=0):
    """
    RECURSIVE

    Parameters
    ==========
        f   :   callable, function to find roots of
        xn  :   x_n (start with x_1)
        xm  :   x_n-1 (start with x_0)
        tol :   tolerance
        n   :   number of steps

    Returns
    =======
        xn  :   root
        n   :   steps until result

    """
    return (xn, n) if abs(xn-xm) < tol or abs(f(xn)) < tol \
            else sekant(f, xn - (xn - xm) / (f(xn) - f(xm)) * f(xn), xn, tol, n+1)


def regula_falsi(f, a, b, tol=1e-12, n=0):
    """
    RECURSIVE, not opimized

    Parameters
    ==========
        f   :   callable, function to find roots of
        a, b:   intervall in which to search
        tol :   tolerance
        n   :   number of steps

    Returns
    =======
        a   :   left border of interval (approx. of root)
        n   :   number of steps

    """
    if abs(a-b) < tol or abs(f(a)) < tol:
        return a, n
    c = (a*f(b) - b*f(a)) / (f(b) - f(a))
    if np.sign(f(a)) == np.sign(f(c)):
        return regula_falsi(f, c, b, tol, n+1)
    else:
        return regula_falsi(f, a, c, tol, n+1)


def newton(f, df, x, tol=1e-12, n=0):
    """
    RECURSIVE

    Parameters
    ==========
        f   :   callable, function
        df  :   callable, derivative of f
        x   :   initial value
        tol :   tolerance
        n   :   number of steps

    Returns
    =======
        x   :   root of f
        n   :   number of steps

    """
    xn = x - f(x) / df(x)
    return (xn, n+1) if abs(x-xn) < tol or abs(f(xn)) < tol \
            else newton(f, df, xn, tol, n+1)


#########
# TEST  #
#########

x = np.linspace(0, 4)
FUNCS = {
    r'$cos(x)$'             :   (np.cos, lambda x: -np.sin(x)),
    r'$-\frac{x^2}{2}+3$'   :   (lambda x: -x**2/2+3, lambda x: -x),
    r'$-x^2 e^{-x}+0.3$'    :   (lambda x: -x**2*np.exp(-x)+.3, lambda x: \
                                 (x**2-2*x)*np.exp(-x)),
    r'$sin(x)^3$'           :   (lambda x: np.sin(x)**3, lambda x: \
                                 3*np.sin(x)**2*np.cos(x)),
}

RESULTS = [
    lambda func: '\tbisek: \t\t %f, %d' % bisek(func, 0.5, 3.5),
    lambda func: '\tsekant: \t %f, %d' % sekant(func, 0.5, 3.5),
    lambda func: '\tregula falsi: \t %f, %d' % regula_falsi(func, 0.5, 3.5),
    lambda func, df: '\tnewton: \t %f, %d' % newton(func, df, 1),
]

def test():
    [plt.plot(x, f[0](x), label=n) for n, f in FUNCS.items()]
    plt.legend()

    print('roots, number of steps:')
    for name, funcs in FUNCS.items():
        func, df = funcs
        print(name)
        for res in RESULTS[:-1]:
            try:
                print(res(func))
            except RecursionError:
                print('\nNO CONVERGENCE!\n')
        try:
            print(RESULTS[-1](func, df))
        except RecursionError:
            print('\nNO CONVERGENCE!\n')
        print('\n')


test()
