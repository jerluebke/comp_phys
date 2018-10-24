# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def rectangle(f, a, b, N):
    h = (b-a) / N
    d = np.arange(a, b, h) + h/2
    return h * np.sum(f(d))

def trapez(f, a, b, N):
    h = (b-a) / N
    y = f(np.arange(a, b+h, h))
    y[0] /= 2
    y[-1] /= 2
    return h * np.sum(y)

def simpson(f, a, b, N):
    h = (b-a) / N
    d = np.arange(a, b+h, h)
    ar = f(d[0]) + f(d[-1])
    au = np.sum(f(d[1:-1:2]))
    ag = np.sum(f(d[2:-2:2]))
    return h * (ar + 4*au + 2*ag) / 3

test_functions = {
    #  'x^3 - x^2' : (lambda x: x**3-x**2,
    #                 lambda x: x**4/4-x**3/3),
    #  'sin^2 x'   : (lambda x: np.sin(x)**2,
    #                 lambda x: x/2-np.sin(x)*np.cos(x)/2),
    #  # divide by zero in [0, 1] (log(0), duh...)
    #  #  'x * log(x)': (lambda x: x*np.log(x),
    #  #                 lambda x: x**2*np.log(x)/2-x**2/4),
    #  'x * e^{-x}': (lambda x: x*np.exp(-x),
    #                 lambda x: -(x+1)*np.exp(-x)),
    'fn' : (lambda x: np.exp(-x)*np.sin(x),
            lambda x: -np.exp(-x)*np.sin(x)/2-np.exp(-x)*np.cos(x)/2)
}

#  integrations = [
#      rectangle,
#      trapez,
#      simpson
#  ]
#  
#  for n, (f, F) in test_functions.items():
#      print(n)
#      for i in integrations:
#          print('\t%f' % i(f, 0.1, 1, 10))
#      print('analytic: %f' % (F(1)-F(0)), '\n\n')

A, B = 0., 10.
M = 25 
Ns = np.arange(2, 2*M+2, 2)
L = 3

res = {
    'rectangle' : (rectangle, np.zeros((M, L)), np.zeros(M)),
    'trapez'    : (trapez, np.zeros((M, L)), np.zeros(M)),
    'simpson'   : (simpson, np.zeros((M, L)), np.zeros(M))
}

fig, ax = plt.subplots()
ax.set_yscale('log')

for ni, (integ, resarr, errarr) in res.items():
    for i, (nf, (f, F)) in enumerate(test_functions.items()):
        resarr[:,i] = np.array([integ(f, A, B, n) for n in Ns])
        resarr[:,i] = np.abs(resarr[:,i] - (F(B)-F(A)))
    #  errarr = np.mean(resarr, axis=1)
    ax.plot(Ns, resarr[:,0], label=ni)

plt.legend()
