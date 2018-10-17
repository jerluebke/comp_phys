# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

F = lambda x: np.sin(x)**3
F_a = lambda x: 3*np.sin(x)**2*np.cos(x)
X0 = 1
H = np.array([np.exp(-i) for i in range(100)])

def right_dif(f, x0, h):
    return (f(x0 + h) - f(x0)) / h

def left_dif(f, x0, h):
    return (f(x0) - f(x0 - h)) / h

def mid_dif(f, x0, h):
    return (f(x0 + h) - f(x0 - h)) / (2*h)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set(xscale='log', yscale='log', title='comparision')

ax.plot(H, np.abs(F_a(X0) - right_dif(F, X0, H)), label='right diff.')
ax.plot(H, np.abs(F_a(X0) - left_dif(F, X0, H)), label='left diff.')
ax.plot(H, np.abs(F_a(X0) - mid_dif(F, X0, H)), label='middle diff.')

plt.legend()
