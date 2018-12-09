# -*- coding: utf-8 -*-
"""
Milstein's method on linear SDE

dX = r*X*(K-X) dt + β*X dW, X(0) = X0
used to describe population dynamics
"""

import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

np.random.seed(100)

# problem parameters
r = 4; K = 1; β = .25; X0 = .5
T = 1; N = 2**10; dt = T / N
t = np.arange(0, T+dt, dt)

# Brownian path
dW = np.sqrt(dt) * randn(N)
W = np.hstack(([0], np.cumsum(dW)))

# numerical solution
Xem = np.zeros(N+1)
Xem[0] = X0 
for i in range(1, N+1):
    Xem[i] = Xem[i-1] + r*Xem[i-1]*(K-Xem[i-1])*dt + β*Xem[i-1]*dW[i-1] \
            + .5*β**2*Xem[i-1] * (dW[i-1]**2 - dt)


plt.plot(t, Xem, 'r-')
#  plt.plot(t, Xem)
plt.xlabel('t')
plt.ylabel('X(t)')

