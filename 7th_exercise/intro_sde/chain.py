# -*- coding: utf-8 -*-
"""
stochastic chain rule

    dX = (α-X) dt + β*sqrt(X) dW, X(0) = X0

this equation is integrated using the Euler-Maruyama method, first directly
and a second time using V = sqrt(X), which yields (after applying the
stochastic chain rule):
    dV = ((4*α-β**2)/(8*V) - V/2) dt + β/2 dW
"""

import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

np.random.seed(100)

α = 2; β = 1; X0 = 1
T = 1; N = 200; dt = T / N
t = np.arange(0, T+dt, dt)

dW = np.sqrt(dt)*randn(N)

# direct
X = np.zeros(N+1)
X[0] = X0
# chain rule
V = np.zeros(N+1)
V[0] = X0**2

for i in range(1, N+1):
    X[i] = X[i-1] + (α - X[i-1]) * dt \
            + β * np.sqrt(np.abs(X[i-1])) * dW[i-1]
    V[i] = V[i-1] + ((4*α-β**2) / (8*V[i-1]) - .5*V[i-1]) * dt \
            + .5*β*dW[i-1]


plt.plot(t, np.sqrt(X), 'b-', label='direct')
plt.plot(t, V, 'r:', label='chain rule')
plt.xlabel('t')
plt.ylabel('V(X)')
plt.legend()

print("err: %f" % np.max(np.abs(np.sqrt(X)-V)))
