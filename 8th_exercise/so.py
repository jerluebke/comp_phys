# -*- coding: utf-8 -*-
"""
harmonic oscillator with stochastic noise
    ddX + ω**2 X dt = μ dW

as a 2nd order system:
    dX = V dt
    dV = μ dW - ω**2 * X dt

note: euler's method causes the system to 'gain' energy (growing amplitudes)
"""
#  TODO:
#      check/discuss results
#      compare with undisturbed and exact analytic solution

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(100)

X0 = 1; V0 = 1;
ω = 10; μ = .2;
T = 1; N = 2**10; dt = T / N;
t = np.arange(0, T+dt, dt)
dW = np.random.randn(N)

X = np.zeros(N+1)
V = np.zeros(N+1)
X[0] = X0
V[0] = V0
for i in range(1, N+1):
    X[i] = X[i-1] + V[i-1] * dt
    V[i] = V[i-1] - ω**2 * X[i-1] * dt + μ * dW[i-1] 


plt.plot(t, X, 'r-', label='X(t)')
plt.plot(t, V, 'g--', label='V(t)')
plt.xlabel('t')
plt.ylabel('Amplitude')
plt.legend()
