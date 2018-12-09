# -*- coding: utf-8 -*-
"""
Euler-Maruyama method on linear SDE

dX = λ*X dt + μ*X dW, X(0) = X0

Brownian path uses dt = 2**(-8)
Euler-Maruyama timestep uses Dt = R * dt
"""

import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

np.random.seed(100)

# problem parameters
λ = 2; μ = 1; X0 = 1
T = 1; N = 2**10; dt = T / N
t = np.arange(0, T+dt, dt)

# Brownian path
dW = np.sqrt(dt) * randn(N)
W = np.hstack(([0], np.cumsum(dW)))

# true solution for comparison
Xtrue = X0 * np.exp((λ-.5*μ**2)*t + μ*W)

# factor and new timestep for em integration
R = 10; Dt = dt*R

# numerical solution
Xem = np.zeros(N//R+1)
Xem[0] = X0 
for i in range(1, N//R+1):
    Winc = np.sum(dW[R*(i-1):R*i])
    Xem[i] = Xem[i-1] + λ*Xem[i-1]*Dt + μ*Xem[i-1]*Winc


plt.plot(t, Xtrue, 'b-')
plt.plot(t[::R], Xem, 'r:')
plt.xlabel('t')
plt.ylabel('X(t)')

print('err: %f' % np.abs(Xtrue[-1]-Xem[-1]))
