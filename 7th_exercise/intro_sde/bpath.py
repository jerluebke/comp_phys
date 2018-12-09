# -*- coding: utf-8 -*-
"""
Brownian path simulation
"""

import numpy as np
# random number from the N(0, 1) distribution
from numpy.random import randn
import matplotlib.pyplot as plt


np.random.seed(100)
T  = 1
N  = 500
dt = T / N
dW = np.zeros(N)
W  = np.zeros(N)

dW = np.sqrt(dt) * randn(N)
W = np.cumsum(dW)

plt.plot(np.arange(0, T+dt, dt), np.hstack(([0], W)))
plt.xlabel('t')
plt.ylabel('W(t)')
