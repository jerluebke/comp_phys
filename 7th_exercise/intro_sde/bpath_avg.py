# -*- coding: utf-8 -*-
"""
function along Browninan path
compute mean of 1000 paths, show 5 individual paths for comparison
"""

import numpy as np
from numpy.random import randn
from numpy.matlib import repmat
import matplotlib.pyplot as plt

np.random.seed(100)
T = 1; N = 500; dt = T / N
t = np.arange(0, 1+dt, dt)

# number of Brownian paths
M = 1000
dW = np.sqrt(dt) * randn(M, N)
W = np.cumsum(dW, 1)
# repmat: replicate t along the rows M times 
U = np.exp(repmat(t[1:], M, 1) + 0.5 * W)
Umean = np.mean(U, 0)

# plot 5 sample paths
# plot(t, [1, U])
plt.plot(repmat(t, 5, 1).T, np.hstack((np.ones((5, 1)), U[0:5,:])).T, 'r:')
# plot of expectation value for comparison
#  plt.plot(t, np.exp(9*t/8), 'g-')
# plot average of all paths
plt.plot(t, np.hstack(([1], Umean)), 'b-')
plt.xlabel('t')
plt.ylabel('U(t)')

print("maximal deviation: %f" % np.max(np.abs(Umean-np.exp(9*t[1:]/8))))
