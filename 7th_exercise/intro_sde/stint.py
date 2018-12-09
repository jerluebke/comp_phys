# -*- coding: utf-8 -*-
"""
Ito and Stratonovich Stochastic Integrals of W dW
"""

import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

#  np.random.seed(100)
T = 1; N = 500; dt = T / N;

dW = np.sqrt(dt) * randn(N)
W = np.hstack(([0], np.cumsum(dW)))

ito = np.sum(W[:-1] * dW)
strat = np.sum((.5 * (W[:-1] + W[1:]) + .5 * np.sqrt(dt) * randn(N)) * dW)

itoerr = np.abs(ito - .5 * (W[-1]**2 - T))
straterr = np.abs(strat - .5 * W[-1]**2)

print("Stochastic Integrals (err)\n"
      "Ito: %f\t(%f)\n"
      "Stratonovich: %f\t(%f)\n"
      % (ito, itoerr, strat, straterr))
