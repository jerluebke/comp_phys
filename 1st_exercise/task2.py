# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('..')
from stencil_diag import stdiag

A, B = 0, 2*np.pi
N = 1000
x = np.linspace(A, B, N)
H = 2*np.pi / N

DD = stdiag(N, [1, -2, 1])
F = np.sin(x)**3
DF = 3*np.sin(x)**2*np.cos(x)

plt.plot(x, DF, label='analytical')
# weird result ...
plt.plot(x, DD@F / H**2, label='discrete operator')
plt.legend()
