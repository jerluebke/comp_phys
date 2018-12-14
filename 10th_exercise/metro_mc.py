# -*- coding: utf-8 -*-
"""
Demonstration of Metropolis-Monte-Carlo algorithm
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(100)

def p(x, σ):
    return np.exp(-(x / σ)**2)

# parameters and initial state
σ = .04
x = .6
px = p(x, σ)
# radius
rd = .5
# sum
s = 0
# steps
N = 10000

# results
xp = np.zeros(N)
pp = np.zeros(N)
accepted = 0

for i in range(N):
    trial = x + rd*(np.random.rand()-.5)
    pt = p(trial, σ)
    r = pt / px
    # accept with probability 1 if pt > px
    # or with probability pt / px
    if r >= 1 or r > np.random.rand():
        px = pt
        x = trial
        accepted += 1
    s += x**2
    xp[i] = x
    pp[i] = px

res = s / N
exact = σ**2 / 2
ratio = accepted / N

print("numerically  : %f\n"
      "exact        : %f\n"
      "ratio        : %f\n"
      % (res, exact, ratio))


xe = np.linspace(-1, 1, 100)
pe = p(xe, σ)

fig, (a, alog) = plt.subplots(1, 2)

a.plot(xe, pe, 'r-', label='exact')
a.plot(xp, pp, 'bo', mfc='none', label='numerically')
a.set_xlabel('x')
a.set_ylabel('p')
alog.semilogy(xe, pe, 'r-', label='exact')
alog.semilogy(xp, pp, 'bo', mfc='none', label='numerically')
alog.set_xlabel('x')
alog.set_ylabel('p')
a.legend()
alog.legend()
