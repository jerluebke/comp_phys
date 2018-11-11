# -*- coding: utf-8 -*-

from cg import *


N = 1000
y = np.zeros(N)
x = np.linspace(-5, 5, N)
b = np.sin(x) * np.exp(-x**2)
L = sparse.diags([1, -2, 1], [-1, 0, 1], (N, N), format='csc')

res, i, err = cg(L, b, y, N+1)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(x, b, label='b')
ax1.plot(x, L@res, label='result')
ax1.legend()

ax2.semilogy(np.arange(i), err[:i], label='error')
ax2.legend()
