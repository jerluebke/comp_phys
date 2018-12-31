# -*- coding: utf-8 -*-

#  To compile cython code:
#      > python setup.py build_ext --inplace

import numpy as np
import matplotlib.pyplot as plt
from visualise import PyQuadtreeEnv

np.random.seed(100)

x = np.linspace(0, 1, 64)
# some random data
d = .4 * np.sin(10*x) + .05 * np.random.randn(64) + .5
d64 = (64 * d).astype(np.uint)
Z = np.zeros((64, 64))
for i in range(64):
    Z[i, d64[i]] = 1
indices = np.ascontiguousarray(np.transpose(Z.nonzero()), dtype=np.uint32)

qtenv = PyQuadtreeEnv(indices)
si = qtenv.get_sorted()

fig, ax = plt.subplots()
ax.grid(False)
ax.tick_params(axis='both', **dict(bottom = False, left = False,
                                   top = False, right = False,
                                   labelbottom = False, labelleft = False))
s = plt.scatter(x, d)
s.set_facecolor(np.array([
    1, 0, 0, 1,
    0, 1, 0, 1,
    0, 0, 1, 1,
    *[1, 0, 0, 1]*61
]).reshape((64, 4)))
ax.set(xlim=(0, 1), ylim=(0, 1))
ax.plot(x[si], d[si], c='gray')

while 1:
    try:
        r = qtenv.insert_next()
        for elem in r:
            x, y, a = elem
            h = plt.axhline(y+a/2, x, x+a)
            v = plt.axvline(x+a/2, y, y+a)
            ax.add_line(h)
            ax.add_line(v)
    except StopIteration:
        break

fig.tight_layout()
plt.savefig('out-plot.svg', dpi=300, format='svg')

print('done.')


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
