# -*- coding: utf-8 -*-

#  To compile cython code:
#      > python setup.py build_ext --inplace

import numpy as np
import matplotlib.pyplot as plt
from visualise import PyQuadtreeEnv

FILENAMETEMPLATE = 'data/plots/out-%03d.png'

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
                                   top = False, right = False))
ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.])
ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.])
ax.set_xticklabels([0, 16, 32, 48, 64])
ax.set_yticklabels([0, 16, 32, 48, 64])
cl = ["gray"]*64
s = ax.scatter(x, d, c='gray', s=200, marker='1', lw=2, zorder=3)
ax.set_title("raw data")
ax.set(xlim=(0, 1), ylim=(0, 1))
plt.savefig(FILENAMETEMPLATE % 0 , dpi=300, format='png')
ax.plot(x[si], d[si], c='gray', lw=1, zorder=1)
ax.set_title("morton curve")
ax.set(xlim=(0, 1), ylim=(0, 1))
plt.savefig(FILENAMETEMPLATE % 1, dpi=300, format='png')

linekwargs = dict(c = 'black', lw = 1, zorder = 2)

for i in range(1, 64):
    try:
        cl[si[i]] = "red"
        cl[si[i-1]] = "green"
        s.set_facecolor(cl)
        ax.set_title("current key: 0x%04X" % qtenv.get_key(i))
        r = qtenv.insert_next()
        for elem in r:
            x, y, a = elem
            h = plt.axhline(y+a/2, x, x+a, **linekwargs)
            v = plt.axvline(x+a/2, y, y+a, **linekwargs)
            ax.add_line(h)
            ax.add_line(v)
        #  fig.draw()
        plt.savefig(FILENAMETEMPLATE % (i+1), dpi=300, format='png')
    except StopIteration:
        break

fig.tight_layout()
#  plt.savefig('out-plot.png', dpi=300, format='png')

print('done.')


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
