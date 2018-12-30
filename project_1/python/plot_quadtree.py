# -*- coding: utf-8 -*-

#  To compile cython code:
#      > python setup.py build_ext --inplace

import numpy as np
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

while 1:
    try:
        r = qtenv.insert_next()
        print(r, '\n')
    except StopIteration:
        break
print('done.')


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
