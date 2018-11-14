# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from timeit import timeit
from fft import direct_ft, recursive_ft, bitrev_ft

a = -np.pi
b = np.pi
N = 512
x = np.linspace(a, b, N)

def func(x, **kwds):
    # note peaks for each mode
    #  return np.sum([np.sin(i*x) for i in (1., 2., 5., 10., 50.)], axis=0)
    # note mirrored peaks
    #  return np.sum(600*x)
    #  return np.piecewise(x, [np.abs(x)<.1], [1.])
    s = kwds.get('s', .1)
    return np.exp(-x**2/(2*s**2)) / np.sqrt(2*np.pi*s**2)

f = func(x)

np.fft.fft.__name__ = 'fft_numpy'
dfts = [direct_ft, recursive_ft, bitrev_ft, np.fft.fft]
for dft_func in dfts:
    print(dft_func.__name__, ':\t',
          timeit('ft(x)', globals={'ft' : dft_func, 'x' : f}, number=1000),
          'ms')


fft = bitrev_ft(f)

fig, (ao, ar, ai) = plt.subplots(3)
fig.tight_layout()
fig.canvas.set_window_title('discrete fourier transform')

ao.set(title='original', xlabel='x', ylabel='A')
ao.plot(x, f)

ar.set(title='ft - real part', xlabel='k', ylabel='Â')
ar.plot(fft.real)

ai.set(title='ft - imaginary part', xlabel='k', ylabel='Â')
ai.plot(fft.imag)
