# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import rfft2, irfft2, fftshift, ifftshift

def rhs(p, ohat):
    ohat[np.where(p.K_sq > p.Nx*p.Ny/9.)] = 0
    o = irfft2(ifftshift(ohat, axes=0))

    uhat = 1j * p.KY * p.ω_hat / p.K_sq
    u = irfft2(ifftshift(uhat, axes=0))
    u *= o
    uhat = fftshift(rfft2(u), axes=0)
    
    res = ohat - 1j * p.KX * uhat * p.dt
    
    uhat = -1j * p.KX * p.ω_hat / p.K_sq
    u = irfft2(ifftshift(uhat, axes=0))
    u *= o
    uhat = fftshift(rfft2(u), axes=0)
    
    res -= 1j * p.KY * uhat * p.dt

    return res
    #  return o, uhat
    #  return o, u

def shu_osher(p):
    # rk3
    ω_1 = p.rhs(p.ω_hat) * p.prop()
    ω_2 = .25 * (p.rhs(ω_1) * p.prop(-.5)
                 + 3. * p.ω_hat * p.prop(.5))
    p.ω_hat = 1./3. * (2. * p.rhs(ω_2) * p.prop(.5)
                          + p.ω_hat * p.prop())
    return ω_2
