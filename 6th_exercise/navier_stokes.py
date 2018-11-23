# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import (rfft2,
                       irfft2,
                       fftshift,
                       ifftshift)
import matplotlib.pyplot as plt
from matplotlib import animation

#  TODO, class PDE
#  
#  __init__: why the choice of kx, ky?
#  rhs: how does the aliasing work?
#  
#  general: why the use of fftshift?


class PDE:
    def __init__(self, params, iv, dt, t0=0.):
        self.__dict__.update(params)

        # k-space
        kx = np.linspace(0, self.Nx//2, self.Nx//2+1)
        ky = np.linspace(-self.Ny//2, self.Ny//2-1, self.Ny)
        self.KX, self.KY = np.meshgrid(kx, ky)
        self.K_sq = np.square(self.KX) + np.square(self.KY)
        self.K_sq[self.Nx//2,0] = 1

        # x-space
        self.dx = (self.xe - self.xb) / self.Nx
        self.x = np.arange(self.xb, self.xe, self.dx)
        self.dy = (self.ye - self.yb) / self.Ny
        self.y = np.arange(self.yb, self.ye, self.dy)
        X, Y = np.meshgrid(self.x, self.y)

        # initial value and its fourier transform
        self.ω = iv(X, Y)
        self.ω_hat = fftshift(rfft2(self.ω), axes=0)
        #  self.t = t0
        self.t = [t0]
        self.dt = dt

        #  self.cfl = 1
        self.cfl = []

        self.scheme = self.shu_osher


    def time_step(self, steps=1):
        # calculate timesteps
        # solution is computed in fourier space, inverse transformation is
        # done for plotting, `steps` gives plotting frequency
        for _ in range(steps):
            self.scheme()
        self.ω = irfft2(ifftshift(self.ω_hat, axes=0))
        ux = irfft2(ifftshift( 1j * self.KY * self.ω_hat / self.K_sq ,
                    axes=0))
        uy = irfft2(ifftshift(-1j * self.KX * self.ω_hat / self.K_sq ,
                    axes=0))

        # check cfl condition, should be < 1
        _cfl = lambda u, dx: np.max(np.abs(u)) * self.dt / dx
        #  self.cfl = np.max(self.u) * self.dt / self.dx
        self.cfl.append(max(_cfl(ux, self.dx), _cfl(uy, self.dy)))

        self.t.append(self.t[-1] + self.dt * steps)
        return self.ω


    def rhs(self, ω_hat):
        # aliasing
        ω_hat[np.where(self.K_sq > self.Nx * self.Ny / 9.)] = 0

        ω = irfft2(ifftshift(ω_hat, axes=0))
        ux = irfft2(ifftshift( 1j * self.KY * self.ω_hat / self.K_sq ,
                    axes=0))
        uy = irfft2(ifftshift(-1j * self.KX * self.ω_hat / self.K_sq ,
                    axes=0))
        tmp = 1j * self.KX * fftshift(rfft2(ux * ω), axes=0) \
            + 1j * self.KY * fftshift(rfft2(uy * ω), axes=0)

        return ω_hat - tmp * self.dt


    def prop(self, delta=1.):
        # propagator
        return np.exp(-self.kappa * self.K_sq * self.dt * delta)

    def euler(self):
        # general method
        self.ω_hat = self.rhs(self.ω_hat) * self.prop()

    def heun(self):
        # rk2
        ω_1 = self.rhs(self.ω_hat) * self.prop()
        self.ω_hat = .5 * (self.rhs(ω_1) + self.ω_hat * self.prop())

    def shu_osher(self):
        # rk3
        ω_1 = self.rhs(self.ω_hat) * self.prop()
        ω_2 = .25 * (self.rhs(ω_1) * self.prop(-.5)
                     + 3. * self.ω_hat * self.prop(.5))
        self.ω_hat = 1./3. * (2. * self.rhs(ω_2) * self.prop(.5)
                              + self.ω_hat * self.prop())


def four_vortices(x, y,
                  x0=np.pi-2, x1=np.pi+2,
                  y0=np.pi-.5, y1=np.pi+.5, 
                  s_sq=4.):
    def _gaussian(xi, yj):
        return np.exp(-s_sq * ((x-xi)**2 + (y-yj)**2))
    return _gaussian(x0, y1) + _gaussian(x1, y0) \
            - _gaussian(x0, y0) - _gaussian(x1, y1)

params = dict(xb=0, xe=2*np.pi,
              yb=0, ye=2*np.pi,
              Nx=64, Ny=64,
              kappa=.0001)
p = PDE(params, four_vortices, dt=.05)

steps = 10
time = 0
tmax = 100
frames = int(tmax / (steps * p.dt))

res = np.array([p.time_step(steps) for _ in range(frames)])


fig = plt.figure()
ax = fig.add_subplot(111)

im = ax.imshow(p.ω, animated=True)

def step(i):
    im.set_array(res[i])
    print('time = %.2f, cfl = %.2f\r' % (p.t[i], p.cfl[i]), end='')
    return im,

def start_anim(fig=fig):
    return animation.FuncAnimation(fig, step, frames=frames, interval=10,
                                   blit=True, repeat=False)

#  anim = start_anim()
