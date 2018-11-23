# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import rfft, irfft
import matplotlib.pyplot as plt
from matplotlib import animation


class PDE:
    def __init__(self, params, iv, dt, t0=0.):
        self.__dict__.update(params)

        # k-space
        self.k = np.linspace(0, self.N//2, self.N//2+1) \
                * 2*np.pi / (self.xe - self.xb)
        self.k_sq = np.square(self.k)

        # x-space, initial value, fourier transform
        self.dx = (self.xe - self.xb) / self.N
        self.x = np.arange(self.xb, self.xe, self.dx)
        self.u = iv(self.x)
        self.û = rfft(self.u)
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
        self.u = irfft(self.û)

        # check cfl condition, should be < 1
        #  self.cfl = np.max(self.u) * self.dt / self.dx
        self.cfl.append(np.max(self.u) * self.dt / self.dx)

        self.t.append(self.t[-1] + self.dt * steps)
        return self.x, self.u

    def prop(self, delta=1.):
        # propagator
        return np.exp(-self.kappa * self.k_sq * self.dt * delta)

    def rhs(self, û):
        # low-pass with estimate for highest frequency
        û[2*self.k.size//3:] = 0
        return û - .5j * self.dt * self.k * rfft(irfft(û)**2)

    def euler(self):
        # general method
        self.û = self.rhs(self.û) * self.prop()

    def heun(self):
        # rk2
        û1 = self.rhs(self.û) * self.prop()
        self.û = .5 * (self.rhs(û1) + self.û * self.prop())

    def shu_osher(self):
        # rk3
        û1 = self.rhs(self.û) * self.prop()
        û2 = .25 * (self.rhs(û1) * self.prop(-.5) + 3. * self.û *
                     self.prop(.5))
        self.û = 1./3. * (2. * self.rhs(û2) * self.prop(.5) + self.û *
                           self.prop())


# PDE
#  f0 = lambda x: np.sin(x)**2
f0 = np.sin
params = dict(xb=0, xe=2*np.pi, N=256, kappa=.02)
params['kappa'] = 0
p = PDE(params, f0, dt=.01)

# params for plotting
steps = 10
tmax = 1000
frames = int(tmax // (steps * p.dt))

# compute complete solution for smoother plotting
res = np.array([p.time_step(steps) for _ in range(frames)])
#  plt.imshow(res[:500,1,:].T)     # x: time, y: x, color: u


###########################################################
# SET UP PLOTTING
###########################################################
fig = plt.figure()
ax = fig.add_subplot(111, xlim=(p.xb, p.xe),
                     #  ylim=(res[:,1].min()-.1, res[:,1].max()+.1),
                     ylim=(-1.1, 1.1),
                     xlabel='$x$', ylabel='$u$')
line, = ax.plot([], [], lw=2, marker='.', mfc='r', ms=10)
#  ttext = ax.text(.02, .95, '', transform=ax.transAxes)   # time
#  ctext = ax.text(.02, .90, '', transform=ax.transAxes)   # cfl

###########################################################
# ANIMATION
###########################################################

def init():
    line.set_data([], [])
    return line,

def step(i):
    line.set_data(res[i,0], res[i,1])
    #  ttext.set_text('time = %.2f' % time)
    #  ctext.set_text('cfl = %.2f' % p.cfl[i])
    print('time = %.2f, cfl = %.2f\r' % (p.t[i], p.cfl[i]), end='')
    return line, # ttext, ctext,

def start_anim():
    return animation.FuncAnimation(fig, step, frames=frames, interval=10,
                                   blit=True, repeat=True, init_func=init)

#  anim = start_anim()
