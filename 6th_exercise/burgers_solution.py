# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

class Parameter:
    def __init__(self, *args, **kwargs):

        if args and isinstance(args[0],Parameter):
            xb = args[0].xb
            xe = args[0].xe
            N = args[0].N
            dt = args[0].dt
            kappa = args[0].kappa

        elif kwargs :
            xb = kwargs['xb']
            xe = kwargs['xe']
            N = kwargs['N']
            dt = kwargs['dt']
            kappa = kwargs['kappa']

        self.xb = xb
        self.xe = xe
        self.N = N
        self.dx = (self.xe-self.xb)/float(self.N)
        self.dt = dt
        self.kappa = kappa
        self.cfl=1

class PDE(Parameter):
    """
    solves: u_t = \kappa \Delta u
    """
    def __init__(self, p ,t0 = 0.0):
        Parameter.__init__(self,p)
        self.kx = np.linspace(0,self.N/2,self.N/2+1)*2*np.pi/(self.xe-self.xb)
        self.kx2 = np.square(self.kx)

        self.x = self.xb + self.dx * np.arange(self.N)
        self.u = np.sin(self.x)
        self.uhat = np.fft.rfft(self.u)
        print(len(self.uhat))
        #self.u = np.piecewise(self.x, [np.abs(self.x-np.pi) < 0.01], [1])
        #self.u = np.sign(np.cos(self.x))
        self.t = t0

        self.scheme = self.heun

        # attributes used for dynamic plotting
        #  self.x_line = None
        #  self.u_line = None

    def time_step(self, Nsteps = 1):
        for i in range(Nsteps):
            self.scheme()
        self.u = np.fft.irfft(self.uhat)
        #print self.t, np.max(-np.fft.irfft(self.kx*1j*self.uhat))
        self.cfl = np.max(self.u)*self.dt/self.dx
        #print "cfl =",self.cfl
        self.t += self.dt * Nsteps

    def rhs(self, uhat):
        # aliasing
        uhat[int(2.*len(self.kx)/3.):len(self.kx)] = 0
        return  uhat - 0.5 * self.kx * 1j * np.fft.rfft(np.fft.irfft(uhat)**2)*self.dt

    def prop(self, delta):
        return np.exp(-self.kappa*self.dt*self.kx2*delta)

    def euler(self):
        self.uhat = self.rhs(self.uhat)*self.prop(1.)

    def heun(self):
        uone = self.rhs(self.uhat)*self.prop(1.)
        self.uhat = 0.5*(self.uhat*self.prop(1.) + self.rhs(uone))

    def shuOsher(self):
        uone =  self.rhs(self.uhat)*self.prop(1.)
        utwo = 3./4.*self.uhat*self.prop(0.5) + 1./4.*self.rhs(uone)*self.prop(-0.5)
        self.uhat = 1./3.*self.uhat*self.prop(1.) + 2./3.*self.rhs(utwo)*self.prop(0.5)

p = Parameter(xb=0,xe=2*np.pi,N=256,dt=0.01,kappa=0.02)

N_steps = 10
t_max = 100
frames = int(t_max / float(N_steps * p.dt))
frames = 500

pde = PDE(p)

######################################################################
# Set up plot
fig = plt.figure()
ax = plt.axes(xlim=(p.xb,p.xe), ylim=(-1.2, 1.2))
u_line, = ax.plot([], [],lw=2,marker='.',markerfacecolor='red', markersize=12)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
cfl_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
#ax.legend(prop=dict(size=14))
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')

######################################################################
# Animate plot
def init():
    u_line.set_data([], [])

    time_text.set_text('')
    cfl_text.set_text('')
    return (u_line, time_text,cfl_text)

def integrate(i):
    pde.time_step(N_steps)
    u_line.set_data(pde.x, pde.u)

    time_text.set_text('time = %.2f' % pde.t)
    cfl_text.set_text('cfl = %.2f' % pde.cfl)
    return (u_line, time_text,cfl_text)

anim = animation.FuncAnimation(fig, integrate, init_func=init, frames=frames,
                               interval=100, blit=True,repeat=False)

plt.show()
