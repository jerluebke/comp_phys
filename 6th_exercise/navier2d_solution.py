#!/usr/bin/env python

import matplotlib
#  matplotlib.use('TkAgg')

import sys, types, pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

class Parameter:
    def __init__(self, *args, **kwargs):

        if args and isinstance(args[0],Parameter):
            xb = args[0].xb
            xe = args[0].xe
            yb = args[0].yb
            ye = args[0].ye
            Nx = args[0].Nx
            Ny = args[0].Ny
            dt = args[0].dt
            kappa = args[0].kappa

        elif kwargs :
            xb = kwargs['xb']
            xe = kwargs['xe']
            yb = kwargs['yb']
            ye = kwargs['ye']
            Nx = kwargs['Nx']
            Ny = kwargs['Ny']
            dt = kwargs['dt']
            kappa = kwargs['kappa']

        self.xb = xb
        self.xe = xe
        self.yb = yb
        self.ye = ye
        self.Nx = Nx
        self.Ny = Ny
        self.dx = (self.xe-self.xb)/float(self.Nx)
        self.dy = (self.ye-self.yb)/float(self.Ny)
        self.dt = dt
        self.kappa = kappa
        self.cfl = 1

class PDE(Parameter):
    """
    solves: u_t = \kappa \Delta u
    """
    def __init__(self, p ,t0 = 0.0):
        Parameter.__init__(self,p)
        Nx = self.Nx
        Ny = self.Ny
        kx = np.linspace(0,Nx//2,Nx//2+1)
        ky = np.linspace(-Ny//2,Ny//2-1,Ny)
        self.KX,self.KY = np.meshgrid(kx, ky)
        self.K2 = np.square(self.KX)+np.square(self.KY)
        self.K2[Nx//2,0] = 1 #fix

        self.x = self.xb + self.dx * np.arange(self.Nx)
        self.y = self.yb + self.dy * np.arange(self.Ny)
        X,Y = np.meshgrid(self.x,self.y)
        self.om = np.sin(2*X)*np.sin(Y)+np.cos(2*Y)
        x0=np.pi-2
        x1=np.pi+2
        y0=np.pi-0.5
        y1=np.pi+0.5
        self.om = -np.exp(-4*((X-x0)**2 + (Y-y0)**2)) + np.exp(-4*((X-x0)**2 + (Y-y1)**2)) \
                +  np.exp(-4*((X-x1)**2 + (Y-y0)**2)) - np.exp(-4*((X-x1)**2 + (Y-y1)**2))
        self.omhat = np.fft.fftshift(np.fft.rfft2(self.om),0)
        print(len(self.omhat))
        self.t = t0

        self.scheme = self.heun

    def time_step(self, Nsteps = 1):
        for i in range(Nsteps):
            self.scheme()
        self.om = np.fft.irfft2(np.fft.ifftshift(self.omhat,0))
        ux = np.fft.irfft2(np.fft.ifftshift( 1j*self.KY/self.K2*self.omhat,0))
        uy = np.fft.irfft2(np.fft.ifftshift(-1j*self.KX/self.K2*self.omhat,0))
        self.cfl = max(np.max(np.abs(ux))*self.dt/self.dx,np.max(np.abs(uy))*self.dt/self.dy)
        #  print("cfl =", self.cfl)
        self.t += self.dt * Nsteps
        return self.om

    def rhs(self, omhat):
        # aliasing
        omhat[np.where(self.K2>self.Nx*self.Ny/9.)]=0

        om = np.fft.irfft2(np.fft.ifftshift(omhat,0))
        ux = np.fft.irfft2(np.fft.ifftshift( 1j*self.KY/self.K2*omhat,0))
        uy = np.fft.irfft2(np.fft.ifftshift(-1j*self.KX/self.K2*omhat,0))
        tmp = 1j*self.KX*np.fft.fftshift(np.fft.rfft2(ux*om),0) \
            + 1j*self.KY*np.fft.fftshift(np.fft.rfft2(uy*om),0)

        return  omhat - tmp*self.dt

    def prop(self, delta):
        return np.exp(-self.kappa*self.dt*self.K2*delta)

    def euler(self):
        self.omhat = self.rhs(self.omhat)*self.prop(1.)

    def heun(self):
        omone = self.rhs(self.omhat)*self.prop(1.)
        self.omhat = 0.5*(self.omhat*self.prop(1.) + self.rhs(omone))

    def shuOsher(self):
        omone =  self.rhs(self.omhat)*self.prop(1.)
        omtwo = 3./4.*self.omhat*self.prop(0.5) + 1./4.*self.rhs(omone)*self.prop(-0.5)
        self.omhat = 1./3.*self.omhat*self.prop(1.) + 2./3.*self.rhs(omtwo)*self.prop(0.5)

p = Parameter(xb=0,xe=2*np.pi,
              yb=0,ye=2*np.pi,
              Nx=64,Ny=64,dt=0.05,kappa=0.0001)

N_steps = 10
t_max = 100
frames = int(t_max / float(N_steps * p.dt))
frames = 500

pde = PDE(p)

######################################################################
# Set up plot
fig = plt.figure()
#time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
#cfl_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)

######################################################################
# Animate plot
im=plt.imshow(pde.om,animated=True)

res = np.array([pde.time_step(N_steps) for _ in range(frames)])

def integrate(i):
    #  pde.time_step(N_steps)
    im.set_array(res[i])
#    time_text.set_text('time = %.2f' % pde.t)
#    cfl_text.set_text('cfl = %.2f' % pde.cfl)
    return im,

anim = animation.FuncAnimation(fig, integrate, frames=frames,
                               interval=100, blit=True,repeat=False)

plt.show()
