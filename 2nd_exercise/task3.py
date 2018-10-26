# -*- coding: utf-8 -*-

import numpy as np
#  import matplotlib.pyplot as plt

def rk4(x, t, F, dt, *fargs, **fkwds):
    """
    classical runge-kutta (4 steps)

    fn = f + a1*k1 + a2*k2 + a3*k3 + a4*k4
    k1 = dt*F(f, t)
    k2 = dt*F(f + nu21*k1, t + nu21*dt)
    k3 = dt*F(f + nu31*k1 + nu32*k2, t + nu31*dt + nu32*dt)
    k4 = dt*F(f + nu41*k1 + nu42*k2 + nu43*k3, t + dt*(nu41 + nu42 + nu43))
    ...
    kn = dt*F(f + sum(nu_ni*ki, i=1,n-1), t + dt*sum(nu_ni, i=1,n-1))
    """
    k1 = dt * F(x, t)
    k2 = dt * F(x + k1/2., t + dt/2.)
    k3 = dt * F(x + k2/2., t + dt/2.)
    k4 = dt * F(x + k3, t + dt)
    xn = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6.
