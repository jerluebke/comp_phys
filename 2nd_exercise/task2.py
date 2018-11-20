# -*- coding: utf-8 -*-

# TODO: Leapfrog

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def euler(f, G, t0, tmax, dt, *gargs, **gkwds):
    """
    Forward Euler method for integrating 1st order ODE systems

    Parameters
    ==========
        f   :   scalar or array, initial value f = f(t0)
        G   :   callable, df = G(f, t)
        t0  :   starting time
        tmax:   end time
        dt  :   time step

    Yields
    ======
        fn  :   n-th time step of f, f_n = f_n-1 + dt * G(f_n-1, t)

    """
    t = t0
    while t < tmax:
        fn = f + dt * G(f, t, *gargs, **gkwds)
        f = fn
        t += dt
        yield fn

def verlet(x0, x1, F, t0, tmax, dt, *fargs, **fkwds):
    """
    Verlet method for explicitly integrating 2nd order ODE systems

    Parameters
    =========
        x0  :   first initial value, x0 = x(t0)
        x1  :   second initial value, x1 = x(t0+dt)
        F   :   callable, ddx = F(x, t), "Force"
        t0  :   starting time
        tmax:   end time
        dt  :   time step

    Yields
    =====
        xn  :   n-th time step of x

    """
    t = t0
    while t < tmax:
        xn = 2*x1 - x0 + dt**2 * F(x1, t, *fargs, **fkwds)
        x0, x1 = x1, xn
        t += dt
        yield xn

def leapfrog_1(x, v, F, t, tmax, dt, *fargs, **fkwds):
    """
    x_h = x_n + dt/2 * v_n
    v_n+1 = v_n + dt * a_h, a_h = F(x_h)
    x_n+1 = x_h + dt/2 * v_n+1
    """
    while t < tmax:
        x_h = x + dt/2. * v
        v = v + dt * F(x_h, t, *fargs, **fkwds)
        x = x_h + dt/2. * v
        t += dt
        yield x

def leapfrog_2(x, v, F, t, tmax, dt, *fargs, **fkwds):
    """
    v_h = v_n + dt/2 * a_n
    x_n+1 = x_n + v_h * dt
    v_n+1 = v_h + dt/2 * a_n+1
    """
    a = F(x, t, *fargs, **fkwds)
    while t < tmax:
        v_h = v + dt/2. * a
        x = x + v_h * dt
        a = F(x, t, *fargs, **fkwds)
        v = v_h + dt/2. * a
        t += dt
        yield x

def leapfrog_3(x, v, F, t, tmax, dt, *fargs, **fkwds):
    """
    x_n+1 = x_n + v_n * dt + dt**2/2 * a_n
    v_n+1 = v_n + dt/2 * (a_n + a_n+1)
    """
    a = F(x, t, *fargs, **fkwds)
    while t < tmax:
        x = x + v * dt + dt**2/2. * a
        a_n = F(x, t, *fargs, **fkwds)
        v = v + dt/2. * (a + a_n)
        a = a_n
        t += dt
        yield x


def G(f, t, *args, **kwds):
    """
    f = [r, v]
    return [dr, dv]
    """
    r = f[0]
    v = f[1]
    dr = v
    dv = -kwds.get('μ', 1) * r / np.sum(r**2)**(3/2)
    return np.array([dr, dv])

def F(x, t, *args, **kwds):
    """
    F(x) = ddx
    """
    return -kwds.get('μ', 1) * x / np.sum(x**2)**(3/2)


def euler_sample():
    f0 = np.array([[1, 1], [1, 0]])
    res_gen = euler(f0, G, 0, 7, .05, μ=2)
    res = np.array([f for f in res_gen])
    rsq = np.sum(res[:,0]**2, axis=1)
    vsq = np.sum(res[:,1]**2, axis=1)

    fig, axs = plt.subplots(2, 2)
    axs = axs.ravel()
    axs[0].plot(rsq)
    axs[0].set_title('distance squared')
    axs[1].plot(vsq)
    axs[1].set_title('speed squared')
    axs[2].plot(res[:,0,0], res[:,0,1])
    axs[2].set_title('radius')
    axs[3].plot(res[:,1,0], res[:,1,1])
    axs[3].set_title('velocity')


def verlet_sample():
    x0 = np.array([1., 1.])
    x1 = np.array([1.05, 1.])
    res_gen = verlet(x0, x1, F, 0, 7, .05, μ=2)
    res = np.array([x for x in res_gen])
    xsq = np.sum(res**2, axis=1)

    fig, axs = plt.subplots(1, 2)
    axs[0].plot(xsq)
    axs[0].set_title('distance squared')
    axs[1].plot(res[:,0], res[:,1])
    axs[1].set_title('radius')


def comparison_verlet_lf():
    x0 = np.array([1., 1.])
    v0 = np.array([1., 0.])
    x1 = np.array([1.05, 1.])   # result of one euler step with x0, v0, dt = 0.05
    lf = np.array([xn for xn in leapfrog_2(x0, v0, F, 0, 7, .05, μ=2)])
    verl = np.array([xn for xn in verlet(x0, x1, F, 0, 7, .05, μ=2)])
    fig, ax = plt.subplots()
    ax.set_title('comparison of trajectories - verlet vs leapfrog, dt=0.05')
    ax.plot(lf[:,0], lf[:,1], label='leapfrog')
    ax.plot(verl[:,0], verl[:,1], label='verlet')
    ax.legend()


def anim(integrator, iv, steps=1000, dt=.05):
    """
    iv : initial value
        verlet iv: x1 = [1.05, 1.]
        leapfrpg iv: v0 = [1., 0.]
    """
    x0 = np.array([1., 1.])
    res_gen = integrator(x0, iv, F, 0, steps, dt, μ=2)

    fa = plt.figure()
    fa.suptitle('animation')
    title = 'integrator: %s, steps: %d, dt: %04.2f' % (integrator.__name__.replace('_', '-'),
                                                       steps, dt)
    ax = fa.add_subplot(111, xlim=(-2, 2), ylim=(-2, 2),
                        title=title)

    c, = ax.plot([0], [0], 'ro', ms=6)
    s, = ax.plot([], [], 'ro', ms=6, zorder=2)
    t, = ax.plot([], [], 'b-', zorder=1)
    traj = np.zeros((steps, 2))

    def step(i):
        traj[i] = next(res_gen)
        t.set_data([traj[:i,0], traj[:i,1]])
        s.set_data(traj[i,0], traj[i,1])
        return c, t, s,

    def init():
        c.set_data([0], [0])
        t.set_data([], [])
        s.set_data([], [])
        return c, t, s,

    anim = animation.FuncAnimation(fa, step, frames=steps, interval=20,
                                   blit=True, init_func=init)

    return fa, anim


fv, av = anim(integrator=verlet, iv=np.array([1.05, 1.]), steps=10000)
fl, al = anim(integrator=leapfrog_1, iv=np.array([1., 0.]), steps=5000)


#  some comments:
#      on my machine verlet is about 10ms faster than leapfrog (260ms vs 270ms)
#      there is no noticable difference between lf_1 and lf_2, lf_2 runs with 295ms
#
#      leapfrog makes twice the amount of timesteps compared to verlet per
#      iteration (makes sense considering the alogrithms themselves...)
#
#      trajectories of leapfrog and verlet show small deviation from one another
#      the lf algos conform with each other, but it the different step sequence
#      between lf_1 and the others is perceptible when looking closely
#
#  timeit code (ipython):
#      >>> %timeit np.array([xn for xn in integrator(x, iv, F, 0, 1000, .05, μ=2)])

