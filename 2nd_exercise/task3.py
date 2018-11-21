# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pyplot as plt
from matplotlib import animation

# implicite procedure - like leapfrog - is not suitable for a problem like this
# TODO: why?

def crank_nicolson(f, D, dt, steps):
    # solve: B @ fn = A @ f (for fn)
    # don't invert B! its not sparse anymore afterwards
    res = np.zeros((f.size, steps))
    res[:,0] = f
    # sparse.linalg.inv is more efficient with csc format
    I = sparse.identity(f.size, format='csc')
    A = I + dt * D
    B = I - dt * D
    B_LU = linalg.splu(B)
    for i in range(1, steps):
        res[:,i] = B_LU.solve(A.dot(res[:,i-1]))
    return res

def fwd_euler(f, D, dt, steps):
    res = np.zeros((f.size, steps))
    res[:,0] = f
    I = sparse.identity(f.size, format='csc')
    for i in range(1, steps):
        res[:,i] = (I + dt * D) @ res[:,i-1]
    return res

def leapfrog(f, f1, D, dt, steps):
    res = np.zeros((f.size, steps))
    res[:,0] = f
    res[:,1] = f1
    for i in range(2, steps):
        # TODO: why?
        res[:,i] = res[:,i-2] + 2. * dt * D @ res[:,i-1]
        #  res[:,i] = 2 * res[:,i-1] - res[:,i-2] + dt**2. * D @ res[:,i-1]
    return res


def test(s=1., N=100, steps=200):
    k = 1
    x = np.linspace(-6, 6, N)
    dx_sq = (x[1] - x[0])**2
    dt = dx_sq / (2. * k)   # cfl condition, nyquist wavelength
    L = sparse.diags([1., -2., 1.], [-1, 0, 1], shape=(N, N), format='csc') / dx_sq
    #  get_f0 = lambda s: np.exp(-x**2/(2*s**2)) / (s*np.sqrt(2*np.pi))
    #  f0 = get_f0(s)
    f0 = np.piecewise(x, [np.abs(x) < 1.], [1.])

    fwd_euler_sol = fwd_euler(f0, k*L, dt/2., steps)

    fig, ax = plt.subplots(3, 1, sharex=True)
    fig.suptitle('heat diffusion in 1D over time')
    ax[0].set(ylabel='x', title='crank-nicolson')
    ax[1].set(xlabel='time', ylabel='x', title='forward euler')
    ax[2].set(xlabel='time', ylabel='x', title='leapfrog')
    im0 = ax[0].imshow(crank_nicolson(f0, L*k/2., dt/2., steps))
    im1 = ax[1].imshow(fwd_euler_sol)
    im2 = ax[2].imshow(leapfrog(f0, fwd_euler_sol[:,1], L*k/2., dt, steps))
    return im0, im1, im2



#  Aim: animate heat equation in two dimensions (x, y) over time
#  ansatz taken from https://en.wikipedia.org/wiki/Discrete_Poisson_equation
#
#  with fwd euler: result diverges...
#  SOLUTION: time step was to large!
#   dt <= dx**2 / (2 * k)
#  floating point error made it go above slightly => EXPLODES!


def fwd_euler_gen_mat(f, D, dt, tmax):
    # input/output shape: (N**2,)
    I = sparse.identity(f.size, format='csc')
    A = I + dt * D
    t = 0
    while t <  tmax:
        f = A @ f
        t += dt
        yield f

def fwd_euler_gen_direct(f, dx, dt, tmax):
    # input/output shape: (N, N)
    # about 0.6 s slower than `fwd_euler_gen_mat` for N = 100
    t = 0
    while t < tmax:
        f = f + dt * (np.roll(f,1,0) + np.roll(f,-1,0) + np.roll(f,1,1) \
                + np.roll(f,-1,1) - 4*f) / dx**2
        t += dt
        yield f


#  N = 100
#  tmax = 100
def test_2d(N=100, tmax=100):
    step = 20 / N
    dt = step**2 / 4.
    x, y = np.mgrid[-10:10:step,-10:10:step]

    # LAPLACIAN
    # bounds are not correct, L2d needs to be a block matrix, not a band matrix
    # but for this problem it is of lower priority
    # variant with correct bounds see below
    L2d = sparse.diags([1, 1, -4, 1, 1], [-N, -1, 0, 1, N], shape=(N**2, N**2),
                       format='csc') / step**2

    # INITIAL CONDITIONS
    # circle, soft
    get_f0 = lambda s, x, y: np.exp(-(x**2+y**2)/(2*s**2)) / (s*np.sqrt(2*np.pi))
    # line
    #  get_f0 = lambda s, x, y: np.exp(-(x**2)/(2*s**2)) / (s*np.sqrt(2*np.pi))
    f0 = get_f0(2, x+4, y+4) + get_f0(2, x-4, y-4) + get_f0(2, x+4, y-4) \
        + get_f0(2, x-4, y+4)
    #  f0 = get_f0(2, x)
    #  f0 *= 1e24
    # circle, sharp
    #  f0 = np.piecewise(x, [np.abs(x**2+y**2)<1.], [1.])
    # rectangle, sharp
    #  f0 = np.zeros(x.shape)
    #  f0[45:55,45:55] = 1.

    fegm = fwd_euler_gen_mat(f0.reshape(N**2), L2d, dt, tmax)
    #  fegd = fwd_euler_gen_direct(f0, step, dt, tmax)

    # PLOTTING
    fig = plt.figure()
    ax = fig.add_subplot(111, xticklabels=[], yticklabels=[],
                         title='2D Heat Equation over time - Forward Euler solver')
    im = ax.imshow(f0, animated=True)

    def step(data):
        im.set_array(data.reshape(N, N))
        return im,

    anim = animation.FuncAnimation(fig, step, frames=fegm, interval=20,
                                   blit=True, repeat=False)
    return fig, anim

#  f, a = test_2d()



# LAPLACIAN WITH CORRECT BOUNDS
#
#  L = sparse.diags([1, -4, 1], [-1, 0, 1], shape=(N, N), format='csc')
#  I = sparse.identity(N, format='csc')
#  mat_grid = np.array([[None]*N]*N)
#  for i in range(1, N-1):
#      mat_grid[i, i] = L
#      mat_grid[i, i+1] = I
#      mat_grid[i, i-1] = I
#  mat_grid[0, 0] = L
#  mat_grid[0, 1] = I
#  mat_grid[N-1, N-1] = L
#  mat_grid[N-1, N-2] = I
#
#  L2d = sparse.bmat(mat_grid, format='csc') / step**2
