# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] ='D:\\source\\Libs\\ffmpeg-20180521-c24d247-win64-static\\bin\\ffmpeg.exe'


def solve_all_two(q0, q1, m, dt, steps, F, *fargs, **fkwds):
    res = np.zeros((2, 2, steps))   # particles, dimensions, steps
    res[:,:,0] = q0
    res[:,:,1] = q1
    for i in range(2, steps):
        r = res[1,:,i-1] - res[0,:,i-1]
        res[0,:,i] = 2. * res[0,:,i-1] - res[0,:,i-2] + dt**2 * F(r, m[1])
        res[1,:,i] = 2. * res[1,:,i-1] - res[1,:,i-2] + dt**2 * F(-r, m[0])
    return res

def F_two(r, m):
    return m * r / np.sum(r**2)**(3./2.)

#  q0 = np.array([[1., 1.], [0., 0.]])
#  q1 = q0 + np.array([[.05, 0.], [0., 0.]])
#  m = np.array([1, 2])
#  res = solve_all_two(q0, q1, m, .05, 1000, F_two)
#  plt.figure()
#  [plt.plot(res[i,0,:], res[i,1,:], label='mass = %f' % m[i]) for i in (0, 1)]



def solve_all_three(q0, q1, m, dt, steps, F, *fargs, **fkwds):
    res = np.zeros((3, 2, steps))   # particles, dimensions, steps
    res[:,:,0] = q0
    res[:,:,1] = q1
    for i in range(2, steps):
        r21 = res[1,:,i-1] - res[0,:,i-1]
        r31 = res[2,:,i-1] - res[0,:,i-1]
        r32 = res[2,:,i-1] - res[1,:,i-1]
        res[0,:,i] = 2. * res[0,:,i-1] - res[0,:,i-2] + dt**2 \
            * F(r21, r31, m[1], m[2])
        res[1,:,i] = 2. * res[1,:,i-1] - res[1,:,i-2] + dt**2 \
                * F(-r21, r32, m[0], m[2])
        res[2,:,i] = 2. * res[2,:,i-1] - res[2,:,i-2] + dt**2 \
                * F(-r31, -r32, m[0], m[1])
    return res

def F_three(rj, rk, mj, mk):
    return mj * rj / np.sum(rj**2)**(3./2.) + mk * rk / np.sum(rk**2)**(3./2.)

#  q0 = np.array([[1., 1.], [0., 0.], [2., 0.]])
#  q1 = q0 + np.array([[.05, 0.], [0., 0.], [0., 0.]])
#  m = np.array([1, 2, 1])
#  res = solve_all_three(q0, q1, m, .05, 60, F_three)
#  plt.figure()
#  [plt.plot(res[i,0,:], res[i,1,:], label='mass = %f' % m[i]) for i in (0, 1, 2)]
#  plt.legend()



def solve_all_n(q0, q1, m, dt, steps, F, *fargs, **fkwds):
    """
    traj[particle, dimension, steps]
    """
    p = q0.shape[0]
    n = q0.shape[1]
    traj = np.zeros((p, n, steps))
    traj[:,:,0] = q0
    traj[:,:,1] = q1
    for i in range(2, steps):
        dist_vec, dist_mag = get_dist(traj[:,:,i-1])
        for j in range(p):
            traj[j,:,i] = 2. * traj[j,:,i-1] - traj[j,:,i-2] + dt**2 \
                    * F(dist_vec, dist_mag, m, j)
    return traj

def get_dist(q):
    """q : positions"""
    p = q.shape[0]  # number of particles
    n = q.shape[1]  # dimensions
    rv = np.zeros((p, p, n))    # distance vectors
    for i in range(p):
        for j in range(i):
            rv[i,j,:] = q[j,:] - q[i,:]     # distances between each particle
    rm = np.sum(rv**2, axis=2)  # distance magnitutes
    #  return np.ma.masked_where(dv==0,dv), np.ma.masked_where(dm==0,dm)
    return rv, rm

def F_n(rv, rm, m, i):
    p = m.size
    return np.sum(
            [m[j]*rv[i,j,:]/rm[i,j]**(3./2.)
            for j in range(i)],
        axis=0) \
        - np.sum(
            [m[j]*rv[j,i,:]/rm[j,i]**(3./2.)
            for j in range(i+1, p)],
        axis=0)


q0 = np.array([
    [1., 1.],
    [-1., 1.],
    [-1., -1.],
    [1., -1.]
])
q1 = q0 + np.array([
    [-.05, 0.],
    [0., -.05],
    [.05, 0.],
    [0., .05]
])
m = np.array([1., 1., 1., 1.])
#  m = np.array([10., 2., 1.])
STEPS = 1000
res = solve_all_n(q0, q1, m, .05, STEPS, F_n)
#  plt.figure()
#  [plt.plot(res[i,0,:], res[i,1,:]) for i in range(m.size)]

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 5), ylim=(-5, 5))
lines = [ax.plot([], [], 'b--', zorder=1)[0] for _ in range(m.size)]
points = [ax.plot([], [], 'ro', zorder=2)[0] for _ in range(m.size)]

def init():
    for i in range(m.size):
        lines[i].set_data([], [])
        points[i].set_data([], [])
    return [*lines, *points]

def step(i):
    for j in range(m.size):
        lines[j].set_data(res[j,0,:i], res[j,1,:i])
        points[j].set_data(res[j,0,i], res[j,1,i])
    return [*lines, *points]


anim = animation.FuncAnimation(fig, step, frames=STEPS, interval=20, blit=True,
                               repeat=False, init_func=init)

FFWriter = animation.FFMpegWriter(fps=60)
anim.save('4_bodies.mp4', writer=FFWriter, dpi=300)
