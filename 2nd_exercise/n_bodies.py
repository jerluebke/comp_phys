# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

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

q0 = np.array([[1., 1.], [0., 0.], [2., 0.]])
q1 = q0 + np.array([[.05, 0.], [0., 0.], [0., 0.]])
m = np.array([1, 2, 1])



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
                    * F(dist_vec[j,:,:], dist_mag[j,:], m, j)
                    #  * F(dist_vec[j,j-1,:], m[j-1])   # for F = F_two
    return traj

def get_dist(q):
    p = q.shape[0]
    n = q.shape[1]
    dv = np.zeros((p, p, n))
    for i in range(p):
        for j in range(p):
            dv[i,j,:] = q[j,:] - q[i,:]
    dm = np.sum(dv**2, axis=2)
    return np.ma.masked_where(dv==0,dv), np.ma.masked_where(dm==0,dm)

def F_n(rv, rm, m, i):
    p = m.size
    #  ddx = np.sum([
    #      m[j] * rv[i,j,:] / rm[i,j]**(3./2.) for j in range(i)
    #  ]) + np.sum([
    #      m[j] * rv[i,j,:] / rm[i,j]**(3./2.) for j in range(i+1,p)
    #  ])
    #  ddx = np.sum([m[j] * rv[j,:] / rm[j]**(3./2.) for j in range(p)])
    ddx = 0
    for j in range(p):
        if i == j:
            continue
        ddx += m[j] * rv[j,:] / rm[j]**(3./2.)
    return ddx
