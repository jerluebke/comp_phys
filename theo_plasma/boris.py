# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mayavi import mlab


def E(r, t):
    return np.array([0., 0., .001])

def B(r, t):
    return np.array([0., 0., 1.])

def cross(a, b):
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ])

def boris_gen(r, v, E, B, qm, dt, t, tmax):
    dtqm = qm * dt
    while t < tmax:
        t += .5 * dt                                # t_n+1/2
        r = r + .5 * dt * v                         # r_n+1/2
        B_rt = B(r, t)
        E_rt = E(r, t)
        p = .5 * dtqm * B_rt
        a_sq = .25 * dtqm**2 * np.dot(B_rt, B_rt)
        v = v + .5 * dtqm * E_rt                    # v-
        v_prime = v + cross(v, p)                   # v'
        v = v + 2 * cross(v_prime, p) / (1 + a_sq)  # v+
        v = v + .5 * dtqm * E_rt                    # v_n+1
        r = r + .5 * dt * v                         # r_n+1
        t += .5 * dt
        yield r


r0 = np.array([0., 0., 0.])
v0 = np.array([.5, 0., 0.])
qm = 1
dt = .25
tmax = 50
steps = int(tmax/dt)

boris = boris_gen(r0, v0, E, B, qm, dt, 0, tmax)

#  fig, axs = plt.subplots(1, 3)

data = np.zeros((steps, 3))
data[0] = r0

#  lines = [ax.plot([], [])[0] for ax in axs]

#  def update(i):
#      pass
#  traj, = ax.plot(*[[r0[i]] for i in (0, 1, 2)], 'b-')
#  part, = ax.plot(*[[r0[i]] for i in (0, 1, 2)], 'ro')
#
#  def update(i):
#      data[i] = next(boris)
#      traj.set_data(data[:i,0], data[:i,1])
#      traj.set_3d_properties(data[:i, 2])
#      part.set_data(data[i,0], data[i,1])
#      part.set_3d_properties(data[i,2])
#      #  xmin = data[:i,0].min()
#      #  xmax = data[:i,0].max()
#      #  ymin = data[:i,1].min()
#      #  ymax = data[:i,1].max()
#      #  zmin = data[:i,2].min()
#      #  zmax = data[:i,2].max()
#      #  ax.set_xlim3d(xmin-.1, xmax+.1)
#      #  ax.set_ylim3d(ymin-.1, ymax+.1)
#      #  ax.set_zlim3d(zmin-.1, zmax+.1)
#      return traj, part
#
#  anim = animation.FuncAnimation(fig, update, range(1, steps), interval=50,
#                                 blit=True, repeat=False)

sc = np.zeros(steps)
fig = mlab.figure()
traj = mlab.plot3d(data[:1,0], data[:1,1], data[:1,2],
                   tube_radius=None)
part = mlab.points3d(data[:1,0], data[:1,1], data[:1,2],
                     scale_factor=.1, color=(1., 0., 0.))
ts = traj.mlab_source
ps = part.mlab_source


@mlab.animate(delay=100, ui=False)
def update():
    i = 1
    while i < steps:
        data[i] = next(boris)
        ts.reset(points=data[:i], scalars=sc[:i])
        ps.reset(points=np.reshape(data[i], (1,3)))
        fig.scene.reset_zoom()
        i += 1
        yield

anim = update()
