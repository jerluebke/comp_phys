# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from mayavi import mlab

plt.rcParams['animation.ffmpeg_path'] = 'D:\\source\\Libs\\ffmpeg-20180521-c24d247-win64-static\\bin\\ffmpeg.exe'


def E(r, t):
    #  return np.array([0., 0., .001])
    return np.array([0., 0., 0.])

def B(r, t):
    #  return np.array([0., 0., 1.])
    d = np.sqrt(r[0]**2 + r[1]**2)
    b0 = 1/d
    return np.array([-b0*r[1]/d, b0*r[0]/d, .5])

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

def b_line(r0, n=100, ds=-.1):
    # a3.plot(l[:,0], l[:,1], l[:,2], 'r-')
    l = np.zeros((n,3))
    l[0] = r0
    for i in range(1, n):
        r = r0 + .5 * ds * B(r0, 0.)
        r = r0 + ds * B(r, 0.)
        l[i] = r
        r0 = r
    return l


r0 = np.array([1., 0., 0.])
v0 = np.array([0., 1., 0.])
qm = 1
dt = .1
tmax = 1000
steps = int(tmax/dt)

#  boris = boris_gen(r0, v0, E, B, qm, dt, 0, tmax)


def setup(anim=False):
    boris = boris_gen(r0, v0, E, B, qm, dt, 0, tmax)
    data = np.array([e for e in boris])
    time = np.arange(0, tmax, dt)
    if time.size == data[:,0].size-1:
        data = data[:-1]

    fig = plt.figure(figsize=(12,9))
    fig.tight_layout()
    #  plt.get_current_fig_manager().window.state('zoomed')    # fullscreen
    a3 = fig.add_subplot(221, projection='3d', title='trajectory',
                         xlabel='x', ylabel='y', zlabel='z')
    axs = [fig.add_subplot(2, 2, i, title=label, xlabel='time',
                           ylabel='%s-compenent' % label)
           for i, label in zip((2, 3, 4), ('x', 'y', 'z'))]

    if anim:
        traj, = a3.plot([data[0,0]], [data[0,1]], [data[0,2]], 'b-')
        part, = a3.plot([data[0,0]], [data[0,1]], [data[0,2]], 'yo')
        lines = [ax.plot([time[0]], [data[0,i]], 'b-')[0]
                 for i, ax in enumerate(axs)]

        return data, time, fig, (a3, axs), (traj, part, lines)

    else:
        a3.plot(data[:,0], data[:,1], data[:,2], 'b-')
        a3.plot([data[-1,0]], [data[-1,1]], [data[-1,2]], 'yo')
        for i, ax in enumerate(axs):
            ax.plot(time, data[:,i], 'b-')
        return (a3, axs)


def update(i, data, time, a3, axs, traj, part, lines):
    traj.set_data(data[:i,0], data[:i,1])
    traj.set_3d_properties(data[:i, 2])
    part.set_data(data[i,0], data[i,1])
    part.set_3d_properties(data[i,2])
    xmin = data[:i,0].min(); xmax = data[:i,0].max()
    ymin = data[:i,1].min(); ymax = data[:i,1].max()
    zmin = data[:i,2].min(); zmax = data[:i,2].max()
    a3.set_xlim3d(xmin-.1, xmax+.1)
    a3.set_ylim3d(ymin-.1, ymax+.1)
    a3.set_zlim3d(zmin-.1, zmax+.1)
    for j, line in enumerate(lines):
        line.set_data(time[:i], data[:i,j])
    axs[0].set(xlim=(time[0]-.1, time[i]+.1), ylim=(xmin-.1, xmax+.1))
    axs[1].set(xlim=(time[0]-.1, time[i]+.1), ylim=(ymin-.1, ymax+.1))
    axs[2].set(xlim=(time[0]-.1, time[i]+.1), ylim=(zmin-.1, zmax+.1))
    #  fig.canvas.draw()
    return [traj, part, *lines]


def save(name='boris.mp4'):
    data, time, fig, (a3, axs), (traj, part, lines) = setup(True)
    FFWriter = animation.FFMpegWriter(fps=30)
    animation.FuncAnimation(fig, update, range(1, steps), interval=10,
                            fargs=(data, time, a3, axs, traj, part, lines),
                            blit=True, repeat=False).save(
                                name, writer=FFWriter, dpi=100)


if __name__ == '__main__':
    save()



def _3d_anim(maxsize=steps):
    boris = boris_gen(r0, v0, E, B, qm, dt, 0, tmax)
    data = np.zeros((maxsize, 3))
    data[0] = r0

    sc = np.zeros(steps)
    fig = mlab.figure()
    traj = mlab.plot3d(data[:1,0], data[:1,1], data[:1,2],
                       tube_radius=None)
    part = mlab.points3d(data[:1,0], data[:1,1], data[:1,2],
                         scale_factor=.1, color=(1., 0., 0.))
    ts = traj.mlab_source
    ps = part.mlab_source

    @mlab.animate(delay=100, ui=False)
    def update_mlab():
        i = 1
        while i < maxsize:
            data[i] = next(boris)
            ts.reset(points=data[:i], scalars=sc[:i])
            ps.reset(points=np.reshape(data[i], (1,3)))
            fig.scene.reset_zoom()
            i += 1
            yield
        while i < steps:
            data[0:maxsize-1] = data[1:maxsize]
            data[-1] = next(boris)
            ts.trait_set(points=data, scalars=sc)
            ps.trait_set(points=np.reshape(data[-1], (1,3)))
            fig.scene.reset_zoom()
            i += 1
            yield

    return boris, data, update_mlab()
#  mlab.show()
