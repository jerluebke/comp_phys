#!/usr/bin/env python

import numpy as np
import math
import matplotlib
#  matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation

# runtime parameters
n = 256     # finest resolution
xmin = -5.0  # left border of domain
xmax =  5.0  # right border of domain
iterations = 30 # number of repititions of multigrid algorithm

# initialization for rho
init_rho = lambda x: np.sin(x) * np.exp(-x**2)

class Grid:
    def __init__(self, n, xmin, xmax):
        '''grid-data for one level; n + 1 gridpoints, range from xmin to xmax'''
        self.n = n + 1
        self.f = np.zeros(n + 1)
        self.rho = np.zeros(n + 1)
        self.defect = np.zeros(n + 1)
        self.dx = (xmax - xmin) / n
        self.x = np.linspace(xmin, xmax, n + 1) # g.x = xmin + (0:n)' * g.dx;

    def __repr__(self):
        return "Grid(" + str(self.n + 1) + ", " \
                + str(self.x[0]) + ", " + str(self.x[-1])+ ")"

def defect(grid):
    '''update defect array of grid'''
    n = grid.n
    grid.defect[1:n-1] = grid.rho[1:n-1] - \
                (grid.f[2:n] - 2.*grid.f[1:n-1] + grid.f[0:n-2]) / grid.dx**2

def prolong(coarse, fine):
    '''coarse to fine, fine is modified'''
    nc = coarse.n
    nf = fine.n
    fine.f[2:nf-2:2] += coarse.f[1:nc-1]
    fine.f[1:nf-1:2] += 0.5 * (coarse.f[0:nc-1] + coarse.f[1:nc])

def restrict(fine, coarse):
    '''fine to coarse, coarse is modified'''
    nc = coarse.n
    nf = fine.n
    coarse.rho[0:nc] = 0.
    # injektion
    #coarse.rho[1:nc-1] = fine.defect[2:nf-2:2]
    # full weighting
    coarse.rho[1:nc-1] = 0.5*fine.defect[2:nf-2:2]\
                        + 0.25*(fine.defect[3:nf-1:2] + fine.defect[1:nf-3:2])
    coarse.f[0:nc] = 0.

def smooth(grid):
    n = grid.n
    # Jacobi:
    # grid.f[1:n-1] = 0.5 * (grid.f[0:n-2] + grid.f[2:n] - grid.dx**2 * grid.rho[1:n-1])
    # omega-Jacobi:
    # omega = 0.5
    # grid.f[1:n-1] = 0.5 * omega * ( grid.f[0:n-2] + grid.f[2:n] - grid.dx**2 * grid.rho[1:n-1] ) \
    #              + (1.-omega) * grid.f[1:n-1];
    # Gauss-Seidel:
    #  for i in range(1, n-1):
    #     grid.f[i] = 0.5 * (grid.f[i-1] + grid.f[i+1] - grid.dx**2 * grid.rho[i])
    # Red-Black Gauss-Seidel:
    grid.f[1:n-1:2] = 0.5 * (grid.f[0:n-2:2] + grid.f[2:n:2] - grid.dx**2 * grid.rho[1:n-1:2])
    grid.f[2:n-2:2] = 0.5 * (grid.f[1:n-3:2] + grid.f[3:n-1:2] - grid.dx**2 * grid.rho[2:n-2:2])

def solve_one(grids, level):
    '''one iteration on the specified level, calls itself recursively'''
    smooth(grids[level])
    if level < len(grids) - 1:
        defect(grids[level])
        restrict(grids[level], grids[level+1])
        solve_one(grids, level+1)
        prolong(grids[level+1], grids[level])
    smooth(grids[level])

# initialize grids
eps = 1.e-2 # for circumventing discrete floating point effects
grids = [Grid(n//2**i, xmin, xmax) for i in range(int(math.log(n, 2) + eps) - 1)]
grids[0].rho[1:grids[0].n-1] = init_rho(grids[0].x[1:grids[0].n-1])

# solve without plotting
#for i in range(iterations):
#    solve_one(grids, 0)

# set up plot
est = 1.5 * np.amax(np.abs(grids[0].rho)) # estimate approiate y-range

fig = plt.figure()

line_rho = plt.Line2D([], [], lw=2, marker='.', markerfacecolor='red', markersize=12)
ax_rho = fig.add_subplot(2, 2, 1)
ax_rho.set_xlabel("$x$")
ax_rho.set_ylabel("rho")
ax_rho.set_xlim(xmin, xmax)
ax_rho.set_ylim(-est, est)
ax_rho.add_line(line_rho)
time_text = ax_rho.text(0.02, 0.75, '', transform=ax_rho.transAxes)

line_f = plt.Line2D([], [], lw=2, marker='.', markerfacecolor='red', markersize=12)
ax_f   = fig.add_subplot(2, 2, 2)
ax_f.set_xlabel("$x$")
ax_f.set_ylabel("$f$")
ax_f.set_xlim(xmin, xmax)
ax_f.set_ylim(-est, est)
ax_f.add_line(line_f)

line_err = plt.Line2D([], [], lw=2, marker='.', markerfacecolor='red', markersize=12)
ax_err = fig.add_subplot(2, 2, 3)
ax_err.set_xlabel("$x$")
ax_err.set_ylabel("error")
ax_err.set_xlim(xmin, xmax)
ax_err.set_ylim(-.1*est, .1*est)
ax_err.add_line(line_err)

line_conv = plt.Line2D([], [], lw=2, marker='.', markerfacecolor='red', markersize=12)
ax_conv = fig.add_subplot(2, 2, 4)
ax_conv.set_xlabel("iteration")
ax_conv.set_ylabel("log(max defect)")
ax_conv.set_xlim(0, iterations)
ax_conv.set_ylim(-15,np.log10(2*est))
ax_conv.add_line(line_conv)

defects=[]

def init_plot():
    '''initialize plot'''
    line_rho.set_data(grids[0].x, init_rho(grids[0].x))
    line_f.set_data([], [])
    line_err.set_data([], [])
    line_conv.set_data([], [])
    time_text.set_text('')
    return (line_f, line_err, time_text)

def step(i):
    '''single step of multigrid method, for use with animation'''
    global defects
    solve_one(grids, 0)
    defects.append(np.max(grids[0].defect))
    #  print(defects)
    est_f = 1.2 * np.amax(np.abs(grids[0].f)) # estimate approiate y-range
    est_err = 1.2 * np.amax(np.abs(grids[0].defect)) # estimate approiate y-range
    ax_f.set_ylim(-est_f, est_f)
    ax_err.set_ylim(-est_err, est_err)
    line_f.set_data(grids[0].x, grids[0].f)
    line_err.set_data(grids[0].x, grids[0].defect)
    line_conv.set_data(np.arange(i+1), np.log10(defects))
    time_text.set_text("iteration = {}".format(i))
    return (line_f, line_err, line_conv, time_text)

anim = animation.FuncAnimation(fig, step, init_func=init_plot, frames=iterations,
                               interval=20, blit=False, repeat=False)

plt.show()
