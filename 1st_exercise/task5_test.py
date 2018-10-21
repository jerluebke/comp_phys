# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from simplex import Step_Gen
from simplex_py import step
from test_functions import *


#############
# globals   #
#############
x   = None
func= None
fig = None
ax  = None
s   = None


#############
# functions #
#############

def setup(func_name, xydoms_distinct, numpoints=100, lines=20):
    global x, func, fig, ax

    fig, ax = plt.subplots()
    ax.set_title(func_name)

    # 'bukin-n6' and 'tal' need special treatment, since their x- and y-domains are not the
    # same (unlike the other functions present here)
    if xydoms_distinct:
        func, xbounds, ybounds = test_functions[func_name]
        xl = np.linspace(*xbounds, numpoints)
        yl = np.linspace(*ybounds, numpoints)
        x = np.zeros((3, 2))
        x[:,0] = np.random.uniform(*xbounds, (3,))
        x[:,1] = np.random.uniform(*ybounds, (3,))
    else:
        func, bounds = test_functions[func_name]
        if func_name == 'shekel':
            #  func = func(*bounds, 10)
            func = func(0, 10, 100)
        xl = np.linspace(*bounds, numpoints)
        yl = xl[:]
        x = np.random.uniform(*bounds, (3, 2))

    X, Y = np.meshgrid(xl, yl)
    ax.contour(X, Y, func(X, Y), lines)


def init():
    s.set_data([], [])
    return s,


def anim_py(tol=1e-7, N=100):
    global s

    def anim_func(i):
        global x
        try:
            x = step(x, lambda arg: func(*arg.T), tol)
        except StopIteration:
            print('found minimum at (%f, %f)\nafter %d steps' \
                  % (*x[0], i))
            raise
        if i >= N:
            print('no convergence after %d steps\nlast position: (%f, %f)' \
                  % (i, *x[0]))
            raise StopIteration
        s.set_data(x[:,0], x[:,1])
        return s,

    #  s, = ax.plot([*x[:,0], x[0, 0]], [*x[:,1], x[1, 0]], 'r')
    s, = ax.plot([], [], 'r')
    anim = animation.FuncAnimation(fig, anim_func, frames=100,
                                   interval=500, blit=True, init_func=init)
    return anim


def anim_fortran(tol=1e-7, N=100):
    global s

    def anim_func(i):
        try:
            simplex = next(step_gen)
        except StopIteration:
            if step_gen.converges:
                print('found minimum at (%f, %f)\nafter %d steps' \
                      % (*step_gen.simplex[:,0], step_gen.steps))
            else:
                print('no convergence after %d steps\nlast position: (%f, %f)' \
                      % (step_gen.steps, *step_gen.simplex[:,0]))
            raise
        s.set_data([*simplex[0], simplex[0, 0]],
                   [*simplex[1], simplex[1, 0]])
        return s,

    step_gen = Step_Gen(x, lambda arg: func(*arg), tol=tol, N=N)
    iter(step_gen)
    #  s, = ax.plot([*step_gen.simplex[0], step_gen.simplex[0, 0]],
    #               [*step_gen.simplex[1], step_gen.simplex[1, 0]], 'r')
    s, = ax.plot([], [], 'r')
    anim = animation.FuncAnimation(fig, anim_func, frames=100,
                                   interval=500, blit=True, init_func=init)
    return anim


#  setup('tal', True)
#  a = anim_fortran()
#  a = anim_py()
