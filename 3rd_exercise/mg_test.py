# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from multigrid import *

N = 256
gridkwds = dict(rho_func=lambda x: np.sin(x)*np.exp(-x**2),
                N=N, xmin=-5, xmax=5, levels=int(math.log(N, 2)+1e-2)-1)


def solve(imax=20, tol=1e-7, ax=None, **gridkwds):
    if not ax:
        ax = plt.gca()
    grids = init(**gridkwds)
    err = [np.max(np.abs(grids[0].defect))]
    i = 1
    while i < imax:
        solve_one_v(grids, 'red_black')
        err.append(np.max(np.abs(grids[0].defect)))
        i += 1
    ax.semilogy(np.arange(i), err)
    return grids


g = solve(**gridkwds)
