# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

#  TODO
#  documentation!


#  multigrid
#  =========
#  while defect large
#      smooth psi
#      get defect
#      defect to coarse grid
#      RECURSIVE solve defect eq for error F
#      error to fine grid
#      correct psi
#      smooth psi


class Multigrid:
    def __init__(self, rho_func, xmin, xmax, N, maxlevel, solver='red_black'):
        self._solver_dict = {
            'jacobi'        :   self._jacobi,
            'omega_jacobi'  :   self._omega_jacobi,
            'gauss_seidel'  :   self._gauss_seidel,
            'red_black'     :   self._red_black
        }
        self.x              = []
        self.dx             = []
        self._injection     = []
        self._weighting     = []
        self._interpolation = []
        self.L              = []
        self.psi            = []
        self.rho            = []
        self.defect         = []
        self.n              = []
        self.N              = N
        self.maxlevel       = maxlevel

        for i in range(maxlevel):
            n = N // 2**i + 1
            self.x.append(np.linspace(xmin, xmax, n))
            self.dx.append(self.x[i][1]-self.x[i][0])
            self._injection.append(
                sparse.block_diag([[0.,1.]]*(n//2), format='csc')
            )
            self._weighting.append(
                sparse.block_diag([[.25,.5]]*(n//2), format='csc') + \
                sparse.block_diag([[0.,0.,.25]]+[[0.,.25]]*(n//2-2)+[0.])
            )
            self._interpolation.append(
                sparse.block_diag([[[.5],[1.]]]*(n//2), format='csc') + \
                sparse.block_diag([[[0.],[0.],[.5]]] \
                                  + [[[0.],[.5]]]*(n//2-2) \
                                  + [0.], format='csc')
            )
            self.L.append(
                sparse.diags([1, -2, 1], [-1, 0, 1], shape=(n-2, n-2),
                             format='csc') / self.dx[i]**2
            )
            self.psi.append(np.zeros(n))
            self.rho.append(np.zeros(n))
            self.defect.append(np.zeros(n))
            self.n.append(n)

        self.rho[0]     = rho_func(self.x[0])
        self.restrict   = self._weighting
        self.prolong    = self._interpolation
        self._solver    = self._solver_dict[solver]


    def solve(self, solver, tol=1e-7, maxiter=30, ax=None):
        self._solver = self._solver_dict[solver]

        # set up plotting
        if not ax:
            ax = plt.gca()
        ax.set(xlabel='iterations', ylabel='error mean',
               title='Multigrid - datapoints = %d' % self.N)

        # reset results
        for i in range(self.maxlevel):
            self.psi[i][:] = 0

        self.set_defect()
        err = [np.max(np.abs(self.defect[0]))]
        i = 1
        while i < maxiter:
            self.solve_one_v()
            self.set_defect()
            err.append(np.max(np.abs(self.defect[0])))
            i += 1

        ax.semilogy(np.arange(i), err)


    def solve_one_v(self, level=0):
        """
        RECURSIVE
        """
        self.smooth(level)
        if level < self.maxlevel-1:
            self.set_defect(level)
            self.rho[level+1] = self.restrict[level] @ self.defect[level]
            self.psi[level+1][:] = 0
            self.solve_one_v(level+1)
            self.psi[level] += self.prolong[level] @ self.psi[level+1]
        self.smooth(level)

    def set_defect(self, level=0):
        n = self.n[level]
        self.defect[level][1:n-1] = self.rho[level][1:n-1] \
                - self.L[level] @ self.psi[level][1:n-1]

    def smooth(self, level, **kwds):
        self._solver(self.psi[level],
                     self.rho[level],
                     self.dx[level],
                     self.n[level],
                     **kwds)
        self.psi[level][0] = self.psi[level][-1] = 0

    def _jacobi(self, psi, rho, dx, n, **kwds):
        psi[1:n-1] = .5 * (psi[0:n-2] + psi[2:n] - dx**2 * rho[1:n-1])

    def _omega_jacobi(self, psi, rho, dx, n, **kwds):
        omega = kwds.get('omega', 1./2.)
        psi[1:n-1] = .5 * omega * (psi[0:n-2] + psi[2:n] - dx**2 * rho[1:n-1]) \
                + (1. - omega) * psi[1:n-1]

    def _gauss_seidel(self, psi, rho, dx, n, **kwds):
        psi[1:n-1] = np.array([.5 * (psi[i-1] + psi[i+1] - dx**2 * rho[i])
                               for i in range(1, n-1)])

    def _red_black(self, psi, rho, dx, n, **kwds):
        psi[1:n-1:2] = .5 * (psi[0:n-2:2] + psi[2:n:2] - dx**2 * rho[1:n-1:2])
        psi[2:n-2:2] = .5 * (psi[1:n-3:2] + psi[3:n-1:2] - dx**2 * rho[2:n-2:2])



def rho_func(x):
    return np.sin(x)*np.exp(-x**2)


N = 256
xmin = 0.
xmax = 1.
levels = int(math.log(N, 2)+1e-2)-1
mg = Multigrid(rho_func, xmin, xmax, N, levels)
mg.solve('red_black', tol=1e-3, maxiter=30, ax=plt.gca())
