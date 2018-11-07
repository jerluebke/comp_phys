# -*- coding: utf-8 -*-

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
    def __init__(self, rho_func, xmin, xmax, N, maxlevel):
        self._solver_dict = {
            'jacobi'        :   self._jacobi,
            'omegajacobi'   :   self._omega_jacobi,
            'gaussseidel'   :   self._gauss_seidel
        }
        self.x              = []
        self.dx             = []
        self._injection     = []
        self._weighting     = []
        self._interpolation = []
        self.L              = []
        self.psi            = []
        self.rho            = []
        #  self.defect         = []
        self.N              = N
        self.maxlevel       = maxlevel
        self.iterations     = 0
        self.conv_psi       = None

        for i in range(maxlevel):
            n = N // 2**i
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
                sparse.diags([1, -2, 1], [-1, 0, 1], shape=(n, n),
                             format='csc') / self.dx[i]**2
            )
            self.psi.append(np.zeros(n))
            self.rho.append(np.zeros(n))
            #  self.defect.append(np.zeros(n))

        self.rho[0]     = rho_func(self.x[0])
        self.restrict   = self._weighting
        self.prolong    = self._interpolation
        self._solver    = self._solver_dict['jacobi']


    def solve(self, solver, steps, scheme,
              tol=1e-7, maxiter=1000,
              ax=None, label=''):
        self._solver = self._solver_dict[solver]

        # set up plotting
        if not ax:
            _, ax = plt.subplots()
        ax.set(xlabel='iterations', ylabel='error mean',
               title='Multigrid - datapoints = %d' % self.N)

        # reset results
        self.iterations = 0
        for i in range(self.maxlevel):
            self.psi[i] = np.zeros(self.N//2**i)

        iter_arr    = np.zeros(maxiter)
        err_arr     = np.zeros(maxiter)
        conv        = False
        for i in range(maxiter):
            self.solve_one(0, steps, scheme)
            iter_arr[i] = self.iterations
            err_arr[i]  = np.mean(np.abs(
                self.rho[0] - self.L[0] @ self.psi[0])[1:-1])
            if not conv and err_arr[i] < tol:
                self.conv_psi = self.psi[0].copy()
                ax.plot([iter_arr[i]], [err_arr[i]], 'ro', zorder=2)
                conv = True
                print('convergence after %d iterations' % self.iterations)

        ax.plot(iter_arr, err_arr,
                label='%s - %s - %s' % (label, str(steps), solver),
                zorder=1)
        ax.plot([0, iter_arr[-1]], [tol, tol], 'k:')
        ax.legend()


    def solve_one(self, level, steps, scheme):
        """
        RECURSIVE
        """
        i = 0
        while i < scheme[level]:
            self.smooth(level, steps[level])
            defect = self.rho[level] - self.L[level] @ self.psi[level]
            if level < self.maxlevel-1:
                self.rho[level+1] = self.restrict[level] @ defect
                self.solve_one(level+1, steps, scheme)
                self.psi[level] += self.prolong[level] @ self.psi[level+1]
            self.smooth(level, steps[level])
            i += 1


    def smooth(self, level, steps):
        i = 0
        while i < steps:
            self.psi[level] = self._solver(self.psi[level],
                                           self.rho[level],
                                           self.dx[level])
            self.psi[level][0] = self.psi[level][-1] = 0
            i += 1
        self.iterations += steps


    def _jacobi(self, psi, rho, dx):
        return (np.roll(psi, 1) + np.roll(psi, -1) - dx**2 * rho) / 2.


    def _omega_jacobi(self, psi, rho, dx, omega=2./3.):
        return (1-omega)*psi + omega*(
            np.roll(psi, 1) + np.roll(psi, -1) - dx**2*rho) / 2.


    def _gauss_seidel(self, psi, rho, dx):
        res = np.zeros(psi.size)
        for i in range(1, psi.size-1):
            res[i] = (psi[i+1] + res[i-1] - rho[i]*dx**2) / 2.
        return res



def rho_func(x):
    return np.sin(x)*np.exp(-x**2)

#  mg  = Multigrid(-5, 5, 128, 3)
#  x   = mg.x[0]
#  rho = rho_func(x)
#  L   = mg.L[0]

#  sol = mg.solve(rho, steps=[10, 0])          # 10740 iterations
#  sol = mg.solve(rho, steps=[100, 21, 0])     # 17182 iterations
#  sol = mg.solve(rho, steps=[160, 15, 120])   # 21240 iterations

#  fig, ax = plt.subplots()
#  ax.set_title('Comparison Multigrid - poisson equation')
#  ax.plot(x, rho, label=r'$\rho$')
#  ax.plot(x, sol, label=r'$\psi$')
#  ax.plot(x, L@sol, label=r'$\Delta\psi$')
#
#  plt.legend()

mg = Multigrid(rho_func, 0, 1, 128, 3)
mg.solve('jacobi', [2, 50, 0], [1, 1, 0], tol=1e-3, maxiter=100)
