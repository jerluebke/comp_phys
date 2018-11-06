# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

#  TODO
#  test other functions
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
    def __init__(self, xmin, xmax, N, maxlevel):
        self.x              = []
        self.dx             = []
        self._injection     = []
        self._weighting     = []
        self._interpolation = []
        self.L              = []
        self.psi            = []
        self.defect         = None
        self.N              = N
        self.maxlevel       = maxlevel
        self.iterations     = 0

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
            #  self.defect.append(np.zeros(n))

        self.restrict   = self._weighting
        self.prolong    = self._interpolation
        self._solver    = self._jacobi


    def solve(self, rho, steps=[0], tol=1e-7, maxiter=1000):
        if callable(rho):
            rho = rho(self.x[0])
        # reset psi
        for i in range(self.maxlevel):
            self.psi[i] = np.zeros(self.N//2**i)
        self.iterations = 0

        for i in range(maxiter):
            self.solve_one(rho, 0, steps)
            self.defect = rho - self.L[0] @ self.psi[0]
            if np.mean(np.abs(self.defect[1:-1])) < tol:
                print('convergence after %d iterations' % self.iterations)
                return self.psi[0]
        return self.psi[0]


    def solve_one(self, rho, level=0, steps=[0]):
        """
        RECURSIVE
        """
        self.smooth(rho, level, steps[level])
        defect = rho - self.L[level] @ self.psi[level]

        if level < self.maxlevel-1 and steps[level+1]:
            self.solve_one(self.restrict[level] @ defect,
                           level+1, steps)
            error = self.prolong[level] @ self.psi[level+1]
            self.psi[level] += error

        self.smooth(rho, level, steps[level])


    def smooth(self, rho, level, steps=1):
        for _ in range(steps):
            self.psi[level] = self._solver(rho,
                                           self.psi[level],
                                           self.dx[level])
            self.psi[level][0] = self.psi[level][-1] = 0
        self.iterations += steps


    def _jacobi(self, rho, psi, dx):
        return (np.roll(psi, 1) + np.roll(psi, -1) - dx**2 * rho) / 2.


    def _omega_jacobi(self, rho, psi, dx, omega):
        return (1-omega)*psi + omega*(
            np.roll(psi, 1) + np.roll(psi, -1) - dx**2*rho) / 2.


    def _gaus_seidel(self, rho, psi, dx):
        res = np.zeros(psi.size)
        for i in range(1, psi.size-1):
            res[i] = (psi[i+1] + res[i-1] - rho[i]*dx**2) / 2.
        return res



def rho_func(x):
    return np.sin(x)*np.exp(-x**2)

mg  = Multigrid(-5, 5, 128, 3)
x   = mg.x[0]
rho = rho_func(x)
L   = mg.L[0]

sol = mg.solve(rho, steps=[10, 0])          # 10740 iterations
sol = mg.solve(rho, steps=[100, 21, 0])     # 17182 iterations
sol = mg.solve(rho, steps=[160, 15, 120])   # 21240 iterations

fig, ax = plt.subplots()
ax.set_title('Comparison Multigrid - poisson equation')
ax.plot(x, rho, label=r'$\rho$')
ax.plot(x, sol, label=r'$\psi$')
ax.plot(x, L@sol, label=r'$\Delta\psi$')

plt.legend()
