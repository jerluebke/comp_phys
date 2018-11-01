# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse


class Multigrid:
    """
    multigrid
    =========
    while defect large
        smooth psi
        get defect
        defect to coarse grid
        RECURSIVE solve defect eq for error F
        error to fine grid
        correct psi
        smooth psi

    """
    def __init__(self, xmin, xmax, N, maxlevel, rpl, spl):
        """
        Parameter
        =========
        xmin, xmax  :   domain boundaries
        N           :   numpoints, N = 2**i, i < maxlevel
        maxlevel    :   deepest level
        rpl         :   repeats per level, list with len = maxlevel-1
        spl         :   smoothing steps per level, list with len = maxlevel

        Initiate domain, solution vectors and transformation matrices

        How to define a symmetric scheme according to which to step through
        the grids:
            >>> rpl = [1, 2, 1]     # 4 levels
            1 -> 2 -> 3 -> 4 -> 3 -> 2 -> 3 -> 4 -> 3 -> 2 -> 1

        """
        self.x              = []
        self.dx             = []
        self._injection     = []
        self._weighting     = []
        self._interpolation = []
        self.L              = []    # Laplacian
        self.psi            = []    # solution vector
        self.repeat         = rpl
        self.steps          = spl
        self.converges      = False
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
        self.N          = N
        self.maxlevel   = maxlevel
        self.restrict   = self._injection
        self.prolong    = self._interpolation
        self._solver    = self._jacobi

    def solve_all(self, rho, eps=1e-7, maxiter=1000):
        # reset solution
        for i in range(self.maxlevel):
            self.psi[i] = np.zeros(self.N // 2**i)
        self.converges  = False
        self.iterations = 0     # is updated with every call of smooth
        if callable(rho):
            rho = rho(self.x[0])
        for i in range(maxiter):
            psi = self.solve_one(rho)
            defect = np.abs(rho - self.L[0] @ psi)
            if np.mean(defect) < eps:
                self.converges = True
                return psi
        return psi

    def solve_one(self, rho, level=0):
        """
        RECURSIVE
        """
        for i in range(self.repeat[level]):
            #  self.psi[level] = self.smooth(rho, level)
            self.smooth(rho, level)
            if level < self.maxlevel-1:
                defect = rho - self.L[level] @ self.psi[level]
                error = self.solve_one(self.restrict[level] @ defect, level+1)
                self.psi[level] -= self.prolong[level] @ error
        # post-smoothing
        #  self.psi[level] = self.smooth(rho, level)
        self.smooth(rho, level)
        return self.psi[level]

    def smooth(self, rho, level):
        for _ in range(self.steps[level]):
            self.psi[level] = self._solver(rho, self.psi[level], self.dx[level])
            self.psi[level][0] = self.psi[level][-1] = 0
        self.iterations += self.steps[level]
        return self.psi[level]

    def _jacobi(self, rho, psi, dx):
        return (np.roll(psi, 1) + np.roll(psi, -1) - dx**2 * rho) / 2.
        #  return (np.roll(psi,1)+np.roll(psi,-1)-dx**2*rho)/4.+psi/2.

    def _gaus_seidel(self, rho, psi, dx):
        pass
