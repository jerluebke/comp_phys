# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse


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
        self.defect         = []
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
            self.defect.append([np.zeros(n),    # before
                                np.zeros(n)])   # after
        self.maxlevel   = maxlevel
        self.restrict   = self._injection
        self.prolong    = self._interpolation
        self._solver    = self._jacobi

    def solve(self, rho, i=0, steps=1):
        """
        RECURSIVE

        V scheme
        """
        if callable(rho):
            rho = rho(self.x[i])
        self.psi[i] = self.smooth(rho, i, steps)
        self.defect[i][0] = rho - self.L[i] @ self.psi[i]
        if i < self.maxlevel-1:
            error = self.solve(self.restrict[i] @ self.defect[i][0], i+1)
            self.psi[i] -= self.prolong[i] @ error
        self.psi[i] = self.smooth(rho, i, steps)
        self.defect[i][1] = rho - self.L[i] @ self.psi[i]
        return self.psi[i]

    def smooth(self, rho, i, j=1):
        """
        perform j solution steps for the poisson equation on gridlevel i
        """
        psi = self.psi[i]
        for k in range(j):
            psi = self._solver(rho, psi, self.dx[i])
            psi[0] = psi[-1] = 0
        return psi

    def _jacobi(self, rho, psi, dx):
        return (np.roll(psi, 1) + np.roll(psi, -1) - dx**2 * rho) / 2.

    def _gaus_seidel(self, rho, psi, dx):
        pass
