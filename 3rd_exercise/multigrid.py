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
    _injection      = []
    _weighting      = []
    _interpolation  = []
    L               = []
    psi             = []
    defect          = []

    def __init__(self, N, maxlevel):
        for i in range(maxlevel):
            n = N // 2**i
            self._injection.append(
                sparse.block_diag([[0.,1.]]*n//2, format='csc')
            )
            self._weighting.append(
                sparse.block_diag([[.25,.5]]*n//2, format='csc') + \
                sparse.block_diag([[0.,0.,.25]]+[[0.,.25]]*(n//2-2)+[0.])
            )
            self._interpolation.append(
                sparse.block_diag([[[.5],[1.]]]*n//2, format='csc') + \
                sparse.block_diag([[[0.],[0.],[.5]]] \
                                  + [[[0.],[.5]]]*(n//2-2) \
                                  + [0.], format='csc')
            )
            self.L.append(
                sparse.diags([1, -2, 1], [-1, 0, 1], shape=(n, n), format='csc')
                # TODO: dx
            )
            self.psi.append(np.zeros(n))
            self.defect.append((np.zeros(n),    # before
                                np.zeros(n)))   # after
        self.maxlevel   = maxlevel
        self.restrict   = self._injection
        self.prolong    = self._interpolation
        self._solver    = self._gaus_seidel

    def solve(self, rho, i):
        """
        RECURSIVE

        V step
        """
        self.psi[i] = self.smooth(self.psi[i], rho)
        self.defect[i][0] = rho - self.L[i] @ self.psi[i]
        if i < self.maxlevel:
            error = self.solve(self.restrict[i] @ self.defect[i], i+1)
            self.psi[i] -= self.prolong[i] @ error
        self.psi[i] = self.smooth(self.psi[i], rho)
        self.defect[i][1] = rho - self.L[i] @ self.psi[i]
        return self.psi[i]

    def smooth(self, psi, rho):
        psi = self._solver(psi, rho)
        psi[0] = psi[-1] = 0
        return psi

    def _jacobi(self, psi, rho, dx):
        return (np.roll(psi, 1) + np.roll(psi, -1) - dx**2 * rho) / 2.

    def _gaus_seidel(self, psi, rho):
        pass
