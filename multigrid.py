# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt



multigrid
=========
while defect large
    smooth psi
    get defect
    defect to coarse grid
    RECURSIVE solve defect eq for error F
    error to fine grid
    correct psi
    smooth pi


class Multigrid:
    def __init__(self, N, maxlevel):
        for i in range(N):
            self.injection.append(sparse.block_diag([inj]*2))

    def solve(self, rho, i)
        """
        RECURSIVE

        V step
        """
        self.psi[i] = self.smooth(self.psi[i], rho)
        defect[i] = rho - self.L[i] @ self.psi[i]
        if i < self.maxlevel:
            error = self.solve(self.restrict[i] @ rho, i+1)
            psi[i] -= self.prolong[i] @ error
        self.psi[i] = self.smooth(psi[i], rho)
        return self.psi[i]

    def smooth(psi, rho):
        pass
