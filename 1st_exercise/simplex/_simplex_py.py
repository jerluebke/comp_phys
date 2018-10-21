# -*- coding: utf-8 -*-

from numpy import std
import _simplex

_step_fortran = _simplex.simplex.step


class Step_Gen:
    """
    Wrapper for Nelder-Mead method written in Fortran

    Usage
    =====
        >>> sg = Step_Gen(simplex, func)
        >>> for step in sg:
        ...     # update

    or
        >>> sg = iter(Step_Gen(simplex, func))
        >>> current_simplex = next(sg)

    """
    required_shape = (2, 3)

    def __init__(self, simplex, func, tol=1e-19, N=1000):
        """
        Parameter
        =========
        simplex :   ndarray, ndim=2, shape=(2, 3), Fortran-contiguous
        func    :   callback function to optimize
                    f(x) = z, x like simplex, z scalar
        tol     :   convergence tolerance (default 1e-19)
        N       :   max. steps (default 1000)
        """
        self.simplex = simplex if simplex.flags.f_contiguous else simplex.T
        if self.simplex.shape != self.required_shape:
            raise ValueError('Required shape: (%d, %d); Received: (%d, %d)'
                             % (*self.required_shape, *self.simplex.shape))
        self.func = func
        self.tol  = tol
        self.N    = N

    def __iter__(self):
        self.steps      = 0
        #  self.converges  = False
        return self

    def __next__(self):
        if self.steps >= self.N or self.converges:
            raise StopIteration
        self.steps += 1
        _step_fortran(self.simplex, self.func)
        return self.simplex

    @property
    def converges(self):
        return all(std(self.simplex, axis=1) < self.tol)
