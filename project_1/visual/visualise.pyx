# -*- coding: utf-8 -*-

from cvisualise cimport *
import numpy as np
cimport numpy as np

cdef class PyQuadtreeEnv:
    cdef QuadtreeEnv *this

    def __cinit__(self, np.ndarray[np.uint, ndim=2, mode='c'] data not None):
        self.is_last = 0
        a, b = data.shape
        self.this = qtenv_setup(data.reshape((a*b,)), a)
        if self.thisptr is NULL:
            raise MemoryError

    def __dealloc__(self):
        if self.this is not NULL:
            qtenv_f(self.this)

    cpdef np.ndarray insert_next(self):
        if self.is_last:
            raise StopIteration
        res = np.empty((maxlvl*dim+1), dtype=np.uint, order='c')
        nl = qtenv_insert(self.this, res)
        self.is_last = qtenv_is_last(self.this);
        return res.reshape((maxlvl, dim+1))[:maxlvl-nl,:]


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
