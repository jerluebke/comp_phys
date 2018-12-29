# -*- coding: utf-8 -*-

from cvisualise cimport *
import cython
import numpy as np
cimport numpy as np

maxlvl = 8
dim = 2

cdef class PyQuadtreeEnv:
    cdef QuadtreeEnv *this
    cdef int is_last

    def __cinit__(self, np.ndarray[np.uint_t, ndim=2, mode='c'] data not None):
        self.is_last = 0
        a, b = data.shape[0], data.shape[1]
        cdef lvl_t[::1] data_memview = data.reshape(a*b)
        self.this = qtenv_setup(&data_memview[0], a)
        if self.this is NULL:
            raise MemoryError

    def __dealloc__(self):
        if self.this is not NULL:
            qtenv_free(self.this)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def insert_next(self):
        if self.is_last:
            raise StopIteration
        res = np.empty((maxlvl*dim+1), dtype=np.uint, order='c')
        cdef lvl_t[::1] res_mv = res
        nl = qtenv_insert(self.this, &res_mv[0])
        self.is_last = qtenv_is_last(self.this);
        return res.reshape((maxlvl, dim+1))[:maxlvl-nl,:]


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
