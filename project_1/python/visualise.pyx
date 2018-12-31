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
    cdef np.ndarray sorted_indices

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, np.ndarray[np.uint32_t, ndim=2, mode='c'] data not None):
        self.is_last = 0
        a, b = data.shape[0], data.shape[1]
        self.sorted_indices = np.empty(a, dtype=np.uint32, order='c')
        cdef unsigned int[::1] data_memview = data.reshape(a*b)
        cdef unsigned int[::1] sorted_memview = self.sorted_indices
        self.this = qtenv_setup(&data_memview[0], a, &sorted_memview[0])
        if self.this is NULL:
            raise MemoryError

    def __dealloc__(self):
        if self.this is not NULL:
            qtenv_free(self.this)

    def get_sorted(self):
        return self.sorted_indices

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def insert_next(self):
        if self.is_last:
            raise StopIteration
        res = np.empty(maxlvl*(dim+1), dtype=np.double, order='c')
        cdef double[::1] res_mv = res
        nl = qtenv_insert(self.this, &res_mv[0])
        self.is_last = qtenv_is_last(self.this);
        #  res = np.asarray(res_mv)
        return (res.reshape((maxlvl, dim+1)))[:nl]

#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
