import cython
import numpy as np
cimport numpy as np

cdef extern from "csearch.h":
    void search_naive( const unsigned int *, size_t, double )
    void search_fast( const unsigned int *, size_t, double )
    void search_fastfast( const unsigned int *, size_t, double )


def search_naive_py(np.ndarray[np.uint32_t, ndim=2, mode='c'] data not None,
                    double r_sq):
    a, b = data.shape[0], data.shape[1]
    cdef unsigned int [::1] data_mv = data.reshape(a*b)
    search_naive(&data_mv[0], a, r_sq)

def search_fast_py(np.ndarray[np.uint32_t, ndim=2, mode='c'] data not None,
                    double r_sq):
    a, b = data.shape[0], data.shape[1]
    cdef unsigned int [::1] data_mv = data.reshape(a*b)
    search_fast(&data_mv[0], a, r_sq)

def search_fastfast_py(np.ndarray[np.uint32_t, ndim=2, mode='c'] data not None,
                    double r_sq):
    a, b = data.shape[0], data.shape[1]
    cdef unsigned int [::1] data_mv = data.reshape(a*b)
    search_fast(&data_mv[0], a, r_sq)
