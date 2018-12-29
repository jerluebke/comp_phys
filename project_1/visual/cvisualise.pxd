# -*- coding: utf-8 -*-

cdef extern from "cvisualise.h" nogil:
    ctypedef struct QuadtreeEnv:
        pass

    unsigned int qtenv_insert(QuadtreeEnv *this, unsigned int* res)
    int qtenv_is_last(QuadtreeEnv *this)
    QuadtreeEnv *qtenv_setup(const unsigned int *vals, size_t size)
    void qtenv_free(QuadtreeEnv *this)

#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
