# -*- coding: utf-8 -*-

cdef extern from "visualise.h" nogil:
    ctypedef unsigned int lvl_t
    ctypedef struct QuadtreeEnv:
        pass

    lvl_t qtenv_insert(QuadtreeEnv *this, lvl_t* res)
    int qtenv_is_last(QuadtreeEnv *this)
    QuadtreeEnv *qtenv_setup(const lvl_t *vals, size_t size)
    void qtenv_free(QuadtreeEnv *this)

#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
