# -*- coding: utf-8 -*-

import stencil_diag_lib

_stdiag = stencil_diag_lib.diag_from_stencil_module.stdiag

def stdiag(n, s):
    """
    diagonal matrix from stencil

    Wrapper for Fortran function ``diag_from_stencil``

    Parameters
    ----------
    n : input int, size of output matrix
    s : input rank-1 array, stencil

    Returns
    -------
    sqdiag : rank-2 array with bounds (n, n)

    """
    res, err = _stdiag(n, s)
    if err == -1:
        raise ValueError("size(s) = %d needs to be odd" % s.shape)
    elif err == -2:
        raise ValueError("n <= size(s)/2")
    else:
        return res
