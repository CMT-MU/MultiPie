"""
This file provides utility functions for multipole basis set.
"""
import sympy as sp
from gcoreutils.nsarray import NSArray


# ==================================================
def matrix_sum(ml):
    """
    sum of matrices.

    Args:
        ml (list): matrix [Matrix].

    Returns:
        Matrix: sum of matrices
    """
    if len(ml) == 0:
        return 0
    else:
        return sum(ml, NSArray.zeros(shape=ml[0].shape, style="matrix", fmt="sympy"))


# ==================================================
def matrix_to_dict(m):
    """
    convert matrix to list of matrix elements.

    Args:
        m (Matrix): matrix.

    Returns:
        tuple: tuple of matrix shape and elements, (shape, [(i,j,v)]).
    """
    m = NSArray(m, style="matrix", fmt="sympy")
    dim_r, dim_c = m.shape
    mat = []
    for i in range(dim_r):
        for j in range(dim_c):
            v = m[i, j]
            if v != 0:
                mat.append((i, j, str(v)))

    shape = (dim_r, dim_c)
    return shape, mat


# ==================================================
def dict_to_matrix(shape, elements):
    """
    convert dict to full matrix.

    Args:
        shape (tuple): matrix shape.
        elements (list): matrix elements, [(i,j,v)].

    Returns:
        NSArray: full matrix.
    """
    m = sp.zeros(shape[0], shape[1])
    for i, j, c in elements:
        m[i, j] = sp.sympify(c)

    fm = NSArray(str(m.tolist()), "matrix")

    return fm
