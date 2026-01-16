"""
For symmetry operation (sympy).
"""

import sympy as sp
import numpy as np

from multipie.util.util import str_to_sympy
from multipie.util.util_tag import TagSymmetryOperation

from others.tools.utils.util_symmetry_operation_spherical import representation_matrix


# ==================================================
def symmetry_operation_matrix(tag, hexagonal=False, axial=False):
    """
    Symmetry operation matrix (cartesian coordinate, conventional cell).

    Args:
        tag (str): symmetry operation tag.
        hexagonal (bool, optional): hexagonal systems ?
        axial (bool, optional): SO matrix for axial quantity ?

    Returns:
        - (ndarray) -- symmetry operation matrix in cartesian basis (3x3 or 4x4).
    """
    d = TagSymmetryOperation.parse(tag, hexagonal)
    g = representation_matrix(1, d["n"], str_to_sympy(d["axis"]), d["inversion"], d["mirror"], axial, tesseral=True)

    if not d["point_group"]:
        t = str_to_sympy(d["partial_trans"])
        g4 = np.full((4, 4), sp.S(0))
        g4[0:3, 0:3] = g
        g4[0:3, 3] = t
        g4[3, 3] = sp.S(1)
        g = g4

    return g


# ==================================================
def symmetry_operation_matrix_multipole(s, tag, hexagonal=False, axial=False, tesseral=False):
    """
    Symmetry operation matrix for multipole (lm or tesseral basis).

    Args:
        s (int): internal rank.
        tag (str): symmetry operation tag.
        hexagonal (bool, optional): hexagonal systems ?
        axial (bool, optional): SO matrix for axial quantity ?
        tesseral (bool, optional): tesseral basis ?

    Returns:
        - (ndarray) -- symmetry operation matrix in lm or tesseral basis, (2s+1 x 2s+1).

    Note:
        - ignore translational part for space group.
    """
    if s == 0:
        d = TagSymmetryOperation.parse(tag, False)
        g = sp.S(-1) if axial and (d["mirror"] or d["inversion"]) else sp.S(1)
        g = np.array(g, ndmin=2)
        return g

    d = TagSymmetryOperation.parse(tag, hexagonal)
    g = representation_matrix(s, d["n"], str_to_sympy(d["axis"]), d["inversion"], d["mirror"], axial, tesseral)
    g = np.array(g, ndmin=2)

    return g
