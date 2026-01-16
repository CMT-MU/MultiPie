"""
For symmetry operation for spherical system (sympy).
"""

import sympy as sp
import numpy as np
from sympy.physics.quantum.spin import Rotation

from others.tools.utils.util_spherical_harmonics import _create_unitary_matrix_tesseral


# ==================================================
def _rotation_matrix(theta, axis):
    """
    Rotation matrix (cartesian coordinate).

    Parameters:
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (ndarray) -- 3x3 rotation matrix, [[sympy]].
    """
    n = np.asarray([sp.S(axis[0]), sp.S(axis[1]), sp.S(axis[2])])
    norm = sp.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    n /= norm

    K = np.array([[sp.S(0), -n[2], n[1]], [n[2], sp.S(0), -n[0]], [-n[1], n[0], sp.S(0)]])
    I = np.full((3, 3), sp.S(0))
    I[0, 0] = I[1, 1] = I[2, 2] = sp.S(1)
    R = I + sp.sin(theta) * K + (1 - sp.cos(theta)) * K @ K

    return R


# ==================================================
def _axis_to_euler_zyz(theta, axis):
    """
    Convert (theta, axis) to zyz-Euler angle, (alpha, beta, gamma).

    Args:
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (tuple) -- angles, ( alpha(sympy), beta(sympy), gamma(sympy) ).
    """
    R = _rotation_matrix(theta, axis)

    if R[2, 2] < 1 and R[2, 2] > -1:  # |Rzz| ne 1.
        alpha = sp.atan2(R[1, 2], R[0, 2])
        beta = sp.acos(R[2, 2])
        gamma = sp.atan2(R[2, 1], -R[2, 0])
    else:  # gimbal lock.
        alpha = sp.atan2(-R[0, 1], R[1, 1])
        beta = sp.S(0) if R[2, 2] == sp.S(1) else sp.pi
        gamma = sp.S(0)

    return alpha, beta, gamma


# ==================================================
def _wigner_D_matrix(j, theta, axis):
    """
    Wigner D matrix.

    Args:
        j (sympy): angule momentum, (integer or half integer).
        theta (sympy): rotation angle (radian).
        axis (array-like): rotation axis (cartesian coordinate, no need for normalization).

    Returns:
        - (ndarray) -- rotation matrix, <m1|D|m2>, (m1,m2=l,l-1,...,-l) [O' = D * O].
    """

    alpha, beta, gamma = _axis_to_euler_zyz(theta, axis)
    rng = [sp.S(m) / 2 for m in range(2 * j, -2 * j - 1, -2)]
    pairs = [(m1, m2) for m1 in rng for m2 in rng]
    D = np.asarray([Rotation.D(j, m1, m2, alpha, beta, gamma).simplify() for m1, m2 in pairs])
    D = D.reshape(2 * j + 1, 2 * j + 1)

    return D


# ==================================================
def representation_matrix(l, n, axis, inversion=False, mirror=False, axial=False, tesseral=False):
    """
    Representation matrix of symmetry operation for spherical harmonics, l.

    Args:
        l (int): rank.
        n (int): n-fold rotation.
        axis (list): rotation or mirror axis (cartesian coordinate, no need for normalization).
        inversion (bool, optional): with inversion ?
        mirror (bool, optional): with mirror ?
        axial (bool, optional): for axial quantity ?
        tesseral (bool, optional): tesseral basis ?

    Returns:
        - (ndarray) -- representation matrix, <lm|g|lm'> or <alpha|g|beta> (tesseral=True) [alpha, beta=cl,sl,cl-1, ..., c0].
    """
    if mirror:
        n = 2
        inversion = True

    angle = 2 * sp.pi / n
    g = _wigner_D_matrix(l, angle, axis)

    p = l
    if axial:
        p += 1
    if inversion:
        g = (-1) ** p * g

    if tesseral:
        U = _create_unitary_matrix_tesseral(l)
        g = np.vectorize(sp.expand_complex)((U.conjugate().T) @ g @ U)

    return g
