"""
For Gram-Schmidt orthonormalization.
"""

import numpy as np
import sympy as sp


# ==================================================
def gram_schmidt(vec, n_max=None):
    """
    Gram-Schmidt orthogonalization (sympy).

    Args:
        vec (array-like): list of vectors to be orthogonalized, [[sympy]].
        n_max (int, optional): max. of nonzero basis.

    Returns:
        - (ndarray) -- list of nonzero orthogonalized vectors, [[sympy]].
        - (ndarray) -- indices of nonzero vectors, [int].
    """
    ip = lambda x, y: np.vectorize(sp.expand)(np.vdot(x, y))
    norm = lambda x: sp.sqrt(ip(x, x))

    vec = np.asarray(vec, dtype=object)
    vec = np.vectorize(sp.expand)(vec)

    if vec.ndim < 2:
        nv = norm(vec)
        if nv != 0:
            vec /= nv
        vec = np.vectorize(sp.radsimp)(vec)
        return vec, np.array([0])

    if n_max is None:
        n_max = len(vec)

    ortho_vec = []
    nonzero = []
    for no, v in enumerate(vec):
        for u in ortho_vec:
            v -= sp.expand(sp.radsimp(ip(v, u) / ip(u, u))) * u
        v = np.vectorize(sp.simplify)(v)
        if not (v == 0).all():
            ortho_vec.append(v)
            nonzero.append(no)
        if len(ortho_vec) == n_max:
            break

    nonzero = np.asarray(nonzero)
    ret = []
    for v in ortho_vec:
        nv = norm(v)
        if nv != 0:
            v /= nv
        v = np.vectorize(sp.radsimp)(v)
        ret.append(v)
    ortho_vec = np.asarray(ret, dtype=object)

    return ortho_vec, nonzero


# ==================================================
def gram_schmidt_float(vec, n_max=None, TOL=1e-8):
    """
    Gram-Schmidt orthogonalization for real vectors by converting to float.

    Args:
        vec (array-like): list of vectors to be orthogonalized, [[sympy]].
        n_max (int, optional): max. of nonzero basis.
        TOL (float, optional): absolute tolerance.

    Returns:
        - (ndarray) -- list of nonzero orthogonalized vectors, [[float]].
        - (ndarray) -- indices of nonzero vectors, [int].
    """
    norm = lambda x: np.sqrt(np.dot(x, x))

    vec = np.asarray(vec, dtype=float)

    if vec.ndim < 2:
        nv = norm(vec)
        if nv > TOL:
            vec /= nv
        return vec, np.array([0])

    if n_max is None:
        n_max = len(vec)

    ortho_vec = []
    nonzero = []
    for no, v in enumerate(vec):
        for u in ortho_vec:
            v -= np.dot(v, u) / np.dot(u, u) * u
        nv = norm(v)
        if nv > TOL:
            ortho_vec.append(v / nv)
            nonzero.append(no)
        if len(ortho_vec) == n_max:
            break

    ortho_vec = np.asarray(ortho_vec)
    nonzero = np.asarray(nonzero)

    return ortho_vec, nonzero


# ==================================================
def gram_schmidt_complex(vec, n_max=None, TOL=1e-8):
    """
    Gram-Schmidt orthogonalization for complex vectors by converting to complex.

    Args:
        vec (array-like): list of vectors to be orthogonalized, [[sympy]].
        n_max (int, optional): max. of nonzero basis.
        TOL (float, optional): absolute tolerance.

    Returns:
        - (ndarray) -- list of nonzero orthogonalized vectors, [[complex]].
        - (ndarray) -- indices of nonzero vectors, [int].
    """
    norm = lambda x: np.sqrt(np.vdot(x, x))

    vec = np.asarray(vec, dtype=complex)

    if vec.ndim < 2:
        nv = norm(vec)
        if nv > TOL:
            vec /= nv
        return vec, np.array([0])

    if n_max is None:
        n_max = len(vec)

    ortho_vec = []
    nonzero = []
    for no, v in enumerate(vec):
        for u in ortho_vec:
            v -= np.vdot(v, u) / np.vdot(u, u) * u
        nv = norm(v)
        if nv > TOL:
            ortho_vec.append(v / nv)
            nonzero.append(no)
        if len(ortho_vec) == n_max:
            break

    ortho_vec = np.asarray(ortho_vec)
    nonzero = np.asarray(nonzero)

    return ortho_vec, nonzero
