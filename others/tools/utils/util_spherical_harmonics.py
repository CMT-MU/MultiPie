"""
For multipolar spherical harmonics (sympy).
"""

import numpy as np
import sympy as sp
from sympy.physics.wigner import wigner_3j

from others.tools.data.data_definition import max_s_rank

TOL = 1e-11


# ==================================================
def _create_unitary_matrix_tesseral(l):
    """
    Create unitary matrix to internal basis.

    Args:
        l (int): rank.

    Returns:
        - (ndarray) -- unitary matrix, [[sympy]] (2l+1,2l+1), <m|gamma> (m=l, l-1, ..., -l).
    """
    U = np.full((2 * l + 1, 2 * l + 1), sp.S(0))
    for m in range(l, 0, -1):
        U[2 * (l - m), l - m] = (-1) ** m / sp.sqrt(2)
        U[2 * (l - m), l + m] = 1 / sp.sqrt(2)
        U[2 * (l - m) + 1, l - m] = -sp.I * ((-1) ** m) / sp.sqrt(2)
        U[2 * (l - m) + 1, l + m] = sp.I / sp.sqrt(2)
    U[2 * l, l] = sp.S(1)

    return U.T


# ==================================================
# U-matrix for internal (tesseral) basis <m|alpha> in each rank.
# alpha in order of [(l,c),(l,s),(l-1,c),...,(0,c)].
_internal_U = [_create_unitary_matrix_tesseral(s) for s in range(max_s_rank + 1)]


# ==================================================
def is_transform_real(U):
    """
    Check if transformed harmonics is real.

    Args:
        U (array-like): coefficient vector, <m|gamma>, for gamma-harmonics (m=l,l-1,...,-l).

    Returns:
        - (bool) -- is real ?
    """
    U = np.asarray(U)
    l = (U.shape[0] - 1) // 2
    if U.dtype == object:
        sgn = np.asarray([sp.S(1) if m % 2 == 0 else sp.S(-1) for m in range(l, -l - 1, -1)])
    else:
        sgn = np.asarray([1.0 if m % 2 == 0 else -1.0 for m in range(l, -l - 1, -1)])
    Ub = sgn * U[::-1].conjugate()

    if U.dtype == object:
        return np.array_equal(U, Ub)
    else:
        return np.allclose(U, Ub, TOL, 0.01 * TOL)


# ==================================================
def _Olm(l, m, rv=None):
    """
    Spherical harmonics (scaled), sqrt(4pi/(2l+1)) r^l Ylm.

    Args:
        l (int): rank.
        m (int): component (-l<= m <=l).
        rv (array-like, optional): (x,y,z) symbol/str/sympy/float values.

    Returns:
        - (sympy) -- spherical harmonics.
    """
    x, y, z = sp.symbols("x y z", real=True)
    r = sp.sqrt(x**2 + y**2 + z**2)
    if m == 0:
        olm = r**l * sp.legendre(l, z / r)
    else:
        rp = sp.sqrt(x**2 + y**2)
        s = m > 0
        m = sp.Abs(m)
        olm = sp.sqrt(sp.factorial(l - m) / sp.factorial(l + m)) * r**l * sp.assoc_legendre(l, m, z / r) / rp ** (m)
        if s:
            olm *= (x + sp.I * y) ** m
        else:
            olm *= (-1) ** m * (x - sp.I * y) ** m

    olm = sp.factor(olm)

    if rv is not None:
        olm = olm.subs({x: rv[0], y: rv[1], z: rv[2]})

    olm = sp.radsimp(olm)

    return olm


# ==================================================
def Olm(l, m, s=0, k=0, rv=None, factor=False):
    """
    Multipolar spherical harmonics (scaled), Olm^{(s)}(k).

    Args:
        l (int): rank.
        m (int): component (-l<= m <= l).
        s (int, optional): internal rank.
        k (int, optional): internal component (-s<= k <=s).
        rv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        factor (bool, optional): factorize resultant expression ? (sometime slow)

    Returns:
        - (sympy or complex) -- multipolar spherical harmonics.

    Note:
        - internal basis (s>0) in order of [(s,"c"),(s,"s"),(s-1,"c"),...,(0,"c")].
        - for s > 0, return ndarray of [sympy(2s+1)].
    """
    assert s < len(_internal_U)

    if rv is not None:
        rv = np.asarray(rv)

    if s == 0:
        molm = _Olm(l, m, rv)
    else:
        ev = _internal_U[s].conjugate()
        fv = [sp.I**i for i in range(2 * s + 1)]  # k=-s, -s+1, ..., s.

        molm = np.full(2 * s + 1, sp.S(0))
        n1 = max(m - l - k, -s)
        n2 = min(m + l + k, s)
        for n in range(n1, n2 + 1):
            w3j = wigner_3j(l + k, l, s, m - n, -m, n)
            olm = _Olm(l + k, m - n, rv)
            molm += w3j * olm * ev[s - n]
        molm *= fv[k + s] * (-1) ** (l + m) * sp.sqrt(2 * l + 1)

    if any(bool(i.free_symbols) for i in np.atleast_1d(molm)) or rv is None or rv.dtype == object:
        if factor:
            molm = np.vectorize(lambda i: sp.factor(sp.radsimp(i)))(molm)
        else:
            molm = np.vectorize(lambda i: sp.expand(sp.radsimp(i)))(molm)
    else:
        molm = np.vectorize(complex)(molm)

    if molm.size == 1:
        molm = molm.item()

    return molm


# ==================================================
def create_harmonics_set(l, s=0, k=0, rv=None, sv=None, factor=False):
    """
    Create harmonics set (m=l, l-1, ..., -l).

    Args:
        l (int): rank.
        s (int, optional): internal rank.
        k (int, optional): internal component (-s<= k <=s).
        rv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        sv (array-like, optional): internal variable.
        factor (bool, optional): factorize resultant expression ? (sometime slow)

    Returns:
        - (ndarray) -- multipolar spherical harmonics set.
    """
    olm = np.asarray([Olm(l, m, s, k, rv, factor) for m in range(l, -l - 1, -1)])

    if sv is not None and s > 0:
        basis = create_internal_basis(s, sv, factor=factor)
        olm = olm @ basis

    return olm


# ==================================================
def transform_harmonics(Olm_set, U):
    """
    Transform harmonics by O_{lg} = sum_m U_{m,g}O_{lm}.

    Args:
        Olm_set (array-like): Olm [m, point, internal comp.].
        U (array-like): U [m,g].

    Returns:
        - (ndarray) -- transformed harmonics [g, point, internal comp.].
    """
    Olm_set = np.asarray(Olm_set)
    U = np.asarray(U)

    Olg = np.vectorize(lambda i: sp.factor(sp.expand(i)))(U.T @ Olm_set)

    if Olg.size == 1:
        Olg = Olg.item()

    return Olg


# ==================================================
def create_internal_basis(s, sv=None, factor=False):
    """
    Create internal basis.

    Args:
        s (int, optional): internal rank.
        sv (array-like, optional): (x,y,z) symbol/str/sympy/float values.
        factor (bool, optional): factorize resultant expression ? (sometime slow)

    Returns:
        - (ndarray) -- internal basis.
    """
    return transform_harmonics(create_harmonics_set(s, rv=sv, factor=factor), _internal_U[s])
