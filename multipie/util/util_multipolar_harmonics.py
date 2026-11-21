# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Multipolar spherical harmonics with internal degree of freedom (sympy).

The definition is that of the PG multipolar harmonics of MultiPieGen (util_multipolar_harmonics.Olm):

    O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} sqrt(2l+1) sum_n (l+k, l, s; m-n, -m, n) O_{l+k,m-n}(r) e_n,

where O_{lm}(r) = sqrt(4pi/(2l+1)) r^l Y_lm(r) is ALWAYS A POLAR FUNCTION of r, and e_n is the internal basis of rank s.
"""

import functools

import numpy as np
import sympy as sp
from sympy.physics.wigner import wigner_3j

# ==================================================
max_s_rank = 3  # max. internal rank (the internal basis symbols are defined for s <= 3).
internal_crystal = ["cubic", "hexagonal"]  # internal basis.
_x, _y, _z = sp.symbols("x y z", real=True)


# ==================================================
@functools.lru_cache(maxsize=None)
def solid_harmonics(l, m):
    """
    Spherical harmonics (Racah normalization, Condon-Shortley phase), O_lm = sqrt(4pi/(2l+1)) r^l Y_lm.

    Args:
        l (int): rank.
        m (int): component, -l <= m <= l.

    Returns:
        - (sympy) -- homogeneous polynomial in x, y, z (real symbols).
    """
    c = sp.sqrt(sp.factorial(l + m) * sp.factorial(l - m))
    a, b = -(_x + sp.I * _y) / 2, (_x - sp.I * _y) / 2
    ex = sp.S(0)
    for q in range(max(0, -m), (l - m) // 2 + 1):
        p = q + m
        s = l - p - q
        ex += a**p * b**q * _z**s / (sp.factorial(p) * sp.factorial(q) * sp.factorial(s))

    return sp.expand(c * ex)


# ==================================================
def internal_type(X, s, k):
    """
    Internal type q (polar) or g (axial) allowed for the multipolar harmonics (X, s, k), whose spatial part is polar.

    Args:
        X (str): Q/G/T/M.
        s (int): internal rank.
        k (int): internal component.

    Returns:
        - (str) -- "q" or "g".
    """
    even = (s + k) % 2 == 0
    return "q" if (X in ["Q", "T"]) == even else "g"


# ==================================================
@functools.lru_cache(maxsize=None)
def _coefficient(l, m, s, k, U):
    """
    Coefficient (polynomial in x, y, z) of each internal basis function.

    Args:
        l (int): rank.
        m (int): component.
        s (int): internal rank.
        k (int): internal component.
        U (tuple): <n|alpha> of the internal basis (row: n = s, ..., -s), as a tuple of tuples.

    Returns:
        - (tuple) -- coefficient of each internal basis function alpha.
    """
    n_s = 2 * s + 1
    acc = [sp.S(0)] * n_s
    for n in range(max(m - l - k, -s), min(m + l + k, s) + 1):
        w3j = wigner_3j(l + k, l, s, m - n, -m, n)
        if w3j == 0:
            continue
        wo = w3j * solid_harmonics(l + k, m - n)
        for a in range(n_s):
            u = U[s - n][a]
            if u != 0:
                acc[a] += wo * sp.conjugate(u)
    pre = sp.I ** (k + s) * (-1) ** (l + m) * sp.sqrt(2 * l + 1)

    return tuple(sp.expand(sp.radsimp(pre * e)) for e in acc)


# ==================================================
def multipolar_harmonics_expression(X, l, m, s, k, x, symbols, U):
    """
    Expression of multipolar spherical harmonics with internal degree of freedom.

    Args:
        X (str): Q/G/T/M.
        l (int): rank.
        m (int): component, -l <= m <= l.
        s (int): internal rank, 0 <= s <= 3.
        k (int): internal component, -s <= k <= s, l + k >= 0.
        x (str): internal type, q (polar) or g (axial).
        symbols (list): internal basis symbols of (s, x) (Q/G + orbital name).
        U (ndarray): <n|alpha> of the internal basis (row: n = s, ..., -s, column: alpha).

    Returns:
        - (sympy) -- expression, 0 for the combination of X and x of the wrong parity.

    Notes:
        - the spatial part is always polar, so that x is fixed by X and the parity of s+k.
        - for X = T/M, the internal symbols Q/G are renamed T/M.
        - for s = 0, the internal symbol is multiplied only for G_{0} and M_{0}.
    """
    if X not in ["Q", "G", "T", "M"]:
        raise ValueError(f"X must be Q, G, T or M, '{X}' is given.")
    if x not in ["q", "g"]:
        raise ValueError(f"x must be q or g, '{x}' is given.")
    for name, v in [("l", l), ("m", m), ("s", s), ("k", k)]:
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)):
            raise ValueError(f"{name} must be an integer, '{v}' is given.")
    if not 0 <= s <= max_s_rank:
        raise ValueError(f"s must be in 0..{max_s_rank}, '{s}' is given.")
    if l < 0 or abs(m) > l:
        raise ValueError(f"l >= 0 and -l <= m <= l are required, (l, m) = ({l}, {m}) is given.")
    if abs(k) > s or l + k < 0:
        raise ValueError(f"-s <= k <= s and l + k >= 0 are required, (l, s, k) = ({l}, {s}, {k}) is given.")

    if internal_type(X, s, k) != x:
        return sp.S(0)

    U = tuple(tuple(sp.S(v) for v in row) for row in np.asarray(U, dtype=object))
    coeff = _coefficient(int(l), int(m), int(s), int(k), U)

    symbols = list(symbols)
    if X in ["T", "M"]:
        rename = {"Q": "T", "G": "M"}
        symbols = [sp.Symbol(rename[str(v)[0]] + str(v)[1:], real=True) for v in symbols]
    if s == 0 and x == "q":  # Q_{0} and T_{0} are not multiplied.
        symbols = [sp.S(1)]

    return sp.expand(sum(c * v for c, v in zip(coeff, symbols)))
