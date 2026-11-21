# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Operator expression of atomic multipoles, following H. Kusunose, R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. 89, 104704 (2020).

Spinless (Eqs. 1-4), where l = -i (r x nabla) and nabla acts only on O_lm(r) = sqrt(4pi/(2l+1)) r^l Y_lm:

    Q_lm = O_lm,
    M_lm = 1/2 * 2/(l+1) [(nabla O_lm).l + l.(nabla O_lm)],
    T_lm = 1/2 * 2/((l+1)(l+2)) [(nabla O_lm).(r x l) - (l x r).(nabla O_lm)],
    G_lm = 1/2 * 4i/((l+1)^2 (l+2)) [(nabla O_lm).l l^2 - l^2 l.(nabla O_lm)].

Spinful (Eq. 13 and Table I):

    X_lm^(s)(k) = i^(s+k) (-1)^(l+m) sqrt(2l+1) sum_n (l+k, l, s; m-n, -m, n) X'_(l+k,m-n) sigma_(s,n),

with sigma_(0,0) = 1, sigma_(1,0) = sigma_z, sigma_(1,+-1) = -+(sigma_x +- i sigma_y)/sqrt(2), and X' of Table I.

The operators are noncommutative symbols: the position x, y, z, the orbital angular momentum l_x, l_y, l_z
(dimensionless), and the Pauli matrices sigma_x, sigma_y, sigma_z (the spin is sigma/2), with
[l_i, l_j] = i eps_ijk l_k, [sigma_i, sigma_j] = 2i eps_ijk sigma_k, [l_i, x_j] = i eps_ijk x_k, [x_i, x_j] = 0, and
sigma commuting with l and the position.
"""

import functools

import numpy as np
import sympy as sp
from sympy.physics.wigner import wigner_3j

from multipie.util.util_atomic_multipole import qmtg_map_spinful

# ==================================================
# commutative x, y, z for the polynomials O_lm.
_x, _y, _z = sp.symbols("x y z", real=True)

# noncommutative operators.
op_position = sp.symbols("x y z", commutative=False)
op_angular_momentum = sp.symbols("l_x l_y l_z", commutative=False)
op_pauli = sp.symbols("sigma_x sigma_y sigma_z", commutative=False)
_to_operator = dict(zip((_x, _y, _z), op_position))
_order = {v: no for no, v in enumerate(op_position)}


# ==================================================
@functools.lru_cache(maxsize=None)
def _solid_harmonics(l, m):
    """
    Spherical harmonics (Racah normalization, Condon-Shortley phase), O_lm = sqrt(4pi/(2l+1)) r^l Y_lm.

    Args:
        l (int): rank.
        m (int): component, -l <= m <= l.

    Returns:
        - (sympy) -- homogeneous polynomial in x, y, z (commutative, real symbols).
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
def _as_operator(poly):
    """
    Polynomial in commutative x, y, z as that in the position operators (the order in a monomial is x, y, z).

    Args:
        poly (sympy): polynomial in x, y, z.

    Returns:
        - (sympy) -- noncommutative expression.
    """
    poly = sp.Poly(sp.expand(poly), _x, _y, _z)
    X, Y, Z = op_position
    return sp.Add(*[c * X**a * Y**b * Z**d for (a, b, d), c in poly.terms()])


# ==================================================
def _dot(a, b):
    return sp.Add(*[ai * bi for ai, bi in zip(a, b)])


# ==================================================
def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


# ==================================================
def canonical_form(ex):
    """
    Canonical form of an operator expression, by the relations that do not change the operator.

    Args:
        ex (sympy): expression in the noncommutative operators.

    Returns:
        - (sympy) -- expression, where a run of consecutive position operators in each term is sorted as x, y, z
          (they commute), and the Pauli matrices are moved to the right end (they commute with the others).

    Note:
        - no other relation (such as the commutation relations of l and the position) is used.
    """
    terms = []
    for term in sp.Add.make_args(sp.expand(ex)):
        c, nc = term.args_cnc()
        factors = []
        for f in nc:
            base, n = f.as_base_exp()
            factors += [base] * int(n) if n.is_Integer and n > 0 else [f]
        pauli = [f for f in factors if f in op_pauli]
        other = [f for f in factors if f not in op_pauli]
        out, run = [], []
        for f in other + [None]:
            if f in _order:
                run.append(f)
                continue
            out += sorted(run, key=_order.get)
            run = []
            if f is not None:
                out.append(f)
        terms.append(sp.Mul(*c) * sp.Mul(*(out + pauli)))

    return sp.expand(sp.Add(*terms))


# ==================================================
@functools.lru_cache(maxsize=None)
def orbital_multipole_operator(X, l, m):
    """
    Spinless (orbital) multipole operator, X_lm^(orb) of Eqs. (1)-(4).

    Args:
        X (str): Q/M/T/G.
        l (int): rank.
        m (int): component.

    Returns:
        - (sympy) -- operator expression (canonical form).

    Note:
        - G_lm with l < 2 is 0 (the operator vanishes since l^2 commutes with l).
    """
    O = _solid_harmonics(l, m)
    if X == "Q":
        return canonical_form(_as_operator(O))

    grad = [_as_operator(sp.diff(O, v)) for v in (_x, _y, _z)]
    L = list(op_angular_momentum)
    r = list(op_position)
    if X == "M":
        ex = sp.Rational(1, l + 1) * (_dot(grad, L) + _dot(L, grad))
    elif X == "T":
        ex = sp.Rational(1, (l + 1) * (l + 2)) * (_dot(grad, _cross(r, L)) - _dot(_cross(L, r), grad))
    else:  # "G"
        if l < 2:
            return sp.S(0)
        L2 = _dot(L, L)
        ex = 2 * sp.I / ((l + 1) ** 2 * (l + 2)) * (_dot(grad, L) * L2 - L2 * _dot(L, grad))

    return canonical_form(ex)


# ==================================================
def _pauli(s, n):
    """
    sigma_(s,n): sigma_(0,0) = 1, sigma_(1,0) = sigma_z, sigma_(1,+-1) = -+(sigma_x +- i sigma_y)/sqrt(2).
    """
    sx, sy, sz = op_pauli
    if s == 0:
        return sp.S(1)
    if n == 0:
        return sz
    return -n * (sx + n * sp.I * sy) / sp.sqrt(2)


# ==================================================
@functools.lru_cache(maxsize=None, typed=True)
def atomic_multipole_operator(X, l, m, s=0, k=0):
    """
    Operator expression of atomic multipole X_lm^(s)(k) (JPSJ 89, 104704 (2020)).

    Args:
        X (str): Q/G/T/M.
        l (int): rank.
        m (int): component, -l <= m <= l.
        s (int, optional): 0 (spinless) or 1 (spinful).
        k (int, optional): -s <= k <= s, l + k >= 0.

    Returns:
        - (sympy) -- operator expression in the noncommutative symbols x, y, z, l_x, l_y, l_z, sigma_x, sigma_y, sigma_z
          (canonical form, see canonical_form).
    """
    if X not in ["Q", "G", "T", "M"]:
        raise ValueError(f"X must be Q, G, T or M, '{X}' is given.")
    for name, v in [("l", l), ("m", m), ("s", s), ("k", k)]:
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)):
            raise ValueError(f"{name} must be an integer, '{v}' is given.")
    l, m, s, k = int(l), int(m), int(s), int(k)
    if s not in [0, 1]:
        raise ValueError(f"s must be 0 or 1, '{s}' is given.")
    if l < 0 or abs(m) > l:
        raise ValueError(f"l >= 0 and -l <= m <= l are required, (l, m) = ({l}, {m}) is given.")
    if abs(k) > s or l + k < 0:
        raise ValueError(f"-s <= k <= s and l + k >= 0 are required, (l, s, k) = ({l}, {s}, {k}) is given.")

    if s == 0:
        return orbital_multipole_operator(X, l, m)

    Xo = qmtg_map_spinful(X, s, k)
    ex = sp.S(0)
    for n in range(-s, s + 1):
        if abs(m - n) > l + k:
            continue
        w3j = wigner_3j(l + k, l, s, m - n, -m, n)
        if w3j == 0:
            continue
        ex += w3j * orbital_multipole_operator(Xo, l + k, m - n) * _pauli(s, n)
    ex *= sp.I ** (s + k) * (-1) ** (l + m) * sp.sqrt(2 * l + 1)

    return canonical_form(ex)
