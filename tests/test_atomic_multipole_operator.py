# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for Group.atomic_multipole_operator: the matrix elements of the operators agree with those of
atomic_multipole_matrix (JPSJ 89, 104704 (2020), with the radial factor <r^k> = 1).

The orbital |L M> is the homogeneous polynomial r^L Y_LM, on which the position acts as the multiplication and
l = -i (r x nabla) as the differential operator; the matrix element is the integral over the unit sphere.
"""

import functools
import math
import os
import runpy

import numpy as np
import pytest
import sympy as sp

from multipie import Group
from multipie.util.util_atomic_multipole import atomic_multipole_matrix
from multipie.util.util_atomic_multipole_operator import (
    op_position,
    op_angular_momentum,
    op_pauli,
    canonical_form,
    _solid_harmonics,
)

x, y, z = sp.symbols("x y z", real=True)
X_, Y_, Z_ = op_position
LX, LY, LZ = op_angular_momentum
SX, SY, SZ = op_pauli
PAULI = {None: np.eye(2), SX: np.array([[0, 1], [1, 0]]), SY: np.array([[0, -1j], [1j, 0]]), SZ: np.diag([1, -1])}
HALF = sp.Rational(1, 2)


# ==================================================
# a polynomial in x, y, z is a dict {(a, b, c): complex coefficient of x^a y^b z^c} (fast, instead of sympy).
def _add(f, g, c=1.0):
    h = dict(f)
    for mon, v in g.items():
        h[mon] = h.get(mon, 0.0) + c * v
    return h


def _mul(f, i):
    # multiply by x_i.
    return {tuple(e + (n == i) for n, e in enumerate(mon)): v for mon, v in f.items()}


def _diff(f, i):
    # derivative with respect to x_i.
    return {tuple(e - (n == i) for n, e in enumerate(mon)): v * mon[i] for mon, v in f.items() if mon[i] > 0}


def act(op, f):
    # act with a single operator on a polynomial f: the position multiplies, and l = -i (r x nabla).
    if op in (X_, Y_, Z_):
        return _mul(f, (X_, Y_, Z_).index(op))
    i = (LX, LY, LZ).index(op)
    j, k = (i + 1) % 3, (i + 2) % 3
    return {mon: -1j * v for mon, v in _add(_mul(_diff(f, k), j), _mul(_diff(f, j), k), -1.0).items()}


@functools.lru_cache(maxsize=None)
def polynomial(L, M):
    # O_LM as a dict.
    poly = sp.Poly(_solid_harmonics(L, M), x, y, z)
    return {mon: complex(c) for mon, c in poly.terms()}


# ==================================================
def split(ex):
    # {pauli (None for identity): [(coefficient, [orbital operators])]}.
    res = {}
    for term in sp.Add.make_args(sp.expand(ex)):
        if term == 0:
            continue
        c, nc = term.args_cnc()
        factors = []
        for f in nc:
            base, n = f.as_base_exp()
            factors += [base] * int(n)
        pauli = [f for f in factors if f in op_pauli]
        assert len(pauli) <= 1 and (not pauli or factors[-1] == pauli[0])
        res.setdefault(pauli[0] if pauli else None, []).append((sp.Mul(*c), [f for f in factors if f not in op_pauli]))
    return res


# ==================================================
@functools.lru_cache(maxsize=None)
def sphere(a, b, c):
    # integral of x^a y^b z^c over the unit sphere.
    if a % 2 or b % 2 or c % 2:
        return 0.0
    g = math.gamma
    return 2 * g((a + 1) / 2) * g((b + 1) / 2) * g((c + 1) / 2) / g((a + b + c + 3) / 2)


# ==================================================
def orbital_matrix(terms, L1, L2):
    # <L1 M1| A |L2 M2> with Y_LM = sqrt((2L+1)/4pi) O_LM on the unit sphere (M = L, ..., -L).
    norm = math.sqrt((2 * L1 + 1) * (2 * L2 + 1)) / (4 * math.pi)
    terms = [(complex(c), tuple(ops)) for c, ops in terms]
    mat = np.zeros((2 * L1 + 1, 2 * L2 + 1), dtype=complex)
    for j, M2 in enumerate(range(L2, -L2 - 1, -1)):
        memo = {(): polynomial(L2, M2)}  # suffix of the operators -> the result on |L2 M2>.

        def apply(ops):
            if ops not in memo:
                memo[ops] = act(ops[0], apply(ops[1:]))
            return memo[ops]

        f = {}
        for c, ops in terms:
            f = _add(f, apply(ops), c)
        for i, M1 in enumerate(range(L1, -L1 - 1, -1)):
            o1 = polynomial(L1, M1)
            mat[i, j] = norm * sum(
                np.conj(v1) * v * sphere(a1 + a, b1 + b, c1 + c) for (a1, b1, c1), v1 in o1.items() for (a, b, c), v in f.items()
            )
    return mat


# ==================================================
def operator_matrix(ex, L1, L2):
    # matrix in the lms basis, (L1, M1, s1) x (L2, M2, s2) with s = 1/2, -1/2 for each M.
    mat = np.zeros((2 * (2 * L1 + 1), 2 * (2 * L2 + 1)), dtype=complex)
    for pauli, terms in split(ex).items():
        mat += np.kron(orbital_matrix(terms, L1, L2), PAULI[pauli])
    return mat


# ==================================================
def reference_matrix(X, l, m, s, k, L1, L2):
    basis1 = [(L1, M, S) for M in range(L1, -L1 - 1, -1) for S in (HALF, -HALF)]
    basis2 = [(L2, M, S) for M in range(L2, -L2 - 1, -1) for S in (HALF, -HALF)]
    return np.array(
        [[complex(sp.N(atomic_multipole_matrix(X, l, m, s, k, b1, b2, "lms"), 20)) for b2 in basis2] for b1 in basis1]
    )


# ==================================================
BLOCKS = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3)]


# ==================================================
@pytest.mark.parametrize("L1, L2", BLOCKS)
@pytest.mark.parametrize("X", ["Q", "M", "T", "G"])
def test_spinless(X, L1, L2):
    # all the ranks and components that can be nonzero in the block (and one more rank, which vanishes).
    for l in range(L1 + L2 + 2):
        for m in range(-l, l + 1):
            ex = Group.atomic_multipole_operator(X, l, m)
            mat = operator_matrix(ex, L1, L2)
            assert np.allclose(mat, reference_matrix(X, l, m, 0, 0, L1, L2), atol=1e-12), (X, l, m, L1, L2)


# ==================================================
@pytest.mark.parametrize("L1, L2", [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)])
@pytest.mark.parametrize("X", ["Q", "M", "T", "G"])
def test_spinful(X, L1, L2):
    for k in (-1, 0, 1):
        for l in range(max(0, -k), L1 + L2 + 2):
            for m in sorted({l, 0, -l, min(1, l)}):
                ex = Group.atomic_multipole_operator(X, l, m, 1, k)
                mat = operator_matrix(ex, L1, L2)
                assert np.allclose(mat, reference_matrix(X, l, m, 1, k, L1, L2), atol=1e-12), (X, l, m, k, L1, L2)


# ==================================================
def test_examples():
    # Sect. 2.4 of the paper.
    f = Group.atomic_multipole_operator
    s3 = sp.sqrt(3)
    assert f("Q", 0, 0) == 1
    assert f("Q", 0, 0, 1, 1) == canonical_form((LX * SX + LY * SY + LZ * SZ) / s3)  # spin-orbit coupling.
    assert f("M", 0, 0, 1, 1) == canonical_form((X_ * SX + Y_ * SY + Z_ * SZ) / s3)  # M monopole r.sigma.
    assert f("M", 1, 0) == LZ  # orbital angular momentum.
    assert f("M", 1, 0, 1, -1) == SZ  # spin.
    assert f("Q", 1, 0) == Z_
    assert f("Q", 2, 0) == canonical_form((2 * Z_**2 - X_**2 - Y_**2) / 2)
    assert f("G", 1, 0) == 0 and f("T", 0, 0) == 0 and f("M", 0, 0) == 0
    # T dipole t = [(r x l) - (l x r)]/6 (Eq. 27): T^(0)_z = t_z.
    t_z = ((X_ * LY - Y_ * LX) - (LX * Y_ - LY * X_)) / 6
    assert f("T", 1, 0) == canonical_form(t_z)
    # G monopole (t.sigma)/sqrt(3) (Eq. 27).
    t = [((Y_ * LZ - Z_ * LY) - (LY * Z_ - LZ * Y_)) / 6, ((Z_ * LX - X_ * LZ) - (LZ * X_ - LX * Z_)) / 6, t_z]
    assert f("G", 0, 0, 1, 1) == canonical_form((t[0] * SX + t[1] * SY + t[2] * SZ) / s3)


# ==================================================
def test_canonical_form():
    # position operators commute with each other, and sigma with the others; l and the position do not.
    assert canonical_form(Z_ * X_ * LX * Y_ * X_) == X_ * Z_ * LX * X_ * Y_
    assert canonical_form(SX * X_ * LX + X_ * LX * SX) == 2 * X_ * LX * SX
    assert canonical_form(SX * LX * X_) == LX * X_ * SX
    assert canonical_form(LX * X_ - X_ * LX) != 0
    assert canonical_form(Y_**2 * X_) == X_ * Y_**2


# ==================================================
@pytest.mark.parametrize(
    "args, match",
    [
        (("A", 1, 0), "X must be Q, G, T or M"),
        (("Q", 1, 2), r"-l <= m <= l"),
        (("Q", -1, 0), r"l >= 0"),
        (("Q", 1, 0, 2, 0), r"s must be 0 or 1"),
        (("Q", 1, 0, 0, 1), r"-s <= k <= s"),
        (("Q", 0, 0, 1, -1), r"l \+ k >= 0"),
        (("Q", 1.0, 0), "l must be an integer"),
        (("Q", 1, 0, True), "s must be an integer"),
    ],
)
def test_error(args, match):
    with pytest.raises(ValueError, match=match):
        Group.atomic_multipole_operator(*args)


# ==================================================
def test_docs_example(capsys):
    example = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples", "atomic_multipole_operator.py")
    runpy.run_path(example, run_name="__main__")
    assert "spin-orbit" in capsys.readouterr().out
