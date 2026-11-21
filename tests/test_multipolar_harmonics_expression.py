# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for Group.multipolar_harmonics_expression, against the multipolar harmonics and the harmonics of the point groups.
"""

import os
import runpy

import numpy as np
import pytest
import sympy as sp

from multipie import Group
from multipie.util.util_multipolar_harmonics import solid_harmonics

# the function does not depend on the group, which gives only U of the symmetry adaptation: representative groups,
# cubic and hexagonal families, Q and G in different irreps (C4v, D2d), another setting (D3h-1), complex groups, and a
# low symmetry (each component in its own irrep).
PG_ID = ["PG:32", "PG:27", "PG:31", "PG:13", "PG:14", "PG:20", "PG:37", "PG:41", "PG:46", "PG:7"]
x, y, z = sp.symbols("x y z", real=True)
FULL_CHECK = ["PG:32", "PG:27"]  # also T/M and the other internal type (Oh, D6h).


# ==================================================
def same(a, b):
    return sp.expand(sp.S(a) - sp.S(b)) == 0


# ==================================================
def internal(group):
    return "hexagonal" if group.info.hexagonal_g else "cubic"


# ==================================================
def symmetry_adapted(group, X, l, Gamma, n, s, k, xi):
    # sum_m U_{m gamma} O_{lm}^{(s)}(k) with U of the harmonics (X, l, Gamma, n) (T -> Q, M -> G).
    U = group.harmonics[({"T": "Q", "M": "G"}.get(X, X), l, Gamma, n, -1, 0, 0, "q")][1]
    ex = np.array([Group.multipolar_harmonics_expression(X, l, m, s, k, xi, internal(group)) for m in range(l, -l - 1, -1)])
    return U.T @ ex


# ==================================================
@pytest.mark.parametrize("pg", PG_ID)
def test_multipolar_harmonics(pg):
    # agrees with the multipolar harmonics of the point group (s = 1, 2, 3, l <= 3).
    group = Group(pg)
    mh = group.multipolar_harmonics
    assert mh
    for (X, l, Gamma, n, p, s, k, xi), (_, expression) in mh.items():
        ex = symmetry_adapted(group, X, l, Gamma, n, s, k, xi)
        for a, b in zip(ex, expression):
            assert same(a, b), (pg, X, l, Gamma, n, p, s, k, xi)

        if pg not in FULL_CHECK:
            continue

        # T/M: the same with the internal symbols renamed.
        XM = {"Q": "T", "G": "M"}[X]
        rename = {
            v: sp.Symbol({"Q": "T", "G": "M"}[str(v)[0]] + str(v)[1:], real=True)
            for v in set().union(*(sp.S(e).free_symbols for e in expression))
            if str(v)[0] in "QG"
        }
        ex = symmetry_adapted(group, XM, l, Gamma, n, s, k, xi)
        for a, b in zip(ex, expression):
            assert same(a, b.xreplace(rename))

        # the other internal type gives 0.
        other = {"q": "g", "g": "q"}[xi]
        for Xi in [X, XM]:
            assert all(
                Group.multipolar_harmonics_expression(Xi, l, m, s, k, other, internal(group)) == 0 for m in range(-l, l + 1)
            )


# ==================================================
@pytest.mark.parametrize("pg", ["PG:32", "PG:27", "PG:1", "PG:24", "PG:15"])
def test_harmonics_s0(pg):
    # s = 0: the harmonics of the point group (l <= 6), with G_{0} or M_{0} for the axial ones.
    group = Group(pg)
    G0 = sp.Symbol("G_{0}", real=True)
    M0 = sp.Symbol("M_{0}", real=True)
    for (X, l, Gamma, n, _, _, _, _), (expression, _) in group.harmonics.items():
        if l > 6:
            continue
        xi = "q" if X == "Q" else "g"
        for Xi, factor in ([("Q", 1), ("T", 1)] if X == "Q" else [("G", G0), ("M", M0)]):
            ex = symmetry_adapted(group, Xi, l, Gamma, n, 0, 0, xi)
            for a, b in zip(ex, expression):
                assert same(a, factor * b), (pg, X, l, Gamma, n)


# ==================================================
def test_solid_harmonics():
    # O_lm = sqrt(4pi/(2l+1)) r^l Y_lm (Racah normalization, Condon-Shortley phase).
    assert solid_harmonics(0, 0) == 1
    assert solid_harmonics(1, 0) == z
    assert sp.expand(solid_harmonics(1, 1) + (x + sp.I * y) / sp.sqrt(2)) == 0
    assert sp.expand(solid_harmonics(1, -1) - (x - sp.I * y) / sp.sqrt(2)) == 0
    assert sp.expand(solid_harmonics(2, 0) - (2 * z**2 - x**2 - y**2) / 2) == 0
    r2 = x**2 + y**2 + z**2
    for l in range(6):
        for m in range(-l, l + 1):
            o = solid_harmonics(l, m)
            # O_{l,-m} = (-1)^m O_lm^*, and harmonic (Laplacian = 0).
            assert sp.expand(solid_harmonics(l, -m) - (-1) ** m * sp.conjugate(o)) == 0
            assert sp.expand(sum(sp.diff(o, v, 2) for v in (x, y, z))) == 0
        # addition theorem: sum_m |O_lm|^2 = r^(2l).
        assert (
            sp.expand(
                sum(sp.expand(solid_harmonics(l, m) * sp.conjugate(solid_harmonics(l, m))) for m in range(-l, l + 1)) - r2**l
            )
            == 0
        )


# ==================================================
def test_examples():
    # the examples in the docs (Oh, s = 1): r.Q, Q x r, and r.G.
    f = Group.multipolar_harmonics_expression
    Qx, Qy, Qz, Gx, Gy, Gz = sp.symbols("Q_{x} Q_{y} Q_{z} G_{x} G_{y} G_{z}", real=True)
    assert sp.expand(f("Q", 0, 0, 1, 1, "q") - (x * Qx + y * Qy + z * Qz) / sp.sqrt(3)) == 0
    assert sp.expand(f("G", 1, 0, 1, 0, "q") - (y * Qx - x * Qy) / sp.sqrt(2)) == 0
    assert sp.expand(f("G", 0, 0, 1, 1, "g") - (x * Gx + y * Gy + z * Gz) / sp.sqrt(3)) == 0
    assert sp.expand(f("Q", 1, 0, 1, -1, "q") - Qz) == 0  # the polar vector itself (l + k = 0).
    # the wrong parity gives 0, and s = 0.
    assert f("Q", 0, 0, 1, 1, "g") == 0 and f("G", 1, 0, 1, 0, "g") == 0
    assert f("Q", 2, 0) == f("T", 2, 0) == solid_harmonics(2, 0)
    assert f("Q", 1, 0, x="g") == 0 and f("G", 1, 0) == 0
    assert f("G", 1, 0, x="g") == sp.Symbol("G_{0}", real=True) * z
    # hexagonal: the same for s = 1, and other symbols for s = 2.
    assert f("Q", 0, 0, 1, 1, "q", "hexagonal") == f("Q", 0, 0, 1, 1, "q")
    assert sp.Symbol("Q_{u}", real=True) in f("Q", 0, 0, 2, 2, "q", "hexagonal").free_symbols


# ==================================================
@pytest.mark.parametrize(
    "args, match",
    [
        (("A", 1, 0), "X must be Q, G, T or M"),
        (("Q", 1, 0, 1, 0, "p"), "x must be q or g"),
        (("Q", 1, 2), r"-l <= m <= l"),
        (("Q", -1, 0), r"l >= 0"),
        (("Q", 1, 0, 4, 0), r"s must be in 0..3"),
        (("Q", 1, 0, 1, 2), r"-s <= k <= s"),
        (("Q", 0, 0, 1, -1), r"l \+ k >= 0"),
        (("Q", 1.0, 0), "l must be an integer"),
        (("Q", 1, 0, True), "s must be an integer"),
        (("Q", 1, 0, 1, 0, "q", "tesseral"), "internal must be cubic or hexagonal"),
    ],
)
def test_error(args, match):
    with pytest.raises(ValueError, match=match):
        Group.multipolar_harmonics_expression(*args)


# ==================================================
def test_docs_example(capsys):
    # the example of the docs runs, and its symmetry adaptation agrees with the multipolar harmonics.
    example = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples", "multipolar_harmonics_expression.py")
    runpy.run_path(example, run_name="__main__")
    out = dict(line.split(":", 1) for line in capsys.readouterr().out.splitlines())
    assert out["wrong parity      "].strip() == "0"
    assert out["Rashba (C4v)      "].strip() == out["multipolar_harmonics"].strip()
