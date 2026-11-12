# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for the definitions of the local operators (Sx/Sy/Sz, Lx/Ly/Lz, Qu/Qv/Qyz/Qxz/Qxy), as documented in post_processing.md.
"""

import numpy as np
import pytest
import sympy as sp

from multipie.util.util_model_analyzer import create_all_local_operator

# index of the orbital in the tesseral basis of each rank.
ORBITAL = {
    1: ["px", "py", "pz"],
    2: ["dv", "dxy", "dxz", "dyz", "du"],
    3: ["fc", "fs", "fbz", "f3", "fx", "fy", "faz"],
}


# ==================================================
@pytest.fixture(scope="module")
def local_operator():
    d = create_all_local_operator()
    return {rank: {name: np.array(sp.Matrix(m).evalf(), dtype=complex) for name, m in dd.items()} for rank, dd in d.items()}


# ==================================================
@pytest.mark.parametrize("rank", [1, 2, 3])
def test_angular_momentum(local_operator, rank):
    op = local_operator[rank]
    Lx, Ly, Lz = op["Lx"], op["Ly"], op["Lz"]
    dim = 2 * rank + 1

    for A in (Lx, Ly, Lz):
        assert np.allclose(A, A.conj().T)
    assert np.allclose(Lx @ Ly - Ly @ Lx, 1j * Lz)
    assert np.allclose(Ly @ Lz - Lz @ Ly, 1j * Lx)
    assert np.allclose(Lz @ Lx - Lx @ Lz, 1j * Ly)
    assert np.allclose(Lx @ Lx + Ly @ Ly + Lz @ Lz, rank * (rank + 1) * np.eye(dim))
    assert np.allclose(np.linalg.eigvalsh(Lz), np.arange(-rank, rank + 1))


# ==================================================
@pytest.mark.parametrize("rank, orbital, value", [(1, "pz", 2 / 5), (2, "du", 2 / 7), (3, "faz", 4 / 15)])
def test_quadrupole_m0(local_operator, rank, orbital, value):
    # <l 0|(3z^2-r^2)/(2r^2)|l 0> = l(l+1)/((2l-1)(2l+3)).
    i = ORBITAL[rank].index(orbital)
    assert local_operator[rank]["Qu"][i, i] == pytest.approx(value)


# ==================================================
def test_quadrupole_p(local_operator):
    # for p orbitals (sqrt(3) x/r, ...), <x^2/r^2> = 3/5, <y^2/r^2> = 1/5, and 3<x^2 y^2/r^4> = 1/5.
    op = local_operator[1]
    px, py, pz = 0, 1, 2
    c = np.sqrt(3) / 5
    assert op["Qu"][px, px] == pytest.approx(-1 / 5)
    assert op["Qv"][px, px] == pytest.approx(c)
    assert op["Qv"][py, py] == pytest.approx(-c)
    assert op["Qxy"][px, py] == pytest.approx(c)
    assert op["Qxz"][px, pz] == pytest.approx(c)
    assert op["Qyz"][py, pz] == pytest.approx(c)


# ==================================================
@pytest.mark.parametrize("rank", [1, 2, 3])
def test_quadrupole_angular_function(local_operator, rank):
    # Q(r) are quadratic in r/|r|, so that they equal the operator equivalents up to a common factor (Wigner-Eckart).
    op = local_operator[rank]
    Lx, Ly, Lz = op["Lx"], op["Ly"], op["Lz"]
    L2 = Lx @ Lx + Ly @ Ly + Lz @ Lz
    s3 = np.sqrt(3) / 2
    equivalent = {
        "Qu": (3 * Lz @ Lz - L2) / 2,
        "Qv": s3 * (Lx @ Lx - Ly @ Ly),
        "Qyz": s3 * (Ly @ Lz + Lz @ Ly),
        "Qxz": s3 * (Lx @ Lz + Lz @ Lx),
        "Qxy": s3 * (Lx @ Ly + Ly @ Lx),
    }
    factor = -2 / ((2 * rank - 1) * (2 * rank + 3))
    for name, A in equivalent.items():
        assert np.allclose(op[name], factor * A), name


# ==================================================
def test_rank0_and_spin(local_operator):
    for name, m in local_operator[0].items():
        assert np.allclose(m, 0), name

    S = local_operator["s"]
    sigma = {
        "Sx": np.array([[0, 1], [1, 0]]),
        "Sy": np.array([[0, -1j], [1j, 0]]),
        "Sz": np.array([[1, 0], [0, -1]]),
    }
    for name, s in sigma.items():
        assert np.allclose(S[name], s / 2), name
