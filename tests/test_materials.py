# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tight-binding models of materials (tests/data/materials), checked by properties that do not depend on the choice of
the SAMB basis: the number of SAMB, and the symmetry of the band energies of a Hamiltonian with random weights.

- SrVO3: cubic perovskite (Pm-3m), V t2g and O p.
- MoS2: monolayer 2H-MoS2 (P-6m2, D3h), Mo d and S p.
- Te: trigonal Te (P3_121, chiral D3), spinful p.

These checks hold for any (complete) SAMB basis of the same model, so that they stay the same when the group database
changes the numbering, signs, or the basis of SAMB.
"""

import os

import numpy as np
import pytest

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict
from multipie.util.util_model_analyzer import fourier_r_to_k

MATERIAL_DIR = os.path.join(os.path.dirname(__file__), "data", "materials")

# number of SAMB of the identity irrep (Q and G, for the Hamiltonian).
N_SAMB = {"SrVO3": 10, "MoS2": 30, "Te": 66}

# degeneracies of the energies at the Γ point (by the irreps of the point group, and the Kramers degeneracy).
GAMMA_DEGENERACY = {
    "SrVO3": [3, 3, 3, 3],  # T2g (V t2g), 2 T1u + T2u (O p).
    "MoS2": [1, 1, 1, 2, 2, 2, 2],  # 2 A1' + A2'' + 2 E' + 2 E''.
    "Te": [2] * 9,  # Kramers pairs.
}


# ==================================================
@pytest.fixture(scope="module", params=list(N_SAMB))
def material(request, tmp_path_factory):
    """
    Create a material model, and analyze it with random SAMB weights (fixed seed).
    """
    name = request.param
    model = read_dict(os.path.join(MATERIAL_DIR, f"{name}_in.py"))
    model["pdf"] = {"create": False}
    model["qtdraw"] = {"create": False}

    topdir = str(tmp_path_factory.mktemp(name))
    create_model(model, topdir=topdir)

    rng = np.random.default_rng(1)
    parameter = {f"z{j}": float(rng.uniform(-1, 1)) for j in range(1, N_SAMB[name] + 1)}
    control = {"samb": {"model": name, "parameter": parameter}, "grid": (4, 4, 4)}
    cwd = os.getcwd()
    try:
        os.chdir(topdir)
        ma = ModelAnalyzer(topdir)
        ma.analyze(control)
    finally:
        os.chdir(cwd)

    return name, ma


def energies(ma, kv):
    """
    Band energies at k points (fractional), Ek[k,n].
    """
    return np.concatenate([Ek for Ek, _ in ma.eigen_chunks(np.asarray(kv, dtype=float), vectors=False)])


def degeneracy(E, tol=1e-8):
    """
    Degeneracies of sorted energies.
    """
    deg, n = [], 1
    for e1, e2 in zip(E[:-1], E[1:]):
        if abs(e2 - e1) < tol:
            n += 1
        else:
            deg.append(n)
            n = 1
    deg.append(n)
    return sorted(deg)


# ==================================================
def test_number_of_samb(material):
    name, ma = material
    combined_id = ma.model["combined_id"]  # z_j: (symbol, UniqueSAMBType, (X, l, Gamma, n, p, s, k, x), component).
    assert len(combined_id) == N_SAMB[name]
    irrep = list(ma.model.group.character["table"].keys())[0]  # identity irrep.
    assert all(tag[2] == irrep and tag[0] in "QG" for _, _, tag, _ in combined_id.values())


# ==================================================
def test_hamiltonian_is_hermitian(material):
    _, ma = material
    HR, atom, tb_gauge = ma.get_hk_data()
    kv = np.random.default_rng(2).uniform(-0.5, 0.5, (4, 3))
    Hk = fourier_r_to_k(HR, atom, kv, tb_gauge)
    assert np.abs(Hk).max() > 1e-3
    assert np.allclose(Hk, Hk.conj().transpose(0, 2, 1))


# ==================================================
def test_energies_have_the_symmetry_of_the_group(material):
    # E(k) = E(kW) for every rotation W of the space group, and E(k) = E(-k) by time reversal.
    _, ma = material
    so = ma.model.group.symmetry_operation["fractional"]
    W = [np.asarray(s, dtype=float)[:3, :3] for s in so]
    for k in np.random.default_rng(3).uniform(-0.5, 0.5, (3, 3)):
        E = energies(ma, [k])[0]
        for w in W:
            assert np.allclose(energies(ma, [k @ w])[0], E, atol=1e-8)
        assert np.allclose(energies(ma, [-k])[0], E, atol=1e-8)


# ==================================================
def test_degeneracy_at_gamma(material):
    name, ma = material
    E = energies(ma, [[0, 0, 0]])[0]
    assert degeneracy(np.sort(E)) == GAMMA_DEGENERACY[name]


# ==================================================
def test_generic_k_point_has_no_accidental_degeneracy(material):
    # at a generic k point, only the Kramers degeneracy (spinful with inversion) or none remains.
    name, ma = material
    E = energies(ma, [[0.123, 0.0456, 0.0789]])[0]
    assert set(degeneracy(np.sort(E))) == {1}
