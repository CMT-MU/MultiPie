# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Coset of a model (model input "coset", util_coset): the SAMB moved by the representative g of the first coset in
the direct parent, and the weights of the eigenvalues of g on the parent structure.

- Te (P3_121, No.152): P-3m1, g = -1, n = 2; the image is P3_221 (No.154), every SAMB with +1.
- No.198 (P2_13): Pa-3, g = -1, n = 2; the image is No.198 itself, the s-p SAMB with -1.
- D2h: D4h, g = 4+[001], n = 2 (rotation).
- D4h: Oh, g = 3+[111], n = 3 (rotation).

The inversion model ("create_inversion_model"): g Z_j g^-1 of a model, with the atoms matched by g, is compared with
the SAMB Z'_k of the inversion model built by MultiPie; <Z'_k, g Z_j g^-1> must be the eps/U of "coset".
"""

import contextlib
import io

import numpy as np
import pytest

from multipie import Group
from multipie.core.material_model import MaterialModel
from multipie.util.util_coset import group_id

MODEL = {
    "TeR": {
        "group": "SG:152",
        "cell": {"a": 4.457, "c": 5.929},
        "site": {"Te": ("[0.274,0,1/3]", ["p"])},
        "bond": [("Te", "Te", [1, 2])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A1"]},
        "atomic_select": {"X": ["Q"]},
    },
    "N198": {
        "group": "SG:198",
        "cell": {"a": 4.5},
        "site": {"Co": ("[0.14,0.14,0.14]", ["s", "p"])},
        "bond": [("Co", "Co", [1])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A"]},
    },
    "D2h": {
        "group": "PG:8",
        "site": {"A": ("[0.30,0.25,0]", ["s"])},
        "bond": [("A", "A", [1, 2])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["Ag"]},
    },
    "D4h": {
        "group": "PG:15",
        "site": {"A": ("[0.30,0.30,0.25]", ["s"])},
        "bond": [("A", "A", [1, 2])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A1g"]},
    },
}

h, t = [0.5, 0.5], [2 / 3, 1 / 3]
# (parent, representative, order, kind, image group, SAMB with -1, weights).
EXPECTED = {
    "TeR": ("SG:164", "-1:[0,0,0]", 2, "inversion", "SG:154", [], [[1, 0]] * 2 + [h] * 2 + [[1, 0]] * 2 + [h] * 2 + [t] * 6),
    "N198": (
        "SG:205",
        "-1:[0,0,0]",
        2,
        "inversion",
        "SG:198",
        ["z2", "z6", "z7", "z8", "z9", "z10", "z11"],
        [[1, 0], [0, 1], [1, 0], [1, 0]] + [h] * 16,
    ),
    "D2h": ("PG:15", "4+[001]", 2, "rotation", "PG:8", [], [[1, 0], h, h]),
    "D4h": ("PG:32", "3+[111]", 3, "rotation", "PG:15", [], [[1, 0, 0], [1 / 3] * 3, [2 / 3, 1 / 6, 1 / 6]]),
}


def _model(name, **kw):
    model = {"model": name, "spinful": False, "qtdraw": {"create": False}, "pdf": {"create": False}}
    model.update(MODEL[name])
    model.update(kw)
    mm = MaterialModel()
    with contextlib.redirect_stdout(io.StringIO()):
        mm.analyze(model)
    return mm


# ==================================================
@pytest.mark.parametrize("name", list(MODEL))
def test_first_coset(name):
    parent, rep, n, kind, image, minus, weight = EXPECTED[name]
    mm = _model(name)
    (r,) = mm["coset"]
    assert (r["parent"], r["rep"], r["order"], r["kind"], r["image_group"]) == (parent, rep, n, kind, image)
    zs = list(mm["combined_id"])
    assert len(zs) == len(weight)
    for z, w in zip(zs, weight):
        # g Z_j = eps_j Z'_j (the SAMB of the image system with the same label).
        assert r["samb"][z] == [(z, -1.0 if z in minus else 1.0)], z
        assert np.allclose(r["weight"][z], w, atol=1e-8), z
    if kind == "inversion":
        site = r["site"][next(iter(mm["site"]["representative"]))]
        assert site["wyckoff"] == site["image"]


# ==================================================
def test_options():
    assert _model("D2h", coset=False)["coset"] == []
    r = _model("D2h", coset_weight=False)["coset"][0]
    assert r["weight"] is None and r["parent_structure"] is None
    # another coset of the same parent, and a coset of the parent of the parent (no weights).
    rs = _model("D2h", coset=[{"parent": "PG:15", "rep": "2[110]"}, {"parent": "PG:32", "rep": "3+[111]"}])["coset"]
    assert [r["kind"] for r in rs] == ["rotation", "rotation"]
    assert rs[0]["weight"] is not None and rs[1]["weight"] is None
    # improper g of an achiral group: g = g' m^-1.
    r = _model("D2h", coset={"parent": "PG:15", "rep": "-4+[001]"})["coset"][0]
    assert r["kind"] == "improper"
    assert all(v == [(z, 1.0)] for z, v in r["samb"].items())  # identity irrep.
    with pytest.raises(ValueError):
        _model("D2h", coset={"parent": "PG:15", "rep": "2[001]"})  # an operation of H.
    with pytest.raises(ValueError):
        _model("D2h", coset={"parent": "PG:27", "rep": "6+[001]"})  # not a parent.


# ==================================================
def test_group_data():
    assert Group("Oh").coset is None and Group("SG:221").coset is None
    info = Group.global_info()
    assert info["parent"]["SG:152"] == ["SG:164", "SG:191"]
    assert info["coset"]["SG:152"]["SG:164"][0][0] == "-1:[0,0,0]"


# ==================================================
def test_no_parent_structure():
    # S at 2h of P-6m2 (MoS2): no value of z makes the structure invariant under the direct parent (P6/mmm).
    model = {"model": "S2h", "group": "SG:187", "cell": {"a": 3.2, "c": 12}, "site": {"S": ("[1/3,2/3,0.25]", ["s"])}}
    model.update({"spinful": False, "qtdraw": {"create": False}, "pdf": {"create": False}})
    mm = MaterialModel()
    with contextlib.redirect_stdout(io.StringIO()):
        mm.analyze(model)
    (r,) = mm["coset"]
    assert r["weight"] is None and r["parent_structure"] is None and "no parent structure" in r["weight_note"]
    assert r["samb"]


# ==================================================
PAR = {"s": 1, "px": -1, "py": -1, "pz": -1}
INVERSION = {
    "TeR": MODEL["TeR"],
    "N198": MODEL["N198"],
    # a hybrid of s and p (norm 1/2 in the label space).
    "P150": {
        "group": "SG:150",
        "cell": {"a": 4.0, "c": 5.0},
        "site": {"A": ("[1/3,2/3,0.2]", ["s", "p"])},
        "bond": [("A", "A", [1, 3])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["E"]},
    },
    # 12f@6e and 12f@6f with U.
    "P208": {
        "group": "SG:208",
        "cell": {"a": 5.0},
        "site": {"A": ("[1/4,1/4,1/4]", ["s"])},
        "bond": [("A", "A", [5])],
        "SAMB_select": {"X": ["Q", "G"]},
    },
    # point group.
    "D3": {
        "group": "PG:18",
        "site": {"A": ("[0.31,0.12,0.2]", ["s", "p"])},
        "bond": [("A", "A", [1])],
        "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A1", "A2"]},
    },
}


def _inv_model(name, **kw):
    model = {"model": name, "spinful": False, "qtdraw": {"create": False}, "pdf": {"create": False}}
    model.update(INVERSION[name])
    model.update(kw)
    mm = MaterialModel()
    with contextlib.redirect_stdout(io.StringIO()):
        mm.analyze(model)
    return mm


def _matrix(mm):
    kets = mm["full_matrix"]["ket"]
    site_of = [(k[0], k[1]) for k in kets]
    orb = [k[4] for k in kets]
    pos = {
        (n, v.sublattice): np.array(v.position, dtype=float) for n, c in mm["site"]["cell"].items() for v in c if v.plus_set == 1
    }
    d = mm.get_combined_samb_matrix(fmt="value", bond=False)
    Z = {z: {(a, b, (n1, n2, n3)): complex(v) for (n1, n2, n3, a, b), v in dd.items()} for z, dd in d.items()}
    return site_of, orb, pos, Z


# ==================================================
@pytest.mark.parametrize("name", list(INVERSION))
def test_inversion_model(name):
    mm = _inv_model(name)
    mi = _inv_model(name, create_inversion_model=True)
    r = mm["coset"][0]
    info = mi["inversion_model"]
    lattice = not mi.group.is_point_group
    assert r["kind"] == "inversion" and mi["model"] == name + "_inv" and info["model"] == name
    assert group_id(mi.group) == r["image_group"] == info["image_group"]
    for site, d in r["site"].items():
        assert np.allclose(mi["site"]["representative"][site].position, d["position"], atol=1e-6)
    for bond, d in r["bond"].items():
        rep = mi["bond"]["representative"][bond]
        assert np.allclose(rep.vector, d["vector"], atol=1e-6) and np.allclose(
            np.mod(rep.center, 1) if lattice else rep.center, d["center"], atol=1e-6
        )

    M = np.array(info["map"])
    W, w = M[:3, :3], M[:3, 3]
    assert np.allclose(W, -np.eye(3))  # p orbitals: parity -1.
    site_of, orb, pos, Z = _matrix(mm)
    site_of2, orb2, pos2, Z2 = _matrix(mi)
    img = {}
    for s, p in pos.items():
        for s2, p2 in pos2.items():
            L = W @ p + w - p2
            if s2[0] == s[0] and np.allclose(L, np.rint(L), atol=1e-6):
                img[s] = (s2, np.rint(L))
    assert len(img) == len(pos)
    idx = {(site_of2[a], orb2[a]): a for a in range(len(site_of2))}

    def act(O):
        out = {}
        for (a, b, n), v in O.items():
            (sa, La), (sb, Lb) = img[site_of[a]], img[site_of[b]]
            key = (idx[(sa, orb[a])], idx[(sb, orb[b])], tuple(int(x) for x in Lb + W @ np.array(n) - La))
            out[key] = out.get(key, 0) + PAR[orb[a]] * PAR[orb[b]] * v
        return out

    assert list(Z) == list(Z2)
    for zj in Z:
        g = act(Z[zj])
        A = {zk: sum(np.conj(v) * g.get(k, 0) for k, v in Z2[zk].items()) for zk in Z2}
        A = {zk: v for zk, v in A.items() if abs(v) > 1e-8}
        expected = dict(r["samb"][zj])
        assert set(A) == set(expected), zj
        assert all(abs(A[zk] - expected[zk]) < 1e-6 for zk in A), zj
        assert np.isclose(sum(abs(c) ** 2 for c in expected.values()), 1), zj  # g is unitary.


# ==================================================
def test_inversion_model_ignored():
    model = {"model": "D2h", "spinful": False, "qtdraw": {"create": False}, "pdf": {"create": False}}
    model.update(MODEL["D2h"])
    model["create_inversion_model"] = True
    mm = MaterialModel()
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        mm.analyze(model)
    assert "is ignored" in out.getvalue()
    assert mm["model"] == "D2h" and mm["inversion_model"] is None and mm["group"] == str(Group("PG:8"))
