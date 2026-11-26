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
"""

import contextlib
import io

import numpy as np
import pytest

from multipie import Group
from multipie.core.material_model import MaterialModel

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
