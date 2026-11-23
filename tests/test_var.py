# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests of the relation between zj and the atomic variables (info/<model>_var.py, ModelAnalyzer.get_var).
"""

import os

import pytest
import sympy as sp

from multipie import ModelAnalyzer
from multipie.core.cmd import create_model
from multipie.util.util import read_dict
from multipie.util.util_model_analyzer import short_orbital_name

OPT = {"SAMB_select": {"Gamma": "IR"}, "pdf": {"create": False}, "qtdraw": {"create": False}}
MODEL = {
    "graphene": {"group": 191, "cell": {"c": 4}, "site": {"C": ("[1/3,2/3,0]", "pz")}, "bond": [("C", "C", [1, 2])]},
    "lgs": {"group": 221, "site": {"A": ("[0,0,0]", "p")}, "bond": [("A", "A", [1])], "spinful": True},
    "jml": {"group": 221, "site": {"A": ("[0,0,0]", ["(3/2,p)", "(3/2,d)"])}, "bond": [("A", "A", [1])], "spinful": True},
    "AB": {
        "group": 221,
        "site": {"A": ("[0,0,0]", "s"), "B": ("[1/2,1/2,1/2]", ["px", "py", "pz"])},
        "bond": [("A", "B", [1])],
    },
}


# ==================================================
@pytest.fixture
def analyzer(tmp_path):
    """
    Create a model and analyze it (the var file is written), returning (ModelAnalyzer, var file).
    """

    def run(name):
        create_model({"model": name} | MODEL[name] | OPT, topdir=str(tmp_path), verbose=False)
        ma = ModelAnalyzer(str(tmp_path))
        ma.analyze({"samb": {"model": name}, "output": {"dispersion": {"k_path": None}}})
        var = read_dict(os.path.join(str(tmp_path), name, "info", f"{name}_var.py"))
        return ma, var

    return run


# ==================================================
def names(ma):
    # the names of the variables in get_var.
    return {
        str(s)
        for d in ma.get_var(ma["samb"]["matrix_info"]).values()
        for zj, ex in d.items()
        if zj != "bond"
        for s in ex.free_symbols
    }


# ==================================================
def test_short_orbital_name():
    assert short_orbital_name("pz", 1, "lg") == "pz"
    assert short_orbital_name("(px,u)", 1, "lgs") == "pxu" and short_orbital_name("(dxy,d)", 2, "lgs") == "dxyd"
    assert short_orbital_name("(3/2,1/2)", 1, "jml") == "p(3/2,1/2)"
    assert short_orbital_name("(3/2,-1/2)", 2, "jml") == "d(3/2,-1/2)"


# ==================================================
def test_graphene(analyzer):
    ma, var = analyzer("graphene")
    assert var == {
        "C;C_001_1": {"bond": ("C(2)", "C(1)"), "z2": "sqrt(6)*g_{pz,pz}"},
        "C;C_002_1": {"bond": ("C(1)", "C(1)"), "z3": "2*sqrt(3)*g_{pz,pz}"},
    }
    assert list(var["C;C_001_1"])[0] == "bond"
    full = ma.get_var(ma["samb"]["matrix_info"])
    assert full["C;C_001_1"]["bond"] == ("C(2)", "C(1)")


# ==================================================
def test_head_tail(analyzer):
    # a bond between different atoms (the bond B;A: tail B, head A): (head, tail) = (A, B), and g_{head orbital, tail
    # orbital}.
    ma, var = analyzer("AB")
    assert list(var) == ["B;A_001_1"]
    assert var["B;A_001_1"]["bond"] == ("A(1)", "B(1)")
    assert names(ma) == {f"g_{{s,{o}}}" for o in ["px", "py", "pz"]}


# ==================================================
def test_lgs(analyzer):
    ma, var = analyzer("lgs")
    assert var["A;A_001_1"]["bond"] == ("A(1)", "A(1)")
    orbital = {f"p{c}{s}" for c in "xyz" for s in "ud"}
    for n in names(ma):
        a, b = n[3:-1].split(",")
        assert a in orbital and b in orbital, n
    assert "g_{pxu,pxu}" in names(ma) and "," not in "".join(orbital)


# ==================================================
def test_jml(analyzer):
    # p and d with the same (j, m) are distinguished.
    ma, var = analyzer("jml")
    n = names(ma)
    prefixes = {x[0] for s in n for x in s[3:-1].replace("),", ")|").split("|")}
    assert prefixes == {"p", "d"}
    assert any("d(3/2," in s for s in n) and any("p(3/2," in s for s in n)
