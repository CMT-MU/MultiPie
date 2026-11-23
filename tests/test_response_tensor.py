# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for the active multipole tags of response tensors.
"""

import pytest
import sympy as sp

from multipie import Group
from multipie.util.util_response_tensor import RANK_LETTER, convert_mp_alias, mp_alias, multipole_rank


# ==================================================
@pytest.mark.parametrize(
    "tag, rank",
    [
        ("Qs", 0),
        ("Gpx", 1),
        ("Mpz", 1),
        ("Tdu", 2),
        ("Qdxy", 2),
        ("Gf3", 3),
        ("Qfs", 3),
        ("Tfc", 3),
        ("Mfbz", 3),
        ("Tg", 4),
        ("Qgbx", 4),
        ("Ggc", 4),
        ("Mgs", 4),
        ("Qhu", 5),
    ],
)
def test_multipole_rank(tag, rank):
    assert multipole_rank(tag) == rank


# ==================================================
def test_multipole_rank_all_symbols():
    # every symbol in mp_alias starts with the letter of its rank, which is also the number in its alias (if any).
    for X in "QGMT":
        for symbol, (sub, sup) in mp_alias.items():
            rank = multipole_rank(X + symbol)
            assert rank == RANK_LETTER.index(symbol[0])
            if sub[0].isdigit() and rank > 3:
                assert sub[0] == str(rank)


# ==================================================
@pytest.mark.parametrize(
    "tag, alias",
    [("Qs", "Q0"), ("Gfax", "Gx^α"), ("Qfc", "Q3c"), ("Tfy", "T3y"), ("Mgs", "M4s"), ("Qgc4", "Q4c4"), ("Ggav", "G4v^α")],
)
def test_alias_from_database(tag, alias):
    # the alias is made from the LaTeX names of the orbitals in the group database.
    assert convert_mp_alias(tag) == alias


# ==================================================
def test_response_tensor_all_from_database():
    # the tensors in the database agree with those given by response_tensor, component by component.
    g = Group("MPG:32.4.121")
    d = g.response_tensor_all("M")
    assert [k for k in d] == [k for k in d if k[0] == "M"] and len(d) == 16
    m, ex = d[("M", 1, "")]
    for i, symbol in enumerate(m[0]):
        assert sp.simplify(ex.get(symbol, sp.S(0)) - g.response_tensor("M", 1, i + 1)) == 0
    m, ex = d[("M", 2, "a")]
    for i, j in [(1, 2), (2, 3), (3, 1)]:
        assert sp.simplify(ex.get(m[i - 1, j - 1], m[i - 1, j - 1]) - g.response_tensor("M", 2, (i, j), "a")) == 0
    with pytest.raises(ValueError, match="unknown tensor type"):
        g.response_tensor_all("X")


# ==================================================
def test_response_tensor_all_of_all_types():
    # without X, the tensors of all the types; the type is the first element of the key.
    g = Group("MPG:32.4.121")
    d = g.response_tensor_all()
    assert {k[0] for k in d} == {"Q", "G", "T", "M"}
    for X in "QGTM":
        assert {k: v for k, v in d.items() if k[0] == X}.keys() == g.response_tensor_all(X).keys()
    m1, ex1 = d[("M", 2, "a")]
    m2, ex2 = g.response_tensor_all("M")[("M", 2, "a")]
    assert (m1 == m2).all() and ex1 == ex2


# ==================================================
def test_multipole_rank_active_multipole():
    # all active multipoles in the database have a rank, and also an alias.
    tags = {tag for tp in Group.global_info()["id_set"]["MPG"]["all"] for tag in Group(tp).active_multipole}
    assert tags
    for tag in tags:
        assert 0 <= multipole_rank(tag) <= 4
        convert_mp_alias(tag)


# ==================================================
@pytest.mark.parametrize("tag", ["", "Q", "s", "Xs", "Qh", "Qdz", "qs", "Qf1", "Gga", None, 1])  # f1, ga: Ver.2.
def test_multipole_rank_unknown(tag):
    with pytest.raises(ValueError, match="unknown active multipole tag"):
        multipole_rank(tag)


# ==================================================
@pytest.mark.parametrize("X", ["Q", "G", "T", "M"])
def test_rank4_t_part(X):
    # the t part of rank 4, X_{abc;d}, is fully symmetric in (a, b, c) and has 30 independent components (C1).
    import itertools

    g = Group("MPG:1.1.1")
    T = {c: sp.expand(g.response_tensor(X, 4, c, "t")) for c in itertools.product(range(1, 4), repeat=4)}
    for c in T:
        for p in itertools.permutations(c[:3]):
            assert sp.expand(T[c] - T[p + (c[3],)]) == 0, (c, p)
    symbols = sorted(set().union(*[e.free_symbols for e in T.values()]), key=str)
    assert len(symbols) == 30
    assert sp.Matrix([[sp.diff(T[c], s) for s in symbols] for c in T]).rank() == 30
