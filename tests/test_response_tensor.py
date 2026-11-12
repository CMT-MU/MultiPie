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

from multipie import Group
from multipie.util.util_response_tensor import convert_mp_alias, mp_alias, multipole_rank


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
        ("Qf1", 3),
        ("Mfbz", 3),
        ("Tg", 4),
        ("Qgbx", 4),
        ("Gga", 4),
    ],
)
def test_multipole_rank(tag, rank):
    assert multipole_rank(tag) == rank


# ==================================================
def test_multipole_rank_all_symbols():
    # every symbol in mp_alias starts with the letter of its rank.
    for X in "QGMT":
        for symbol in mp_alias:
            assert multipole_rank(X + symbol) == "spdfg".index(symbol[0])


# ==================================================
def test_multipole_rank_active_multipole():
    # all active multipoles in the database have a rank, and also an alias.
    tags = {tag for tp in Group.global_info()["id_set"]["MPG"]["all"] for tag in Group(tp).active_multipole}
    assert tags
    for tag in tags:
        assert 0 <= multipole_rank(tag) <= 4
        convert_mp_alias(tag)


# ==================================================
@pytest.mark.parametrize("tag", ["", "Q", "s", "Xs", "Qh", "Qdz", "qs", None, 1])
def test_multipole_rank_unknown(tag):
    with pytest.raises(ValueError, match="unknown active multipole tag"):
        multipole_rank(tag)
