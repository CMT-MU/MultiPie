# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Regression tests for parallel orthogonalization of SAMBs (parallel_map).
"""

from collections.abc import Mapping

import numpy as np

import multipie.util.util as util
from multipie.core.material_model import MaterialModel

MODEL = dict(
    model="D3h_sp",
    group="D3h",
    site={"A": ("[1,0,0]", ["s", "p"])},
    bond=[("A", "A", [1])],
    pdf={"create": False},
    qtdraw={"create": False},
)


# ==================================================
def same(x, y):
    """
    Check that two objects are identical (keys in the same order, and equal elements).
    """
    if type(x) is not type(y):
        return False
    if isinstance(x, Mapping):
        return list(x.keys()) == list(y.keys()) and all(same(x[k], y[k]) for k in x.keys())
    if isinstance(x, (list, tuple)):
        return len(x) == len(y) and all(same(i, j) for i, j in zip(x, y))
    if isinstance(x, np.ndarray):
        return x.shape == y.shape and same(list(x.ravel()), list(y.ravel()))
    return x == y


# ==================================================
def create_samb():
    mm = MaterialModel()
    mm.analyze(MODEL)
    return mm["atomic_samb"], mm["combined_samb"], mm["combined_id"]


# ==================================================
def test_parallel_samb_same_as_serial(monkeypatch):
    monkeypatch.setenv("MULTIPIE_N_JOBS", "1")
    serial = create_samb()

    # run in parallel even for this small model, and record the functions run in parallel.
    monkeypatch.setenv("MULTIPIE_N_JOBS", "2")
    monkeypatch.setattr(util, "PARALLEL_MIN_SIZE", 0)
    funcs = set()
    delayed = util.delayed

    def record(func):
        funcs.add(func.__name__)
        return delayed(func)

    monkeypatch.setattr(util, "delayed", record)
    parallel = create_samb()

    assert funcs == {"orthogonalize_samb", "gram_schmidt"}
    for s, p in zip(serial, parallel):
        assert same(s, p)
