# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Regression tests for model creation (MaterialModel) and analysis (ModelAnalyzer), using graphene.
"""

import os

import numpy as np
import pytest

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples")


# ==================================================
def read_example(filename):
    """
    Read example input in docs/src/examples.
    """
    return read_dict(os.path.abspath(os.path.join(EXAMPLE_DIR, filename)))


# ==================================================
@pytest.fixture(scope="module")
def graphene_dir(tmp_path_factory):
    """
    Create and analyze the graphene example in docs/src/examples (once per module).
    """
    model = read_example("graphene_in.py")
    model["pdf"] = {"create": False}
    model["qtdraw"] = {"create": False}
    control = read_example("graphene_ctrl.py")
    control["samb"]["samb_figure"] = False
    control["grid"] = (6, 6, 1)

    topdir = str(tmp_path_factory.mktemp("graphene"))
    cwd = os.getcwd()
    try:
        create_model(model, topdir=topdir)
        os.chdir(topdir)
        ma = ModelAnalyzer(topdir)
        ma.analyze(control)
    finally:
        os.chdir(cwd)

    return topdir, ma


# ==================================================
def test_model_files(graphene_dir):
    topdir, _ = graphene_dir
    assert os.path.isfile(os.path.join(topdir, "graphene", "graphene.pkl"))


# ==================================================
def test_number_of_samb(graphene_dir):
    _, ma = graphene_dir
    assert len(ma.model["combined_id"]) == 20


# ==================================================
def test_dispersion(graphene_dir):
    """
    Nearest-neighbor graphene (z2=1, hopping t=1/sqrt(6)) along Gamma-M-K-Gamma:
    E = ±3t at Gamma, ±t at M, and both bands touch at E = 0 at K.
    """
    topdir, ma = graphene_dir
    assert ma.HR is not None

    data = np.loadtxt(os.path.join(topdir, "graphene", "output", "graphene_dispersion.txt"))
    # rows are ordered as band 1 (3 segments), then band 2 (3 segments); each segment includes both end points.
    n_band = 2
    assert len(data) % (3 * n_band) == 0
    energy = data[:, 1].reshape(n_band, -1)
    n_seg = energy.shape[1] // 3
    idx = {"Gamma": 0, "M": n_seg - 1, "K": 2 * n_seg - 1}

    t = 1 / np.sqrt(6)
    expected = {"Gamma": [-3 * t, 3 * t], "M": [-t, t], "K": [0.0, 0.0]}
    for point, i in idx.items():
        assert np.sort(energy[:, i]) == pytest.approx(expected[point], abs=1e-8), point


# ==================================================
_CLUSTER_ORDER_SCRIPT = """
from multipie import MaterialModel
mm = MaterialModel()
mm.analyze({"model": "si", "group": 227, "cell": {"a": 5.43}, "site": {"Si": ("[1/8,1/8,1/8]", "s")},
            "bond": [("Si", "Si", [1, 2, 3])], "pdf": {"create": False}, "qtdraw": {"create": False}})
print(list(mm["cluster_samb"].keys()))
"""


def test_cluster_order_does_not_depend_on_hash_seed():
    # clusters of the same multiplicity (48a@48f and 48b@16d for Si) were ordered as in a set, which depends on
    # PYTHONHASHSEED, so that the numbering y# changed from run to run.
    import subprocess
    import sys

    import multipie

    parent = os.path.dirname(os.path.dirname(os.path.abspath(multipie.__file__)))
    orders = set()
    for seed in ["1", "4"]:
        env = os.environ | {"PYTHONHASHSEED": seed, "PYTHONPATH": parent}
        result = subprocess.run([sys.executable, "-c", _CLUSTER_ORDER_SCRIPT], capture_output=True, text=True, env=env)
        assert result.returncode == 0, result.stderr
        orders.add(result.stdout.strip().splitlines()[-1])
    assert orders == {"['8a', '16a@16c', '48a@48f', '48b@16d']"}


# ==================================================
@pytest.mark.parametrize("version, ok", [("2.5.0", False), ("2.3.4", False), ("3.0.0.dev0", True), ("3.1.0", True)])
def test_model_of_older_version(tmp_path, version, ok):
    # a model created by Ver.2 (another group database) must be created again.
    from multipie import MaterialModel

    mm = MaterialModel(topdir=str(tmp_path))
    mm.analyze(read_example("graphene_in.py") | {"pdf": {"create": False}, "qtdraw": {"create": False}})
    mm["version"] = version
    mm.save()
    mm = MaterialModel(topdir=str(tmp_path))
    if ok:
        mm.load("graphene")
        assert mm["version"] == version
    else:
        with pytest.raises(ValueError, match=f"created by MultiPie {version}.*Create the model again"):
            mm.load("graphene")


# ==================================================
def _chain_bonds(neighbor):
    from multipie import MaterialModel

    mm = MaterialModel()
    mm.analyze(
        {
            "model": "chain",
            "group": 123,
            "cell": {"a": 5.0, "c": 0.98},
            "site": {"H": ("[0,0,0]", "s")},
            "bond": [("H", "H", neighbor)],
            "pdf": {"create": False},
            "qtdraw": {"create": False},
        }
    )
    return [(v.neighbor, round(v.distance, 4)) for v in mm["bond"]["representative"].values()]


def test_neighbor_shells_beyond_initial_search_range():
    # a chain with c = 0.98 A: the initial search range (+-2 cells) reaches only 1.96 A along c, so the shells at
    # 2.94, 3.92 and 4.90 A were missing and neighbor 3 was the 5.0 A in-plane bond. The search grid is now
    # enlarged per bond until it contains every site up to the requested shell.
    expected = [(1, 0.98), (2, 1.96), (3, 2.94), (4, 3.92), (5, 4.9), (6, 5.0), (7, 5.0951), (8, 5.3704)]
    assert _chain_bonds(8) == expected


def test_neighbor_shells_of_a_slab_with_a_large_c():
    # graphene with c = 1000 A: the initial search range holds 15 in-plane shells, and the next one found is the next
    # layer at 1000 A. The grid grows step by step to the in-plane shells, not at once to 1000 A (millions of sites).
    from multipie import MaterialModel

    mm = MaterialModel()
    model = {"model": "g", "group": 191, "cell": {"a": 2.435, "c": 1000}, "site": {"C": ("[2/3,1/3,0]", "pz")}}
    mm.analyze(model | {"bond": [("C", "C", list(range(1, 21)))], "pdf": {"create": False}, "qtdraw": {"create": False}})
    dist = [v.distance for v in mm["bond"]["representative"].values()]
    assert len(set(round(d, 4) for d in dist)) == 20 and max(dist) < 10


def test_neighbor_from_bond_input():
    # the number of neighbors is given by "bond" only (no upper limit such as the former "max_neighbor" = 10).
    bonds = _chain_bonds(12)
    assert [n for n, _ in bonds] == list(range(1, 13))


@pytest.mark.parametrize("key, value", [("max_neighbor", 10), ("search_cell_range", (-2, 3, -2, 3, -2, 3))])
def test_removed_model_keys(key, value):
    from multipie import MaterialModel

    model = {"model": "chain", "group": 123, "site": {"H": ("[0,0,0]", "s")}, "bond": [("H", "H", 1)], key: value}
    with pytest.raises(ValueError, match=f"'{key}' is removed in MultiPie Ver.3"):
        MaterialModel().analyze(model)
