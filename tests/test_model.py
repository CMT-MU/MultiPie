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
