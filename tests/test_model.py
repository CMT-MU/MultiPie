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
    Nearest-neighbor graphene: E = ±sqrt(3/2) at Gamma (z2=1), and the bands touch at K.
    """
    topdir, ma = graphene_dir
    assert ma.HR is not None

    data = np.loadtxt(os.path.join(topdir, "graphene", "output", "graphene_dispersion.txt"))
    energy = data[:, 1]
    assert energy.max() == pytest.approx(np.sqrt(1.5), abs=1e-8)
    assert energy.min() == pytest.approx(-np.sqrt(1.5), abs=1e-8)
    assert np.abs(energy).min() == pytest.approx(0.0, abs=1e-8)
