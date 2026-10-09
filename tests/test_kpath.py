"""
Regression tests for the primitive cell of centred lattices and the default k path.
"""

import os
import warnings

import numpy as np
import pytest
import seekpath
from scipy.spatial.transform import Rotation

from multipie import Group, MaterialModel
from multipie.core.cmd import create_model
from multipie.core import model_analyzer
from multipie.core.model_analyzer import ModelAnalyzer, _join_kpath, _kpath_structure
from multipie.util.util_crystal import convert_to_primitive, convert_to_primitive_vector, get_cell_info

# name: (group, cell, site position).
CASES = {
    "bcc": (229, {"a": 3.0}, "[0,0,0]"),
    "bct": (139, {"a": 3.0, "c": 8.0}, "[0,0,1/3]"),
    "fcc": (225, {"a": 4.0}, "[1/4,1/4,1/4]"),
    "c_centred": (65, {"a": 3.0, "b": 4.5, "c": 6.0}, "[0,1/4,1/2]"),
    "a_centred": (38, {"a": 3.0, "b": 4.5, "c": 6.0}, "[0,0,0.1]"),
    "rhombohedral": (166, {"a": 3.0, "c": 20.0}, "[0,0,0.1]"),
    "hexagonal": (191, {"a": 2.5, "c": 6.0}, "[1/3,2/3,0]"),
    "cubic": (221, {"a": 3.0}, "[0,0,0]"),
    "monoclinic": (10, {"a": 5.3, "b": 4.1, "c": 3.0, "beta": 100.0}, "[0,0,0]"),  # seekpath changes the basis.
    "triclinic": (2, {"a": 3.0, "b": 4.0, "c": 5.0, "alpha": 80.0, "beta": 100.0, "gamma": 110.0}, "[0,0,0]"),
}


# generic rotation of the Cartesian frame.
ROTATION = Rotation.from_euler("zyx", [23, 41, -17], degrees=True).as_matrix()


def cell_vectors(a, b, c, alpha, beta, gamma):
    """
    Lattice vectors [a1,a2,a3] (rows) of a general cell.
    """
    ca, cb, cc = np.cos(np.radians([alpha, beta, gamma]))
    sc = np.sin(np.radians(gamma))
    s = np.sqrt(1 - ca * ca - cb * cb - cc * cc + 2 * ca * cb * cc)
    return np.array([[a, 0, 0], [b * cc, b * sc, 0], [c * cb, c * (ca - cb * cc) / sc, c * s / sc]])


def check_kpath(k_path, k_point, B, structure):
    """
    Check k path and k points (fractional in B) against the standardised cell of seekpath.
    """
    info = seekpath.get_path(structure)
    # Cartesian k in the original frame = (Cartesian k in the standardised frame) @ rotation_matrix.
    B_std = np.asarray(info["reciprocal_primitive_lattice"]) @ np.asarray(info["rotation_matrix"])
    assert k_path == _join_kpath(info["path"])
    assert set(k_point) == {"Γ" if k == "GAMMA" else k for k in info["point_coords"]}
    for label, k in k_point.items():
        k_std = np.asarray(info["point_coords"]["GAMMA" if label == "Γ" else label])
        assert np.allclose(np.asarray(k, dtype=float) @ B, k_std @ B_std, atol=1e-6), label
    return info


# ==================================================
@pytest.fixture(autouse=True)
def unchanged_cwd():
    """
    Check that the working directory is not changed by MultiPie (restored for the following tests).
    """
    cwd = os.getcwd()
    yield
    changed = os.getcwd()
    os.chdir(cwd)
    assert changed == cwd


@pytest.fixture(scope="module")
def models(tmp_path_factory):
    topdir = str(tmp_path_factory.mktemp("kpath"))
    cwd = os.getcwd()
    try:
        for name, (group, cell, pos) in CASES.items():
            model = {"model": name, "group": group, "cell": cell, "site": {"X": (pos, "s")}, "bond": [("X", "X", [1])]}
            create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir=topdir)
    finally:
        changed = os.getcwd()
        os.chdir(cwd)
    assert changed == cwd  # the working directory is not changed (checked here, before the autouse fixture).
    return topdir


# ==================================================
@pytest.mark.parametrize("name", CASES)
def test_unit_vector_primitive(models, name):
    # primitive fractional coordinates of the sites with the primitive vectors give the Cartesian positions,
    # modulo primitive lattice vectors (position_primitive is in the home cell).
    mm = MaterialModel(topdir=models)
    mm.load(name)
    A = np.asarray(mm["unit_vector"], dtype=float)
    Ap = np.asarray(mm["unit_vector_primitive"], dtype=float)
    assert np.allclose(Ap, convert_to_primitive_vector(mm.group.info.lattice, A))
    assert abs(np.linalg.det(Ap)) * len(mm.group.symmetry_operation.get("plus_set", [0])) == pytest.approx(abs(np.linalg.det(A)))
    for sites in mm["site"]["cell"].values():
        for site in sites:
            d = (
                np.asarray(site.position_primitive, dtype=float) @ Ap - np.asarray(site.position, dtype=float) @ A
            ) @ np.linalg.inv(Ap)
            assert np.allclose(d, np.rint(d), atol=1e-6)


# ==================================================
@pytest.mark.parametrize("name", CASES)
def test_default_kpath(models, name):
    ma = ModelAnalyzer(models)
    ma.analyze({"samb": {"model": name, "parameter": {"z1": 1.0}}, "grid": (10, 10, 10)})
    disp = ma["output"]["dispersion"]
    assert os.path.isfile(os.path.join(models, name, "output", f"{name}_dispersion.txt"))

    # the same path and Cartesian k points as for the standardised cell of seekpath.
    k_point = {k: [float(i) for i in v.strip("[]").split(",")] for k, v in disp["k_point"].items()}
    structure = _kpath_structure(ma.model.group, ma["info"]["A"])
    info = check_kpath(disp["k_path"], k_point, np.asarray(ma["info"]["B"]), structure)
    assert info["spacegroup_number"] == CASES[name][0]


# ==================================================
@pytest.mark.parametrize("no", range(1, 231))
def test_kpath_structure_all_space_groups(no):
    # for oblique and rotated cells of all space groups, the structure for seekpath has the space group of the group
    # in its primitive cell (not a supercell), and the k points are given in the reciprocal basis of the primitive cell.
    group = Group(no)
    if group.info.crystal == "triclinic":
        A = cell_vectors(3.0, 4.1, 5.3, 80.0, 100.0, 110.0)
    elif group.info.crystal == "monoclinic":
        A = cell_vectors(3.0, 4.1, 5.3, 90.0, 100.0, 90.0)
    else:
        A = get_cell_info(group.info.crystal, {"a": 3.0, "b": 4.1, "c": 5.3})["A"][0:3, 0:3].T
    Ap = convert_to_primitive_vector(group.info.lattice, A) @ ROTATION.T
    structure = _kpath_structure(group, Ap)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        info = seekpath.get_path_orig_cell(structure)
        assert (info["spacegroup_number"], info["is_supercell"]) == (no, False)
        B = 2 * np.pi * np.linalg.inv(Ap).T
        k_point = {"Γ" if k == "GAMMA" else k: v for k, v in info["point_coords"].items()}
        check_kpath(_join_kpath(info["path"]), k_point, B, structure)


# ==================================================
@pytest.mark.parametrize(
    "path, k_path",
    [
        ([("GAMMA", "X")], "Γ-X"),
        ([("GAMMA", "X"), ("X", "M"), ("M", "GAMMA")], "Γ-X-M-Γ"),
        ([("GAMMA", "X"), ("X", "M"), ("Z", "R"), ("R", "A")], "Γ-X-M|Z-R-A"),
        ([("X", "R"), ("M", "GAMMA"), ("GAMMA", "R")], "X-R|M-Γ-R"),
    ],
)
def test_join_kpath(path, k_path):
    assert _join_kpath(path) == k_path


# ==================================================
@pytest.mark.parametrize("lattice", ["A", "B", "C", "I", "F", "R"])
def test_primitive_vector_and_coordinates(lattice):
    # x_p @ A_p = x_c @ A for the primitive fractional coordinates x_p of a conventional point x_c.
    A = cell_vectors(3.0, 4.1, 5.3, 80.0, 100.0, 110.0)
    Ap = convert_to_primitive_vector(lattice, A)
    for x_c in np.random.default_rng(0).random((5, 3)):
        x_p = np.asarray(convert_to_primitive(lattice, x_c, shift=False), dtype=float)
        assert np.allclose(x_p @ Ap, x_c @ A)


# ==================================================
def test_legacy_unit_vector_primitive(models):
    # a model saved with the former unit_vector_primitive (A Pi^T): ModelAnalyzer uses P^T A.
    mm = MaterialModel(topdir=models)
    mm.load("bct")
    A = np.asarray(mm["unit_vector"], dtype=float)
    mm["model"] = "bct_legacy"
    mm["unit_vector_primitive"] = np.asarray(convert_to_primitive("I", A, shift=False), dtype=float).tolist()
    mm.save()
    ma = ModelAnalyzer(models)
    ma.analyze({"samb": {"model": "bct_legacy", "parameter": {"z1": 1.0}}, "grid": (10, 10, 10)})
    assert np.allclose(ma["info"]["A"], convert_to_primitive_vector("I", A))
    assert ma["output"]["dispersion"]["k_path"].startswith("Γ-")


# ==================================================
def test_kpath_space_group_mismatch(models, monkeypatch):
    # a structure of another space group is reported with a clear message.
    structure = lambda group, A: (A, [[0, 0, 0], [0.1, 0.2, 0.3]], [1, 2])  # P1.
    monkeypatch.setattr(model_analyzer, "_kpath_structure", structure)
    ma = ModelAnalyzer(models)
    with pytest.raises(ValueError, match="differs from that of the model"):
        ma.analyze({"samb": {"model": "bct", "parameter": {"z1": 1.0}}, "grid": (10, 10, 10)})


# ==================================================
@pytest.mark.parametrize("group", [1, 2])
def test_triclinic_cell(tmp_path, group):
    # the given cell is used for triclinic groups, also for the neighbour order of bonds.
    cell = CASES["triclinic"][1]
    model = {"model": "tri", "group": group, "cell": cell, "site": {"X": ("[0,0,0]", "s")}, "bond": [("X", "X", [1, 2, 3])]}
    create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir=str(tmp_path))
    mm = MaterialModel(topdir=str(tmp_path))
    mm.load("tri")
    assert mm["cell_info"]["cell"] == cell
    assert np.allclose(np.asarray(mm["unit_vector"], dtype=float), cell_vectors(**cell))
    # neighbours: a1 (3.0), a2 (4.0), a1+a2 (4.098, shorter than a3 = 5.0 since gamma > 90).
    distance = [round(b.distance, 3) for b in mm["bond"]["representative"].values()]
    assert distance == [3.0, 4.0, 4.098]


# ==================================================
@pytest.mark.parametrize(
    "cell, match",
    [
        ({"a": 0.0}, "lattice constants must be positive"),
        ({"b": float("nan")}, "lattice constants must be positive"),
        ({"c": float("inf")}, "lattice constants must be positive"),
        ({"beta": float("nan")}, r"must be in \(0, 180\)"),
        ({"gamma": 180.0}, r"must be in \(0, 180\)"),
        ({"alpha": 120.0, "beta": 120.0, "gamma": 120.0}, "do not form a cell"),  # zero volume.
        ({"alpha": 10.0, "beta": 10.0, "gamma": 170.0}, "do not form a cell"),  # negative (V/abc)^2.
    ],
)
def test_invalid_cell(cell, match):
    with pytest.raises(ValueError, match=match):
        get_cell_info("triclinic", cell)


# ==================================================
@pytest.mark.parametrize(
    "crystal, gamma",
    [
        ("triclinic", 90),
        ("monoclinic", 90),
        ("orthorhombic", 90),
        ("tetragonal", 90),
        ("trigonal", 120),
        ("hexagonal", 120),
        ("cubic", 90),
    ],
)
def test_default_cell(crystal, gamma):
    # without cell, a = b = c = 1 and right angles (except for gamma of trigonal and hexagonal).
    info = get_cell_info(crystal, {})
    assert info["cell"] == {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": float(gamma)}


def test_nearly_flat_monoclinic_cell():
    # a valid but nearly flat cell is accepted.
    info = get_cell_info("monoclinic", {"a": 3.0, "b": 4.0, "c": 5.0, "beta": 179.999})
    assert info["volume"] == pytest.approx(3.0 * 4.0 * 5.0 * np.sin(np.radians(179.999)))
