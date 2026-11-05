# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for DOS, fermi energy from the number of electrons, orbital weight (fat band), saved eigen system,
and user-defined callbacks in ModelAnalyzer.
"""

import os
import pickle
import re
import shutil

import numpy as np
import pytest

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util_model_analyzer import local_expectation, output_dispersion, plot_save_dispersion, power_energy
from multipie.util.util import read_dict

# spinful s and p orbitals at A, and s orbital at B (simple cubic, CsCl type), with SAMBs of identity irrep.
MODEL = {
    "model": "sp",
    "group": 221,
    "cell": {"a": 1.0},
    "site": {"A": ("[0,0,0]", ["s", "p"]), "B": ("[1/2,1/2,1/2]", "s")},
    "bond": [("A", "B", 1), ("A", "A", 1)],
    "spinful": True,
    "pdf": {"create": False},
    "qtdraw": {"create": False},
}
SELECT = {"Gamma": ["A1g"]}

CALLBACK = """
import os
import pickle


def save(key, Ek, Uk, info):
    with open(os.path.join(info["_topdir"], key + ".pkl"), "wb") as f:
        pickle.dump((Ek, Uk, info), f)


def callback_k_path(Ek, Uk, info):
    save("path", Ek, Uk, info)


def callback_k_full(Ek, Uk, info):
    save("full", Ek, Uk, info)
"""


# ==================================================
@pytest.fixture
def in_tmp(tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
        assert os.getcwd() == str(tmp_path)  # MultiPie does not change the working directory.
    finally:
        os.chdir(cwd)


# ==================================================
@pytest.fixture
def sp_model(in_tmp):
    create_model(dict(MODEL), topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "sp", "select": SELECT}, "output": {"dispersion": {"k_path": None}}})
    parameter = {z: 0.1 * (no + 1) * (-1) ** no for no, z in enumerate(ma["samb"]["matrix_info"]["matrix"])}

    return in_tmp, parameter


# ==================================================
def run(parameter, output, grid=(4, 4, 4)):
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "sp", "select": SELECT, "parameter": parameter}, "grid": grid, "output": output})
    data = np.loadtxt(os.path.join("out", "sp", "output", "sp_dos.txt")) if ma.output["dos"] else None

    return ma, data


# ==================================================
def load(topdir, key):
    with open(topdir / "out" / f"{key}.pkl", "rb") as f:
        return pickle.load(f)


# ==================================================
def test_dos_auto(sp_model):
    # partial DOS for each site cluster and orbital (summed over spin), total DOS, and integrated total DOS.
    _, parameter = sp_model
    ma, data = run(parameter, {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 500})
    assert ma["output"]["dos"]["column"] == ["E", "total", "N(E)", "s@A", "px@A", "py@A", "pz@A", "s@B"]
    assert data.shape == (500, 8)

    dE = data[1, 0] - data[0, 0]
    assert data[:, 1].sum() * dE == pytest.approx(10.0)  # number of bands.
    assert data[:, 3:].sum(axis=0) * dE == pytest.approx([2.0] * 5)  # two spins for each orbital.
    assert np.allclose(data[:, 3:].sum(axis=1), data[:, 1])

    # N(E) increases from 0 to the number of bands.
    assert data[0, 2] == 0.0
    assert data[-1, 2] == pytest.approx(10.0)
    assert np.all(np.diff(data[:, 2]) >= -1e-12)

    # the energy range includes all eigen values.
    e_min, e_max, width = (ma["output"]["dos"][k] for k in ["e_min", "e_max", "width"])
    assert ma["output"]["dos"]["method"] == "gaussian"
    assert data[0, 0] == pytest.approx(e_min - 5 * width)
    assert data[-1, 0] == pytest.approx(e_max + 5 * width)
    assert width == pytest.approx(0.01 * (e_max - e_min))

    with open(os.path.join("out", "sp", "output", "sp_dos.txt")) as f:
        comment = [line for line in f if line.startswith("#")]
    assert "# column 3: N(E), integrated total DOS\n" in comment
    assert "# column 4: s@A = (s,u)@A(1) + (s,d)@A(1)\n" in comment
    assert "# column 8: s@B = (s,u)@B(1) + (s,d)@B(1)\n" in comment

    # plot and gnuplot.
    for f in ["sp_dos.pdf", "sp_dos.eps", "plot_dos.gnu"]:
        assert os.path.isfile(os.path.join("out", "sp", "output", f))

    # True is the same as "auto".
    _, data_true = run(parameter, {"dispersion": {"k_path": None}, "dos": True, "dos_e_num": 500})
    assert np.allclose(data_true, data)


# ==================================================
def test_dos_list(sp_model):
    # sum of partial DOS of the kets in each list.
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "dos": [["(s,u)@A(1)"], ["(s,u)@B(1)", "(s,d)@B(1)"]], "dos_e_num": 400}
    ma, data = run(parameter, output)
    assert ma["output"]["dos"]["column"] == ["E", "total", "N(E)", "(s,u)@A(1)", "(s,u)@B(1)+(s,d)@B(1)"]
    dE = data[1, 0] - data[0, 0]
    assert data[:, 3:].sum(axis=0) * dE == pytest.approx([1.0, 2.0])

    _, auto = run(parameter, {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 400})
    assert np.allclose(data[:, 4], auto[:, 7])  # s@B.


# ==================================================
def test_dos_spin(sp_model):
    # spin-resolved partial DOS for "auto".
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 400}
    ma, data = run(parameter, output | {"dos_spin": True})
    column = ["(s,u)@A", "(s,d)@A", "(px,u)@A", "(px,d)@A", "(py,u)@A", "(py,d)@A", "(pz,u)@A", "(pz,d)@A"]
    assert ma["output"]["dos"]["column"][3:] == column + ["(s,u)@B", "(s,d)@B"]
    dE = data[1, 0] - data[0, 0]
    assert data[:, 3:].sum(axis=0) * dE == pytest.approx([1.0] * 10)

    _, auto = run(parameter, output)
    assert np.allclose(data[:, 3::2] + data[:, 4::2], auto[:, 3:])


# ==================================================
def test_dos_gaussian(sp_model):
    # DOS agrees with the direct sum of Gaussians, (1/Nk) sum_kn w_kn g(E-E_kn).
    _, parameter = sp_model
    width = 0.05
    ma, data = run(parameter, {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 4000, "dos_width": width})
    kv, Ek, Uk = ma.get_eigen_full()
    E = data[:, 0][:, None, None]
    g = np.exp(-((E - Ek[None]) ** 2) / (2 * width**2)) / (np.sqrt(2 * np.pi) * width)  # [E, k, n].
    total = g.sum(axis=(1, 2)) / len(kv)
    w = (np.abs(Uk[:, 8, :]) ** 2 + np.abs(Uk[:, 9, :]) ** 2)[None]  # s@B.
    partial = (g * w).sum(axis=(1, 2)) / len(kv)
    assert np.max(np.abs(data[:, 1] - total)) < 1e-2 * np.max(total)
    assert np.max(np.abs(data[:, 7] - partial)) < 1e-2 * np.max(partial)


# ==================================================
def test_dos_tetrahedron(sp_model):
    # optimized tetrahedron method by libtetrabz agrees with Gaussian broadening for fine grid and small width.
    pytest.importorskip("libtetrabz")
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 600}
    grid = (12, 10, 8)  # different sizes, to check the order of k points.
    ma, tetra = run(parameter, output | {"dos_method": "tetrahedron"}, grid=grid)
    assert ma["output"]["dos"]["method"] == "tetrahedron" and ma["output"]["dos"]["width"] is None
    dE = tetra[1, 0] - tetra[0, 0]
    assert tetra[:, 1].sum() * dE == pytest.approx(10.0, rel=1e-2)
    assert np.allclose(tetra[:, 3:].sum(axis=1), tetra[:, 1])
    e_min, e_max = ma["output"]["dos"]["e_min"], ma["output"]["dos"]["e_max"]
    assert tetra[0, 0] == pytest.approx(e_min - 0.05 * (e_max - e_min))

    with open(os.path.join("out", "sp", "output", "sp_dos.txt")) as f:
        assert "# optimized tetrahedron method (libtetrabz)\n" in f.readlines()

    # the k points on the full grid are in the order of eig[i1,i2,i3] for libtetrabz.
    kv, _ = ma.get_eigenvalue_full()
    idx = np.stack(np.meshgrid(*[np.arange(n) for n in grid], indexing="ij"), axis=-1)
    assert np.allclose(kv.reshape(*grid, 3), idx / np.asarray(grid))

    # N(E) of both methods are close (DOS itself differs by the broadening).
    _, gauss = run(parameter, output | {"dos_width": 0.005, "dos_e_num": 4000}, grid=(24, 24, 24))
    N_gauss = np.interp(tetra[:, 0], gauss[:, 0], gauss[:, 2])
    assert np.max(np.abs(tetra[:, 2] - N_gauss)) < 0.1


# ==================================================
def test_dos_tetrahedron_without_libtetrabz(sp_model, monkeypatch):
    _, parameter = sp_model
    monkeypatch.setitem(__import__("sys").modules, "libtetrabz", None)
    with pytest.raises(ImportError, match="libtetrabz is required"):
        run(parameter, {"dispersion": {"k_path": None}, "dos": "auto", "dos_method": "tetrahedron"})


# ==================================================
def test_dos_fermi_energy(sp_model):
    # the energy is shifted by the Fermi energy, as in the dispersion.
    _, parameter = sp_model
    ma = ModelAnalyzer("out")
    control = {"samb": {"model": "sp", "select": SELECT, "parameter": parameter}, "grid": (3, 3, 3), "output": {"dos": "auto"}}
    ma.reset(control)
    ma.exec_samb()
    ma.set_fermi_energy(0.3)
    ma.compute_physical_quantity()
    data = np.loadtxt(os.path.join("out", "sp", "output", "sp_dos.txt"))
    assert data[0, 0] == pytest.approx(ma["output"]["dos"]["e_min"] - 5 * ma["output"]["dos"]["width"] - 0.3)


# ==================================================
@pytest.mark.parametrize("n_electron", [0, 4, 4.5, 10])
def test_n_electron(sp_model, n_electron):
    # fermi energy from the number of occupied states per primitive cell (one per band).
    _, parameter = sp_model
    output = {"dispersion": {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}}, "n_electron": n_electron}
    ma, _ = run(parameter, output, grid=(3, 3, 3))
    ef = ma["info"]["fermi_energy"]
    assert len(ma.get_k_full()) == 27  # n_electron x Nk is not integer for 4.5.
    _, Ek, _ = ma.get_eigen_full()
    e = np.sort(Ek.ravel())
    # Ek is computed again here, so that it is compared with tolerance (states are doubly degenerate by spin).
    if n_electron == 0:
        assert ef == pytest.approx(e[0])
    elif n_electron == 10:
        assert ef == pytest.approx(e[-1])
    elif n_electron == 4:
        assert ef == pytest.approx(0.5 * (e[4 * len(Ek) - 1] + e[4 * len(Ek)]))
    else:
        assert ef == pytest.approx(e[int(4.5 * len(Ek))])
    assert ma["info"]["n_electron"] == n_electron

    # dispersion is shifted by the fermi energy.
    with open(os.path.join("out", "sp", "output", "sp_dispersion.txt")) as f:
        assert f"# shifted by fermi energy = {ef}\n" in f.readlines()


# ==================================================
def test_n_electron_dos(sp_model):
    # N(E_F) is the number of occupied states.
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 2000, "n_electron": 4}
    _, data = run(parameter, output, grid=(8, 8, 8))
    assert np.interp(0.0, data[:, 0], data[:, 2]) == pytest.approx(4.0, abs=0.1)

    pytest.importorskip("libtetrabz")
    ma, data = run(parameter, output | {"dos_method": "tetrahedron"}, grid=(8, 8, 8))
    assert np.interp(0.0, data[:, 0], data[:, 2]) == pytest.approx(4.0, abs=0.05)


# ==================================================
@pytest.mark.parametrize(
    "output, match",
    [
        ({"dos": [["s@A(1)"]]}, r"unknown ket\(s\) in 'output/dos': \['s@A\(1\)'\]"),
        ({"dos": "all"}, r"'output/dos' must be False, 'auto', or a list of lists"),
        ({"dos": ["(s,u)@A(1)"]}, r"'output/dos' must be False, 'auto', or a list of lists"),
        ({"dos": [[]]}, r"'output/dos' must be False, 'auto', or a list of lists"),
        ({"dos": "auto", "dos_e_num": 1}, r"'output/dos_e_num' must be integer"),
        ({"dos": "auto", "dos_e_num": 11}, r"'output/dos_e_num' must be at least 12 for 'gaussian'"),
        ({"dos": "auto", "dos_e_num": 100.0}, r"'output/dos_e_num' must be integer"),
        ({"dos": "auto", "dos_width": 0.0}, r"'output/dos_width' must be positive"),
        ({"dos": "auto", "dos_method": "tetra"}, r"'output/dos_method' must be 'gaussian' or 'tetrahedron'"),
        ({"n_electron": -1}, r"'output/n_electron' must be a number in \[0, 10\]"),
        ({"n_electron": 11}, r"'output/n_electron' must be a number in \[0, 10\]"),
        ({"n_electron": True}, r"'output/n_electron' must be a number"),
        ({"save_k_path": True}, r"'output/save_k_path' requires the dispersion"),
    ],
)
def test_dos_error(sp_model, output, match):
    _, parameter = sp_model
    with pytest.raises(ValueError, match=match):
        run(parameter, {"dispersion": {"k_path": None}} | output)


# ==================================================
def test_dos_point_group(in_tmp):
    # a model of a point group (molecule) has Γ point only.
    model = {"model": "mol", "group": "D3h", "site": {"A": ("[1,0,0]", "s")}, "bond": [("A", "A", [1])]}
    create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "mol", "parameter": {"z1": 1.0}}, "output": {"dos": "auto"}})
    data = np.loadtxt(os.path.join("out", "mol", "output", "mol_dos.txt"))
    dE = data[1, 0] - data[0, 0]
    assert data[:, 1].sum() * dE == pytest.approx(len(ma.get_ket()))  # three sites (with spin if spinful).
    assert ma["output"]["dos"]["column"][3:] == ["s@A"]

    with pytest.raises(ValueError, match="'tetrahedron' cannot be used for a model of point group"):
        ma.analyze({"samb": {"model": "mol", "parameter": {"z1": 1.0}}, "output": {"dos": "auto", "dos_method": "tetrahedron"}})


# ==================================================
def test_weight(sp_model):
    # orbital weight along k path (fat band).
    _, parameter = sp_model
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"}
    output = {"dispersion": {"k_path": "Γ-X-M", "k_point": k_point, "weight": "auto"}, "save_k_path": True}
    run(parameter, output)
    path = os.path.join("out", "sp", "output")
    for f in ["sp_weight.txt", "sp_weight.pdf", "plot_weight.gnu"]:
        assert os.path.isfile(os.path.join(path, f))
    with open(os.path.join(path, "sp_weight.txt")) as f:
        assert f.readline() == "# k Energy s@A px@A py@A pz@A s@B\n"
    data = np.loadtxt(os.path.join(path, "sp_weight.txt"))
    assert data.shape[1] == 7
    assert np.allclose(data[:, 2:].sum(axis=1), 1.0)  # sum of weights of all kets is 1.

    # weight of the kets in each list.
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    run(parameter, output | {"dispersion": output["dispersion"] | {"weight": [["(s,u)@B(1)", "(s,d)@B(1)"]]}})
    data = np.loadtxt(os.path.join(path, "sp_weight.txt"))
    w = np.abs(z["Uk"][:, 8, :]) ** 2 + np.abs(z["Uk"][:, 9, :]) ** 2
    assert np.allclose(np.sort(data[:, 2]), np.sort(w.ravel()))

    with pytest.raises(ValueError, match=r"unknown ket\(s\) in 'output/dispersion/weight'"):
        run(parameter, output | {"dispersion": output["dispersion"] | {"weight": [["s@B"]]}})


# ==================================================
def test_save_eigen(sp_model):
    # eigen system along k path and on full k grid in .npz.
    _, parameter = sp_model
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}
    output = {"dispersion": {"k_path": "Γ-X", "k_point": k_point}, "save_k_path": True, "save_k_full": True}
    ma, _ = run(parameter, output, grid=(4, 3, 2))
    path = os.path.join("out", "sp", "output")

    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    assert z["Ek"].shape == (len(z["k_grid"]), 10) and z["Uk"].shape == (len(z["k_grid"]), 10, 10)
    assert len(z["k_linear"]) == len(z["k_grid"])
    Hk = np.einsum("kmi,ki,kni->kmn", z["Uk"], z["Ek"], z["Uk"].conj())
    assert np.allclose(np.linalg.eigvalsh(Hk), z["Ek"])

    z = np.load(os.path.join(path, "sp_eigen_k_full.npz"))
    kv, Ek, Uk = ma.get_eigen_full()
    assert np.allclose(z["k_grid"], kv) and np.allclose(z["Ek"], Ek)
    assert z["grid"].tolist() == [4, 3, 2]
    assert z["ket"].tolist() == ma.get_ket()
    assert float(z["fermi_energy"]) == 0.0


# ==================================================
def test_callback(sp_model):
    topdir, parameter = sp_model
    (topdir / "out" / "my_callback.py").write_text(CALLBACK)
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"}
    output = {
        "dispersion": {"k_path": "Γ-X-M", "k_point": k_point, "power": 0.5},
        "callback_k_path": "my_callback.py",
        "callback_k_full": "my_callback.py",
    }
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "sp", "select": SELECT, "parameter": parameter}, "grid": (4, 3, 2), "output": output})

    # along k path: Ek is not raised to the power of the dispersion.
    Ek, Uk, info = load(topdir, "path")
    kv = info["k_grid"]
    assert Ek.shape == (len(kv), 10) and Uk.shape == (len(kv), 10, 10)
    assert np.min(Ek) < 0
    Hk = np.einsum("kmi,ki,kni->kmn", Uk, Ek, Uk.conj())
    assert np.allclose(np.linalg.eigvalsh(Hk), Ek)
    assert info["k_path"] == "Γ-X-M"
    assert info["k_point"] == {"Γ": [0.0, 0.0, 0.0], "X": [0.5, 0.0, 0.0], "M": [0.5, 0.5, 0.0]}
    assert np.allclose(kv[0], [0, 0, 0]) and np.allclose(kv[-1], [0.5, 0.5, 0])
    assert len(info["k_linear"]) == len(kv)
    assert info["ket"] == ma.get_ket() and info["ket"][0] == "(s,u)@A(1)"
    assert info["fermi_energy"] == 0.0
    assert info["name"] == "sp"
    assert info["_topdir"] == str(topdir / "out")
    assert info["_verbose"] is False

    # on full grid.
    Ek, Uk, info = load(topdir, "full")
    assert info["grid"] == [4, 3, 2]
    assert info["k_grid"].shape == (24, 3)
    assert Ek.shape == (24, 10) and Uk.shape == (24, 10, 10)
    _, Ek_ref, _ = ma.get_eigen_full()
    assert np.allclose(Ek, Ek_ref)
    assert info["k_path"] == "Γ-X-M"


# ==================================================
def test_callback_with_dos(sp_model):
    # DOS computed from the eigen system for callback_k_full is the same as that computed in chunks.
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    output = {"dispersion": {"k_path": None}, "dos": "auto"}
    _, data = run(parameter, output)
    _, data_cb = run(parameter, output | {"callback_k_full": "cb.py"})
    assert np.allclose(data, data_cb)


# ==================================================
@pytest.mark.parametrize(
    "output, error, match",
    [
        ({"callback_k_full": "cb"}, ValueError, r"'output/callback_k_full' must be '\.py' file"),
        ({"callback_k_full": "none.py"}, FileNotFoundError, r"callback file '.*none\.py' is not found"),
        ({"callback_k_path": "other.py"}, ValueError, r"function 'callback_k_path' is not defined"),
        ({"callback_k_path": "cb.py", "dispersion": {"k_path": None}}, ValueError, r"'output/callback_k_path' requires"),
    ],
)
def test_callback_error(sp_model, output, error, match):
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    (topdir / "out" / "other.py").write_text("def callback_k_full(Ek, Uk, info):\n    pass\n")
    output = {"dispersion": {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}}} | output
    with pytest.raises(error, match=match):
        run(parameter, output, grid=(2, 2, 2))


# ==================================================
@pytest.mark.parametrize("method", ["gaussian", "tetrahedron"])
def test_flat_band_insulator(in_tmp, method):
    # only on-site terms: all bands are flat. the weight of flat bands is kept, and fermi energy is at the middle of the gap.
    if method == "tetrahedron":
        pytest.importorskip("libtetrabz")
    model = {"model": "flat", "group": 221, "cell": {"a": 1.0}, "site": {"A": ("[0,0,0]", ["s", "p"])}, "bond": [("A", "A", 1)]}
    create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir="out")
    ma = ModelAnalyzer("out")
    output = {"dispersion": {"k_path": None}, "dos": "auto", "dos_method": method, "n_electron": 1}
    ma.analyze({"samb": {"model": "flat", "parameter": {"z1": -1.0, "z2": 1.0}}, "grid": (4, 4, 4), "output": output})
    _, Ek = ma.get_eigenvalue_full()
    assert np.allclose(np.ptp(Ek, axis=0), 0.0)
    e = np.unique(np.round(Ek, 10))
    assert len(e) == 2  # s and p levels.
    assert ma["info"]["fermi_energy"] == pytest.approx(0.5 * (e[0] + e[1]))

    data = np.loadtxt(os.path.join("out", "flat", "output", "flat_dos.txt"))
    dE = data[1, 0] - data[0, 0]
    assert data[:, 1].sum() * dE == pytest.approx(4.0)
    assert data[-1, 2] == pytest.approx(4.0, rel=1e-3)
    assert data[:, 3:].sum(axis=0) * dE == pytest.approx([1.0, 1.0, 1.0, 1.0])
    assert np.interp(0.0, data[:, 0], data[:, 2]) == pytest.approx(1.0, abs=1e-3)  # N(E_F) in the gap.


# ==================================================
def test_dos_narrow_width(sp_model):
    # the Gaussian width must not be smaller than the energy step, and the DOS agrees with the direct sum at width = step.
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "dos": "auto", "dos_e_num": 200}
    with pytest.raises(ValueError, match="Gaussian width of DOS, .* is smaller than the energy step"):
        run(parameter, output | {"dos_width": 0.001})

    run(parameter, output | {"dos_e_num": 12})  # the smallest number of points with default width.
    ma, data = run(parameter, output)  # default width is not smaller than the energy step.
    dE = data[1, 0] - data[0, 0]
    width = ma["output"]["dos"]["width"]
    assert width >= dE * (1 - 1e-10)

    _, Ek, _ = ma.get_eigen_full()
    E = data[:, 0][:, None, None]
    total = (np.exp(-((E - Ek[None]) ** 2) / (2 * width**2)) / (np.sqrt(2 * np.pi) * width)).sum(axis=(1, 2)) / len(Ek)
    assert np.max(np.abs(data[:, 1] - total)) < 0.05 * np.max(total)
    assert data[:, 1].sum() * dE == pytest.approx(10.0)
    assert np.sum(data[:, 0] * data[:, 1]) * dE / 10.0 == pytest.approx(Ek.mean(), abs=1e-8)  # mean energy is kept.


# ==================================================
@pytest.mark.parametrize("grid", [(0, 4, 4), (4, -1, 4), (4, 4, 2.0), (True, 4, 4), (4, 4), (4, 4, 4, 4), 4])
def test_invalid_grid(sp_model, grid):
    _, parameter = sp_model
    with pytest.raises(ValueError, match="'grid' must be three positive integers"):
        run(parameter, {"dispersion": {"k_path": None}, "dos": "auto"}, grid=grid)


# ==================================================
def test_callback_same_file(sp_model):
    # a file given for both callbacks is loaded once.
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(
        "import os\n"
        "count = 0\n"
        "def callback_k_path(Ek, Uk, info):\n"
        "    global count\n"
        "    count += 1\n"
        "def callback_k_full(Ek, Uk, info):\n"
        "    with open(os.path.join(info['_topdir'], 'count.txt'), 'w') as f:\n"
        "        f.write(str(count))\n"
    )
    output = {
        "dispersion": {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}},
        "callback_k_path": "cb.py",
        "callback_k_full": "cb.py",
    }
    run(parameter, output, grid=(2, 2, 2))
    assert (topdir / "out" / "count.txt").read_text() == "1"


# ==================================================
def test_docs_post_processing_example(in_tmp):
    # the example of docs/src/post_processing.md works as documented.
    example = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples")
    for f in ["graphene_in.py", "graphene_post_ctrl.py", "graphene_callback.py"]:
        shutil.copy(os.path.join(example, f), in_tmp)
    model = read_dict("graphene_in.py")
    model["pdf"] = {"create": False}
    model["qtdraw"] = {"create": False}
    create_model(model, topdir=".")
    control = read_dict("graphene_post_ctrl.py")
    ma = ModelAnalyzer(".")
    ma.analyze(control)

    assert ma["info"]["fermi_energy"] == pytest.approx(0.0, abs=1e-12)
    out = in_tmp / "graphene" / "output"
    for f in ["graphene_dos.txt", "graphene_weight.pdf", "graphene_eigen_k_path.npz", "my_band.txt", "plot_dos.gnu"]:
        assert os.path.isfile(out / f)
    with open(out / "graphene_dos.txt") as f:
        assert "# column 4: pz@C = pz@C(1) + pz@C(2)\n" in f.readlines()
    occupation = dict(line.split() for line in (out / "my_occupation.txt").read_text().splitlines())
    assert {k: float(v) for k, v in occupation.items()} == pytest.approx({"pz@C(1)": 0.5, "pz@C(2)": 0.5})

    # weight of pz@C(1) is 1/2 by the symmetry between the sublattices.
    z = np.load(out / "graphene_eigen_k_path.npz")
    ket = list(z["ket"])
    assert np.allclose(np.abs(z["Uk"][:, ket.index("pz@C(1)"), :]) ** 2, 0.5)


# ==================================================
@pytest.mark.parametrize(
    "power, e, expected",
    [
        (0.5, [-4.0, 0.0, 4.0], [-2.0, 0.0, 2.0]),  # negative Ek.
        (1.5, [-4.0, 4.0], [-8.0, 8.0]),  # power > 1.
        (2, [-3.0, 3.0], [9.0, 9.0]),  # Ek^power is computed for integer power.
        (3, [-2.0, 2.0], [-8.0, 8.0]),
        (-1, [-2.0, 2.0], [-0.5, 0.5]),
        (-0.5, [-4.0, 4.0], [-0.5, 0.5]),
        (-2, [-2.0, 2.0], [0.25, 0.25]),
        (0, [-2.0, 0.0, 2.0], [1.0, 1.0, 1.0]),  # nothing diverges at power = 0.
    ],
)
def test_power_energy(power, e, expected):
    assert power_energy(np.array(e), power) == pytest.approx(expected)


# ==================================================
def test_power_energy_divergence():
    # Ek = 0 with negative power cannot be plotted.
    assert np.array_equal(power_energy(np.array([0.0, 4.0]), -0.5), [np.nan, 0.5], equal_nan=True)
    # also within round-off error.
    assert np.array_equal(power_energy(np.array([1e-16, -1e-16, 4.0]), -0.5), [np.nan, np.nan, 0.5], equal_nan=True)
    # the threshold is relative to the largest value.
    assert power_energy(np.array([1e-11, 1e-3]), -1) == pytest.approx([1e11, 1e3])
    assert np.array_equal(power_energy(np.array([1e-7, 1e4]), -1), [np.nan, 1e-4], equal_nan=True)
    assert np.array_equal(power_energy(np.array([5e-11, 1.0]), -1), [np.nan, 1.0], equal_nan=True)
    assert power_energy(np.array([2e-10, 1.0]), -1) == pytest.approx([5e9, 1.0])
    # the largest value is taken over all the bands and k points: a k point with small values only is also skipped.
    assert np.array_equal(
        power_energy(np.array([[1e-12, 1e-12], [1.0, 2.0]]), -1), [[np.nan, np.nan], [1.0, 0.5]], equal_nan=True
    )
    # small values are kept for positive power.
    assert power_energy(np.array([1e-16, 4.0]), 0.5) == pytest.approx([1e-8, 2.0])


# ==================================================
def test_dispersion_power(sp_model):
    # negative Ek are not skipped: all values are written as -|Ek|^power, for each band and k point (power > 0 keeps the order).
    _, parameter = sp_model
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}
    output = {"dispersion": {"k_path": "Γ-X", "k_point": k_point}, "save_k_path": True}
    run(parameter, output)
    path = os.path.join("out", "sp", "output")
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    Ek = z["Ek"]
    assert Ek.min() < 0 < Ek.max()

    run(parameter, output | {"dispersion": output["dispersion"] | {"power": 0.5}})
    data = np.loadtxt(os.path.join(path, "sp_dispersion.txt"))
    assert data.shape[0] == Ek.size
    n_k, n_band = Ek.shape
    written = data[:, 1].reshape(n_band, n_k).T  # the file is written for each band.
    assert written == pytest.approx(np.sign(Ek) * np.abs(Ek) ** 0.5)
    assert np.all(np.diff(written, axis=1) >= -1e-12)
    assert data[:, 0].reshape(n_band, n_k)[0] == pytest.approx(z["k_linear"])


# ==================================================
def test_dispersion_power_fermi_energy(sp_model):
    # the power is taken from the fermi energy, also for the orbital weight.
    _, parameter = sp_model
    dispersion = {"k_path": "Γ-X-M", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"}, "weight": "auto"}
    output = {"n_electron": 4, "save_k_path": True}
    path = os.path.join("out", "sp", "output")

    ma, _ = run(parameter, output | {"dispersion": dispersion}, grid=(3, 3, 3))
    ef = ma["info"]["fermi_energy"]
    assert ef != 0.0
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    n_k, n_band = z["Ek"].shape
    reference = np.loadtxt(os.path.join(path, "sp_dispersion.txt"))
    assert reference[:, 1].reshape(n_band, n_k).T == pytest.approx(z["Ek"] - ef)

    # power = 1 is the same as without power.
    run(parameter, output | {"dispersion": dispersion | {"power": 1}}, grid=(3, 3, 3))
    assert np.loadtxt(os.path.join(path, "sp_dispersion.txt")) == pytest.approx(reference)

    run(parameter, output | {"dispersion": dispersion | {"power": 0.5}}, grid=(3, 3, 3))
    x = z["Ek"] - ef
    expected = np.sign(x) * np.abs(x) ** 0.5
    for f in ["sp_dispersion.txt", "sp_weight.txt"]:
        written = np.loadtxt(os.path.join(path, f))[:, 1].reshape(n_band, n_k).T
        assert written == pytest.approx(expected)
    assert np.array_equal(written < 0, z["Ek"] < ef)  # states below the fermi energy remain below it.
    with open(os.path.join(path, "sp_dispersion.txt")) as f:
        assert f"# shifted by fermi energy = {ef}\n" in f.readlines()


# ==================================================
def test_dispersion_power_all_divergent(sp_model):
    _, parameter = sp_model
    parameter = {z: 0.0 for z in parameter}
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}
    with pytest.raises(ValueError, match=r"no value of \(Ek-Ef\)"):
        run(parameter, {"dispersion": {"k_path": "Γ-X", "k_point": k_point, "power": -1}})


# ==================================================
def test_dispersion_power_partial_divergence(in_tmp):
    # Ek = 0 with negative power is not written, and the line is broken there. the plot can be created.
    k = np.array([0.0, 1.0, 2.0, 3.0])
    Ek = power_energy(np.array([[-1.0, 1.0], [0.0, 2.0], [-4.0, 4.0], [-2.0, 2.0]]), -1)
    assert np.isnan(Ek[1, 0]) and np.count_nonzero(np.isnan(Ek)) == 1

    fname = str(in_tmp / "disp.txt")
    output_dispersion(fname, k, 0.0, Ek)
    with open(fname) as f:
        body = "".join(line for line in f if not line.startswith("#")).strip()
    block = [b.split("\n") for b in re.split(r"\n\s*\n+", body)]
    assert [len(b) for b in block] == [1, 2, 4]  # the first band is broken at the second k point.
    assert np.loadtxt(fname) == pytest.approx(
        np.array([[0, -1], [2, -0.25], [3, -0.5], [0, 1], [1, 0.5], [2, 0.25], [3, 0.5]], dtype=float)
    )
    plot_save_dispersion(fname, {0.0: "G", 3.0: "X"})
    assert os.path.isfile(in_tmp / "disp.pdf")


# ==================================================
def test_local_expectation():
    # known spin states: Uk[k,:,n], in the (u, d) basis.
    s = 1 / np.sqrt(2)
    Uk = np.array(
        [
            [[1, 0], [0, 1]],  # Sz = 1/2, -1/2.
            [[s, s], [s, -s]],  # Sx = 1/2, -1/2.
            [[s, s], [1j * s, -1j * s]],  # Sy = 1/2, -1/2.
        ]
    )
    operator = {
        "Sz": ([[0.5, 0], [0, -0.5]], [[0.5, -0.5], [0, 0], [0, 0]]),
        "Sx": ([[0, 0.5], [0.5, 0]], [[0, 0], [0.5, -0.5], [0, 0]]),
        "Sy": ([[0, -0.5j], [0.5j, 0]], [[0, 0], [0, 0], [0.5, -0.5]]),
    }
    for op, expected in operator.values():
        assert local_expectation(Uk, np.array(op)) == pytest.approx(np.array(expected))


# ==================================================
def test_dispersion_local(sp_model):
    # columns of local operators in the dispersion file, for each band and k point.
    _, parameter = sp_model
    keywords = ["Sx", "Sz", "Lz"]
    dispersion = {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}, "local": keywords}
    ma, _ = run(parameter, {"dispersion": dispersion, "save_k_path": True})
    path = os.path.join("out", "sp", "output")
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    data = np.loadtxt(os.path.join(path, "sp_dispersion.txt"))
    n_k, n_band = z["Ek"].shape
    assert data.shape == (n_k * n_band, 2 + len(keywords))
    for no, tag in enumerate(keywords):
        column = data[:, 2 + no].reshape(n_band, n_k).T
        expected = np.einsum("kmi,mn,kni->ki", z["Uk"].conj(), ma.get_local_operator(tag), z["Uk"]).real
        assert column == pytest.approx(expected, abs=1e-10)
        assert column.sum(axis=1) == pytest.approx(np.zeros(n_k), abs=1e-10)  # the operators are traceless.


# ==================================================
@pytest.mark.parametrize(
    "site",
    [
        {"A": ("[0,0,0]", ["s", "p", "d"])},
        {"A": ("[0,0,0]", ["s", "p"]), "B": ("[1/2,1/2,1/2]", "s")},  # two sites.
        {"A": ("[1/2,0,0]", ["s", "p"])},  # a site with three sublattices.
        {"A": ("[0,0,0]", ["dv", "dxy"])},  # a subset of orbitals.
        {"A": ("[0,0,0]", "f")},
    ],
)
def test_spin_operator(in_tmp, site):
    # spin operators in the (u, d) basis, the commutation relation of S and L, and S^2 = 3/4.
    model = dict(MODEL) | {"model": "spd", "site": site, "bond": []}
    create_model(model, topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "spd", "select": {}}, "output": {"dispersion": {"k_path": None}}})

    ket = ma.get_ket()
    orbital = [k.split(",")[0] + k.split("@")[1] for k in ket]  # orbital and site.
    S2 = sum(np.linalg.matrix_power(ma.get_local_operator(t), 2) for t in ["Sx", "Sy", "Sz"])
    assert np.allclose(S2, 0.75 * np.eye(len(ket)))
    if "d" in str(site["A"][1]) and "s" in site["A"][1]:
        # s, p and d orbitals with spin: Lz = 0 (x 2), -1, 0, 1 (x 2), and -2, ..., 2 (x 2).
        assert np.linalg.eigvalsh(ma.get_local_operator("Lz")) == pytest.approx(
            sorted(2 * [0] + 2 * [-1, 0, 1] + 2 * [-2, -1, 0, 1, 2]), abs=1e-10
        )
    pauli = {"Sx": [[0, 0.5], [0.5, 0]], "Sy": [[0, -0.5j], [0.5j, 0]], "Sz": [[0.5, 0], [0, -0.5]]}
    for tag, s in pauli.items():
        op = ma.get_local_operator(tag)
        for i, ki in enumerate(ket):
            for j, kj in enumerate(ket):
                expected = s[int(",d)" in ki)][int(",d)" in kj)] if orbital[i] == orbital[j] else 0.0
                assert op[i, j] == pytest.approx(expected)

    for a, b, c in [("x", "y", "z"), ("y", "z", "x"), ("z", "x", "y")]:
        # L of a subset of orbitals in a shell does not close.
        for t in ["S"] if "dv" in site["A"][1] else ["S", "L"]:
            A, B, C = (ma.get_local_operator(t + i) for i in (a, b, c))
            assert np.allclose(A @ B - B @ A, 1j * C)
        assert np.allclose(ma.get_local_operator("S" + a) @ ma.get_local_operator("L" + b),
                           ma.get_local_operator("L" + b) @ ma.get_local_operator("S" + a))  # fmt: skip


# ==================================================
def test_callback_info(sp_model):
    # expectation values on the full k grid are passed as info[keyword].
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    keywords = ["Sx", "Sy", "Sz", "Lx", "Ly", "Lz", "Qu", "Qv", "Qyz", "Qxz", "Qxy"]
    output = {
        "dispersion": {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}},
        "callback_k_path": "cb.py",
        "callback_k_full": "cb.py",
        "callback_info": keywords,
    }
    ma, _ = run(parameter, output, grid=(4, 3, 2))
    Ek, Uk, info = load(topdir, "full")
    for tag in keywords:
        expected = np.einsum("kmi,mn,kni->ki", Uk.conj(), ma.get_local_operator(tag), Uk).real
        assert info[tag].shape == Ek.shape
        assert np.allclose(info[tag], expected)
    # the band sum of S is the trace of the operator, which is zero.
    for tag in ["Sx", "Sy", "Sz"]:
        assert np.allclose(info[tag].sum(axis=1), 0.0)
    # |<S>| <= 1/2 for each state.
    assert np.all(info["Sx"] ** 2 + info["Sy"] ** 2 + info["Sz"] ** 2 <= 0.25 + 1e-10)
    # the sum over degenerate bands of a time-odd operator vanishes by inversion and time-reversal symmetries.
    for k in range(len(Ek)):
        start = np.flatnonzero(np.r_[True, np.diff(Ek[k]) > 1e-8])
        for tag in ["Sx", "Sy", "Sz", "Lz"]:
            assert np.allclose(np.add.reduceat(info[tag][k], start), 0.0, atol=1e-8)

    _, _, info_path = load(topdir, "path")
    assert not set(keywords) & set(info_path)

    # no quantity is added by default, nor by an empty list.
    default = {key: value for key, value in output.items() if key != "callback_info"}
    for out in [default, default | {"callback_info": []}]:
        run(parameter, out, grid=(2, 2, 2))
        assert not set(keywords) & set(load(topdir, "full")[2])


# ==================================================
@pytest.mark.parametrize(
    "output, match",
    [
        ({"callback_info": ["Sw"]}, r"'output/callback_info' must be a list of"),
        ({"callback_info": "Sx"}, r"'output/callback_info' must be a list of"),
        ({"callback_info": [1]}, r"'output/callback_info' must be a list of"),
        ({"callback_info": [None]}, r"'output/callback_info' must be a list of"),
        ({"callback_info": [np.array(["Sx"])]}, r"'output/callback_info' must be a list of"),
        ({"callback_info": ["Sx"], "callback_k_full": None}, r"'output/callback_info' requires 'output/callback_k_full'"),
    ],
)
def test_callback_info_error(sp_model, output, match):
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    with pytest.raises(ValueError, match=match):
        run(parameter, {"dispersion": {"k_path": None}, "callback_k_full": "cb.py"} | output, grid=(2, 2, 2))


# ==================================================
def test_callback_info_basis(in_tmp):
    # spin operators for a spinless model are rejected, and orbital operators can be used.
    model = dict(MODEL) | {"model": "lg", "spinful": False, "site": {"A": ("[0,0,0]", "p")}, "bond": [("A", "A", 1)]}
    create_model(model, topdir="out")
    (in_tmp / "out" / "cb.py").write_text(CALLBACK)
    ma = ModelAnalyzer("out")
    control = {"samb": {"model": "lg", "select": SELECT, "parameter": {"z1": 1.0}}, "grid": (2, 2, 2)}
    output = {"dispersion": {"k_path": None}, "callback_k_full": "cb.py"}
    with pytest.raises(ValueError, match="spin operator 'Sx' requires a spinful model"):
        ma.analyze(control | {"output": output | {"callback_info": ["Sx"]}})
    ma.analyze(control | {"output": output | {"callback_info": ["Lz"]}})
    assert load(in_tmp, "full")[2]["Lz"].shape == (8, 3)

    # no basis (wannier mode without a model) or jml basis.
    for basis_type in [None, "jml"]:
        ma.set_basis_type(basis_type)
        with pytest.raises(ValueError, match="requires a model with lg or lgs basis"):
            ma.get_callback_info_operator(True)


# ==================================================
def test_power_keeps_saved_eigen_system(sp_model):
    # the eigen system passed to the callback and saved in .npz is not raised to the power.
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    dispersion = {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}}
    filename = os.path.join("out", "sp", "output", "sp_eigen_k_path.npz")
    output = {"save_k_path": True, "n_electron": 4}
    ma, _ = run(parameter, output | {"dispersion": dispersion})
    assert ma["info"]["fermi_energy"] != 0.0
    raw = np.load(filename)
    Ek, Uk = raw["Ek"], raw["Uk"]
    os.remove(filename)
    run(parameter, output | {"dispersion": dispersion | {"power": 0.5}, "callback_k_path": "cb.py"})
    z = np.load(filename)
    assert z["Ek"] == pytest.approx(Ek) and z["Uk"] == pytest.approx(Uk)
    Ek_callback, Uk_callback, _ = load(topdir, "path")
    assert Ek_callback == pytest.approx(Ek) and Uk_callback == pytest.approx(Uk)


# ==================================================
def test_power_keeps_full_grid(sp_model):
    # the DOS, the eigen system on the full k grid, and the callback are not raised to the power.
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    dispersion = {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}}
    filename = os.path.join("out", "sp", "output", "sp_eigen_k_full.npz")
    output = {"dos": "auto", "save_k_full": True, "n_electron": 4, "callback_k_full": "cb.py"}
    ma, dos = run(parameter, output | {"dispersion": dispersion})
    assert ma["info"]["fermi_energy"] != 0.0
    raw = np.load(filename)
    Ek, Uk = raw["Ek"], raw["Uk"]
    assert load(topdir, "full")[0] == pytest.approx(Ek)
    _, dos_power = run(parameter, output | {"dispersion": dispersion | {"power": 0.5}})
    assert dos_power == pytest.approx(dos)
    z = np.load(filename)
    assert z["Ek"] == pytest.approx(Ek) and z["Uk"] == pytest.approx(Uk)
    Ek_callback, Uk_callback, _ = load(topdir, "full")
    assert Ek_callback == pytest.approx(Ek) and Uk_callback == pytest.approx(Uk)


# ==================================================
def test_dispersion_power_plot_range(sp_model):
    # the range of the gnuplot script is that of the raised energies, with a margin of 10 %, at a non-zero Fermi energy.
    _, parameter = sp_model
    dispersion = {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}}
    output = {"save_k_path": True, "n_electron": 4}
    path = os.path.join("out", "sp", "output")
    ma, _ = run(parameter, output | {"dispersion": dispersion})
    ef = ma["info"]["fermi_energy"]
    assert ef != 0.0
    x = np.load(os.path.join(path, "sp_eigen_k_path.npz"))["Ek"] - ef
    assert x.min() < 0.0 < x.max()

    run(parameter, output | {"dispersion": dispersion | {"power": -1}})
    # a Kramers doublet at the Fermi energy within round-off error is skipped, also at a non-zero Fermi energy.
    cutoff = 1e-10 * np.abs(x).max()
    assert np.count_nonzero(np.abs(x) <= cutoff) == 2
    y = 1.0 / x[np.abs(x) > cutoff]
    assert np.sort(np.loadtxt(os.path.join(path, "sp_dispersion.txt"))[:, 1]) == pytest.approx(np.sort(y))
    with open(os.path.join(path, "plot_band.gnu")) as f:
        value = dict(line.split("=") for line in f if line.startswith(("ymin =", "ymax =")))
    margin = 0.1 * (y.max() - y.min())
    assert float(value["ymin "]) == pytest.approx(y.min() - margin)
    assert float(value["ymax "]) == pytest.approx(y.max() + margin)


# ==================================================
def test_dispersion_power_skip(in_tmp):
    # flat s band at zero within round-off error with negative power: only the s band is skipped.
    model = {"model": "flat", "group": 221, "cell": {"a": 1.0}, "site": {"A": ("[0,0,0]", ["s", "p"])}, "bond": [("A", "A", 1)]}
    create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir="out")
    dispersion = {"k_path": "Γ-X", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]"}, "power": -1}
    control = {"samb": {"model": "flat", "parameter": {"z1": 1e-13, "z2": 1.0}}, "grid": (4, 4, 4)}
    ModelAnalyzer("out").analyze(control | {"output": {"dispersion": dispersion, "save_k_path": True}})
    path = os.path.join("out", "flat", "output")
    z = np.load(os.path.join(path, "flat_eigen_k_path.npz"))
    data = np.loadtxt(os.path.join(path, "flat_dispersion.txt"))
    n_k, n_band = z["Ek"].shape
    assert np.max(np.abs(z["Ek"][:, 0])) < 1e-12 and 0.0 < np.min(np.abs(z["Ek"][:, 0]))
    assert data.shape[0] == n_k * (n_band - 1)
    assert data[:, 1] == pytest.approx(np.full(len(data), 1 / z["Ek"].max()))


# ==================================================
def test_weight_order(sp_model):
    # orbital weight is written for each band and k point, with the energy raised to the power.
    _, parameter = sp_model
    k_point = {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"}
    dispersion = {"k_path": "Γ-X-M", "k_point": k_point, "weight": [["(s,u)@B(1)", "(s,d)@B(1)"]], "power": 0.5}
    ma, _ = run(parameter, {"dispersion": dispersion, "save_k_path": True, "n_electron": 4}, grid=(3, 3, 3))
    path = os.path.join("out", "sp", "output")
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    n_k, n_band = z["Ek"].shape
    ket = list(z["ket"])
    weight = np.abs(z["Uk"][:, ket.index("(s,u)@B(1)"), :]) ** 2 + np.abs(z["Uk"][:, ket.index("(s,d)@B(1)"), :]) ** 2
    data = np.loadtxt(os.path.join(path, "sp_weight.txt"))
    x = z["Ek"] - ma["info"]["fermi_energy"]
    assert data[:, 1].reshape(n_band, n_k).T == pytest.approx(np.sign(x) * np.abs(x) ** 0.5)
    assert data[:, 2].reshape(n_band, n_k).T == pytest.approx(weight, abs=1e-12)


# ==================================================
@pytest.mark.parametrize("spinful", [False, True])
def test_orbital_operator(in_tmp, spinful):
    # (L_a)_bc = -i eps_abc in (px, py, pz), and Q are the quadrupoles of L with common factors in each shell.
    model = dict(MODEL) | {"model": "pd", "spinful": spinful, "site": {"A": ("[0,0,0]", ["p", "d"])}, "bond": []}
    create_model(model, topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "pd", "select": {}}, "output": {"dispersion": {"k_path": None}}})
    ket = ma.get_ket()

    def block(tag, shell):
        idx = [i for i, k in enumerate(ket) if k.replace("(", "").startswith(shell) and (not spinful or ",u" in k)]
        return ma.get_local_operator(tag)[np.ix_(idx, idx)]

    # orbital operators are diagonal in the spin, the same for u and d, and do not mix the shells.
    p = [i for i, k in enumerate(ket) if k.replace("(", "").startswith("p")]
    d = [i for i, k in enumerate(ket) if k.replace("(", "").startswith("d")]
    assert len(p) == 3 * (1 + spinful) and len(d) == 5 * (1 + spinful)
    for tag in ["Lx", "Ly", "Lz", "Qu", "Qv", "Qyz", "Qxz", "Qxy"]:
        op = ma.get_local_operator(tag)
        assert np.allclose(op[np.ix_(p, d)], 0.0)
        if spinful:
            u = [i for i, k in enumerate(ket) if ",u" in k]
            dn = [i for i, k in enumerate(ket) if ",d" in k]
            assert np.allclose(op[np.ix_(u, u)], op[np.ix_(dn, dn)])
            assert np.allclose(op[np.ix_(u, dn)], 0.0)

    eps = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[a, b, c], eps[a, c, b] = 1, -1
    for a, tag in enumerate(["Lx", "Ly", "Lz"]):
        assert block(tag, "p") == pytest.approx(-1j * eps[a])

    for shell in ["p", "d"]:
        x, y, z = (block(tag, shell) for tag in ["Lx", "Ly", "Lz"])
        L2 = x @ x + y @ y + z @ z
        Qxy = x @ y + y @ x
        factor = np.sum(block("Qxy", shell) * Qxy.conj()).real / np.sum(np.abs(Qxy) ** 2)
        assert factor != 0.0
        assert np.allclose(block("Qxy", shell), factor * Qxy)
        assert np.allclose(block("Qyz", shell), factor * (y @ z + z @ y))
        assert np.allclose(block("Qxz", shell), factor * (x @ z + z @ x))
        assert np.allclose(block("Qv", shell), factor * (x @ x - y @ y))
        assert np.allclose(block("Qu", shell), factor * (3 * z @ z - L2) / np.sqrt(3))
