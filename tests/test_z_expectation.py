# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for SAMB along the k path (output/dispersion/z), in callback_info, and their expectation values (output/z_expectation).
"""

import os
import pickle

import numpy as np
import pytest

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict
from multipie.util.util_model_analyzer import density_matrix_r, fourier_r_to_k, operator_expectation

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
K_PATH = {"k_path": "Γ-X-M", "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"}}

CALLBACK = """
import os
import pickle


def callback_k_full(Ek, Uk, info):
    with open(os.path.join(info["_topdir"], "full.pkl"), "wb") as f:
        pickle.dump((Ek, Uk, info), f)
"""


# ==================================================
@pytest.fixture
def sp_model(tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        create_model(dict(MODEL), topdir="out")
        ma = ModelAnalyzer("out")
        ma.analyze({"samb": {"model": "sp", "select": SELECT}, "output": {"dispersion": {"k_path": None}}})
        names = list(ma["samb"]["matrix_info"]["matrix"])
        parameter = {z: 0.1 * (no + 1) * (-1) ** no for no, z in enumerate(names)}
        yield tmp_path, parameter
    finally:
        os.chdir(cwd)


# ==================================================
def run(parameter, output, grid=(4, 3, 2), select=SELECT):
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "sp", "select": select, "parameter": parameter}, "grid": grid, "output": output})
    return ma


# ==================================================
def occupation(ma, Ek):
    # step function at T = 0, where the states at the fermi energy share the remaining electrons.
    x = Ek - ma["info"]["fermi_energy"]
    return np.where(x < -1e-8, 1.0, np.where(x > 1e-8, 0.0, ma.get_occupation_at_fermi_energy()))


# ==================================================
def brute_force(ma, tags):
    # (1/N_k) sum_(k,n) f <u_kn|Z_j(k)|u_kn> with Z_j(k) on the full k grid.
    kv, Ek, Uk = ma.get_eigen_full()
    f = occupation(ma, Ek)
    return {tag: float(np.sum(f * operator_expectation(Uk, ma.get_samb_operator_k(tag, kv)))) / len(kv) for tag in tags}


# ==================================================
def test_samb_orthonormal(sp_model):
    # the SAMBs are Hermitian and orthonormal, Tr[Z_i Z_j] = delta_ij, which is used for <Z_j> = Tr[rho Z_j].
    _, parameter = sp_model
    ma = run(parameter, {"dispersion": {"k_path": None}})
    Z = {tag: ma.get_samb_operator(tag) for tag in parameter}
    for a, Za in Z.items():
        for (R, m, n), v in Za.items():
            assert Za[((-R[0], -R[1], -R[2]), n, m)] == pytest.approx(np.conj(v))
        for b, Zb in Z.items():
            tr = sum(v * Zb.get(((-R[0], -R[1], -R[2]), n, m), 0.0) for (R, m, n), v in Za.items())
            assert tr == pytest.approx(float(a == b), abs=1e-12)


# ==================================================
@pytest.mark.parametrize("n_electron", [None, 3, 4.5, 5])
@pytest.mark.parametrize("tb_gauge", [True, False])
def test_z_expectation(sp_model, n_electron, tb_gauge):
    # <Z_j> agrees with the direct sum over k and occupied states, and sum_j z_j <Z_j> is the band energy.
    topdir, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "z_expectation": "all", "n_electron": n_electron}
    ma = run(parameter, output | {"fourier": {"tb_gauge": tb_gauge}})
    d = ma["output"]["z_expectation"]
    assert list(d["expectation"]) == list(parameter)
    assert d["expectation"] == pytest.approx(brute_force(ma, parameter), abs=1e-12)

    kv, Ek, _ = ma.get_eigen_full()
    f = occupation(ma, Ek)
    assert d["band_energy"] == pytest.approx(np.sum(f * Ek) / len(kv))
    assert d["z_sum"] == pytest.approx(d["band_energy"], abs=1e-12)
    assert d["n_k"] == 24 and d["grid"] == [4, 3, 2]
    assert d["fermi_energy"] == ma["info"]["fermi_energy"] and d["n_electron"] == n_electron
    if n_electron is not None:
        # the number of electrons is n_electron even with degenerate states at the fermi energy.
        assert np.sum(f) / len(kv) == pytest.approx(n_electron)

    # independent of the gauge.
    other = run(parameter, output | {"fourier": {"tb_gauge": not tb_gauge}})
    assert other["output"]["z_expectation"]["expectation"] == pytest.approx(d["expectation"], abs=1e-12)

    # output file.
    saved = read_dict(str(topdir / "out" / "sp" / "output" / "sp_z_expectation.py"))
    assert saved["expectation"] == pytest.approx(d["expectation"])
    assert saved["z_sum"] == pytest.approx(d["z_sum"])


# ==================================================
def test_z_expectation_list(sp_model):
    # a list of SAMBs (duplicates removed), with the eigen system in chunks or on the full grid (save_k_full).
    _, parameter = sp_model
    output = {"dispersion": {"k_path": None}, "z_expectation": ["z3", "z1", "z3"], "n_electron": 4}
    ma = run(parameter, output)
    d = ma["output"]["z_expectation"]
    assert list(d["expectation"]) == ["z3", "z1"]
    assert d["z_sum"] is None  # not all z_j of H are computed.
    full = run(parameter, output | {"save_k_full": True})["output"]["z_expectation"]
    assert full["expectation"] == pytest.approx(d["expectation"], abs=1e-12)
    assert d["expectation"] == pytest.approx(brute_force(ma, ["z3", "z1"]), abs=1e-12)

    # the number of electrons, from the SAMB proportional to the identity in the s orbital at B.
    ma = run(parameter, {"dispersion": {"k_path": None}, "z_expectation": "all", "n_electron": 10})
    kv, Ek, Uk = ma.get_eigen_full()
    assert np.sum(occupation(ma, Ek)) / len(kv) == pytest.approx(10.0)  # all bands are occupied.
    for tag, value in ma["output"]["z_expectation"]["expectation"].items():
        Z = ma.get_samb_operator(tag)
        trace = sum(v for (R, m, n), v in Z.items() if R == (0, 0, 0) and m == n)
        assert value == pytest.approx(trace.real, abs=1e-12)  # <Z_j> = Tr[Z_j] for the fully occupied bands.

    # nothing is computed by default.
    ma = run(parameter, {"dispersion": {"k_path": None}})
    assert "z_expectation" not in ma["output"]


# ==================================================
def test_density_matrix_r():
    # sum_k Tr[P(k) O(k)] = sum_(R,m,n) O(R)_mn A[R,m,n].
    rng = np.random.default_rng(0)
    d, nk = 3, 5
    atom = rng.random((d, 3))
    kv = rng.random((nk, 3))
    R = [(0, 0, 0), (1, 0, 0), (0, -1, 1)]
    O = {(r, m, n): complex(*rng.normal(size=2)) for r in R for m in range(d) for n in range(d)}
    Uk = np.linalg.qr(rng.normal(size=(nk, d, d)) + 1j * rng.normal(size=(nk, d, d)))[0]
    f = rng.random((nk, d))
    for s in [True, False]:
        Ok = fourier_r_to_k(O, atom, kv, s)
        expected = np.sum(f * np.einsum("kmi,kmn,kni->ki", Uk.conj(), Ok, Uk))
        A = density_matrix_r(Uk, f, atom, kv, R, s)
        assert sum(v * A[R.index(r), m, n] for (r, m, n), v in O.items()) == pytest.approx(expected)


# ==================================================
def test_dispersion_z(sp_model):
    # columns of SAMBs after those of local operators, for each band and k point.
    _, parameter = sp_model
    ma = run(parameter, {"dispersion": K_PATH | {"local": ["Sz"], "z": ["z2", "z5"]}, "save_k_path": True})
    path = os.path.join("out", "sp", "output")
    with open(os.path.join(path, "sp_dispersion.txt")) as f:
        assert f.readline() == "# k Energy Sz z2 z5\n"
    data = np.loadtxt(os.path.join(path, "sp_dispersion.txt"))
    z = np.load(os.path.join(path, "sp_eigen_k_path.npz"))
    n_k, n_band = z["Ek"].shape
    assert data.shape == (n_k * n_band, 5)
    kv = ma._k_path_info["k_grid"]
    for no, tag in enumerate(["z2", "z5"]):
        column = data[:, 3 + no].reshape(n_band, n_k).T
        expected = operator_expectation(z["Uk"], ma.get_samb_operator_k(tag, kv))
        assert column == pytest.approx(expected, abs=1e-10)
        # the band sum is Tr[Z_j(k)].
        assert column.sum(axis=1) == pytest.approx(np.trace(ma.get_samb_operator_k(tag, kv), axis1=1, axis2=2).real)

    # SAMBs only.
    run(parameter, {"dispersion": K_PATH | {"z": ["z1"]}})
    with open(os.path.join(path, "sp_dispersion.txt")) as f:
        assert f.readline() == "# k Energy z1\n"


# ==================================================
def test_callback_info_z(sp_model):
    # SAMBs in callback_info, together with local operators.
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    output = {"dispersion": {"k_path": None}, "callback_k_full": "cb.py", "callback_info": ["Sz", "z4"]}
    ma = run(parameter, output)
    with open(topdir / "out" / "full.pkl", "rb") as f:
        Ek, Uk, info = pickle.load(f)
    kv = info["k_grid"]
    assert info["z4"].shape == Ek.shape
    assert np.allclose(info["z4"], operator_expectation(Uk, ma.get_samb_operator_k("z4", kv)))
    assert np.allclose(info["Sz"], np.einsum("kmi,mn,kni->ki", Uk.conj(), ma.get_local_operator("Sz"), Uk).real)


# ==================================================
@pytest.mark.parametrize(
    "output, match",
    [
        ({"z_expectation": ["z99"]}, r"'output/z_expectation': \['z99'\] are not in the selected SAMBs"),
        ({"z_expectation": "z1"}, r"'output/z_expectation' must be a list of SAMB names .* or \"all\""),
        ({"z_expectation": [1]}, r"'output/z_expectation' must be a list of SAMB names"),
        ({"dispersion": {"k_path": None, "z": "all"}}, r"'output/dispersion/z' must be a list of SAMB names"),
        ({"dispersion": {"k_path": None, "z": ["z1", "z0"]}}, r"\['z0'\] are not in the selected SAMBs"),
        ({"callback_info": ["z99"], "callback_k_full": "cb.py"}, r"'output/callback_info': \['z99'\] are not in"),
        ({"callback_info": ["zx"], "callback_k_full": "cb.py"}, r"'output/callback_info' must be a list of"),
        ({"callback_info": [""], "callback_k_full": "cb.py"}, r"'output/callback_info' must be a list of"),
    ],
)
def test_z_error(sp_model, output, match):
    topdir, parameter = sp_model
    (topdir / "out" / "cb.py").write_text(CALLBACK)
    with pytest.raises(ValueError, match=match):
        run(parameter, {"dispersion": {"k_path": None}} | output)


# ==================================================
def test_z_not_selected(sp_model):
    # only the selected SAMBs can be used, whose names are those of all SAMBs of the model.
    _, _ = sp_model
    select = {"X": ["Q"], "Gamma": ["A1g", "Eg"]}
    ma = ModelAnalyzer("out")
    control = {"samb": {"model": "sp", "select": select, "parameter": {"z1": 1.0}}, "grid": (2, 2, 2)}
    ma.analyze(control | {"output": {"dispersion": {"k_path": None}, "z_expectation": ["z14", "z38"]}})
    assert list(ma["output"]["z_expectation"]["expectation"]) == ["z14", "z38"]
    assert "z15" not in ma["samb"]["matrix_info"]["matrix"]
    with pytest.raises(ValueError, match=r"\['z15'\] are not in the selected SAMBs"):
        ma.analyze(control | {"output": {"dispersion": {"k_path": None}, "z_expectation": ["z15"]}})


# ==================================================
def test_z_point_group(tmp_path):
    # model of a point group: Γ point only.
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        model = {"model": "mol", "group": "D3h", "site": {"A": ("[1,0,0]", "s")}, "bond": [("A", "A", [1])]}
        create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir="out")
        ma = ModelAnalyzer("out")
        ma.analyze({"samb": {"model": "mol"}, "output": {"dispersion": {"k_path": None}}})
        parameter = {z: 0.3 * (no + 1) for no, z in enumerate(ma["samb"]["matrix_info"]["matrix"])}
        ma.analyze({"samb": {"model": "mol", "parameter": parameter}, "output": {"z_expectation": "all", "n_electron": 1}})
        d = ma["output"]["z_expectation"]
        assert d["n_k"] == 1 and d["grid"] is None
        assert d["z_sum"] == pytest.approx(d["band_energy"])
        # one electron in the lowest level (shared by the degenerate states).
        H = sum(c * ma.get_samb_operator_k(z, np.zeros((1, 3)))[0] for z, c in parameter.items())
        assert d["band_energy"] == pytest.approx(np.min(np.linalg.eigvalsh(H)))
    finally:
        os.chdir(cwd)
