# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests of the conversion of a model of Ver.2 into Ver.3 (convert_model, mp_convert).

The models of Ver.2 are in tests/data/convert_v2 (see its README.md).
"""

import os
import shutil

import numpy as np
import pytest
from click.testing import CliRunner

from multipie import MaterialModel, convert_model
from multipie.scripts import mp_convert
from multipie.util.util import read_dict, str_to_sympy
from multipie.util import util_convert
from multipie.util.util_convert import model_version

DATA = os.path.join(os.path.dirname(__file__), "data", "convert_v2")
OPT = {"pdf": False, "qtdraw": False}


# ==================================================
@pytest.fixture
def v2(tmp_path, monkeypatch):
    """
    Copy of the models and controls of Ver.2 in a temporary directory (cwd).
    """
    shutil.copytree(DATA, tmp_path, dirs_exist_ok=True)
    monkeypatch.chdir(tmp_path)
    return tmp_path


# ==================================================
def hr(data, parameter):
    """
    H(R) = sum_j z_j Z_j of a model (data of Ver.2 or Ver.3), dict[(n1, n2, n3, m, n), value].
    """
    mm = MaterialModel()
    mm.clear()
    mm.update(data)
    Z = mm._combined_samb_matrix(mm["combined_id"], fmt="value", bond=False)
    H = {}
    for z, v in parameter.items():
        v = float(str_to_sympy(v, rational=False)) if isinstance(v, str) else v
        for Rmn, x in Z[z].items():
            H[Rmn] = H.get(Rmn, 0) + v * x
    return H


# ==================================================
def max_diff(H1, H2):
    return max(abs(H1.get(k, 0) - H2.get(k, 0)) for k in set(H1) | set(H2))


# ==================================================
@pytest.mark.parametrize(
    "name, spinful, n_samb",
    [("graphene", False, 4), ("ch4", True, 120), ("pt_jml", True, 67), ("ce_f", False, 28)],
)
def test_convert_model(v2, name, spinful, n_samb):
    version, data2 = model_version(f"{name}/{name}.pkl")
    report = convert_model(name, **OPT)

    # the model of Ver.2 is kept, and the model of Ver.3 is created from the rebuilt input.
    assert os.path.isfile(f"{name}/{name}_v2.pkl")
    assert model_version(f"{name}/{name}_v2.pkl")[0] == version
    assert model_version(f"{name}/{name}.pkl")[0].startswith("3.")
    model_in = read_dict(f"{name}_in.py")
    assert model_in["model"] == name and model_in["spinful"] == spinful
    assert "max_neighbor" not in model_in and "search_cell_range" not in model_in

    mm = MaterialModel()
    mm.load(name)
    assert mm["full_matrix"]["ket"] == [k[:4] + [util_convert.renamed_orbital(k[4])] for k in data2["full_matrix"]["ket"]]
    assert report["SAMB"] == (n_samb, n_samb)
    assert report["span_residual"] < 1e-10
    assert report["neighbor"] == {}
    assert report["removed_key"].keys() == {"max_neighbor", "search_cell_range"}


# ==================================================
def test_renamed_orbital(v2):
    report = convert_model("ce_f", **OPT)
    assert report["renamed_orbital"] == {"Ce": [("f2", "fc"), ("f1", "fs"), ("f3x", "fx"), ("f3y", "fy")]}
    assert read_dict("ce_f_in.py")["site"]["Ce"] == ("[0,0,0]", ["fc", "fs", "fbz", "f3", "fx", "fy", "faz"])


# ==================================================
def test_samb_not_orthonormal(v2):
    # a part of the f orbitals: the SAMB are not orthonormal (in Ver.2 and Ver.3), and the weights are fitted.
    _, data2 = model_version("ce_f_part/ce_f_part.pkl")
    parameter2 = {"z1": 1.0, "z2": -0.4, "z4": 0.3, "z12": 0.2}
    (v2 / "ce_f_part_ctrl.py").write_text(f"ce_f_part_ctrl = {{'samb': {{'model': 'ce_f_part', 'parameter': {parameter2!r}}}}}\n")
    report = convert_model("ce_f_part", "ce_f_part_ctrl", **OPT)
    assert report["span_residual"] < 1e-10
    assert report["control"]["ce_f_part_ctrl.py"]["residual"] < 1e-10

    mm = MaterialModel()
    mm.load("ce_f_part")
    assert not util_convert.ModelConverter(data2, mm).orthonormal

    _, data3 = model_version("ce_f_part/ce_f_part.pkl")
    parameter3 = read_dict("ce_f_part_ctrl_v3.py")["samb"]["parameter"]
    assert max_diff(hr(data2, parameter2), hr(data3, parameter3)) < 1e-10


# ==================================================
def test_control(v2):
    _, data2 = model_version("graphene/graphene.pkl")
    parameter2 = read_dict("graphene_ctrl.py")["samb"]["parameter"]
    report = convert_model("graphene", "graphene_ctrl.py", **OPT)

    r = report["control"]["graphene_ctrl.py"]
    assert r["output"] == os.path.join(str(v2), "graphene_ctrl_v3.py")
    assert r["residual"] < 1e-10
    ctrl3 = read_dict("graphene_ctrl_v3.py")
    assert ctrl3["samb"]["model"] == "graphene"
    parameter3 = ctrl3["samb"]["parameter"]
    assert all(isinstance(v, float) for v in parameter3.values())

    # H = sum_j z_j Z_j is the same with the SAMB and weights of Ver.2 and of Ver.3.
    _, data3 = model_version("graphene/graphene.pkl")
    assert max_diff(hr(data2, parameter2), hr(data3, parameter3)) < 1e-10

    # the original control is not changed.
    assert read_dict("graphene_ctrl.py")["samb"]["parameter"] == parameter2


# ==================================================
def test_neighbor_missed_by_ver2(v2):
    # Ver.2.3.4 missed a part of the shells: its 9th neighbor is the 11th of Ver.3.
    _, data2 = model_version("cubic_s/cubic_s.pkl")
    assert sorted({b.neighbor for b in data2["bond"]["representative"].values()}) == list(range(1, 10))
    parameter2 = read_dict("cubic_s_ctrl.py")["samb"]["parameter"]

    report = convert_model("cubic_s", ["cubic_s_ctrl"], **OPT)
    assert report["neighbor"] == {("A", "A"): (list(range(1, 10)), [1, 2, 3, 4, 5, 6, 7, 8, 11])}
    assert read_dict("cubic_s_in.py")["bond"] == [("A", "A", [1, 2, 3, 4, 5, 6, 7, 8, 11])]
    n2, n3 = report["SAMB"]
    assert n3 > n2  # the 8th shell of Ver.3 has the bonds of a Wyckoff position that Ver.2 missed.
    assert report["span_residual"] < 1e-10

    _, data3 = model_version("cubic_s/cubic_s.pkl")
    parameter3 = read_dict("cubic_s_ctrl_v3.py")["samb"]["parameter"]
    assert max_diff(hr(data2, parameter2), hr(data3, parameter3)) < 1e-10


# ==================================================
def test_parameter(v2):
    _, data2 = model_version("graphene/graphene.pkl")
    parameter2 = {"z1": 0.1, "z2": -1.0, "z3": "1/5"}
    report = convert_model("graphene", parameter=parameter2, **OPT)
    assert report["parameter_residual"] < 1e-10
    _, data3 = model_version("graphene/graphene.pkl")
    assert max_diff(hr(data2, parameter2), hr(data3, report["parameter"])) < 1e-10

    with pytest.raises(ValueError, match="not SAMB of the model of Ver.2"):
        convert_model("graphene", parameter={"z100": 1.0}, **OPT)


# ==================================================
def test_control_of_converted_model(v2):
    convert_model("graphene", **OPT)
    created = os.path.getmtime("graphene/graphene.pkl")

    # only the control is converted for a converted model.
    report = convert_model("graphene", "graphene_ctrl", **OPT)
    assert report["input"] is None
    assert os.path.getmtime("graphene/graphene.pkl") == created
    assert report["control"]["graphene_ctrl.py"]["residual"] < 1e-10


# ==================================================
def test_refuse(v2):
    # an existing input file is not overwritten, and nothing is changed.
    (v2 / "graphene_in.py").write_text("x = {}\n")
    with pytest.raises(FileExistsError, match="graphene_in.py"):
        convert_model("graphene", **OPT)
    assert not os.path.exists("graphene/graphene_v2.pkl")

    # a control of another model.
    with pytest.raises(ValueError, match="samb/model"):
        convert_model("ch4", "graphene_ctrl", **OPT)
    assert not os.path.exists("ch4/ch4_v2.pkl")

    # a model of Ver.3 without the model of Ver.2.
    convert_model("ch4", **OPT)
    os.remove("ch4/ch4_v2.pkl")
    with pytest.raises(ValueError, match="no model of Ver.2"):
        convert_model("ch4", **OPT)

    with pytest.raises(FileNotFoundError):
        convert_model("unknown", **OPT)


# ==================================================
def test_restore_on_failure(v2, monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("the SAMB of Ver.2 are not in the space of those of Ver.3")

    monkeypatch.setattr(util_convert.ModelConverter, "__init__", fail)
    with pytest.raises(ValueError, match="not in the space") as e:
        convert_model("graphene", **OPT)
    assert any("restored" in note for note in e.value.__notes__)
    # the model of Ver.2 is restored.
    assert model_version("graphene/graphene.pkl")[0] == "2.5.2"
    assert not os.path.exists("graphene/graphene_v2.pkl")


# ==================================================
def test_converter_refuses_other_space(v2):
    # the SAMB of the model of Ver.2 are compared with a model of Ver.3 that lacks the 3rd neighbor.
    _, data2 = model_version("graphene/graphene.pkl")
    model_in, _ = util_convert.rebuild_model_input(data2)
    model_in["bond"] = [("C", "C", [1, 2])]
    model_in |= {"pdf": {"create": False}, "qtdraw": {"create": False}, "model": "graphene12"}
    from multipie import create_model

    create_model(model_in)
    mm = MaterialModel()
    mm.load("graphene12")
    with pytest.raises(ValueError, match="not in the space"):
        util_convert.ModelConverter(data2, mm)


# ==================================================
def test_mp_convert(v2):
    result = CliRunner().invoke(mp_convert.cmd, ["--no-pdf", "--no-qtdraw", "cubic_s", "cubic_s_ctrl"])
    assert result.exit_code == 0, result.output
    assert "neighbors of A;A" in result.output
    assert "removed key: max_neighbor = 10" in result.output
    assert os.path.isfile("cubic_s_ctrl_v3.py")

    result = CliRunner().invoke(mp_convert.cmd, ["unknown"])
    assert result.exit_code == 1
    assert "FileNotFoundError" in result.output and "Traceback" not in result.output
