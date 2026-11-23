# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests of the order of the combined SAMB (from 3.2.0), and of the renumbering of a model of Ver.3 before 3.2.0
(convert_model, mp_convert).

The models of MultiPie 3.1.3 are in tests/data/convert_v3 (see its README.md).
"""

import os
import shutil

import pytest
from click.testing import CliRunner

from multipie import MaterialModel, Group, __version__, convert_model
from multipie.scripts import mp_convert
from multipie.util.util import read_dict, str_to_sympy
from multipie.util.util_convert import model_version, needs_renumber, version_tuple
from multipie.util.util_dict import Dict
from multipie.util.util_material_model import parse_samb_select

DATA = os.path.join(os.path.dirname(__file__), "data", "convert_v3")
EM = {"Q": 0, "G": 0, "T": 1, "M": 1}
XO = {"Q": 0, "G": 1, "T": 2, "M": 3}


# ==================================================
@pytest.fixture
def v3(tmp_path, monkeypatch):
    """
    Copy of the models and controls of 3.1.3 in a temporary directory (cwd).
    """
    shutil.copytree(DATA, tmp_path, dirs_exist_ok=True)
    monkeypatch.chdir(tmp_path)
    return tmp_path


# ==================================================
def load(topdir, name):
    mm = MaterialModel(str(topdir))
    mm.load(name)
    return mm


# ==================================================
def data_model(data):
    mm = MaterialModel()
    mm.clear()
    mm.update(data)
    return mm


# ==================================================
def hr(mm, parameter):
    """
    H(R) = sum_j z_j Z_j, dict[(n1, n2, n3, m, n), value].
    """
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
def order_key(mm, v):
    # the order of the combined SAMB: Gamma (character table), site/bond block, electric/magnetic, s, k, l, X, n, p,
    # component.
    irreps = list(mm.group.character["table"].keys())
    _, u, key, comp = v
    c = u.samb_type
    block = (0 if u.neighbor == 0 else 1, c.tail, c.head, u.neighbor, u.n, c)
    X, l, Gamma, n, p, s, k, _ = key
    return (irreps.index(Gamma), block, EM[X], s, k, l, XO[X], n, p, comp)


# ==================================================
@pytest.mark.parametrize("name", ["graphene", "cubic_pd"])
def test_order(name):
    # z_j of a model created now are in the order of the rule.
    mm = MaterialModel(verbose=False)
    mm.analyze(read_dict(f"{name}_in.py", DATA))
    cid = mm["combined_id"]
    keys = [order_key(mm, v) for v in cid.values()]
    assert keys == sorted(keys)
    assert list(cid) == [f"z{i}" for i in range(1, len(cid) + 1)]


# ==================================================
def test_order_irreps():
    # the irreps of the selection and of the cluster SAMB follow the character table (not the alphabetical order).
    irreps = list(Group("PG:32").character["table"].keys())
    assert parse_samb_select({}, irreps)["Gamma"] == irreps
    assert parse_samb_select({"Gamma": ["T1u", "A1u", "A1g"]}, irreps)["Gamma"] == ["A1g", "A1u", "T1u"]
    assert parse_samb_select({"Gamma": "IR"}, irreps)["Gamma"] == ["A1g"]

    mm = MaterialModel(verbose=False)
    mm.analyze(read_dict("graphene_in.py", DATA))
    irreps = list(mm.group.character["table"].keys())
    for wp in {wp for wp, _, _ in mm["cluster_id"].values()}:
        g = [idx[2] for w, idx, _ in mm["cluster_id"].values() if w == wp]
        assert [irreps.index(i) for i in g] == sorted(irreps.index(i) for i in g)


# ==================================================
def test_dict_sort_rank():
    # a dict as the custom order: values of the same rank are not distinguished.
    from collections import namedtuple

    T = namedtuple("T", ["X", "l"])
    d = Dict(T, {T("M", 0): 1, T("Q", 2): 2, T("T", 1): 3, T("G", 1): 4, T("A", 0): 5})
    keys = list(d.sort(("X", {"Q": 0, "G": 0, "T": 1, "M": 1}), "l").keys())
    assert keys == [T("G", 1), T("Q", 2), T("M", 0), T("T", 1), T("A", 0)]


# ==================================================
def test_version():
    assert version_tuple("3.1.3") == (3, 1) and version_tuple("10.0") == (10, 0) and version_tuple("x") == (0, 0)
    assert needs_renumber("3.0.0") and needs_renumber("3.1.3")
    assert not needs_renumber("2.5.2") and not needs_renumber("3.2.0") and not needs_renumber(__version__)


# ==================================================
@pytest.mark.parametrize("name, n_samb", [("graphene", 20), ("cubic_pd", 256)])
@pytest.mark.parametrize("with_input", [True, False])
def test_renumber(v3, name, n_samb, with_input):
    # the model is created again (from model_in.py, or from the input rebuilt from the old model), the old model is
    # kept, and the renumbered weights give the same Hamiltonian.
    if not with_input:
        os.remove(v3 / f"{name}_in.py")
    _, old = model_version(v3 / name / f"{name}.pkl")
    report = convert_model(name, [f"{name}_ctrl.py"])

    assert report["kind"] == "renumber"
    assert report["SAMB"] == (n_samb, n_samb) and report["span_residual"] < 1e-12
    assert os.path.isfile(v3 / name / f"{name}_v3.1.pkl") and os.path.isfile(v3 / f"{name}_in.py")
    assert model_version(v3 / name / f"{name}_v3.1.pkl")[0] == "3.1.3"
    assert model_version(v3 / name / f"{name}.pkl")[0] == __version__
    assert sorted(report["renumber"].values(), key=lambda z: int(z[1:])) == [f"z{i}" for i in range(1, n_samb + 1)]
    assert any(a != b for a, b in report["renumber"].items())

    new = load(v3, name)
    old = data_model(old)
    for z_old, z_new in report["renumber"].items():
        a, b = old["combined_id"][z_old], new["combined_id"][z_new]
        assert a[1:] == b[1:]

    ctrl = read_dict(f"{name}_ctrl.py", str(v3))
    out = report["control"][f"{name}_ctrl.py"]["output"]
    assert out == str(v3 / f"{name}_ctrl_v3.2.py")
    ctrl3 = read_dict(out, str(v3))
    assert ctrl3["samb"]["model"] == name
    p_old, p_new = ctrl["samb"]["parameter"], ctrl3["samb"]["parameter"]
    assert len(p_old) == len(p_new)
    assert max_diff(hr(old, p_old), hr(new, p_new)) < 1e-12
    assert list(p_new) == sorted(p_new, key=lambda z: int(z[1:]))

    # the other SAMB in the control.
    rn = report["renumber"]
    o_old, o_new = ctrl.get("output", {}), ctrl3.get("output", {})
    if "callback_info" in o_old:
        assert o_new["callback_info"] == [rn.get(i, i) for i in o_old["callback_info"]]
        assert o_new["dispersion"]["z"] == [rn[i] for i in o_old["dispersion"]["z"]]
    if isinstance(o_old.get("z_expectation"), list):
        assert o_new["z_expectation"] == [rn[i] for i in o_old["z_expectation"]]
    else:
        assert o_new.get("z_expectation") == o_old.get("z_expectation")


# ==================================================
def test_renumber_graphene_values(v3):
    # the renumbering of graphene: the irreps in the order of the character table (E2g before the u irreps).
    report = convert_model("graphene", ["graphene_ctrl.py"])
    rn = report["renumber"]
    assert [rn[f"z{i}"] for i in range(1, 5)] == ["z1", "z2", "z3", "z4"]
    assert [rn[f"z{i}"] for i in range(15, 21)] == [f"z{i}" for i in range(5, 11)]  # E2g.
    assert [rn[f"z{i}"] for i in range(5, 15)] == [f"z{i}" for i in range(11, 21)]  # B1u, B2u, E1u.
    ctrl3 = read_dict("graphene_ctrl_v3.2.py", str(v3))
    assert ctrl3["output"]["callback_info"] == ["Sz", "z5"]


# ==================================================
def test_parameter_and_converted_model(v3):
    # the weights of parameter, and the conversion of the controls of a model already converted (only the weights).
    report = convert_model("graphene", parameter={"z15": 1.0, "z5": "sqrt(2)"})
    assert report["parameter"] == {"z5": 1.0, "z11": float(str_to_sympy("sqrt(2)", rational=False))}
    assert report["parameter_residual"] == 0.0
    report2 = convert_model("graphene", ["graphene_ctrl.py"])
    assert report2["kind"] == "renumber" and report2["input"] is None
    assert report2["renumber"] == report["renumber"]
    assert os.path.isfile(v3 / "graphene_ctrl_v3.2.py")


# ==================================================
def test_refuse(v3, monkeypatch):
    # an input that gives other SAMB is refused, and the old model is restored.
    model_in = read_dict("graphene_in.py", str(v3))
    model_in["bond"] = [("C", "C", [1])]
    from multipie.util.util import write_dict

    write_dict(model_in, "graphene_in.py", var="graphene", w_dir=str(v3))
    with pytest.raises(ValueError, match="do not agree"):
        convert_model("graphene", ["graphene_ctrl.py"])
    assert model_version(v3 / "graphene" / "graphene.pkl")[0] == "3.1.3"
    assert not os.path.exists(v3 / "graphene" / "graphene_v3.1.pkl")
    assert not os.path.exists(v3 / "graphene_ctrl_v3.2.py")

    # an unknown z_j.
    os.remove(v3 / "graphene_in.py")
    with pytest.raises(ValueError, match="not SAMB of the old model"):
        convert_model("graphene", parameter={"z99": 1.0})

    # an old model kept already.
    shutil.copy(v3 / "cubic_pd" / "cubic_pd.pkl", v3 / "cubic_pd" / "cubic_pd_v3.1.pkl")
    with pytest.raises(FileExistsError):
        convert_model("cubic_pd")

    # both an old model of Ver.2 and of Ver.3 for a converted model.
    shutil.copy(v3 / "graphene" / "graphene_v3.1.pkl", v3 / "graphene" / "graphene_v2.pkl")
    with pytest.raises(ValueError, match="both"):
        convert_model("graphene")


# ==================================================
def test_mp_convert(v3):
    result = CliRunner().invoke(mp_convert.cmd, ["graphene", "graphene_ctrl.py"])
    assert result.exit_code == 0, result.output
    assert "renumbered (" in result.output and "old model kept as:" in result.output
    assert os.path.isfile(v3 / "graphene_ctrl_v3.2.py")
