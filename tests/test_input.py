"""
Tests for input keys, parameter files, PDF failures, and the number of parallel jobs.
"""

import os
import shutil

import pytest

from multipie import MaterialModel
from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import get_n_jobs, write_dict

MODEL = {
    "model": "chain",
    "group": 1,
    "cell": {"a": 1.0, "b": 5.0, "c": 5.0},
    "site": {"Fe-1": ("[0,0,0]", "s")},  # site names are free (except for '_' and ';').
    "bond": [("Fe-1", "Fe-1", [1])],
    "pdf": {"create": False},
    "qtdraw": {"create": False},
}


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
@pytest.mark.parametrize(
    "model, match",
    [
        ({"spinfull": True}, r"'spinfull' \(did you mean 'spinful'\?\)"),
        ({"cell": {"a": 1.0, "alpha_": 90.0}}, r"'cell/alpha_'"),
        ({"SAMB_select": {"Gama": ["A"]}}, r"'SAMB_select/Gama' \(did you mean 'Gamma'\?\)"),
        ({"qtdraw": {"creat": False}}, r"'qtdraw/creat' \(did you mean 'create'\?\)"),
    ],
)
def test_unknown_model_key(model, match):
    with pytest.raises(ValueError, match=match):
        MaterialModel().analyze(MODEL | model)


# ==================================================
@pytest.mark.parametrize(
    "name, bond, match",
    [
        ("A;B", [], "without ';'"),
        ("", [], "without ';'"),
        (None, [], "without ';'"),
        (1, [], "without ';'"),
        ("my_site", [("my_site", "my_site", [1])], "used in bond must not contain '_'"),
    ],
)
def test_invalid_site_name(name, bond, match):
    with pytest.raises(ValueError, match=match):
        MaterialModel().analyze(MODEL | {"site": {name: ("[0,0,0]", "s")}, "bond": bond})


# ==================================================
def test_site_name_with_underscore_without_bond():
    # '_' is allowed for a site not used in bond.
    mm = MaterialModel()
    mm.analyze(MODEL | {"site": {"Fe_1": ("[0,0,0]", "s")}, "bond": []})
    assert "Fe_1" in mm["site"]["representative"]


# ==================================================
@pytest.mark.parametrize(
    "control, match",
    [
        ({"mdoe": "samb"}, r"'mdoe' \(did you mean 'mode'\?\)"),
        ({"samb": {"model": "chain", "selct": {}}}, r"'samb/selct' \(did you mean 'select'\?\)"),
        ({"samb": {"model": "chain", "select": {"gamma": ["A"]}}}, r"'samb/select/gamma'"),
        ({"output": {"dispersion": {"kpath": ""}}}, r"'output/dispersion/kpath' \(did you mean 'k_path'\?\)"),
    ],
)
def test_unknown_control_key(in_tmp, control, match):
    create_model(dict(MODEL), topdir="out")
    with pytest.raises(ValueError, match=match):
        ModelAnalyzer("out").analyze(control)


# ==================================================
def test_free_names(in_tmp):
    # site names, SAMB names in samb/parameter, and k-point labels are chosen by the user.
    create_model(dict(MODEL), topdir="out")
    control = {
        "samb": {"model": "chain", "parameter": {"z1": 1.0}},
        "output": {"dispersion": {"k_path": "Gamma0-Xpoint", "k_point": {"Gamma0": "[0,0,0]", "Xpoint": "[1/2,0,0]"}}},
    }
    ma = ModelAnalyzer("out")
    ma.analyze(control)
    assert ma["output"]["dispersion"]["k_path"] == "Gamma0-Xpoint"


# ==================================================
def test_parameter_file(in_tmp):
    # samb/parameter given as a file name under topdir/model/info.
    create_model(dict(MODEL), topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "chain", "parameter": {"z1": 1.0, "z2": -0.5}}, "output": {"dispersion": {"k_path": None}}})
    expected = dict(ma.HR)

    write_dict({"z1": 1.0, "z2": "-1/2"}, "my_z.py", w_dir=str(in_tmp / "out" / "chain" / "info"))
    ma.analyze({"samb": {"model": "chain", "parameter": "my_z.py"}, "output": {"dispersion": {"k_path": None}}})
    assert ma.parameter == {"z1": 1.0, "z2": -0.5}
    assert ma.HR.keys() == expected.keys()
    assert all(complex(ma.HR[k]) == pytest.approx(complex(v)) for k, v in expected.items())


# ==================================================
def test_pdf_failure(in_tmp, monkeypatch, capsys):
    # LaTeX is found, but ptex2pdf is not: the .tex and .pkl files are written, and only the PDF is skipped.
    monkeypatch.setattr("multipie.core.material_model.check_latex", lambda: True)
    which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda cmd: None if cmd == "ptex2pdf" else which(cmd))
    create_model(MODEL | {"pdf": {"create": True}}, topdir="out")  # warning via logging (set up by create_model).
    path = in_tmp / "out" / "chain"
    assert os.path.isfile(path / "chain.pkl")
    assert os.path.isfile(path / "chain.tex")
    assert not os.path.exists(path / "chain.pdf")
    out = capsys.readouterr()
    assert "skip PDF creation: ptex2pdf is not found" in out.out + out.err


# ==================================================
@pytest.mark.parametrize("value, n", [(None, -1), ("", -1), (" 4 ", 4), ("1", 1), ("-2", -2)])
def test_n_jobs(monkeypatch, value, n):
    if value is None:
        monkeypatch.delenv("MULTIPIE_N_JOBS", raising=False)
    else:
        monkeypatch.setenv("MULTIPIE_N_JOBS", value)
    assert get_n_jobs() == n


@pytest.mark.parametrize("value, match", [("0", "must not be 0"), ("two", "must be integer")])
def test_n_jobs_invalid(monkeypatch, value, match):
    monkeypatch.setenv("MULTIPIE_N_JOBS", value)
    with pytest.raises(ValueError, match=match):
        get_n_jobs()
