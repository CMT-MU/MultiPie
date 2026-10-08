"""
Regression tests for the command-line tools and input error messages.
"""

import os
import subprocess
import sys

import pytest
from click.testing import CliRunner

import multipie
from multipie import MaterialModel
from multipie.core.cmd import analyze_model, create_model
from multipie.scripts import mp_analyze, mp_create
from multipie.util.util import read_dict, str_to_sympy

# run subprocesses with the same multipie as this test (installed or source checkout).
PACKAGE_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(multipie.__file__)))

MODEL = {
    "model": "chain",
    "group": 1,
    "cell": {"a": 1.0, "b": 5.0, "c": 5.0},
    "site": {"A": ("[0,0,0]", "s")},
    "bond": [("A", "A", [1])],
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
    finally:
        os.chdir(cwd)


def model(**kwargs):
    return dict(MODEL) | kwargs


def read_dict_text(d, text):
    (d / "template.py").write_text(text)
    return read_dict("template.py")


# ==================================================
@pytest.mark.parametrize("script", [mp_create.cmd, mp_analyze.cmd])
def test_missing_argument_exit_status(script):
    result = CliRunner().invoke(script, [])
    assert result.exit_code == 2
    assert "Missing argument" in result.output


# ==================================================
@pytest.mark.parametrize("script", [mp_create.cmd, mp_analyze.cmd])
def test_input_template_is_valid_input(in_tmp, script):
    result = CliRunner().invoke(script, ["-i"])
    assert result.exit_code == 0
    (in_tmp / "template.py").write_text(result.output)
    assert isinstance(read_dict("template.py"), dict)


# ==================================================
def test_template_creates_model(in_tmp):
    result = CliRunner().invoke(mp_create.cmd, ["-i"])
    dic = read_dict_text(in_tmp, result.output)
    dic.update(model=MODEL["model"], group=1, cell=MODEL["cell"], site=MODEL["site"], bond=MODEL["bond"])
    dic["pdf"]["create"] = False
    dic["qtdraw"]["create"] = False
    create_model(dic)
    assert os.path.isfile(in_tmp / "chain" / "chain.pkl")


# ==================================================
def test_error_is_reported_in_one_message(in_tmp):
    (in_tmp / "bad.py").write_text(f"bad = {model(site={'A': ('[0,0,0]', 'pzz')})!r}\n")
    result = CliRunner().invoke(mp_create.cmd, ["bad.py"])
    assert result.exit_code == 1
    assert "Traceback" not in result.output
    assert "unknown orbital 'pzz'" in result.output
    assert "in site 'A'" in result.output
    assert "while creating model 'chain'" in result.output
    assert "bad.py" in result.output


# ==================================================
def test_verbose_shows_traceback(in_tmp):
    (in_tmp / "bad.py").write_text(f"bad = {model(site={'A': ('[0,0,0]', 'pzz')})!r}\n")
    result = CliRunner().invoke(mp_create.cmd, ["-v", "bad.py"])
    assert result.exit_code == 1
    assert isinstance(result.exception, ValueError)


# ==================================================
def test_analyze_error_exit_status(in_tmp):
    (in_tmp / "ctrl.py").write_text("ctrl = {'samb': {'model': 'none'}}\n")
    result = CliRunner().invoke(mp_analyze.cmd, ["ctrl.py"])
    assert result.exit_code == 1
    assert "Traceback" not in result.output
    assert "while analyzing model 'none'" in result.output


# ==================================================
def test_create_model_traceback_once(in_tmp, capsys):
    with pytest.raises(ValueError) as e:
        create_model(model(site={"A": ("[0,0,0]", "pzz")}))
    assert "while creating model 'chain'." in e.value.__notes__
    out = capsys.readouterr()
    assert "Traceback" not in out.out + out.err


# ==================================================
@pytest.mark.parametrize(
    "site, match",
    [
        ({"A": ("[0,0,0]", "pzz")}, r"unknown orbital 'pzz', acceptable: p, px, py, pz\."),
        ({"A": ("[0,0,0]", "g")}, "invalid orbital"),
        ({"A": ("[0,0,0]", "(5/2,p)")}, r"acceptable: \(1/2,p\), \(3/2,p\)"),
        ({"A": ("[0,0,0", "s")}, "invalid string"),
        ({}, "no site"),
        ({"A": ("[0,0,0]",)}, "must be \\(position, orbital\\)"),
    ],
)
def test_invalid_site(site, match):
    with pytest.raises(ValueError, match=match):
        MaterialModel().analyze(model(site=site, bond=[]))


# ==================================================
def test_str_to_sympy_unbalanced_bracket():
    with pytest.raises(ValueError, match="invalid string"):
        str_to_sympy("[0,0,0")


# ==================================================
def test_non_literal_input(in_tmp):
    (in_tmp / "m.py").write_text('m = {\n    "a": 1,\n    "b": [1, 2*2],\n}\n')
    with pytest.raises(ValueError, match=r"\(line 3\): '2 \* 2' is not a literal"):
        read_dict("m.py")


# ==================================================
def test_samb_qtdraw_without_qtdraw(in_tmp, monkeypatch, capsys):
    create_model(model(), topdir="out")
    mm = MaterialModel("out", verbose=True)
    mm.load("chain")
    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: False)
    mm.save_samb_qtdraw()
    out = capsys.readouterr()
    assert "QtDraw not found" in out.out + out.err
    assert "save SAMB QtDraw files" not in out.out
    assert not os.path.exists(in_tmp / "out" / "chain" / "samb")


# ==================================================
def test_missing_tools_notice(in_tmp, monkeypatch, capsys):
    create_model(model(), topdir="out")
    mm = MaterialModel("out", verbose=True)
    mm.load("chain")
    mm["qtdraw_prop"]["create"] = True
    mm["pdf_ctrl"]["create"] = True
    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: False)
    monkeypatch.setattr("multipie.core.material_model.check_latex", lambda: False)
    mm.save_view()
    mm.save_pdf()
    out = capsys.readouterr().out
    assert "QtDraw not found" in out
    assert "LaTeX not found" in out


# ==================================================
@pytest.mark.parametrize("module", ["multipie.scripts.mp_create", "multipie.scripts.mp_analyze"])
def test_python_m(module):
    env = os.environ | {"PYTHONPATH": PACKAGE_PARENT}
    result = subprocess.run([sys.executable, "-m", module, "-i"], capture_output=True, text=True, env=env)
    assert result.returncode == 0
    assert result.stdout.startswith("default_")


# ==================================================
@pytest.mark.parametrize(
    "orbital, basis_type, expected",
    [
        ("p", "lg", [[], ["px", "py", "pz"], [], []]),
        (["dxy", "s", "pz"], "lg", [["s"], ["pz"], ["dxy"], []]),
        ("P X", "lg", [[], ["px"], [], []]),
        ("px", "lgs", [[], ["(px,u)", "(px,d)"], [], []]),
        ("s", "lgs", [["(s,u)", "(s,d)"], [], [], []]),
        ("(3/2,p)", "jml", [[], ["(3/2,3/2)", "(3/2,1/2)", "(3/2,-1/2)", "(3/2,-3/2)"], [], []]),
        (["(1/2,s)", "( 1/2 , P )"], "jml", [["(1/2,1/2)", "(1/2,-1/2)"], ["(1/2,1/2)", "(1/2,-1/2)"], [], []]),
    ],
)
def test_valid_orbital(orbital, basis_type, expected):
    from multipie import Group
    from multipie.util.util_material_model import parse_orbital

    basis_info = {k: Group(1).atomic_basis(k) for k in ["jml", "lgs", "lg"]}
    assert parse_orbital(orbital, basis_type, basis_info) == expected


# ==================================================
def test_non_literal_after_valid_call(in_tmp):
    (in_tmp / "m.py").write_text('m = {\n    "ok": set(),\n    "bad": [[2*2]],\n}\n')
    with pytest.raises(ValueError, match=r"\(line 3\): '2 \* 2' is not a literal"):
        read_dict("m.py")


# ==================================================
def test_all_files_are_read_first(in_tmp):
    (in_tmp / "good.py").write_text(f"good = {model()!r}\n")
    (in_tmp / "bad.py").write_text("bad = {'model': 2*2}\n")
    result = CliRunner().invoke(mp_create.cmd, ["good.py", "bad.py"])
    assert result.exit_code == 1
    assert "bad.py" in result.output
    assert not os.path.exists(in_tmp / "chain")

    result = CliRunner().invoke(mp_create.cmd, ["good.py", "missing.py"])
    assert result.exit_code == 1
    assert "missing.py" in result.output
    assert not os.path.exists(in_tmp / "chain")


# ==================================================
def test_python_m_error(in_tmp):
    (in_tmp / "bad.py").write_text(f"bad = {model(site={'A': ('[0,0,0]', 'pzz')})!r}\n")
    cmd = [sys.executable, "-m", "multipie.scripts.mp_create"]
    env = os.environ | {"PYTHONPATH": PACKAGE_PARENT}
    result = subprocess.run(cmd + ["bad.py"], capture_output=True, text=True, env=env)
    assert result.returncode == 1
    assert "Traceback" not in result.stderr
    assert "unknown orbital 'pzz'" in result.stderr

    result = subprocess.run(cmd + ["-v", "bad.py"], capture_output=True, text=True, env=env)
    assert result.returncode == 1
    assert "Traceback" in result.stderr
    assert "while creating model 'chain'." in result.stderr


# ==================================================
def test_samb_qtdraw_create_false(in_tmp, monkeypatch, capsys):
    create_model(model(), topdir="out")
    mm = MaterialModel("out")
    mm.load("chain")
    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: True)
    mm.save_samb_qtdraw()
    out = capsys.readouterr()
    assert "'qtdraw/create' of the model is False" in out.out + out.err
    assert not os.path.exists(in_tmp / "out" / "chain" / "samb")


# ==================================================
def test_multiple_files(in_tmp):
    (in_tmp / "a.py").write_text(f"a = {model(model='a')!r}\n")
    (in_tmp / "b.py").write_text(f"b1 = {model(model='b1')!r}\nb2 = {model(model='b2')!r}\n")
    result = CliRunner().invoke(mp_create.cmd, ["a", "b.py"])
    assert result.exit_code == 0, result.output
    for name in ["a", "b1", "b2"]:
        assert os.path.isfile(in_tmp / name / f"{name}.pkl")


# ==================================================
def test_analyze_without_model_name(in_tmp):
    with pytest.raises(ValueError, match="no model is specified"):
        analyze_model({"samb": {}, "wannier": {}})
