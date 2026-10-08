"""
Regression tests for reading input files, the working directory, and missing files.
"""

import os
import sys
import types

import pytest

from multipie import MaterialModel
from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict, read_dict_file, write_dict
from multipie.util.util_binary import BinaryManager

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
    """
    Run in tmp_path, and check that the working directory is not changed by the test.
    """
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
        assert os.getcwd() == str(tmp_path)
    finally:
        os.chdir(cwd)


# ==================================================
def test_read_dict_with_equal_in_comment(in_tmp):
    (in_tmp / "ctrl.py").write_text('"""\ncontrol, a = b.\n"""\nctrl = {"mode": "samb"}  # s=0\n')
    assert read_dict("ctrl.py") == {"mode": "samb"}
    assert read_dict_file("ctrl.py") == {"ctrl": {"mode": "samb"}}


# ==================================================
def test_read_dict_requires_one_dict(in_tmp):
    (in_tmp / "two.py").write_text("a = {'x': 1}\nb = {'y': 2}\n")
    with pytest.raises(ValueError, match="exactly one dict"):
        read_dict("two.py")
    assert read_dict_file("two.py") == {"a": {"x": 1}, "b": {"y": 2}}


# ==================================================
def test_read_dict_file_does_not_change_cwd(in_tmp):
    sub = in_tmp / "sub"
    sub.mkdir()
    (sub / "m.py").write_text("m = {'model': 'x'}\n")

    assert read_dict_file({"a": 1}, topdir="sub") == {"dict1": {"a": 1}}
    assert read_dict_file("m.py", topdir="sub") == {"m": {"model": "x"}}
    with pytest.raises(FileNotFoundError):
        read_dict_file("missing.py", topdir="sub")


# ==================================================
def test_create_model_relative_topdir(in_tmp):
    create_model(dict(MODEL), topdir="out")
    assert os.path.isfile(in_tmp / "out" / "chain" / "chain.pkl")
    assert not os.path.exists(in_tmp / "out" / "out")


# ==================================================
def test_model_analyzer_relative_topdir(in_tmp):
    create_model(dict(MODEL), topdir="out")
    (in_tmp / "out" / "ctrl.py").write_text(
        "ctrl = {'samb': {'model': 'chain', 'parameter': {'z1': 1.0}}, 'output': {'dispersion': {'k_path': None}}}\n"
    )

    ma = ModelAnalyzer("out")
    ma.analyze("ctrl.py")  # relative to topdir.
    assert os.path.isfile(in_tmp / "out" / "chain" / "info" / "chain_hr.dat")


# ==================================================
def test_missing_binary(in_tmp):
    with pytest.raises(FileNotFoundError):
        BinaryManager("missing", topdir=str(in_tmp))
    with pytest.raises(FileNotFoundError):
        MaterialModel(str(in_tmp)).load("missing")


# ==================================================
def test_save_samb_qtdraw_restores_cwd(in_tmp, monkeypatch):
    create_model(dict(MODEL), topdir="out")
    mm = MaterialModel("out")
    mm.load("chain")

    def fail(*args, **kwargs):
        raise RuntimeError("failure while writing")

    mm["qtdraw_prop"]["create"] = True
    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: True)
    monkeypatch.setattr(mm, "save_atomic_samb", fail)
    with pytest.raises(RuntimeError, match="failure while writing"):
        mm.save_samb_qtdraw()


# ==================================================
def test_save_view_restores_cwd(in_tmp, monkeypatch):
    create_model(dict(MODEL), topdir="out")
    mm = MaterialModel("out")
    mm.load("chain")
    mm["qtdraw_prop"]["create"] = True

    def fail(*args, **kwargs):
        raise RuntimeError("failure while writing")

    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: True)
    monkeypatch.setitem(sys.modules, "qtdraw", types.SimpleNamespace(create_qtdraw_file=fail))
    with pytest.raises(RuntimeError, match="failure while writing"):
        mm.save_view()


# ==================================================
def test_read_dict_rejects_other_statements(in_tmp):
    (in_tmp / "z.py").write_text("z = {'z1': 1.0}\nz.update({'z1': 2.0})\n")
    with pytest.raises(ValueError, match="only a dict"):
        read_dict("z.py")


# ==================================================
def test_analyze_graphene_from_other_cwd(in_tmp):
    # create and analyze in "out" while the current directory is elsewhere, then read the generated z file back.
    example = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples")
    model = read_dict(os.path.abspath(os.path.join(example, "graphene_in.py")))
    model["pdf"] = {"create": False}
    model["qtdraw"] = {"create": False}
    control = read_dict(os.path.abspath(os.path.join(example, "graphene_ctrl.py")))
    control["samb"]["samb_figure"] = False
    control["grid"] = (6, 6, 1)

    create_model(model, topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze(control)
    disp = ma["output"]["dispersion"]
    out = in_tmp / "out" / "graphene"
    assert os.path.isfile(out / "info" / "graphene_z.py")
    assert os.path.isfile(out / "output" / "graphene_dispersion.txt")

    control["samb"]["parameter"] = ""  # read info/graphene_z.py.
    ma = ModelAnalyzer("out")
    ma.analyze(control)
    assert ma["output"]["dispersion"]["e_max"] == pytest.approx(disp["e_max"])
    assert ma["output"]["dispersion"]["e_min"] == pytest.approx(disp["e_min"])


# ==================================================
@pytest.mark.parametrize("stem", ["plain_z", "my-model_z", "1abc_z"])
def test_write_dict_read_back(in_tmp, stem):
    # write_dict uses the file name as the variable name, which need not be an identifier.
    dic = {"z1": 1.0, "s": "a = b"}
    write_dict(dic, f"{stem}.py", comment="comment, x = y.\n", w_dir=str(in_tmp))
    assert read_dict(f"{stem}.py") == dic
    assert read_dict_file(f"{stem}.py") == {stem: dic}


# ==================================================
def test_save_view_output(in_tmp, monkeypatch):
    create_model(dict(MODEL), topdir="out")
    mm = MaterialModel("out")
    mm.load("chain")
    mm["qtdraw_prop"]["create"] = True

    def write(filename, callback):
        with open(filename, "w") as f:
            f.write("qtdraw")

    monkeypatch.setattr("multipie.core.material_model.check_qtdraw", lambda: True)
    monkeypatch.setitem(sys.modules, "qtdraw", types.SimpleNamespace(create_qtdraw_file=write))
    mm.save_view()
    assert os.path.isfile(in_tmp / "out" / "chain" / "chain.qtdw")


# ==================================================
def test_dict_file_non_identifier_edge_cases(in_tmp):
    # multiline strings are not changed.
    (in_tmp / "a.py").write_text('my-model_z = {"s": """x\nmy-model = {1: 2}\n"""}\n')
    assert read_dict("a.py") == {"s": "x\nmy-model = {1: 2}\n"}

    # an existing name like the temporary one does not collide.
    (in_tmp / "b.py").write_text('_multipie_var0 = {"a": 1}\nmy-model = {"b": 2}\n')
    assert read_dict_file("b.py") == {"_multipie_var0": {"a": 1}, "my-model": {"b": 2}}

    # CRLF line endings.
    (in_tmp / "c.py").write_bytes(b'"""\r\ncomment\r\n"""\r\nmy-model_z = {\r\n    "z1": 1.0,\r\n}\r\n')
    assert read_dict("c.py") == {"z1": 1.0}

    # valid targets are not renamed, even if another target is not an identifier.
    (in_tmp / "e.py").write_text('my-model = {"a": 1}\nobj.attr = {"b": 2}\n')
    assert read_dict_file("e.py") == {"my-model": {"a": 1}}

    # line separators other than newline (U+2028) in a string.
    (in_tmp / "f.py").write_text('"""doc\u2028more"""\nmy-model_z = {"z1": 1.0}\n', encoding="utf-8")
    assert read_dict("f.py") == {"z1": 1.0}

    # valid python statements keep their meaning: rejected by read_dict, ignored by read_dict_file.
    for src in ['x = {"a": 1}\nx += {"b": 2}\n', 'x = {"a": 1}\nobj.attr = {"b": 2}\n']:
        (in_tmp / "d.py").write_text(src)
        with pytest.raises(ValueError, match="only a dict"):
            read_dict("d.py")
        assert read_dict_file("d.py") == {"x": {"a": 1}}
