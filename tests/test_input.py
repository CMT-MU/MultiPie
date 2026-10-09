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
from multipie.util.util_pdf_latex import LaTeXError, PDFviaLaTeX, tex_text

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
def test_parameter_str_value(in_tmp):
    # values in samb/parameter given as str are evaluated by SymPy, as in a parameter file.
    create_model(dict(MODEL), topdir="out")
    ma = ModelAnalyzer("out")
    ma.analyze({"samb": {"model": "chain", "parameter": {"z1": 1.0, "z2": 2**0.5}}, "output": {"dispersion": {"k_path": None}}})
    expected = dict(ma.HR)

    ma.analyze(
        {"samb": {"model": "chain", "parameter": {"z1": "1", "z2": "sqrt(2)"}}, "output": {"dispersion": {"k_path": None}}}
    )
    assert ma.parameter == pytest.approx({"z1": 1.0, "z2": 2**0.5})
    assert ma.HR.keys() == expected.keys()
    assert all(complex(ma.HR[k]) == pytest.approx(complex(v)) for k, v in expected.items())


# ==================================================
@pytest.mark.parametrize(
    "control, error, match",
    [
        ({"parameter": "my_z"}, ValueError, r"must be '\.py' file"),
        ({"parameter": "info/my_z.py"}, FileNotFoundError, r"info/info/my_z\.py' is not found \(samb/parameter is relative to"),
        ({"parameter": os.path.abspath("my_z.py")}, ValueError, r"must be relative to '.*out/chain/info'"),
    ],
)
def test_parameter_file_error(in_tmp, control, error, match):
    create_model(dict(MODEL), topdir="out")
    os.makedirs(in_tmp / "out" / "chain" / "info")
    write_dict({"z1": 1.0}, "my_z.py", w_dir=str(in_tmp / "out" / "chain" / "info"))
    with pytest.raises(error, match=match):
        ModelAnalyzer("out").analyze({"samb": {"model": "chain"} | control, "output": {"dispersion": {"k_path": None}}})


# ==================================================
def test_pdf_failure(in_tmp, monkeypatch, capsys):
    # LaTeX is found, but ptex2pdf is not: the .tex and .pkl files are written, and only the PDF is skipped.
    monkeypatch.setattr("multipie.core.material_model.check_latex", lambda: True)
    which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda cmd: None if cmd == "ptex2pdf" else which(cmd))
    path = in_tmp / "out" / "chain"
    os.makedirs(path)
    (path / "chain.pdf").write_text("PDF of an earlier run")
    create_model(MODEL | {"pdf": {"create": True}}, topdir="out")  # warning via logging (set up by create_model).
    assert os.path.isfile(path / "chain.pkl")
    assert os.path.isfile(path / "chain.tex")
    assert not os.path.exists(path / "chain.pdf")  # the PDF of an earlier run is removed.
    out = capsys.readouterr()
    assert "skip PDF creation: ptex2pdf is not found" in out.out + out.err


# ==================================================
def test_tex_text():
    assert tex_text("Fe1") == "Fe1"
    assert tex_text(r"a_b#c%d&e$f{g}h") == r"a\texttt{\symbol{95}}b\#c\%d\&e\$f\{g\}h"
    assert tex_text("a\\b^c~d") == r"a\textbackslash{}b\^{}c\~{}d"
    assert tex_text("a\\b^c~d_e", math=True) == r"a\backslash{}b\wedge{}c\sim{}d\texttt{\symbol{95}}e"


# ==================================================
def test_pdf_names_escaped(in_tmp, monkeypatch):
    # model and site names with '_' are written with the underscore glyph in the LaTeX source.
    monkeypatch.setattr("multipie.core.material_model.check_latex", lambda: True)
    which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda cmd: None if cmd == "ptex2pdf" else which(cmd))
    model = MODEL | {
        "model": "my_chain",
        "site": {"Fe_1": ("[0,0,0]", "s"), "A": ("[1/2,0,0]", "s")},
        "bond": [("A", "A", [1])],
        "pdf": {"create": True},
    }
    create_model(model, topdir="out")
    tex = (in_tmp / "out" / "my_chain" / "my_chain.tex").read_text()
    u = r"\texttt{\symbol{95}}"
    assert r"\texttt{my" + u + "chain}" in tex
    assert r"\texttt{Fe" + u + "1}" in tex
    assert r"{\rm Fe" + u + "1}" in tex  # basis in the full matrix (math mode).
    assert "Fe_1" not in tex and "my_chain" not in tex


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


# ==================================================
@pytest.mark.parametrize("fail", [None, 1, 2, "interrupt"])
def test_pdf_build_cleanup(tmp_path, monkeypatch, fail):
    # the PDF of an earlier run is removed, and no PDF is left when a pass fails or is interrupted.
    (tmp_path / "x.pdf").write_text("earlier run")
    (tmp_path / "other.pdf").write_text("unrelated")
    n_pass = []

    def run_tex(cmd, timeout):
        n_pass.append(1)
        (tmp_path / "x.pdf").write_text(f"pass {len(n_pass)}")
        if fail == "interrupt" and len(n_pass) == 2:
            raise KeyboardInterrupt
        return 1 if fail == len(n_pass) else 0

    which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda cmd: "/bin/ptex2pdf" if cmd == "ptex2pdf" else which(cmd))
    monkeypatch.setattr("multipie.util.util_pdf_latex._run_tex", run_tex)
    monkeypatch.setattr(PDFviaLaTeX, "_check_package", lambda self: None)
    pdf = PDFviaLaTeX("x", dir=str(tmp_path))
    pdf.table([["a"]], ["1"], ["b"], long=True)  # a long table needs two passes.
    cwd = os.getcwd()
    if fail is None:
        pdf.build()
        assert (tmp_path / "x.pdf").read_text() == "pass 2"
    else:
        with pytest.raises(KeyboardInterrupt if fail == "interrupt" else Exception):
            pdf.build()
        assert not os.path.exists(tmp_path / "x.pdf")
    assert os.getcwd() == cwd
    assert (tmp_path / "other.pdf").read_text() == "unrelated"


# ==================================================
def test_pdf_compile_special_names(in_tmp, capsys):
    # with TeX installed, model and site names with LaTeX special characters are compiled.
    if shutil.which("latex") is None or shutil.which("ptex2pdf") is None:
        pytest.skip("LaTeX is not installed.")
    name = "Fe#1&%$^~{x}\\"
    model = MODEL | {
        "model": "my_chain",
        "site": {name: ("[0,0,0]", ["s", "p"]), "A": ("[1/2,0,0]", "s")},
        "bond": [(name, "A", [1]), ("A", "A", [1])],
        "pdf": {"create": True, "common_samb": True},
    }
    create_model(model, topdir="out")
    out = capsys.readouterr()
    if "not found" in out.out + out.err:
        pytest.skip("a LaTeX package or ptex2pdf is not found.")
    assert "skip PDF creation" not in out.out + out.err
    assert os.path.isfile(in_tmp / "out" / "my_chain" / "my_chain.pdf")


# ==================================================
def test_pdf_build_cleanup_error(tmp_path, monkeypatch):
    # an error in removing the PDF does not hide the compile error, and the working directory is restored.
    (tmp_path / "x.pdf").write_text("earlier run")
    which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda cmd: "/bin/ptex2pdf" if cmd == "ptex2pdf" else which(cmd))
    monkeypatch.setattr("multipie.util.util_pdf_latex._run_tex", lambda cmd, timeout: 1)
    monkeypatch.setattr(PDFviaLaTeX, "_check_package", lambda self: None)
    remove = os.remove

    def remove_no_pdf(path):
        if path.endswith(".pdf"):
            raise PermissionError(path)
        remove(path)

    monkeypatch.setattr(os, "remove", remove_no_pdf)
    cwd = os.getcwd()
    with pytest.raises(LaTeXError, match="compile error"):
        PDFviaLaTeX("x", dir=str(tmp_path)).build()
    assert os.getcwd() == cwd
