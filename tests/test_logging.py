"""
Tests for the logging setup of MultiPie.
"""

import contextlib
import io
import logging
import os

import pytest

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import setup_logging

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
@pytest.fixture(autouse=True)
def fresh_loggers():
    """
    Start without the handlers of the "multipie" logger, and restore the root and "multipie" loggers.
    """
    root, mp = logging.getLogger(), logging.getLogger("multipie")
    saved = [(lg, lg.handlers[:], lg.level, lg.propagate) for lg in (root, mp)]
    mp.handlers.clear()
    mp.setLevel(logging.NOTSET)
    mp.propagate = True
    yield
    for lg, handlers, level, propagate in saved:
        lg.handlers[:] = handlers
        lg.setLevel(level)
        lg.propagate = propagate


@pytest.fixture
def in_tmp(tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(cwd)


# ==================================================
def test_root_logger_unchanged(in_tmp):
    # create_model sets up the "multipie" logger only, not the root logger of the application.
    root = logging.getLogger()
    handler = logging.StreamHandler(io.StringIO())
    root.addHandler(handler)
    root.setLevel(logging.ERROR)
    before = (root.handlers[:], root.level)
    create_model(dict(MODEL), topdir="out", verbose=True)
    assert (root.handlers, root.level) == before
    assert handler.stream.getvalue() == ""  # not propagated to the root logger.


# ==================================================
def test_messages_to_current_stdout(in_tmp):
    # messages are written to sys.stdout at the time of writing, which may be replaced after the setup.
    setup_logging()
    for _ in range(2):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            create_model(dict(MODEL), topdir="out", verbose=True)
        assert "=== (create model='chain') end" in out.getvalue()


# ==================================================
def test_user_handler_kept(in_tmp, capsys):
    # a handler of the "multipie" logger set by the user is kept, and no handler is added.
    stream = io.StringIO()
    mp = logging.getLogger("multipie")
    mp.addHandler(logging.StreamHandler(stream))
    mp.setLevel(logging.INFO)
    create_model(dict(MODEL), topdir="out", verbose=True)
    assert len(mp.handlers) == 1
    assert "=== (create model='chain') end" in stream.getvalue()
    assert "=== (create model='chain')" not in capsys.readouterr().out


# ==================================================
def test_setup_logging_once():
    # the first setup is kept.
    setup_logging(logging.WARNING)
    setup_logging(logging.INFO)
    mp = logging.getLogger("multipie")
    assert len(mp.handlers) == 1
    assert mp.level == logging.WARNING


# ==================================================
@pytest.mark.parametrize("verbose", [False, True])
def test_dos_message_verbose(in_tmp, capsys, verbose):
    # the message is written only with verbose.
    create_model(dict(MODEL), topdir="out")
    capsys.readouterr()
    ModelAnalyzer("out", verbose=verbose).analyze(
        {"samb": {"model": "chain", "parameter": {"z1": 1.0}}, "grid": (4, 4, 4), "output": {"dos": True}}
    )
    assert ("compute and output dos." in capsys.readouterr().out) == verbose
