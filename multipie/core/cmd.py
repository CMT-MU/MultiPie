"""
Create command for model, matrix, and qtdraw.
"""

import os
from multipie.util.util import timer, read_dict_file
from multipie.core.material_model import MaterialModel


# ==================================================
def create_samb(models, topdir=None, verbose=False):
    """
    Create SAMB.

    Args:
        models (str or [str] or dict or [dict]): model file name(s) or dict(s).
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Note:
        - model name directory is created to store results in top directory.
        - if topdir is None, current directory is used.
    """

    def create(mm, model):
        @timer(f"create model='{model["model"]}'", verbose=verbose)
        def create0():
            mm.analyze(model)
            mm.save()
            mm.save_view()
            mm.save_pdf()

        create0()

    if topdir is not None:
        os.makedirs(topdir, exist_ok=True)

    models = list(read_dict_file(models, topdir, verbose).values())

    mm = MaterialModel(topdir, verbose)
    for model in models:
        create(mm, model)


# ==================================================
def create_samb_qtdraw(models, topdir=None, verbose=False):
    """
    Create SAMB QtDraw file.

    Args:
        models (str or [str] or dict or [dict]): model name(s) or dict(s).
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (bool) -- if error occurs, return True, otherwise False.

    Note:
        - samb directory is created to store results in model directory.
        - if topdir is None, current directory is used, which should contain model directory.
    """

    def create(mm):
        @timer(f"create SAMB QtDraw file for model='{mm["model"]}'", verbose=verbose)
        def create0():
            mm.save_samb_qtdraw()

        create0()

    if topdir is None:
        topdir = os.getcwd()

    if not isinstance(models, (list, tuple)):
        models = [models]

    mm = MaterialModel(topdir)
    for model in models:
        if type(model) == dict:
            model = model["model"]
        try:
            mm.load(model)
        except FileNotFoundError:
            return True
        create(mm)

    return False


# ==================================================
def create_samb_matrix(select_inputs, topdir=None, verbose=False):
    """
    Create SAMB matrix (hr).

    Args:
        select_inputs (str or [str] or dict or [dict]): select input file name(s) or dict(s).
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (bool) -- if error occurs, return True, otherwise False.

    Note:
        - matrix (and hr) file is created in model directory.
        - if topdir is None, current directory is used, which should contain model directory.
    """

    def create(mm, select, parameter):
        @timer(f"create matrix (hr) for model='{mm["model"]}'", verbose=verbose)
        def create0():
            mm.save_samb_matrix(select, parameter)

        create0()

    select_inputs = list(read_dict_file(select_inputs, topdir, verbose).values())

    mm = MaterialModel(topdir, verbose)
    for inp in select_inputs:
        model = inp["model"]
        select = inp["select"]
        parameter = inp.get("parameter", None)
        try:
            mm.load(model)
        except FileNotFoundError:
            return True
        create(mm, select, parameter)

    return False
