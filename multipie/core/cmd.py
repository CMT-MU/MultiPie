"""
Create command for model and analyze
"""

import os
import logging
from multipie.util.util import timer, read_dict_file, setup_logging
from multipie.core.material_model import MaterialModel
from multipie.core.model_analyzer import ModelAnalyzer


# ==================================================
def create_model(models, topdir=None, verbose=False):
    """
    Create model.

    Args:
        models (str or [str] or dict or [dict]): model file name(s) or dict(s).
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (bool) -- if no error occurs, return False.

    Note:
        - model name directory is created to store results in top directory.
        - if topdir is None, current directory is used.
    """
    setup_logging()

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
        try:
            create(mm, model)
        except Exception:
            logging.exception("in create_model")
            raise

    return False


# ==================================================
def analyze_model(control, N1=50, N2=50, N3=50, topdir=None, verbose=False):
    """
    Analyze model.

    Args:
        control (str or dict): control file or dict.
        N1 (int, optional): number of divisions in a1.
        N2 (int, optional): number of divisions in a2.
        N3 (int, optional): number of divisions in a3.
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (bool) -- if no error occurs, return False.

    Note:
        - if topdir is None, current directory is used, which should contain model directory.
    """
    setup_logging()

    def create(ma):
        @timer(f"analyze model by '{control}'", verbose=verbose)
        def create0():
            ma.exec(control)

        create0()

    if topdir is None:
        topdir = os.getcwd()

    ma = ModelAnalyzer(N1, N2, N3, topdir, verbose=verbose)
    try:
        create(ma)
    except Exception:
        logging.exception("in analyze_model")
        raise

    return False
