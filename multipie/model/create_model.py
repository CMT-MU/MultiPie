"""
create model view, info, and basis.
"""
import pickle

from multipie.model.material_model import MaterialModel
from multipie.model.multipie_manager import MultiPieManager
from multipie.model.symmetry_adapted_model import SymmetryAdaptedModel
from multipie.model.util.create_pdf import ModelPDF


# ==================================================
def _create_single_model(model_dict, mpm, view_mode):
    """
    create a model view, info, and basis.

    Args:
        model_dict (dict or str): minimal model information.
        mpm (MultiPieManager): multipie manager.
        view_mode (str): mode for QtDraw.
    """
    model_dict = mpm.read(model_dict)

    model_dict = MaterialModel.regularize(model_dict)
    model_name = model_dict["model"]
    mpm.set_group(model_dict["group"])
    mpm.create_subdir(model_dict["option"]["output"])

    mpm.log(f"=== creating files for '{model_name}' in {mpm.dirname} ...", None)

    # create model.
    mpm.log("creating model ... ", None)
    cmodel = MaterialModel(model_dict, mpm)
    mpm.write(model_name + "_model.py", cmodel, cmodel._header(), model_name)

    # create view.
    if mpm.qtdraw:
        mpm.log("creating view ... ", None)
        view = cmodel.create_view(mode=view_mode)
        mpm.write(model_name + "_view.qtdw", view, "\nQtDraw data file in Python dict format.\n")

    # create SAMB.
    mpm.log("creating SAMB ... ", None)
    if cmodel["info"]["generate"]["time_reversal_type"] == "electric":
        head = ["Q", "G"]
    elif cmodel["info"]["generate"]["time_reversal_type"] == "magnetic":
        head = ["T", "M"]
    else:
        head = ["Q", "G", "T", "M"]
    samb = SymmetryAdaptedModel(cmodel, mpm=mpm, head=head)
    samb_dict = samb.create_dict()
    mpm.write(model_name + "_samb.py", samb_dict, SymmetryAdaptedModel._header(), model_name)

    # create matrix (real space)
    mpm.log("creating SAMB matrix (real space) ... ", None)
    samb_matrix_real = samb.create_matrix(fmt="sympy")

    if model_dict["option"]["binary_output"]:
        filename = model_name + "_matrix.pkl"
        full = mpm.filename(filename)
        f = open(full, "wb")
        pickle.dump(samb_matrix_real, f)
        f.close()
        mpm.log(f"  * wrote '{filename}'.", None)
    else:
        mpm.write(model_name + "_matrix.py", samb_matrix_real, SymmetryAdaptedModel._matrix_header(), model_name)

    # create matrix (momentum space)
    if not cmodel["info"]["molecule"] and cmodel["info"]["generate"]["fourier_transform"]:
        mpm.log("creating SAMB matrix (momentum space) ... ", None)
        samb_matrix = samb.create_matrix_k(full=True, fmt="sympy")

        if model_dict["option"]["binary_output"]:
            filename = model_name + "_matrix_k.pkl"
            full = mpm.filename(filename)
            f = open(full, "wb")
            pickle.dump(samb_matrix, f)
            f.close()
            mpm.log(f"  * wrote '{filename}'.", None)
        else:
            mpm.write(model_name + "_matrix_k.py", samb_matrix, SymmetryAdaptedModel._matrix_header_k(), model_name)

    # create LaTeX and PDF.
    if mpm.pdf:
        mpm.log("creating LaTeX and PDF ... ", None)
        ModelPDF(cmodel, samb_dict, mpm)

    mpm.formatter()
    mpm.log("=== total elapsed_time:", "start")


# ==================================================
def create_model(
    model_list,
    topdir=None,
    symbolic=True,
    parallel=True,
    verbose=False,
    pdf=True,
    formatter=True,
    view_mode=None,
    qtdraw=True,
):
    """
    create model view, info, and basis.

    Args:
        model_list (dict or list): minimal model information.
        topdir (str, optional): top directory for output.
        symbolic (bool, optional): output in sympy format ?
        parallel (bool, optional): use parallel code.
        verbose (bool, optional): verbose parallel info.
        pdf (bool, optional): create pdf/tex file.
        formatter (bool, optional): format by using black.
        view_mode (str, optional): mode for QtDraw, standard/arrow/debug. if None, option in the model is used.
        qtdraw (bool, optional): use QtDraw ?
    """
    mpm = MultiPieManager(topdir, verbose, symbolic, formatter, pdf, parallel, qtdraw)

    if type(model_list) == list:
        for model in model_list:
            _create_single_model(model, mpm, view_mode)
    else:
        _create_single_model(model_list, mpm, view_mode)
