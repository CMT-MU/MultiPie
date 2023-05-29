"""
create model view.
"""
from multipie.model.material_model import MaterialModel
from multipie.model.multipie_manager import MultiPieManager


# ==================================================
def _create_single_view(model_dict, mpm, view_mode):
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

    mpm.log(f"=== creating files for '{model_name}' in {mpm.dirname} ...", None)

    # create model.
    mpm.log("creating model ... ", None)
    cmodel = MaterialModel(model_dict, mpm)

    # create view.
    if mpm.qtdraw:
        mpm.log("creating view ... ", None)
        view = cmodel.create_view(mode=view_mode)
        mpm.write(model_name + "_view.qtdw", view, "\nQtDraw data file in Python dict format.\n")

    mpm.formatter()
    mpm.log("=== total elapsed_time:", "start")


# ==================================================
def create_view(model_list, verbose=False, formatter=True, view_mode=None):
    """
    create model view.

    Args:
        model_list (dict or list): minimal model information.
        verbose (bool, optional): verbose parallel info.
        formatter (bool, optional): format by using black.
        view_mode (str, optional): mode for QtDraw, standard/arrow/debug. if None, option in the model is used.
    """
    mpm = MultiPieManager(verbose=verbose, formatter=formatter)

    if type(model_list) == list:
        for model in model_list:
            _create_single_view(model, mpm, view_mode)
    else:
        _create_single_view(model_list, mpm, view_mode)
