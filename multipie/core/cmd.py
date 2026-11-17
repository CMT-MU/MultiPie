# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create command for model and analyze
"""

import os
import copy
import logging
from multipie.util.util import timer, read_dict_file, read_dict, write_dict, str_to_sympy, setup_logging
from multipie.core.material_model import MaterialModel
from multipie.core.model_analyzer import ModelAnalyzer

logger = logging.getLogger(__name__)


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
        @timer(f"create model='{model['model']}'", verbose=verbose)
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
        except Exception as e:
            e.add_note(f"while creating model '{model.get('model')}'.")
            raise

    return False


# ==================================================
def analyze_model(controls, topdir=None, verbose=False):
    """
    Analyze model.

    Args:
        controls (str or [str] or dict or [dict]): control file(s) or dict(s).
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (bool) -- if no error occurs, return False.

    Note:
        - if topdir is None, current directory is used, which should contain model directory.
    """
    setup_logging()

    def create(ma, control, name):
        @timer(f"analyze model by '{name}'", verbose=verbose)
        def create0():
            ma.analyze(control)

        create0()

    if topdir is None:
        topdir = os.getcwd()

    controls = list(read_dict_file(controls, topdir, verbose).values())

    ma = ModelAnalyzer(topdir, verbose=verbose)
    for control in controls:
        name = None
        try:
            name = control.get("samb", {}).get("model")
            if name is None:
                name = control.get("wannier", {}).get("seedname")
                if name is None:
                    raise ValueError("no model is specified in 'samb/model' or 'wannier/seedname'.")
            create(ma, control, name)
        except Exception as e:
            if name is not None:
                e.add_note(f"while analyzing model '{name}'.")
            raise

    return False


# ==================================================
def convert_model(model, controls=None, parameter=None, topdir=None, pdf=None, qtdraw=None, verbose=False):
    """
    Convert a model (.pkl) of MultiPie Ver.2 into Ver.3, and the weights z_j of its control files.

    Args:
        model (str): model name, i.e., "topdir/model/model.pkl".
        controls (str or [str], optional): control file(s) of mp_analyze whose "samb/parameter" are converted.
        parameter (dict, optional): weights of the SAMB of Ver.2 to convert, dict[zj, float or str], e.g., a fit.
        topdir (str, optional): top directory. [default: cwd]
        pdf (bool, optional): create PDF ? [default: as in the model of Ver.2]
        qtdraw (bool, optional): create QtDraw file ? [default: as in the model of Ver.2]
        verbose (bool, optional): verbose ?

    Returns:
        - (dict) -- report of the conversion.

    Note:
        - the model is created again with Ver.3 from the input rebuilt from the model of Ver.2, which is written as
          "topdir/model_in.py". The model of Ver.2 is kept as "topdir/model/model_v2.pkl".
        - the f orbitals are renamed (f3x, f3y, f2, f1 -> fx, fy, fc, fs), and max_neighbor and search_cell_range are
          dropped. The neighbors are chosen by the distances of the bonds of Ver.2 (Ver.2 may have missed a part of
          the neighbor shells).
        - the conversion is refused if the kets differ or the SAMB of Ver.2 are not in the space of those of Ver.3.
        - the weights z_j of Ver.3 are the projections of H = sum_j z_j Z_j of Ver.2 onto the SAMB of Ver.3,
          z_j(Ver.3) = Tr[Z_j(Ver.3) H]. A converted control file is written as "control_v3.py" next to the original.
        - the weights given by parameter are converted into report["parameter"], dict[zj, float].
        - for a model already converted (model.pkl of Ver.3 and model_v2.pkl), only the weights are converted.
    """
    from multipie.util.util_convert import model_version, major_version, rebuild_model_input, remap_neighbor, ModelConverter

    setup_logging()
    if topdir is None:
        topdir = os.getcwd()
    if isinstance(controls, str):
        controls = [controls]
    controls = [] if controls is None else list(controls)
    controls = [i if i.endswith(".py") else i + ".py" for i in controls]

    name = model[:-4] if model.endswith(".pkl") else model
    name = os.path.basename(os.path.normpath(name))
    model_dir = os.path.join(topdir, name)
    pkl = os.path.join(model_dir, name + ".pkl")
    pkl2 = os.path.join(model_dir, name + "_v2.pkl")
    input_file = os.path.join(topdir, name + "_in.py")
    if not os.path.isfile(pkl):
        raise FileNotFoundError(f"model file '{pkl}' is not found.")

    # read the control files first, so that nothing is written when one of them cannot be read.
    ctrl_data = {}
    for c in controls:
        ctrl = read_dict(c, topdir)
        ctrl_name = ctrl.get("samb", {}).get("model")
        if ctrl_name != name:
            raise ValueError(f"'{c}': samb/model is {ctrl_name!r}, not '{name}'.")
        out = os.path.join(topdir, c[:-3] + "_v3.py")
        if os.path.exists(out):
            raise FileExistsError(f"'{out}' exists. Remove it to convert '{c}'.")
        ctrl_data[c] = (ctrl, out)

    report = {"model": name, "input": None, "model_v2": pkl2, "renamed_orbital": {}, "removed_key": {}, "neighbor": {}}

    version, data = model_version(pkl)
    if major_version(version) == 2:
        if os.path.exists(pkl2):
            raise FileExistsError(f"'{pkl2}' exists. Remove it to convert '{pkl}'.")
        if os.path.exists(input_file):
            raise FileExistsError(f"'{input_file}' exists. Remove or rename it to convert '{pkl}'.")
        data2 = data
        model_in, note = rebuild_model_input(data2)
        report["renamed_orbital"] = note["renamed_orbital"]
        report["removed_key"] = note["removed_key"]
        report["neighbor"] = remap_neighbor(model_in, data2)
        if pdf is not None:
            model_in.setdefault("pdf", {})["create"] = pdf
        if qtdraw is not None:
            model_in.setdefault("qtdraw", {})["create"] = qtdraw

        comment = f"model input of '{name}' rebuilt from the model of MultiPie {version} by convert_model (mp_convert).\n"
        write_dict(model_in, os.path.basename(input_file), var=name, comment=comment, w_dir=topdir)
        report["input"] = input_file

        os.replace(pkl, pkl2)
        try:
            create_model(model_in, topdir=topdir, verbose=verbose)
            mm = MaterialModel(topdir)
            mm.load(name)
            converter = ModelConverter(data2, mm)
        except Exception as e:
            os.replace(pkl2, pkl)  # restore the model of Ver.2.
            e.add_note(f"while converting model '{name}' (the model of Ver.2 is restored as '{pkl}').")
            raise
    elif major_version(version) >= 3:
        if not os.path.isfile(pkl2):
            raise ValueError(f"'{pkl}' was created by MultiPie {version}; no model of Ver.2 ('{pkl2}') to convert from.")
        version2, data2 = model_version(pkl2)
        if major_version(version2) != 2:
            raise ValueError(f"'{pkl2}' was created by MultiPie {version2}, not Ver.2.")
        mm = MaterialModel(topdir)
        mm.load(name)
        converter = ModelConverter(data2, mm)
    else:
        raise ValueError(f"'{pkl}' was created by MultiPie {version}, which cannot be converted (only Ver.2).")

    report["SAMB"] = (len(converter.name2), len(converter.name3))
    report["span_residual"] = converter.span_residual

    def to_float(dic):
        return {z: float(str_to_sympy(v, rational=False)) if isinstance(v, str) else v for z, v in dic.items()}

    if parameter is not None:
        report["parameter"], report["parameter_residual"] = converter.convert(to_float(parameter))

    # control files.
    report["control"] = {}
    for c, (ctrl, out) in ctrl_data.items():
        z2 = ctrl["samb"].get("parameter", {})
        if isinstance(z2, str):  # z file relative to topdir/model/info.
            z2 = read_dict(z2, os.path.join(model_dir, "info"))
        z2 = to_float(z2)
        new, residual = converter.convert(z2)
        ctrl3 = copy.deepcopy(ctrl)
        ctrl3["samb"]["parameter"] = new
        comment = f"control converted from '{c}' (SAMB of MultiPie {version if major_version(version) == 2 else version2}) "
        comment += "by convert_model (mp_convert).\n"
        if report["neighbor"] and ctrl["samb"].get("select", {}).get("bond"):
            comment += "samb/select/bond is not converted: the neighbors of Ver.3 differ from those of Ver.2.\n"
            logger.warning(f"'{c}': samb/select/bond is not converted, while the neighbors changed: {report['neighbor']}.")
        write_dict(ctrl3, os.path.basename(out), comment=comment, w_dir=os.path.dirname(out))
        report["control"][c] = {"output": out, "parameter": (len(z2), len(new)), "residual": residual}

    return report
