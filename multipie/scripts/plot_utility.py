"""
utility for model construction and plot.
"""

import math
import os
import numpy as np
from gcoreutils.nsarray import NSArray
from gcoreutils.plot_util import init_plot, plot_dispersion
from multipie.model.construct_model import construct_samb_matrix, convert_samb_to_matrix_set
from multipie.model.multipie_manager import MultiPieManager


# ==================================================
def _solve_model(matrix_dict, param, N1=50):
    """
    solve model.

    Args:
        matrix_dict (dict): matrix dict.
        param (dict): z_j parameter dict (crystal) or z_j parameter list dict + x for molecule.
        N1 (int, optional): number of division in each k-section for crystal.

    Returns: tuple
        - NSArray: eigen energy, [k_index/param_index, basis_index].
        - NSArray: eigenvectors, [k_index/param_index, basis_index, eigenstate_index].
        - NSArray: linear position of k point (crystal).
        - dict: name of k point, {disconnected linear position: label} (crystal).
        - NSArray: linear k grid (crystal).
    """
    if matrix_dict["molecule"]:
        # construct SAMB matrices as a function of z_j parameters.
        Z = construct_samb_matrix(matrix_dict)

        # convert SAMB matrices to the linear combination with z_j parameters for all sets.
        pm = [i for k, i in param.items() if k[0] == "z"]
        pm = np.array(pm).T.tolist()
        M = convert_samb_to_matrix_set(Z, pm)
    else:
        # set k-grid for high-symmetry line.
        k_point = {name: eval(val) for name, val in matrix_dict["k_point"].items()}
        k_path = matrix_dict["k_path"]
        B = NSArray(matrix_dict["A"], fmt="value").inverse().T
        k_grid, k_linear, k_name = NSArray.grid_path(k_point, k_path, N1, B)

        # construct SAMB matrices as a function of z_j parameters for all k grid points.
        Z = construct_samb_matrix(matrix_dict, k_grid)

        # convert SAMB matrices to the linear combination with z_j parameters for all k grid points.
        pm = list(param.values())
        M = convert_samb_to_matrix_set(Z, pm)

    # solve eigensystem.
    E, U = np.linalg.eigh(M)

    if matrix_dict["molecule"]:
        return E, U
    else:
        return E, U, k_linear, k_name, k_grid


# ==================================================
def _plot_model(k_linear, E, k_name=None, title="", xlabel=None, model=None, verbose=False):
    """
    plot a model.

    Args:
        k_linear (NSArray): linear position of k point.
        E (ndarray): eigen energy, [k_index, basis_index].
        k_name (dict, optional): name of k point, {disconnected linear position: label}.
        title (str, optional): title of plot.
        xlabel (str, optional): x label.
        model (str, optional): model name.
        verbose (bool, optional): verbose ?
    """
    # initialize plot.
    plt, figure = init_plot()
    grid = plt.GridSpec(1, 1, hspace=0.1, wspace=0.07, width_ratios=[1], height_ratios=[1])
    ax = figure.add_subplot(grid[0, 0])

    # plot energy.
    rm = max(abs(math.floor(E.min())), abs(math.ceil(E.max())))
    E_range = [-rm, rm]
    plot_dispersion(k_linear, E, k_name, title=title, ax=ax, E_range=E_range, xlabel=xlabel)
    # plt.show()
    fname = model
    fname = os.getcwd() + "/" + fname + ".png"
    plt.savefig(fname, transparent=True)
    if verbose:
        print("  * create '" + fname + "'.")


# ==================================================
def _plot_single_model(topdir, param_file, verbose):
    """
    plot a model with parameters.

    Args:
        topdir (str): top directory.
        param_file (str): parameter file.
        verbose (bool): verbose parallel info.
        no (int, optional): plot no.
    """
    mpm = MultiPieManager(topdir=topdir + "/", verbose=verbose)
    pfile = mpm.read(param_file)
    model = pfile["model"]

    # mpm = MultiPieManager(topdir=topdir + "/" , verbose=verbose)
    d = mpm.read(model + "/" + model + "_matrix.py")

    if "k_point" in pfile.keys():
        d["k_point"] = pfile["k_point"]
    if "k_path" in pfile.keys():
        d["k_path"] = pfile["k_path"]

    n_param = len(d["matrix"])
    param = {f"z_{i+1:03d}": 0.0 for i in range(n_param)}
    for zno, val in pfile["param"].items():
        param[f"z_{zno:03d}"] = val

    if d["molecule"]:
        # construct matrices, and diagonalize them.
        E, _ = _solve_model(d, param)
        # plot energy.
        _plot_model(pfile["x"], E, title=f"energy for {model} molecule", xlabel=pfile["x_name"], model=model, verbose=verbose)
    else:
        # construct matrices, and diagonalize them.
        E, _, k_linear, k_name, _ = _solve_model(d, param)
        # plot dispersion.
        _plot_model(k_linear, E, k_name=k_name, title=f"energy dispersion for {model}", model=model, verbose=verbose)


# ==================================================
def plot_model(model_param_list, topdir=".", verbose=False):
    """
    plot models.

    Args:
        model_param_list (list): list of model name and parameter dict.
        topdir (str, optional): top data directory.
        verbose (bool, optional): verbose parallel info.

    Notes:
        - topdir must be the parent directory of model data.
    """
    if type(model_param_list) == list:
        for model_param in model_param_list:
            _plot_single_model(topdir, model_param, verbose)
    else:
        _plot_single_model(topdir, model_param_list, verbose)
