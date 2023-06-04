import math
import numpy as np
from multipie.model.multipie_manager import MultiPieManager
from multipie.model.construct_model import construct_samb_matrix, convert_samb_to_matrix_set
from gcoreutils.plot_util import init_plot, plot_dispersion


# ==================================================
def plot_model(x, E, title, xlabel, ylabel):
    # initialize plot.
    plt, figure = init_plot()
    grid = plt.GridSpec(1, 1, hspace=0.1, wspace=0.07, width_ratios=[1], height_ratios=[1])
    ax = figure.add_subplot(grid[0, 0])

    # plot parameter dependence.
    rm = max(abs(math.floor(E.min())), abs(math.ceil(E.max())))
    E_range = [-rm, rm]
    plot_dispersion(x, E[:, :], title=title, ax=ax, xlabel=xlabel, ylabel=ylabel, E_range=E_range)
    plt.show()


# ==================================================
def construct_model(matrix_dict, param):
    # construct SAMB matrices as a function of z_j parameters.
    Z = construct_samb_matrix(matrix_dict)

    # convert SAMB matrices to the linear combination with z_j parameters for all sets.
    M = convert_samb_to_matrix_set(Z, param)

    # return matrices.
    return M


# ==================================================
def plot():
    # model name.
    model = "C3v"

    # parameters for model.
    x = np.linspace(0, 1, 100)
    param = [[0.0, 1.0, 0.0, 0.0, 0.0, 0.2, 0.0, xi, 0.0] for xi in x]

    # initialize manager.
    mpm = MultiPieManager(verbose=True)

    # construct matrices, and diagonalize them.
    matrix_dict = mpm.read(f"{model}/{model}_matrix.py")
    M = construct_model(matrix_dict, param)
    E, U = np.linalg.eigh(M)

    # plot results.
    mpm.log(f"  parameters (zj; z8=0) = {param[0]}", None)
    mpm.log(f"  (#set, #energies) =  {E.shape}", None)
    plot_model(x, E, f"Energy for {model} molecule", "$z_8$", "$E_n$")


# ==================================================
plot()
