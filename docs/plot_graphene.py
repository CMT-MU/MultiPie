import numpy as np
from gcoreutils.nsarray import NSArray
from multipie.model.multipie_manager import MultiPieManager
from multipie.model.construct_model import construct_samb_matrix, convert_samb_to_matrix_set
from gcoreutils.plot_util import init_plot, plot_dispersion


# ==================================================
def plot_model(k_linear, E, k_name):
    # initialize plot.
    plt, figure = init_plot()
    grid = plt.GridSpec(1, 1, hspace=0.1, wspace=0.07, width_ratios=[1], height_ratios=[1])
    ax = figure.add_subplot(grid[0, 0])

    # plot dispersion relation.
    plot_dispersion(k_linear, E[:, :], k_name, title="dispersion relation for graphene", ax=ax)
    plt.show()


# ==================================================
def construct_model(matrix_dict, param, N1=50):
    # set k-grid for high-symmetry line.
    k_point = {name: eval(val) for name, val in matrix_dict["k_point"].items()}
    k_path = matrix_dict["k_path"]
    B = NSArray(matrix_dict["A"], fmt="value").inverse().T
    k_grid, k_linear, k_name = NSArray.grid_path(k_point, k_path, N1, B)

    # construct SAMB matrices as a function of z_j parameters for all k grid points.
    Z = construct_samb_matrix(matrix_dict, k_grid)

    # convert SAMB matrices to the linear combination with z_j parameters for all k grid points.
    M = convert_samb_to_matrix_set(Z, param)

    # return matrices, linear position of k, and names of k points.
    return M, k_linear, k_name


# ==================================================
def plot_graphene():
    # initialize manager.
    mpm = MultiPieManager(verbose=True)

    # parameters for graphene.
    param = [-0.163, -7.274, 0.880, -0.693, 0.0761, 0.202, -0.080]

    matrix_dict = mpm.read("graphene/graphene_matrix.py")
    M, k_linear, k_name = construct_model(matrix_dict, param)
    E, U = np.linalg.eigh(M)

    print(f"(#k points, #energy) = {E.shape}")
    plot_model(k_linear, E, k_name)


# ==================================================
plot_graphene()
