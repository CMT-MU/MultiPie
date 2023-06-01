import numpy as np
from gcoreutils.nsarray import NSArray
from multipie.model.multipie_manager import MultiPieManager
from multipie.model.construct_model import construct_samb_matrix, convert_samb_to_matrix_set
from gcoreutils.plot_util import init_plot, plot_dispersion

data_dir = __file__[: __file__.rfind("/")] + "/material_model"


# ==================================================
def plot_model(x, E, molecule, name, k_name=None):
    plt, figure = init_plot()
    grid = plt.GridSpec(1, 1, hspace=0.1, wspace=0.07, width_ratios=[1], height_ratios=[1])
    ax = figure.add_subplot(grid[0, 0])

    if molecule:
        plot_dispersion(x, E[:, :], title=f"Energy for '{name}'", ax=ax, xlabel=r"$x$", ylabel=r"$E_n$")
    else:
        plot_dispersion(x, E[:, :], k_name, title=f"dispersion relation for '{name}'", ax=ax)
    plt.show()


# ==================================================
def construct_model(matrix_dict, param, N1=50):
    molecule = matrix_dict["molecule"]

    if molecule:
        Z = construct_samb_matrix(matrix_dict)
        M = convert_samb_to_matrix_set(Z, param)
        return M

    else:
        k_point = {name: eval(val) for name, val in matrix_dict["k_point"].items()}
        k_path = matrix_dict["k_path"]
        B = NSArray(matrix_dict["A"], fmt="value").inverse().T
        k_grid, k_linear, k_name = NSArray.grid_path(k_point, k_path, N1, B)
        Z = construct_samb_matrix(matrix_dict, k_grid)
        M = convert_samb_to_matrix_set(Z, param)
        return M, k_linear, k_name


# ================================================== test
def test_construct_model():
    mpm = MultiPieManager(data_dir, verbose=True)
    models = ["graphene", "Te", "SrVO3"]
    params = [[0.0, 1.0, 0.0], None, [22.650918659029713, -0.85036445, 0.5324245, -0.15995012, 0.2726717, -0.04474434]]

    # crystals.
    for model, param in zip(models, params):
        matrix_dict = mpm.read(model + "/" + model + "_matrix.py")
        name = matrix_dict["model"]
        print(f"=== {name} ===")
        M, k_linear, k_name = construct_model(matrix_dict, param)
        E, U = np.linalg.eigh(M)
        print(f"shape: {E.shape}")
        plot_model(k_linear, E, False, name, k_name)

    # molecule.
    x = np.linspace(0, 3, 100)
    param = [[0, np.sqrt(6), np.sqrt(2), xi, xi, 0, 0, 0, 0, 0, 0] for xi in x]
    matrix_dict = mpm.read("CH4/CH4_matrix.py")
    name = matrix_dict["model"]
    print(f"=== {name} ===")
    M = construct_model(matrix_dict, param)
    E, U = np.linalg.eigh(M)
    print(f"shape: {E.shape}")
    plot_model(x, E, True, name)


test_construct_model()
