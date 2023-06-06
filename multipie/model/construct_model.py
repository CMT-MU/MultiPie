"""
construct model from matrx dict.
"""
import numpy as np
from gcoreutils.nsarray import NSArray


# ==================================================
def construct_samb_matrix(matrix_dict, k_grid=None):
    """
    construct SAMB matrix from matrix dict.

    Args:
        matrix_dict (dict): matrix dict from _matrix.py file.
        k_grid (NSArray, optional): k_grid points (reduced coordinate).

    Returns:
        ndarray: basis matrix, [#bases, #k, dim, dim].
    """
    z_samb = matrix_dict["matrix"]
    z_mat = {k: str(NSArray(v).evalf()) for k, v in z_samb.items()}

    if matrix_dict["molecule"] or k_grid is None:
        Nk = 1
        bond = [{}]
    else:
        Nk = len(k_grid)
        d = {}
        for name, ex in matrix_dict["bond"].items():
            ex = NSArray(ex, fmt="value")
            ex = np.broadcast_to(ex, (Nk, 3))
            k_grid = np.array(k_grid)
            p = 2 * np.pi * np.sum(ex * np.array(k_grid), axis=1)
            d[name.replace("bond_", "c")] = np.cos(p)
            d[name.replace("bond_", "s")] = np.sin(p)

        bond = [{name: lst[i] for name, lst in d.items()} for i in range(len(k_grid))]

    def multiple_replace(s, d):
        s_ = s
        [s_ := s_.replace(name, f"({val})") for name, val in d.items()]
        return s_

    vmultiple_replace = np.vectorize(multiple_replace)

    Z_list = [Zj.replace("I", "1j") for Zj in z_mat.values()]
    Z = np.array([np.array([np.array(eval(m), dtype="complex64") for m in vmultiple_replace(Zj, bond)]) for Zj in Z_list])

    return Z


# ==================================================
def convert_samb_to_matrix_set(Z, param=None):
    """
    convert SAMB to a set of matrix.

    Args:
        Z (ndarray): samb matrix, [#bases, #k, dim, dim].
        param (list, optional): parameter set(s), [z_j] for crystal or [z_j]/[[z_j]] for molecule.

    Returns:
        ndarray: matrix, [#params/#k, dim, dim].
    """
    if param is None:
        param = [1.0] * Z.shape[0]
    param = np.array(param)
    if param.ndim == 1:
        M = [Z.transpose(1, 2, 3, 0) @ param]
    else:
        M = [Z.transpose(1, 2, 3, 0) @ p for p in param]
    M = np.vstack(M)

    return M
