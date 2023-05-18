from gcoreutils.convert_util import text_to_sympy
from gcoreutils.nsarray import NSArray


def convert_to_matrix(n_bra, n_ket, lst):
    mat = NSArray.zeros((n_bra, n_ket), "matrix")
    for i, j, m in lst:
        mat[i, j] = text_to_sympy(m)
    return mat
