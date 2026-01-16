"""
Create binary data (information for representation matrix in each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp

from multipie import RepMatrixType
from multipie.util.util import timer
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

# ==================================================
h_representation_matrix = """
* Representation matrix for all PG.
- PG_id (str): (dict) representation matrix for each group.
  - "tag" (str): ([str]) list of SO_tag.
  - "matrix" (str): (Dict) representation matrix.
    - Gamma (str): (ndarray(n,dim,dim,sympy)) list of rep. matrix.
NOTE:
  - rep_matrix is in the same order of SO.
"""


# ==================================================
def create_rep_matrix_data(group_pg, harmonics, so_matrix, axial, first=False):
    """
    Create representation matrix.

    Args:
        group_pg (dict): group data for given point group.
        harmonics (dict): harmonics data for given point group.
        so_matrix (dict): SO matrix data for given point group.
        axial (bool): for axial harmonics ?
        first (bool, optional): return first one only ?

    Returns:
        - (dict) -- representation matrix in the same order of symmetry operations, {irrep: [matrix]}.
    """
    head = "G" if axial else "Q"
    so = so_matrix.select(X=head)
    rep_harmonics = harmonics["irrep"].select(X=head)

    # get representative harmonics for each irrep.
    h = {}
    for (_, irrep), lst in rep_harmonics.items():
        if first:
            rank, mul = lst[0]
            h_irrep = (rank, harmonics["harmonics"][(head, rank, irrep, mul, -1, 0, 0, "q")][1])  # u-matrix.
        else:
            h_irrep = {}
            for rank, mul in lst:
                h_irrep[(rank, mul)] = harmonics["harmonics"][(head, rank, irrep, mul, -1, 0, 0, "q")][1]  #  u-matrix.
        h[irrep] = h_irrep

    # get rep matrix.
    so_tag = group_pg["symmetry_operation"]["tag"]
    ch = list(group_pg["character"]["table"].keys())
    rep_mat = {"tag": so_tag}

    dic = Dict(RepMatrixType)
    for irrep, lst in h.items():
        if first:
            rank, harm = lst
            Um = harm
            Umi = np.vectorize(sp.conjugate)(Um.T)
            mat = [Umi @ so[(head, rank)][tag] @ Um for tag in so_tag]
            rep_mat_irrep = np.vectorize(sp.expand_complex)(mat)
        else:
            rep_mat_irrep = {}
            for (rank, mul), harm in lst.items():
                Um = harm
                Umi = np.vectorize(sp.conjugate)(Um.T)
                mat = [Umi @ so[(head, rank)][tag] @ Um for tag in so_tag]
                rep_mat_irrep[(rank, mul)] = np.vectorize(sp.expand_complex)(mat)
        dic[irrep] = rep_mat_irrep
    dic = dic.sort(("Gamma", ch))

    rep_mat["matrix"] = dic

    return rep_mat


# ==================================================
@timer
def create_rep_matrix():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)
    so_matrix = BinaryManager("symmetry_operation_matrix", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)

    rep_matrix = BinaryManager(verbose=True, topdir=BIN_DIR)

    for no in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        g_info = group[no]
        hexagonal = g_info["info"].hexagonal_g
        so = so_matrix["hexagonal"] if hexagonal else so_matrix["cubic"]
        harm = harmonics[no]

        # rep_matrix is common both for polar and axial, thus write polar only.
        rep_matrix[no] = create_rep_matrix_data(g_info, harm, so, axial=False, first=True)

    rep_matrix.add_comment(h_representation_matrix)
    rep_matrix.save_binary("representation_matrix")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_rep_matrix()
