"""
Check cluster SAMB data.

- dimension of irrep.
- orthonormality.
"""

import numpy as np

from multipie.util.util_binary import BinaryManager
from multipie.core.multipie_info import __top_dir__

BIN_DIR = __top_dir__ + "others/tools/binary_data"


# ==================================================
def inner_matrix(samb):
    vec = []
    for val in samb.values():
        for i in val[0]:
            vec.append(i.reshape(-1).tolist())
    vec = np.asarray(vec).astype(float)
    on_mat = vec @ vec.T
    return on_mat


# ==================================================
def is_orthonormal_samb(samb, tol=1e-8):
    on_mat = inner_matrix(samb)
    return np.allclose(on_mat, np.eye(on_mat.shape[0]), atol=tol)


# ==================================================
def check_irrep_dimension(samb):
    dim = {"A": 1, "B": 1, "E": 2, "T": 3}
    lst = []
    for idx, val in samb.items():
        irrep = idx[2]
        if dim[irrep[0]] != len(val[0]):
            lst.append(idx)
    return lst


# ==================================================
def check(samb_set):
    for wp, samb in samb_set.items():
        lst = check_irrep_dimension(samb)
        if lst:
            print(f"   {lst} in {wp} have invalid dimensions.")
        if not is_orthonormal_samb(samb):
            print(f"   {wp} is not orthonormalized.")


# ================================================== main
if __name__ == "__main__":
    info = BinaryManager("info", verbose=False, topdir=BIN_DIR)
    samb_data = BinaryManager("cluster_samb", verbose=False, topdir=BIN_DIR)

    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        group_samb = samb_data[no]
        print(f"checking {group_tag}")
        for samb_type in ["site", "bond_s", "bond_a", "vector"]:
            samb_set = group_samb[samb_type]
            check(samb_set)

    for no in info["id_set"]["SG"]["all"]:
        group_tag = info["tag"][no]
        group_samb = samb_data[no]
        print(f"checking {group_tag}")
        for samb_type in ["site", "bond_s", "bond_a", "vector"]:
            samb_set = group_samb[samb_type]
            check(samb_set)
