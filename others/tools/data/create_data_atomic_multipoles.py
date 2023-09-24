"""
This file functions to create data of atomic multipoles.
"""
import os
import numpy as np
import sympy as sp
from joblib import Parallel, delayed
from util.atomic_multipole_matrix import create_atomic_multipole, create_atomic_multipole_harmonics_basis
from util.output_util import write_data

create_dir = __file__[: __file__.rfind("/")]


# ==================================================
header_lm = """
matrix elements of atomic multipoles in lm basis
    multipole_tag: [[matrix element]]
        multipole_tag = Xa(rank, irrep., multiplicity, component|s,k)
"""

header_jm = """
matrix elements of atomic multipoles in jm basis
    multipole_tag: [[matrix element]]
        multipole_tag = Xa(rank, irrep., multiplicity, component|s,k)
"""

header_cubic = """
matrix elements of atomic multipoles in cubic basis
    multipole_tag: [[matrix element]]
        multipole_tag = Xa(rank, irrep., multiplicity, component|s,k)
"""

header_hexagonal = """
matrix elements of atomic multipoles in hexagonal basis
    multipole_tag: [[matrix element]]
        multipole_tag = Xa(rank, irrep., multiplicity, component|s,k)
"""

header_dict = {"lm": header_lm, "jm": header_jm, "cubic": header_cubic, "hexagonal": header_hexagonal}

var_dict = {
    "lm": "_data_atomic_multipoles_lm",
    "jm": "_data_atomic_multipoles_jm",
    "cubic": "_data_atomic_multipoles_cubic",
    "hexagonal": "_data_atomic_multipoles_hexagonal",
}


# ==================================================
def create_data_atomic_multipoles(output_dir=None):
    if output_dir is not None:
        os.chdir(output_dir)

    ofile = "data_atomic_multipoles.py"

    print("=== create data_atomic_multipoles ===")
    # LMS basis
    print(" write _data_atomic_multipoles_lm ... ", end="")
    Xlmsk_lm = create_atomic_multipole(b_type="lm")
    dic = {str(tag): str(m.tolist()).replace(" ", "") for tag, m in Xlmsk_lm.items()}
    header, var = header_dict["lm"], var_dict["lm"]
    write_data(ofile, dic, header, var, mode="w")
    print("done.")

    # JM basis
    print(" write _data_atomic_multipoles_jm ... ", end="")
    Xlmsk_jm = create_atomic_multipole(b_type="jm")
    dic = {str(tag): str(m.tolist()).replace(" ", "") for tag, m in Xlmsk_jm.items()}
    header, var = header_dict["jm"], var_dict["jm"]
    write_data(ofile, dic, header, var, mode="a")
    print("done.")

    # cubic basis
    print(" write _data_atomic_multipoles_cubic ... ", end="")
    Xlmsk_cubic = create_atomic_multipole_harmonics_basis(Xlmsk_lm, harmonics_type="cubic")
    dic = {str(tag): str(m.tolist()).replace(" ", "") for tag, m in Xlmsk_cubic.items()}
    header, var = header_dict["cubic"], var_dict["cubic"]
    write_data(ofile, dic, header, var, mode="a")
    print("done.")

    # hexagonal basis
    print(" write _data_atomic_multipoles_hexagonal ... ", end="")
    Xlmsk_hexagonal = create_atomic_multipole_harmonics_basis(Xlmsk_lm, harmonics_type="hexagonal")
    dic = {str(tag): str(m.tolist()).replace(" ", "") for tag, m in Xlmsk_hexagonal.items()}
    header, var = header_dict["hexagonal"], var_dict["hexagonal"]
    write_data(ofile, dic, header, var, mode="a")
    print("done.")


# ==================================================
def check_orthogonalization(b_type="lm"):
    if b_type in ("cubic", "hexagonal"):
        Xlmsk_lm = create_atomic_multipole(b_type="lm")
        Xlmsk = create_atomic_multipole_harmonics_basis(Xlmsk_lm, b_type)
    else:
        Xlmsk = create_atomic_multipole(b_type)

    # normalize
    Xlmsk = {tag: M / M.norm() for tag, M in Xlmsk.items()}

    print(f"* check orthogonalization (btype={b_type}) ...")

    def proc(tag1, M1, tag2, M2):
        c = sp.simplify(sp.trace(M1 * M2.adjoint()))
        return tag1, tag2, c

    sub = Parallel(n_jobs=-1, verbose=10)(
        [delayed(proc)(tag1, M1, tag2, M2) for tag1, M1 in Xlmsk.items() for tag2, M2 in Xlmsk.items()]
    )

    for tag1, tag2, c in sub:
        if tag1 != tag2 and c != 0 or tag1 == tag2 and c != 1:
            raise Exception(f"multipoles are not orthogonalized. (tag1 = {tag1}, tag2 = {tag2}, c = {c})")

    print("success !")


# ==================================================
def check_completeness(b_type="lm"):
    if b_type in ("cubic", "hexagonal"):
        Xlmsk_lm = create_atomic_multipole(b_type="lm")
        Xlmsk = create_atomic_multipole_harmonics_basis(Xlmsk_lm, b_type)
    else:
        Xlmsk = create_atomic_multipole(b_type)

    print(f"=== check completeness (btype={b_type}) ===")

    # normalize
    M_list = [np.array(M / M.norm()).reshape(1024, 1) for M in Xlmsk.values()]
    diag_matrix = sum([M @ M.transpose().conjugate() for M in M_list])
    diag_matrix = sp.Matrix(diag_matrix)

    if not diag_matrix.is_Identity:
        # print(diag_matrix)
        raise Exception(f"atomic multipoles are not complete (btype={b_type}).")

    print("success !")


# ================================================== main
create_data_atomic_multipoles(output_dir=create_dir)

# check_orthogonalization(b_type="lm")
# check_orthogonalization(b_type="jm")
# check_orthogonalization(b_type="cubic")
# check_orthogonalization(b_type="hexagonal")

# check_completeness(b_type="lm")
# check_completeness(b_type="jm")
# check_completeness(b_type="cubic")
# check_completeness(b_type="hexagonal")
