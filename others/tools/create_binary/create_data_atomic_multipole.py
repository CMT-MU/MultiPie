"""
Create binary data (matrix of atomic multipoles).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np

from multipie import SphericalMultipoleType
from multipie.util.util import timer
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager
from multipie.util.util_atomic_multipole import create_atomic_multipole_matrix

from others.tools.data.data_definition import braket_block

# ==================================================
h_atomic_multipole = """
* Matrix of atomic multipole.
- "lms/jml/lm" (str): (dict) matrix in LMs/JML/LM basis.
  - (L1,L2) (int/int): (Dict) atomic multipole data for each bra-ket block (x="q").
    - (X,l,s,k,x) (str,int,int,int,str): (ndarray(2l+1,2(2L1+1),2(2L2+1),sympy)) list of matrix.
NOTE:
  - half size of matrix for spinless "lm".
  - multipoles are sorted in order of [Q/G/T/M, s, k, l].
"""


# ==================================================
def create_atomic_multipole_data(atomic_multipole):
    for bt in ["lms", "jml"]:
        print(f"creating for {bt} basis", flush=True)
        Xlsk = create_atomic_multipole_matrix(b_type=bt)

        atomic_multipole[bt] = {}
        for bi, bj in braket_block:
            samb = Dict(SphericalMultipoleType)
            for (X, l, s, k), mat in Xlsk.items():
                m = np.asarray(mat)[:, 2 * bi**2 : 2 * (bi + 1) ** 2, 2 * bj**2 : 2 * (bj + 1) ** 2]
                if not (m == 0).all():
                    samb[(X, l, s, k, "q")] = m
            samb = samb.sort(("X", ["Q", "G", "T", "M"]), "s", "k", "l")
            atomic_multipole[bt][(bi, bj)] = samb

    print(f"creating for lm basis", flush=True)
    atomic_multipole["lm"] = {  # s=0, up-spin only.
        block: Dict(samb.key_type, {idx: mat[:, ::2, ::2] for idx, mat in samb.select(s=0).items()})
        for block, samb in atomic_multipole["lms"].items()
    }

    atomic_multipole.add_comment(h_atomic_multipole)


# ==================================================
@timer
def create_atomic_multipole():
    atomic_multipole = BinaryManager(verbose=True, topdir=BIN_DIR)
    create_atomic_multipole_data(atomic_multipole)
    atomic_multipole.save_binary("atomic_multipole")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_atomic_multipole()
