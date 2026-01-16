"""
Create binary data (point-group symmetry operation for multipole).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging

from multipie import SOMatrixType
from multipie.util.util import timer
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

from others.tools.utils.util_symmetry_operation import symmetry_operation_matrix_multipole
from others.tools.data.data_group_pg import pg_info

# ==================================================
h_symmetry_operation_matrix = """
* SO matrix of multipole in lm basis for cubic and hexagonal subgroups.
- "cubic/hexagonal" (str): (Dict) polar (Q) or axial (G) SO matrix.
  - (X,l) (str,int): (dict) SO matrix data.
    - SO_tag (str): (ndarray(2l+1,2l+1,sympy)) SO matrix.
"""


# ==================================================
def create_so_matrix_data(so_matrix):
    info = BinaryManager("info", topdir=BIN_DIR)
    max_rank = info["harmonics"]["max_rank"]
    Oh = pg_info[info["id"]["Oh"]]
    D6h = pg_info[info["id"]["D6h"]]

    cubic = Dict(SOMatrixType)
    hexagonal = Dict(SOMatrixType)
    for l in range(max_rank + 1):
        print(f"creating rank {l}", flush=True)

        cubic[("Q", l)] = {
            tag: symmetry_operation_matrix_multipole(l, tag, hexagonal=False, axial=False, tesseral=False) for tag in Oh.SO
        }
        cubic[("G", l)] = {
            tag: symmetry_operation_matrix_multipole(l, tag, hexagonal=False, axial=True, tesseral=False) for tag in Oh.SO
        }

        hexagonal[("Q", l)] = {
            tag: symmetry_operation_matrix_multipole(l, tag, hexagonal=True, axial=False, tesseral=False) for tag in D6h.SO
        }
        hexagonal[("G", l)] = {
            tag: symmetry_operation_matrix_multipole(l, tag, hexagonal=True, axial=True, tesseral=False) for tag in D6h.SO
        }

    so_matrix["cubic"] = cubic
    so_matrix["hexagonal"] = hexagonal

    so_matrix.add_comment(h_symmetry_operation_matrix)


# ==================================================
@timer
def create_so_matrix():

    so_matrix = BinaryManager(verbose=True, topdir=BIN_DIR)
    create_so_matrix_data(so_matrix)
    so_matrix.save_binary("symmetry_operation_matrix")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_so_matrix()
