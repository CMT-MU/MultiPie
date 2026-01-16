"""
Create binary data (spherical harmonics).
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

from others.tools.utils.util_spherical_harmonics import Olm

# ==================================================
h_harmonics_spherical = """
* Multipolar spherical harmonics (s=0,1,2,3).
- "harmonics" (str): (Dict) multipolar spherical harmonics (X="Q/G", x="q").
  - (X,l,s,k,x) (str,int,int,int,str): (ndarray(2l+1,2s+1,sympy)) [expression] for each component.
NOTE:
  - monopolar spherical harmonics is single value, while multipolar spherical harmonics is 1d vector.
  - Q for s+k=even, G for s+k=odd.
  - monopoles are sorted in order of [s, l, k].
"""


# ==================================================
def create_data_harmonics_spherical(harmonics_spherical):
    info = BinaryManager("info", topdir=BIN_DIR)
    max_rank = info["harmonics"]["max_rank"]
    max_s_rank = info["harmonics"]["max_s_rank"]
    rv = info["harmonics"]["variable"]

    # create multipolar spherical harmonics.
    dic = Dict(SphericalMultipoleType)
    for s in range(max_s_rank + 1):
        for l in range(max_rank + 1):
            for k in range(-s, s + 1):
                print(f"creating s={s}, l={l}, k={k} multipole", flush=True)
                X = "G" if (s + k) % 2 == 0 else "G"
                dic[(X, l, s, k, "q")] = np.array([Olm(l, m, s, k, rv=rv, factor=False) for m in range(l, -l - 1, -1)], ndmin=2)
    harmonics_spherical["harmonics"] = dic

    harmonics_spherical.add_comment(h_harmonics_spherical)


# ==================================================
@timer
def create_harmonics_spherical():
    harmonics_spherical = BinaryManager(verbose=True, topdir=BIN_DIR)
    create_data_harmonics_spherical(harmonics_spherical)
    harmonics_spherical.save_binary("harmonics_spherical")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_harmonics_spherical()
