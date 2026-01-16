"""
Create binary data (information for multipole harmonics in each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp

from multipie import PGMultipoleType
from multipie.util.util import timer
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

# ==================================================
h_harmonics_multipole = """
* Multipolar harmonics harmonics for all PG.
- PG_id (str): (Dict) multipolar harmonics (s=1,2,3) with polar (q) and axial (g) internal basis.
  - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,sympy), ndarray(dim,2s+1,sympy), ndarray(dim,sympy)) [cartesian ex.], [multipole basis ex.], and [multipole ex.] for each component.
NOTE:
  - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
  - multipoles are sorted in order of [q/g, s, Gamma, l, Q/G, k, n, p].
"""


# ==================================================
def create_harmonics_multipole_data(harmonics, sh, ivar, max_s_rank):
    harmonics_multipole = Dict(PGMultipoleType)

    for s in range(1, max_s_rank + 1):
        dic1 = {}
        for (X, l, irrep, n, p, _, _, x), (ex, U, _) in harmonics["harmonics"].items():
            print(f"\r  s={s}, l={l}, X={X}, {irrep}    ", end="", flush=True)
            harm = sh["harmonics"].select(s=s, l=l)
            for (_, _, _, k, _), v in harm.items():
                mharm = np.vectorize(sp.expand)(U.T @ v)
                if (X == "Q" and (s + k) % 2 == 0) or (X == "G" and (s + k) % 2 == 1):
                    intra = "q"
                elif (X == "Q" and (s + k) % 2 == 1) or (X == "G" and (s + k) % 2 == 0):
                    intra = "g"
                if l < 5:
                    mharm = np.vectorize(sp.factor)(mharm)
                m_ex = mharm @ ivar[(s, intra)]

                if not (m_ex == 0).all():
                    dic1[(intra, l, X, irrep, k)] = dic1.get((intra, l, X, irrep, k), []) + [(ex, mharm, m_ex)]
        print()

        # renumber multiplicity.
        for (intra, l, X, irrep, k), v in dic1.items():
            if len(v) == 1:
                vi = v[0]
                harmonics_multipole[(X, l, irrep, -1, -1, s, k, intra)] = vi
            else:
                for no, vi in enumerate(v):
                    harmonics_multipole[(X, l, irrep, no + 1, -1, s, k, intra)] = vi

    harmonics_multipole = harmonics_multipole.sort(("x", ["q", "g"]), "s", "Gamma", "l", ("X", ["Q", "G"]), "k", "n", "p")

    return harmonics_multipole


# ==================================================
@timer
def create_harmonics_multipole():
    info = BinaryManager("info", topdir=BIN_DIR)
    sh = BinaryManager("harmonics_spherical", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)

    max_s_rank = info["harmonics"]["max_s_rank"]
    ivar = info["harmonics"]["internal_basis"]

    harmonics_multipole = BinaryManager(verbose=True, topdir=BIN_DIR)
    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        harmonics_multipole[no] = create_harmonics_multipole_data(harmonics[no], sh, ivar, max_s_rank)

    harmonics_multipole.add_comment(h_harmonics_multipole)
    harmonics_multipole.save_binary("harmonics_multipole")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_harmonics_multipole()
