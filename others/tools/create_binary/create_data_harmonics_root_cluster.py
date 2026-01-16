"""
Create binary data (cluster harmonics on root cluster).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import sympy as sp
import numpy as np
import logging

from multipie import PGMultipoleType
from multipie.util.util import timer
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

# ==================================================
h_harmonics_root_cluster = """
* Root-cluster harmonics for all real PG.
- PG_id (str): (dict) root-cluster harmonics.
  - "monopole" (str): (Dict) monopole harmonics.
    - (X,l,Gamma,n,p,0,0,"q") (str,int,str,int,int,int,int,str): (ndarray(dim,npts,sympy), ndarray(dim,sympy)) [basis] and [cartesian ex.] for each component.
  - "dipole" (str): (Dict) dipolar harmonics.
    - (X,l,Gamma,n,p,1,k,"q") (str,int,str,int,int,int,int,str): (ndarray(dim,npts,3,sympy), ndarray(dim,sympy), ndarray(dim,sympy)) [basis], [cartesian ex.], and [vector cartesian ex.] for each component.
NOTE:
  - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
  - monopole harmonics are sorted in order of [Gamma, l, Q/G, n, p].
  - dipole harmonics are sorted in order of [Gamma, l, Q/G, k, n, p].
"""


# ==================================================
@timer
def create_harmonics_root_cluster():
    info = BinaryManager("info", topdir=BIN_DIR)
    rc = BinaryManager("root_cluster", topdir=BIN_DIR)
    monopole = BinaryManager("harmonics", topdir=BIN_DIR)
    multipole = BinaryManager("harmonics_multipole", topdir=BIN_DIR)

    rv = sp.symbols("x y z", real=True)
    rc_pos = rc["position"]
    n = len(rc_pos)

    def rc_exp(basis, dipole=False):
        data = np.asarray([np.vectorize(lambda i: i.subs(dict(zip(rv, r))).expand())(basis) for r in rc_pos])
        if dipole:
            data = data.transpose(1, 0, 2).reshape(-1, n, 3)
        else:
            data = data.transpose(1, 0).reshape(-1, n)
        return data

    # create point-group harmonics.
    rc_harmonics = BinaryManager(verbose=True, comment=h_harmonics_root_cluster, topdir=BIN_DIR)
    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        monopole_pg = monopole[no]["harmonics"]
        dipole_pg = multipole[no].select(s=1, x="q")  # polar/axial, internal polar dipole only.

        m_harm = Dict(PGMultipoleType)
        for idx, basis in monopole_pg.items():
            m_harm[idx] = (rc_exp(basis[0], dipole=False), basis[0])
        d_harm = Dict(PGMultipoleType)
        for idx, basis in dipole_pg.items():
            d_harm[idx] = (rc_exp(basis[1], dipole=True), basis[0], basis[2])

        rc_harmonics[no] = {"monopole": m_harm, "dipole": d_harm}

    rc_harmonics.save_binary("harmonics_root_cluster")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_harmonics_root_cluster()
