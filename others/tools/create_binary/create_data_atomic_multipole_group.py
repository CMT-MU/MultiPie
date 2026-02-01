"""
Create binary data (matrix of atomic multipoles for each group).
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
from multipie.util.util_samb import orthogonalize_samb

# ==================================================
h_atomic_multipole_group = """
* Matrix of atomic multipole for all real PG.
- PG_id (str): (dict) matrix data.
  - "jml/lgs/lg" (str): (dict) matrix of atomic multipole in (J,M;L)/(L,gama,s)/(L,gamma) basis.
    - (L1,L2) (int,int): (Dict) atomic multipole data for each bra-ket block (x="q).
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,2(2L1+1),2(2L2+1),sympy), ndarray(dim,sympy)) [matrix] and [cartesian ex.] for each component.
NOTE:
  - half size of matrix for spinless "lg".
  - multipoles are sorted in order of [Q/G/T/M, s, k, l, Gamma, n, p].
"""


# ==================================================
def renumber(d, key):
    """
    Renumber Dict.

    Args:
        d (Dict): Dict.
        key (str): key for renumber.

    Returns:
        - (Dict) -- renumbered Dict.
    """
    dic = Dict(d.key_type)
    for idx, v in d.named_items():
        idx = idx._replace(**{key: -1})
        dic[idx] = dic.get(idx, []) + [v]

    rd = Dict(d.key_type)
    for idx, v in dic.named_items():
        if len(v) == 1:
            rd[idx] = v[0]
        else:
            for no, vi in enumerate(v):
                idx = idx._replace(**{key: no + 1})
                rd[idx] = vi

    return rd


# ==================================================
def orthogonalize_pg_multipole(samb, diagonal_block):
    irreps = set([i.Gamma for i in samb.named_keys()])

    samb_ortho = Dict(PGMultipoleType)
    for X in ["Q", "G", "T", "M"]:
        for Gamma in irreps:
            samb1 = samb.select(X=X, Gamma=Gamma)
            samb_ortho.update(orthogonalize_samb(samb1, diagonal_block))

    samb_ortho = renumber(samb_ortho, "n")
    samb_ortho = samb_ortho.sort(("X", ["Q", "G", "T", "M"]), "s", "k", "l", "Gamma", "n", "p")
    return samb_ortho


# ==================================================
def create_atomic_multipole_pg(atomic_multipole, harmonics, Us):
    # create spinful pg atomic multipole.
    gam = {"lms": {}, "jml": {}}
    for bt in ["lms", "jml"]:
        for (bra, ket), samb in atomic_multipole[bt].items():
            dic1 = {}
            for (X, l, irrep, n, p, _, _, _), (ex, U, _) in harmonics["harmonics"].items():
                for (head, _, s, k, _), mat in samb.select(l=l).items():
                    if (X == "Q" and head in ["Q", "T"]) or (X == "G" and head in ["G", "M"]):
                        gx = np.vectorize(sp.expand)(np.einsum("mg,mij->gij", U, mat))
                        dic1[(head, l, irrep, s, k)] = dic1.get((head, l, irrep, s, k), []) + [(gx, ex)]

            # renumber multiplicity.
            dic = Dict(PGMultipoleType)
            for (head, l, irrep, s, k), v in dic1.items():
                if len(v) == 1:
                    vi = v[0]
                    dic[(head, l, irrep, -1, -1, s, k, "q")] = vi
                else:
                    for no, vi in enumerate(v):
                        dic[(head, l, irrep, no + 1, -1, s, k, "q")] = vi

            gam[bt][(bra, ket)] = dic

    # convert lm and lms to cubic and hexagonal one.
    gam1 = {}

    # spinful.
    dic1 = {}
    for (bra, ket), samb in gam["jml"].items():
        dic1[(bra, ket)] = orthogonalize_pg_multipole(samb, bra == ket)

    gam1["jml"] = dic1

    # spinful.
    dic1 = {}
    for (bra, ket), samb in gam["lms"].items():
        Ubra = Us[bra].conjugate().T
        Uket = Us[ket]

        dic = Dict(PGMultipoleType)
        for idx, (mat, ex) in samb.items():
            matp = np.vectorize(sp.expand)(Ubra @ mat @ Uket)
            dic[tuple(idx)] = (matp, ex)

        dic1[(bra, ket)] = orthogonalize_pg_multipole(dic, bra == ket)

    gam1["lgs"] = dic1

    # spinless.
    gam1["lg"] = {  # s=0, up-spin only.
        block: Dict(samb.key_type, {idx: (mat[:, ::2, ::2] * sp.sqrt(2), ex) for idx, (mat, ex) in samb.select(s=0).items()})
        for block, samb in gam1["lgs"].items()
    }

    return gam1


# ==================================================
def create_atomic_multipole_group_data(atomic_multipole, gam, info, U):
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)

    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        gam[no] = create_atomic_multipole_pg(atomic_multipole, harmonics[no], U)

    gam.add_comment(h_atomic_multipole_group)


# ==================================================
@timer
def create_atomic_multipole_group():
    info = BinaryManager("info", topdir=BIN_DIR)

    # create U-matrix for point-group atomic basis.
    U = []
    for lst in info["harmonics"]["atomic_basis"]["spinless"]["lg"].values():
        UL = []
        for name in lst:
            UL.append(info["harmonics"]["basis_function"][name][1].tolist())
        U1 = np.asarray(UL).T
        sh = U1.shape
        Uf = np.full((2 * sh[0], 2 * sh[1]), sp.S(0))
        Uf[::2, ::2] = U1
        Uf[1::2, 1::2] = U1
        U.append(Uf)

    atomic_multipole = BinaryManager("atomic_multipole", topdir=BIN_DIR)

    atomic_multipole_group = BinaryManager(verbose=True, topdir=BIN_DIR)
    create_atomic_multipole_group_data(atomic_multipole, atomic_multipole_group, info, U)
    atomic_multipole_group.save_binary("atomic_multipole_group")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_atomic_multipole_group()
