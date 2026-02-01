"""
Create binary data (overall information of groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import sympy as sp
import numpy as np

from multipie import InternalBasisType
from multipie.util.util import normalize_vector, timer, str_to_sympy
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

from others.tools.data.data_definition import (
    max_rank,
    max_s_rank,
    rc_representative_site,
    rc_scale,
    rc_representative_vector,
    rc_point_group_set,
    internal_basis,
    atomic_basis,
    rank_name,
    orbital,
    wannier90,
    tesseral,
    orbital_decomposition,
)
from others.tools.data.data_compatibility_relation import compatibility_relation
from others.tools.data.data_group_pg import pg_id_set, pg_tag, pg_id
from others.tools.data.data_group_sg import sg_id_set, sg_tag, sg_id
from others.tools.data.data_group_mpg import mpg_id_set, mpg_tag, mpg_id
from others.tools.data.data_group_msg import msg_id_set, msg_tag, msg_id
from others.tools.data.data_response_tensor import rt_cartesian_multipole, rt_multipole_list

# ==================================================
h_info = """
* Information of all groups.
- "id_set" (str): (dict) info. of group's IDs.
  - "PG/SG/MPG/MSG" (str):  (dict) ID set.
    - PG: all/crystal/complex/irrep.
    - SG: all/crystal/PG.
    - MPG: all/crystal/PG/type.
    - MSG: all/crystal/PG/SG/MPG/type.
- "tag" (str): (dict) from id to tag.
- "id" (str):  (dict) from tag to id.
- "character" (str): (dict) info. of character.
  - "alias" (str): (dict) alias.
    - "wp/wm" (str): (sympy) definition of wp and wm, exp(2 pi i/3) or exp(-2 pi i/3).
  - "compatibility_relation" (str): (dict) compatibility relation.
    - (PG_tag1,PG_tag2) (str,str): ([(str,str)]) list of compatibility relation, (irrep1,irrep2).
- "root_cluster" (str): (dict) info. of root cluster.
  - "point_group" (str): ([str]) list of PG_tag for root cluster.
  - "rep_site_c" (str): (ndarray(3,sympy)) representative site (fractional for cubic, r=1; = cartesian).
  - "rep_site_h" (str): (ndarray(3,sympy)) representative site (fractional for hexagonal).
  - "scaled_rep_site_c" (str): (ndarray(3,sympy)) scaled representative site (fractional for cubic, r=1/5; = cartesian).
  - "scaled_rep_site_h" (str): (ndarray(3,sympy)) scaled representative site (fractional for hexagonal).
  - "rep_vector_c" (str): (ndarray(3,sympy)) representative vector (fractional for cubic, v=1; = cartesian).
  - "rep_vector_h" (str): (ndarray(3,sympy)) representative vector (fractional for hexagonal).
- "harmonics" (str): (dict) info. of harmonics.
  - "max_rank" (str): (int) max rank of harmonics.
  - "max_s_rank" (str): (int) max internal rank.
  - "variable" (str): (ndarray(3,sympy)) coordinate variable, [x,y,z].
  - "internal_basis" (str): (Dict) internal basis.
    - (s,x) (int,str): (ndarray(2s+1,sympy)) list of basis.
  - "atomic_basis" (str): (dict) atomic basis.
    - "spinless" (str): (dict) spinless basis.
      - "lm/lg" (str): (dict) basis.
        - L (int): ([str]) list of basis, M or gamma.
    - "spinful" (str): (dict) spinful basis.
      - "jml/lms/lgs" (str): (dict) basis.
        - L (int): ([str]) list of basis, (J,M) or (M,s) or (gamma,s).
  - "basis_function" (str): (dict) orbital basis function.
    - name (str): ((sympy,ndarray(2L+1,sympy))) cartesian expression and u-matrix, <m|gamma>.T.
  - "wannier90" (str): (dict) alias for Wannier90.
    - tag (str): (str) name corresponding to Wannier90 tag.
  - "tesseral" (str): (dict) tesseral tag.
    - name (str): ((int,int,str)) tesseral tag (L,M,"c/s").
  - "rank_name" (str): (dict) conversion between rank and name.
    - name/rank (str/int): (int/str) rank/name.
  - "orbital_decomposition" (str): (dict) spherical harmonics irrep. decomposition.
    - PG_tag (str): (dict) spherical harmonics irrep. decomposition for each PG.
      - L (int): ([(int,str)]) list of (n,irrep).
  NOTE:
    - see "multipie/data/data_definition.py" for detail.
- "response_tensor" (str): (dict) info. of response tensor.
  - "multipole" (str): (dict) multipole info.
    - "cubic/hexagonal" (str): ([str]) list of multipoles.
  - "cartesian_multipole" (str): (dict) linear combination info.
    - "cubic/hexagonal" (str): { comp (str): ([(int,str)]) }, relation between cartesian component and multipole.
"""


# ==================================================
def create_group_info(info):
    info["id_set"] = {
        "PG": pg_id_set,
        "SG": sg_id_set,
        "MPG": mpg_id_set,
        "MSG": msg_id_set,
    }
    info["tag"] = pg_tag | sg_tag | mpg_tag | msg_tag
    info["id"] = pg_id | sg_id | mpg_id | msg_id


# ==================================================
def create_character_info(info):
    # complex character.
    subs = {"wp": sp.exp(2 * sp.pi * sp.I / 3), "wm": sp.exp(-2 * sp.pi * sp.I / 3)}

    # compatibility table.
    sconv = {1: "", -1: "-"}
    cconv = {-1: "", 0: "(x)", 1: "(y)", 2: "(z)"}

    ct = {}
    for tag in ["Oh", "D6h"]:
        dp = {}
        for pg, (top, dic) in compatibility_relation.items():
            if top != tag:
                continue
            irrep_p = [sconv[vp] + ip + cconv[cp] for ((ip, cp, vp), (ia, ca, va)) in dic.values()]
            dp[pg] = irrep_p
        for pg1, tbl1 in dp.items():
            for pg2, tbl2 in dp.items():
                ct[(pg1, pg2)] = [(i, j) for i, j in zip(tbl1, tbl2)]

    dic = {"alias": subs, "compatibility_relation": ct}
    info["character"] = dic


# ==================================================
def create_root_cluster_info(info):
    # representative site.
    rep_site = str_to_sympy(rc_representative_site)
    rep_site = normalize_vector(rep_site)
    scaled_rep_site = np.vectorize(sp.simplify)(rep_site / rc_scale)
    x, y, z = rep_site[0], rep_site[1], rep_site[2]
    fractional = np.vectorize(sp.simplify)([x + y / sp.sqrt(3), 2 * y / sp.sqrt(3), z])
    sx, sy, sz = scaled_rep_site[0], scaled_rep_site[1], scaled_rep_site[2]
    scaled_fractional = np.vectorize(sp.simplify)([sx + sy / sp.sqrt(3), 2 * sy / sp.sqrt(3), sz])

    # representative bond.
    rep_bond = str_to_sympy(rc_representative_vector)
    rep_bond = normalize_vector(rep_bond)
    vx, vy, vz = rep_bond[0], rep_bond[1], rep_bond[2]
    fractional_bond = np.vectorize(sp.simplify)([vx + vy / sp.sqrt(3), 2 * vy / sp.sqrt(3), vz])

    # point group for root cluster.
    pg_set = rc_point_group_set

    dic = {
        "point_group": pg_set,
        "rep_site_c": rep_site,
        "rep_site_h": fractional,
        "scaled_rep_site_c": scaled_rep_site,
        "scaled_rep_site_h": scaled_fractional,
        "rep_vector_c": rep_bond,
        "rep_vector_h": fractional_bond,
    }

    info["root_cluster"] = dic


# ==================================================
def create_harmonics_info(info):
    # variable.
    r = sp.symbols("x y z", real=True)

    # internal basis.
    bv = Dict(InternalBasisType)
    for (s, x), basis in internal_basis.items():
        bv[(s, x)] = [sp.Symbol(i[0] + "_" + i[1:], real=True) for i in basis]

    # orbital dict.
    orbital_d = {
        name: (str_to_sympy(ex, subs={"x": r[0], "y": r[1], "z": r[2]}), str_to_sympy(u)) for name, (ex, u) in orbital.items()
    }

    # info. dict.
    dic = {
        "max_rank": max_rank,
        "max_s_rank": max_s_rank,
        "variable": r,
        "internal_basis": bv,
        "atomic_basis": atomic_basis,
        "rank_name": rank_name,
        "wannier90": wannier90,
        "tesseral": tesseral,
        "basis_function": orbital_d,
        "orbital_decomposition": orbital_decomposition,
    }

    info["harmonics"] = dic


# ==================================================
def create_response_tensor_info(info):
    info["response_tensor"] = {"multipole": rt_multipole_list, "cartesian_multipole": rt_cartesian_multipole}


# ==================================================
@timer
def create_info():
    info = BinaryManager(verbose=True, topdir=BIN_DIR)

    create_group_info(info)
    create_character_info(info)
    create_root_cluster_info(info)
    create_harmonics_info(info)
    create_response_tensor_info(info)

    info.add_comment(h_info)

    info.save_binary("info")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_info()
