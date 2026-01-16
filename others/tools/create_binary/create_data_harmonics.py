"""
Create binary data (information for harmonics in each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp
from collections import namedtuple

from multipie import PGMultipoleType, IrrepListType
from multipie.util.util import timer, str_to_sympy
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager

from others.tools.data.data_top_harmonics import top_harmonics, top_harmonics_tag
from others.tools.data.data_compatibility_relation import compatibility_relation
from others.tools.utils.util_samb import expand_component

# ==================================================
h_harmonics = """
* Harmonics for all PG.
- PG_id (str): (dict) harmonics data for each PG.
  - "harmonics" (str): (Dict) polar (Q) and axial (G) harmonics (s=k=0,x="q").
    - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,sympy), ndarray(2l+1,dim,sympy), ndarray(dim,sympy)) [cartesian ex.], u-matrix, <m|gamma>, and [tesseral ex.] for each component.
  - "irrep" (str): (Dict) Q and G harmonics grouped by irrep.
    - (X,Gamma) (str,str): ([(int,int)]) list of (l,n).
  - "index" (str): (dict) harmonics index by (X,expression).
    - (X,ex) (str,str): ((X,l,Gamma,n,p,s,k,x),comp).
NOTE:
  - u-matrix is in order of m=[l, l-1, ..., -l][gamma].
  - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
  - for complex expresson, use point group tag with "-c".
  - multipoles are sorted in order of [Gamma, l, Q/G, n, p].
"""


# ==================================================
def create_top_harmonics():
    """
    Create top harmonics (polar up to rank 11 for Oh and D6h).

    Returns:
        - (dict) -- top harmonics data.

    Notes:
        - "Oh/D6h" (dict): data.
          {("Oh/D6h",l,irrep.,multiplicity,component): (tesseral ex., cartesian ex., u-matrix)} ({(str,int,str,int,int): (sympy,sympy,[sympy(2l+1)])}).
        - u-matrix <m|gamma> is in order of m=[l, l-1, ..., -l].
        - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
        - complex expressions (postfixed with a/b in irrep) are added for creating subgroup harmonics.
        - axial irrep is obtained by interchanging "g" and "u".
        - representation matrices of symmetry operations are common upto sign and order for the same irrep.
    """
    info = BinaryManager("info", topdir=BIN_DIR)
    x, y, z = info["harmonics"]["variable"]
    sub = {"x": x, "y": y, "z": z}

    dic_harm = {}
    for tag, (t_ex, c_ex, u_mat) in top_harmonics.items():
        t_ex = str_to_sympy(t_ex, subs=sub)
        c_ex = str_to_sympy(c_ex, subs=sub)
        u_mat = str_to_sympy(u_mat)
        dic_harm[tag] = (t_ex, c_ex, u_mat)
    dic = {"data": dic_harm, "tag": top_harmonics_tag}

    return dic


# ==================================================
def create_group_harmonics(group_tag):
    top, cr = compatibility_relation[group_tag]
    top_harm = create_top_harmonics()
    top_harmonics_pg = top_harm["tag"][top]

    # create mapping.
    mapping = {"Q": {}, "G": {}}
    for rank, irrep, mul, comp in top_harmonics_pg:
        polar, axial = cr[(irrep, comp)]
        if polar[0].count("a") > 0 or polar[0].count("b") > 0:
            if comp == 0:
                irrep_p = irrep + "a"
                comp_p = -1
            elif comp == 1:
                irrep_p = irrep + "b"
                comp_p = -1
            else:
                irrep_p = irrep
                comp_p = comp
        else:
            irrep_p = irrep
            comp_p = comp
        if axial[0].count("a") > 0 or axial[0].count("b") > 0:
            if comp == 0:
                irrep_a = irrep + "a"
                comp_a = -1
            elif comp == 1:
                irrep_a = irrep + "b"
                comp_a = -1
            else:
                irrep_a = irrep
                comp_a = comp
        else:
            irrep_a = irrep
            comp_a = comp
        p_idx = (rank, polar[0], polar[1])
        a_idx = (rank, axial[0], axial[1])
        mapping["Q"][p_idx] = mapping["Q"].get(p_idx, []) + [(polar[2], (top, rank, irrep_p, mul, comp_p))]
        mapping["G"][a_idx] = mapping["G"].get(a_idx, []) + [(axial[2], (top, rank, irrep_a, mul, comp_a))]

    # renumber multiplicity.
    harmonics_data = {"Q": {}, "G": {}}
    for pa, mapping_pa in mapping.items():
        for (rank, irrep, comp), data in mapping_pa.items():
            if len(data) == 1:
                v = data[0]
                harm = top_harm["data"][v[1]]
                harmonics_data[pa][(rank, irrep, -1, comp)] = (v[0] * harm[0], v[0] * harm[1], v[0] * harm[2])
            else:
                for no, v in enumerate(data):
                    harm = top_harm["data"][v[1]]
                    harmonics_data[pa][(rank, irrep, no + 1, comp)] = (v[0] * harm[0], v[0] * harm[1], v[0] * harm[2])

    # irrep list.
    irrep_list = Dict(IrrepListType)
    for pa, data_pa in harmonics_data.items():
        for rank, irrep, mul, comp in data_pa.keys():
            if comp == -1 or comp == 0:
                irrep_list[(pa, irrep)] = irrep_list.get((pa, irrep), []) + [(rank, mul)]
    irrep_list = irrep_list.sort(("X", ["Q", "G"]))

    # reformat.
    tmp = Dict(namedtuple("tmp", ["l", "X", "irrep", "n", "comp"]))
    for X in ["Q", "G"]:
        for (rank, irrep, mul, comp), v in harmonics_data[X].items():
            tmp[(rank, X, irrep, mul, comp)] = v
    tmp = tmp.sort("l", ("X", ["Q", "G"]), "irrep", "n", "comp")

    harmonics1 = Dict(namedtuple("Harmonics", ["l", "X", "irrep", "n"]))
    for (l, X, irrep, mul, comp), (t, c, u) in tmp.items():
        harmonics1[(l, X, irrep, mul)] = harmonics1.get((l, X, irrep, mul), []) + [(t, c, u)]

    harmonics = Dict(PGMultipoleType)
    for (l, X, Gamma, n), v in harmonics1.items():
        tset = []
        cset = []
        uset = []
        for t, c, u in v:
            tset.append(t)
            cset.append(c)
            uset.append(u)
        tset = np.asarray(tset)
        if group_tag.count("-c") == 0:  # to avoid factorization for non-complex point group due to extremely slow.
            cset = np.vectorize(sp.factor)(np.asarray(cset))
        uset = np.asarray(uset).T
        harmonics[(X, l, Gamma, n, -1, 0, 0, "q")] = (cset, uset, tset)
    harmonics = harmonics.sort("Gamma", "l", ("X", ["Q", "G"]))

    # reverse dict.
    index = {}
    for X in ["Q", "G"]:
        q_harm = expand_component(harmonics.select(X=X), vector=False)
        for idx, val in q_harm.items():
            v = str(sp.expand(val[0])).replace(" ", "")
            index[(X, v)] = idx

    harmonics_data = {"harmonics": harmonics, "irrep": irrep_list, "index": index}

    return harmonics_data


# ==================================================
@timer
def create_harmonics():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)

    harmonics = BinaryManager(verbose=True, topdir=BIN_DIR)
    for no in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        g_info = group[no]
        group_tag = g_info["info"].tag
        print("creating", group_tag, flush=True)
        harmonics[no] = create_group_harmonics(group_tag)

    harmonics.add_comment(h_harmonics)
    harmonics.save_binary("harmonics")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_harmonics()
