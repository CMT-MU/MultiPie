"""
Create binary data (information for each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp

from multipie import SOMatrixType
from multipie.util.util_dict import Dict
from multipie.util.util_binary import BinaryManager
from multipie.util.util import timer, str_to_sympy, normalize_vector
from multipie.util.util_tag import TagIrrep
from multipie.util.util_crystal import (
    convert_to_cartesian_hexagonal,
    convert_to_fractional_hexagonal,
    convert_to_fractional_hexagonal_matrix,
    convert_to_primitive_matrix,
    convert_to_primitive,
    t_dict,
)

from others.tools.data.data_group_pg import pg_info
from others.tools.data.data_group_sg import sg_info
from others.tools.data.data_group_mpg import mpg_info
from others.tools.data.data_group_msg import msg_info
from others.tools.data.data_character_table import character, character_plus
from others.tools.data.data_product_so_pg import point_group_product
from others.tools.data.data_wyckoff_site import wyckoff_site_data
from others.tools.data.data_wyckoff_mpg import wyckoff_mpg_data
from others.tools.data.data_wyckoff_msg import wyckoff_msg_data
from others.tools.utils.util_symmetry_operation import symmetry_operation_matrix, symmetry_operation_matrix_multipole
from others.tools.create_binary.create_data_group_mapping import create_group_mapping
from others.tools.create_binary.create_data_group_bond import create_group_bond

# ==================================================
h_group = """
* Information for all PG ("PG").
- PG_id (str): (dict) id, "PG:{ID}".
  - "info" (str): (namedtuple) info.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: = id.
    - SG: associated SG, 1st SG in the same PGs.
    - MPG: associated MPG, 1st MPG with type II in the same PGs.
    - MSG: associated MSG, PG -> MPG -> MSG.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    - "product" (str): (dict) product table.
      - (SO_tag,SO_tag) (str,str): (str) SO_tag of product.
    - "mapping" (str): (dict) mapping of SO to parent PG.
      - SO_tag (str): ([str]) list of SO_tag.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
      - monopole (s=0) and dipole (s=1) basis are [1] and [x,y,z].
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).
  - "character" (str): (dict) character table of group.
    - "conjugacy" (str): ([[str]]) list of SO_tag in each conjugacy class.
    - "dimension" (str): (dict) dimension of irrep.
      - irrep_tag (str): (int) dimension.
    - "table" (str): (dict) character table (first SO in each conjugacy class).
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "table_full" (str): (dict) character table for all SO.
      - irrep_tag (str): (ndarray(n,sympy) list of character.
    - "polar_axial_conversion" (str): (dict) conversion between polar and axial.
      - irrep_tag (str): (str) converted irrep_tag.
    - "symmetric_product" (str): (dict) symmetric product decomposition.
      - (irrep_tag,irrep_tag) (str,str): [(int,str)] list of (n, irrep_tag).
    - "anti_symmetric_product" (str): (dict) anti-symmetric product decomposition.
      - irrep_tag (str): ([(int,str)]) list of (n,irrep_tag).
    NOTE:
      - wp = exp(+2pi i/3), wm = exp(-2pi i/3).
      - for complex character, use PG_tag with "-c", otherwise without "-c", real version is given.

* Information for all SG ("SG").
- SG_id (str): (dict) id, "SG:{ID}".
  - "info" (str): (namedtuple) info.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting comment.
    - PG: associated PG, unique.
    - SG: = id.
    - MPG: associated MPG, SG -> PG -> MPG.
    - MSG: associated MSG, 1st MSG with type II in the same SGs.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "generator" (str): ([str]) list of generator SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "cartesian_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, primitive).
    - "fractional_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, primitive).
    - "plus_set" (str): (ndarray(n,3,sympy)) list of plus set vector (fractional, conventional).
    - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
      - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SOs of SG and PG (without translational part) are in the same order.
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site/bond of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, primitive).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
        - "bond" (str): ([str]) list of bw_tag.
    - "bond" (str): (dict) info. of Wyckoff bond.
      - bw_tag (str): (dict) bond Wyckoff tag.
        - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, conventional, no plus set).
        - "primitive" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, primitive).
        - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
    NOTE:
      - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
      - Wyckoff sites are in the same order of ITA (BCS).

* Information for all MPG ("MPG").
- MPG_id (str): (dict) id, "MPG:{PG}.{no}.{ID}".
  - "info" (str): (namedtuple) info.
    - tag: MPG_id (PG.no.ID).
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, MPG -> MSG -> SG.
    - MPG: = id.
    - MSG: associated MSG, 1st MSG in the same MPGs.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III.
    - SO: symmetry operations.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
    - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
  - "active_multipole" (list): active MPs list, [str].

* Information for all MSG ("MSG").
- MSG_id (str): (dict) id, "MSG:{SG}.{no}".
  - "info" (str): (namedtuple) info.
    - tag: BNS_id (SG.no).
    - BNS: Belov-Neronova-Smirnova (international) notation in LaTeX.
    - OG: Opechowski-Guccione notation in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, unique.
    - MPG: associated MPG, unique.
    - MSG: = id.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III, IV.
    - SO: symmetry operations.
  - "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
    - "tag" (str): ([str]) list of SO tag.
    - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
    - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
    - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
    - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
    NOTE:
      - SOs are in the same order of ITA (BCS).
      - SO matrices are in the same order of "tag".
  - "wyckoff" (str): (dict) Wyckoff site of group.
    - "site" (str): (dict) info. of Wyckoff site.
      - sw_tag (str): (dict) site Wyckoff tag.
        - "symmetry" (str): (str) site symmetry.
        - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
        - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
        - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
"""


# ==================================================
def create_character(table, so, conjugacy):
    """
    Create character for all symmetry operations.

    Args:
        table (dict): character table for irrep and character of conjugacy class.
        so (list): list of symmetry operation tag.
        conjugacy (list): list of symmetry operation in cojugacy class.

    Returns:
        - (dict) -- character table, {irrep, character}.
    """
    ch_all = {}
    for irrep, ch in table.items():
        ch_all_irrep = []
        for s in so:
            pos = [(i, j) for i, row in enumerate(conjugacy) for j, element in enumerate(row) if element == s]
            ch_all_irrep.append(ch[pos[0][0]])
        ch_all_irrep = np.asarray(ch_all_irrep)

        ch_all_irrep = np.vectorize(
            lambda x: x.subs({"wp": -1 / sp.S(2) + sp.sqrt(3) * sp.I / 2, "wm": -1 / sp.S(2) - sp.sqrt(3) * sp.I / 2})
        )(ch_all_irrep)
        ch_all[irrep] = ch_all_irrep

    return ch_all


# ==================================================
def create_group_symmetry_operation(g_info, group_tag):
    dic = {}

    hexagonal = g_info.hexagonal_g
    gen, so = g_info.SO_gen, g_info.SO

    # tag, generator, matrices, plus set.
    dic = {}
    dic["tag"] = so
    dic["generator"] = gen
    som = np.asarray([symmetry_operation_matrix(i, hexagonal, axial=False) for i in so])
    dic["cartesian"] = som
    if hexagonal:
        dic["fractional"] = convert_to_fractional_hexagonal_matrix(som)
    else:
        dic["fractional"] = som
    dic["det"] = np.asarray([(sp.Matrix(i).det()) for i in som])
    if group_tag.count("^") > 0:
        lattice = g_info.lattice
        dic["cartesian_primitive"] = convert_to_primitive_matrix(lattice, dic["cartesian"])
        dic["fractional_primitive"] = convert_to_primitive_matrix(lattice, dic["fractional"])
        dic["plus_set"] = t_dict[lattice]
    else:
        # product table.
        dic["product"] = point_group_product[group_tag.replace("-c", "")]

    so_mat = Dict(SOMatrixType)
    for s in range(2):
        so_mat[("Q", s)] = [symmetry_operation_matrix_multipole(s, i, hexagonal, axial=False, tesseral=True) for i in so]
        so_mat[("G", s)] = [symmetry_operation_matrix_multipole(s, i, hexagonal, axial=True, tesseral=True) for i in so]

    dic["so_matrix"] = so_mat

    return dic


# ==================================================
def create_m_group_symmetry_operation(g_info):
    dic = {}

    hexagonal = g_info.hexagonal_g
    so = g_info.SO

    # tag, generator, matrices.
    dic = {}
    dic["tag"] = so
    dic["tr_sign"] = np.asarray([-sp.S(1) if i.count("'") > 0 else sp.S(1) for i in so])
    som = np.asarray([symmetry_operation_matrix(i.replace("'", ""), hexagonal, axial=False) for i in so])
    dic["cartesian"] = som
    if hexagonal:
        dic["fractional"] = convert_to_fractional_hexagonal_matrix(som)
    else:
        dic["fractional"] = som
    dic["det"] = np.asarray([sp.Matrix(i).det() for i in som])

    return dic


# ==================================================
def create_group_character(group_tag, so):
    val = character[group_tag]
    val_plus = character_plus[group_tag]

    dic = {}

    # table.
    conjugacy = val[0]
    dic["conjugacy"] = conjugacy
    dim = {i: TagIrrep.parse(i)["dimension"] for i in val[1]}
    dic["dimension"] = dim
    table = {irrep: str_to_sympy(v) for irrep, v in val[1].items()}
    dic["table"] = table
    dic["table_full"] = create_character(table, so["tag"], conjugacy)

    # polar-axial conv.
    dic["polar_axial_conversion"] = val_plus["polar_axial"]

    # sym. and anti-sym. product decomposition.
    dic["symmetric_product"] = val_plus["symmetric_product"]
    dic["anti_symmetric_product"] = val_plus["anti_symmetric_product"]

    return dic


# ==================================================
def create_group_wyckoff(g_info, group_tag, wyckoff_data, sym_op, s1c, s1h):
    # create wyckoff mapping.
    point_group = group_tag.count("^") == 0
    hexagonal = g_info.hexagonal_g
    lattice = "0" if point_group else g_info.lattice
    s = s1h if hexagonal else s1c
    so = sym_op["fractional"]
    if not point_group:
        plus_set = sym_op["plus_set"]
    wyckoff_mapping_group = {}
    for wp, (wp_data, sym) in wyckoff_data[group_tag.replace("-c", "")].items():
        # specific wyckoff sites obtained by inserting (x,y,z)-wyckoff expression with representative site, s.
        pos = wp_data.replace("x", f"({s[0]})").replace("y", f"({s[1]})").replace("z", f"({s[2]})")
        pos = str_to_sympy(pos)
        # create mapping by applying SO to first wyckoff site (pos0).
        if point_group:
            pos_lst = np.vectorize(sp.radsimp)(so @ pos[0])
            mapping = [[] for _ in range(len(pos))]
            for i, p in enumerate(pos):
                for no, pw in enumerate(pos_lst):
                    if (np.vectorize(sp.radsimp)(p - pw) == 0).all():
                        mapping[i].append(no + 1)
        else:
            pos_lst = so @ np.asarray([pos[0][0], pos[0][1], pos[0][2], 1])
            pos_lst = np.vectorize(sp.radsimp)(pos_lst[:, 0:3])
            mapping = [[] for _ in range(len(pos))]
            for i, p in enumerate(pos):
                for no, pw in enumerate(pos_lst):
                    for t in plus_set:
                        # thaks to symbolic calc., mod for shift works fine.
                        if (np.vectorize(sp.radsimp)(np.mod(p - pw - t, 1)) == 0).all():
                            mapping[i].append(no + 1)
            pos = np.mod(np.concatenate([pos + t for t in plus_set]), 1)  # shift by mod works fine.

        if hexagonal:
            cart = np.vectorize(sp.radsimp)(convert_to_cartesian_hexagonal(pos))
        else:
            cart = pos

        if point_group:
            cart = np.vectorize(sp.radsimp)(normalize_vector(cart))
            if hexagonal:
                pos = np.vectorize(sp.radsimp)(convert_to_fractional_hexagonal(cart))
            else:
                pos = cart

        pos_ex = str_to_sympy(wp_data)
        if point_group:
            prim = pos_ex
        else:
            prim = convert_to_primitive(lattice, pos_ex, shift=True)

        if point_group:
            wyckoff_mapping_group[wp] = {
                "symmetry": sym,
                "conventional": pos_ex,
                "reference": pos,
                "mapping": mapping,
            }
        else:
            wyckoff_mapping_group[wp] = {
                "symmetry": sym,
                "conventional": pos_ex,
                "primitive": prim,
                "reference": pos,
                "mapping": mapping,
            }

    return wyckoff_mapping_group


# ==================================================
def create_m_group_wyckoff(g_info, group_tag, wyckoff_data, sym_op, s1c, s1h):
    # create wyckoff mapping.
    point_group = group_tag.split(":")[0] == "MPG"
    hexagonal = g_info.hexagonal_g
    s = s1h if hexagonal else s1c
    so = sym_op["fractional"]
    wyckoff_mapping_group = {}
    for wp, (wp_data, sym) in wyckoff_data[group_tag.split(":")[1]].items():
        # specific wyckoff sites obtained by inserting (x,y,z)-wyckoff expression with representative site, s.
        pos = wp_data.replace("x", f"({s[0]})").replace("y", f"({s[1]})").replace("z", f"({s[2]})")
        pos = str_to_sympy(pos)
        # create mapping by applying SO to first wyckoff site (pos0).
        if point_group:
            pos_lst = np.vectorize(sp.radsimp)(so @ pos[0])
            mapping = [[] for _ in range(len(pos))]
            for i, p in enumerate(pos):
                for no, pw in enumerate(pos_lst):
                    if (np.vectorize(sp.radsimp)(p - pw) == 0).all():
                        mapping[i].append(no + 1)
        else:
            pos_lst = so @ np.asarray([pos[0][0], pos[0][1], pos[0][2], 1])
            pos_lst = np.vectorize(sp.radsimp)(pos_lst[:, 0:3])
            mapping = [[] for _ in range(len(pos))]
            for i, p in enumerate(pos):
                for no, pw in enumerate(pos_lst):
                    # thaks to symbolic calc., mod for shift works fine.
                    if (np.vectorize(sp.radsimp)(np.mod(p - pw, 1)) == 0).all():
                        mapping[i].append(no + 1)
            pos = np.mod(pos, 1)  # shift by mod works fine.

        if hexagonal:
            cart = np.vectorize(sp.radsimp)(convert_to_cartesian_hexagonal(pos))
        else:
            cart = pos

        if point_group:
            cart = np.vectorize(sp.radsimp)(normalize_vector(cart))
            if hexagonal:
                pos = np.vectorize(sp.radsimp)(convert_to_fractional_hexagonal(cart))
            else:
                pos = cart

        pos_ex = str_to_sympy(wp_data)

        wyckoff_mapping_group[wp] = {
            "symmetry": sym,
            "conventional": pos_ex,
            "reference": pos,
            "mapping": mapping,
        }

    return wyckoff_mapping_group


# ==================================================
def create_active_multipole(g_info, info, so):
    x, y, z = sp.symbols("x y z", real=True)
    v = np.array([x, y, z])

    hexagonal = g_info.hexagonal_g
    mp_list = np.array(
        info["response_tensor"]["multipole"]["hexagonal"] if hexagonal else info["response_tensor"]["multipole"]["cubic"]
    )
    bf = np.array([info["harmonics"]["basis_function"][i][0] for i in mp_list])

    det = so["det"]
    tr = so["tr_sign"]

    soc = so["cartesian"].transpose(0, 2, 1)  # transpose for each matrix.
    vt = soc @ v

    tQ = [[sp.simplify(bi.subs({x: vi[0], y: vi[1], z: vi[2]}, simultaneous=True)) for vi in vt] for bi in bf]
    tG = [[sp.simplify(d * bi.subs({x: vi[0], y: vi[1], z: vi[2]}, simultaneous=True)) for vi, d in zip(vt, det)] for bi in bf]
    tT = [[sp.simplify(t * bi.subs({x: vi[0], y: vi[1], z: vi[2]}, simultaneous=True)) for vi, t in zip(vt, tr)] for bi in bf]
    tM = [
        [sp.simplify(t * d * bi.subs({x: vi[0], y: vi[1], z: vi[2]}, simultaneous=True)) for vi, t, d in zip(vt, tr, det)]
        for bi in bf
    ]

    iir_Q = ["Q" + i for i in mp_list[[bool(np.all(np.array(i) == i[0])) for i in tQ]]]
    iir_G = ["G" + i for i in mp_list[[bool(np.all(np.array(i) == i[0])) for i in tG]]]
    iir_T = ["T" + i for i in mp_list[[bool(np.all(np.array(i) == i[0])) for i in tT]]]
    iir_M = ["M" + i for i in mp_list[[bool(np.all(np.array(i) == i[0])) for i in tM]]]

    iir = iir_Q + iir_G + iir_T + iir_M

    return iir


# ==================================================
def create_group():
    info = BinaryManager("info", topdir=BIN_DIR)

    # representative site.
    s1c = info["root_cluster"]["scaled_rep_site_c"]
    s1h = info["root_cluster"]["scaled_rep_site_h"]

    group = BinaryManager(verbose=True, topdir=BIN_DIR)
    group.add_comment(h_group)

    for pg_id in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        group_tag = info["tag"][pg_id]
        print("creating", group_tag, flush=True)
        so = create_group_symmetry_operation(pg_info[pg_id], group_tag)
        group[pg_id] = {
            "info": pg_info[pg_id],
            "symmetry_operation": so,
            "wyckoff": {"site": create_group_wyckoff(pg_info[pg_id], group_tag, wyckoff_site_data, so, s1c, s1h)},
            "character": create_group_character(group_tag, so),
        }

    for sg_id in info["id_set"]["SG"]["all"]:
        group_tag = info["tag"][sg_id]
        print("creating", group_tag, flush=True)
        so = create_group_symmetry_operation(sg_info[sg_id], group_tag)
        group[sg_id] = {
            "info": sg_info[sg_id],
            "symmetry_operation": so,
            "wyckoff": {"site": create_group_wyckoff(sg_info[sg_id], group_tag, wyckoff_site_data, so, s1c, s1h)},
        }

    for mpg_id in info["id_set"]["MPG"]["all"]:
        group_tag = info["tag"][mpg_id]
        print("creating", mpg_id, flush=True)
        g_info = mpg_info[mpg_id]
        so = create_m_group_symmetry_operation(g_info)
        active_mp = create_active_multipole(g_info, info, so)
        group[mpg_id] = {
            "info": mpg_info[mpg_id],
            "symmetry_operation": so,
            "wyckoff": {"site": create_m_group_wyckoff(mpg_info[mpg_id], mpg_id, wyckoff_mpg_data, so, s1c, s1h)},
            "active_multipole": active_mp,
        }

    for msg_id in info["id_set"]["MSG"]["all"]:
        group_tag = info["tag"][msg_id]
        print("creating", msg_id, flush=True)
        so = create_m_group_symmetry_operation(msg_info[msg_id])
        group[msg_id] = {
            "info": msg_info[msg_id],
            "symmetry_operation": so,
            "wyckoff": {"site": create_m_group_wyckoff(msg_info[msg_id], msg_id, wyckoff_msg_data, so, s1c, s1h)},
        }

    group.save_binary("group")


# ==================================================
@timer
def create_group_all():
    create_group()  # basic.
    create_group_mapping()  # add point-group mapping.
    create_group_bond()  # add wyckoff bond info.


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_group_all()
