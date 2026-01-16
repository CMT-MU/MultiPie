"""
Check created data.

- compatibility relation.
- LC of Clm, Slm and cartesian expression (Oh, D6h harmonics).
- consistency of polar and axial irrep. (PG harmonics).
- orthonormality and completeness (spherical, PG atomic SAMB).
- consistency of tag between PG atomic SAMB and PG harmonics.
- symmetry operation between PG and SG.
- symmetry operation mapping between parent and child PG.
- find wyckoff site and bond (PG and SG).
- consistency between character and trace of rep. matrix.
- consistency between polar and axial rep. matrix.

Todo: Check the following points.
- compatibility (polar/axial): (Oh,D6h) subgroup l: Gamma,gamma -> (sgn)Gamma'gamma' ?
- atomic MP irrep = product decomp of orbital irrep ?
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import numpy as np
import sympy as sp
from tqdm import tqdm

from multipie.util.util_binary import BinaryManager
from multipie.util.util_wyckoff import find_wyckoff_site, find_wyckoff_bond, shift_site, shift_bond, find_vector

from others.tools.data.data_compatibility_relation import compatibility_relation
from others.tools.utils.util_spherical_harmonics import Olm

from others.tools.create_binary.create_data_harmonics import create_top_harmonics
from others.tools.create_binary.create_data_rep_matrix import create_rep_matrix_data


# ==================================================
def bar(n, desc):
    return tqdm(total=n, desc=desc, unit="step", bar_format="{bar:20} {n_fmt}/{total_fmt} {desc}{postfix}")


# ==================================================
# check atomic multipole.
# ==================================================
def _check_amp_orthonormality(vec, same_block=True):
    n = vec.shape[0]
    if same_block:
        vij = vec @ vec.conjugate().T
    else:
        vij = 2 * np.vectorize(sp.re)((vec @ vec.conjugate().T))
    im = np.eye(n, dtype=object) * sp.S(1)
    return (vij == im).all()


# ==================================================
def _check_amp_completeness(vec, same_block):
    n = vec.shape[1]
    vij = vec.conjugate().T @ vec
    im = np.eye(n, dtype=object) * sp.S(1)
    if same_block:
        return (vij == im).all()
    uij = vec.T @ vec
    zm = np.full((n, n), sp.S(0), dtype=object)
    return (vij == im).all() and (uij == zm).all()


# ==================================================
def check_amp_orthonormality_pg(pg_tag, amg, no):
    # check orthonormality of atomic multipole (# of AMP as well).

    if pg_tag.count("^") > 0 or pg_tag.count("-c") > 0:  # space group or complex point group.
        return []

    pbar = bar(30, f"{pg_tag}: atomic multipole orthonormality")
    err = []
    for b_type in ["jml", "lgs", "lg"]:
        n_mp = 0
        for (bra, ket), data in amg[no][b_type].items():
            pbar.set_postfix_str(f"type={b_type}, bra-ket:({bra},{ket})")
            pbar.update(1)
            vec = np.vstack([mat.reshape(mat.shape[0], -1) for mat, _ in data.values()])
            n_mp += len(vec)
            if not _check_amp_orthonormality(vec, bra == ket):
                err.append({"chk": "amp_ortho", "tag": pg_tag, "b_type": b_type, "bra": bra, "ket": ket, "vec": vec})
        if b_type == "lg":
            if n_mp != 256:
                err.append({"chk": "amp_ortho", "tag": pg_tag, "b_type": b_type, "n_mp": n_mp})
        else:
            if n_mp != 1024:
                err.append({"chk": "amp_ortho", "tag": pg_tag, "b_type": b_type, "n_mp": n_mp})
    pbar.close()

    return err


# ==================================================
def check_amp_completeness_pg(pg_tag, amg, no):
    # check completeness of atomic multipole.

    if pg_tag.count("^") > 0 or pg_tag.count("-c") > 0:  # space group or complex point group.
        return []

    pbar = bar(30, f"{pg_tag}: atomic multipole completeness")
    err = []
    for b_type in ["jml", "lgs", "lg"]:
        for (bra, ket), data in amg[no][b_type].items():
            pbar.set_postfix_str(f"type={b_type}, bra-ket:({bra},{ket})")
            pbar.update(1)
            vec = np.vstack([mat.reshape(mat.shape[0], -1) for mat, _ in data.values()])
            if not _check_amp_completeness(vec, bra == ket):
                err.append({"chk": "amp_comp", "tag": pg_tag, "b_type": b_type, "bra": bra, "ket": ket, "vec": vec})
    pbar.close()

    return err


# ==================================================
def check_amp_harmonics_consistency_pg(pg_tag, amg, harmonics, no):
    # check consistency with harmonics label of atomic multipole.

    if pg_tag.count("^") > 0 or pg_tag.count("-c") > 0:  # space group or complex point group.
        return []

    head = {"Q": "Q", "G": "G", "T": "Q", "M": "G"}

    pbar = bar(30, f"{pg_tag}: atomic multipole harmonics label")
    err = []
    for b_type in ["jml", "lgs", "lg"]:
        for (bra, ket), samb in amg[no][b_type].items():
            pbar.set_postfix_str(f"type={b_type}, bra-ket:({bra},{ket})")
            pbar.update(1)
            for (X, l, Gamma, n, p, s, k, _), (mat, ex) in samb.items():
                c, U, t = harmonics[no]["harmonics"][(head[X], l, Gamma, n, -1, 0, 0, "q")]
                if not (ex == c).all():
                    err.append(
                        {
                            "chk": "amp_label",
                            "tag": pg_tag,
                            "X": X,
                            "l": l,
                            "Gamma": Gamma,
                            "n": n,
                            "p": p,
                            "s": s,
                            "k": k,
                            "ex": ex,
                            "cart": c,
                        }
                    )
    pbar.close()

    return err


# ==================================================
# check group.
# ==================================================
def check_symmetry_operation_sg(sg_tag, group, no):
    # check symmetry operation between space and point group.

    if sg_tag.count("^") == 0:  # point group.
        return []

    sg = group[no]
    pg_no = sg["info"].PG
    pg = group[pg_no]

    so_sg = sg["symmetry_operation"]["fractional"][:, 0:3, 0:3]
    so_pg = pg["symmetry_operation"]["fractional"]

    pbar = bar(1, "symmetry operation consistency")
    err = []
    if (so_sg == so_pg).all():
        pbar.set_postfix_str("ok")
        pbar.update(1)
    else:
        err.append({"chk": "group_so", "tag": sg_tag})

    return err


# ==================================================
def check_so_mapping_pg(pg_tag, group, no):
    # check symmetry operation mapping between parent and child point group.

    if pg_tag.count("^") > 0:  # space group.
        return []

    pg = group[no]
    mp = pg["symmetry_operation"]["mapping"]
    mno = "PG:27" if pg["info"].hexagonal_g else "PG:32"
    mtag = "D6h" if pg["info"].hexagonal_g else "Oh"
    prod = group[mno]["symmetry_operation"]["product"]
    mp_so = mp["1"]

    pbar = bar(len(mp), f"symmetry operation mapping to {mtag}")
    err = []
    for v in mp.values():
        pbar.set_postfix_str(f"SO={v[0]}")
        pbar.update(1)
        lst = [prod[(v[0], i)] for i in mp_so]
        if v != lst:
            err.append({"chk": "group_map", "tag": pg_tag, "v": v})
    pbar.close()

    return err


# ==================================================
def check_find_wyckoff_site(tag, group, no):
    # check find Wyckoff site.

    err = []
    g = group[no]
    wp_data = g["wyckoff"]["site"]
    pset = g["symmetry_operation"].get("plus_set", None)
    pbar = bar(len(wp_data), f"{tag}: Wyckoff position (site)")
    for wp, data in wp_data.items():
        pbar.set_postfix_str(f"{wp}")
        pbar.update(1)
        site = data["reference"]
        site = site[len(site) // 2].astype(float)  # rsite for check.
        if pset is not None:
            site = shift_site(site)
        s_wp, all_site = find_wyckoff_site(g, site)
        if s_wp is None:
            err.append({"chk": "wyckoff_site", "tag": tag, "wp": wp})
            continue
        if wp != s_wp:
            err.append({"chk": "wyckoff_site", "tag": tag, "wp": wp, "s_wp": s_wp})
        idx = find_vector(site, all_site)
        if idx is None:
            err.append({"chk": "wyckoff_site", "tag": tag, "wp": wp, "s_wp": s_wp, "idx": idx})
    pbar.close()

    return err


# ==================================================
def check_find_wyckoff_bond(tag, group, no):
    # check find Wyckoff bond.

    err = []
    g = group[no]
    wp_data = g["wyckoff"]["bond"]
    pset = g["symmetry_operation"].get("plus_set", None)
    pbar = bar(len(wp_data), f"{tag}: Wyckoff position (bond)")
    for wp, data in wp_data.items():
        pbar.set_postfix_str(f"{wp}")
        pbar.update(1)
        bond = data["reference"]
        bond = bond[len(bond) // 2].astype(float)  # bond for check.
        if pset is not None:
            bond = shift_bond(bond)
        b_wp, all_bond = find_wyckoff_bond(g, bond)
        if b_wp is None:
            err.append({"chk": "wyckoff_bond", "tag": tag, "wp": wp})
            continue
        if wp != b_wp:
            err.append({"chk": "wyckoff_bond", "tag": tag, "wp": wp, "b_wp": b_wp})
        idx = find_vector(bond, all_bond)
        if idx is None:
            err.append({"chk": "wyckoff_bond", "tag": tag, "wp": wp, "b_wp": b_wp, "idx": idx})
    pbar.close()

    return err


# ==================================================
def check_character(pg_tag, group, harmonics, so_matrix, no):
    # check if character equals to trace of representation matrix.

    if pg_tag.count("^") > 0:  # space group.
        return []

    pg = group[no]
    ch_all = pg["character"]["table_full"]

    harm = harmonics[no]
    so = so_matrix["hexagonal"] if pg["info"].hexagonal_g else so_matrix["cubic"]

    p = create_rep_matrix_data(pg, harm, so, axial=False)["matrix"]
    a = create_rep_matrix_data(pg, harm, so, axial=True)["matrix"]

    pbar = bar(len(p), f"{pg_tag}: character trace")
    err = []
    for irrep in p.keys():
        pbar.set_postfix_str(f"irrep={irrep}")
        pbar.update(1)
        ch = ch_all[irrep]

        pt = []
        for idx, lst in p[irrep].items():
            pt.append((idx, np.trace(lst, axis1=1, axis2=2)))
        at = []
        for idx, lst in a[irrep].items():
            at.append((idx, np.trace(lst, axis1=1, axis2=2)))

        for idx, t in pt:
            if not ((t - ch) == 0).all():
                err.append({"chk": "char_char", "tag": pg_tag, "polar": (idx, t.tolist())})
        for idx, t in at:
            if not ((t - ch) == 0).all():
                err.append({"chk": "char_char", "tag": pg_tag, "axial": (idx, t.tolist())})
    pbar.close()

    return err


# ==================================================
def check_rep_matrix(pg_tag, group, harmonics, so_matrix, no):
    # check representation matrix (polar/axial); given = U_Q^dagger g^(Q,l) U_Q = U_G^dagger g^(G,l) U_G ?

    if pg_tag.count("^") > 0:  # space group.
        return []

    pg = group[no]
    dim = pg["character"]["dimension"]

    harm = harmonics[no]
    so = so_matrix["hexagonal"] if pg["info"].hexagonal_g else so_matrix["cubic"]

    p = create_rep_matrix_data(pg, harm, so, axial=False)["matrix"]
    a = create_rep_matrix_data(pg, harm, so, axial=True)["matrix"]

    pbar = bar(len(p), f"{pg_tag}: character representation matrix")
    err = []
    for irrep in p.keys():
        pbar.set_postfix_str(f"irrep={irrep}")
        pbar.update(1)
        if dim[irrep] == 1:
            continue
        pm = []
        for idx, lst in p[irrep].items():
            pm.append((idx, lst))
        am = []
        for idx, lst in a[irrep].items():
            am.append((idx, lst))

        _, m0 = pm[0]
        for idx, m in pm:
            if not ((m - m0) == 0).all():
                err.append({"chk": "char_repmat", "tag": pg_tag, "irrep": irrep, "polar": idx})
        for idx, m in am:
            if not ((m - m0) == 0).all():
                err.append({"chk": "char_repmat", "tag": pg_tag, "irrep": irrep, "axial": idx})
    pbar.close()

    return err


# ==================================================
# check info.
# ==================================================
def check_top_harmonics():
    # check Oh, D6h harmonics: compatibility between cartesian expression and linear combination of Clm, Slm.

    info = create_top_harmonics()

    err = []
    rv = sp.symbols("x y z", real=True)
    rvs = {"x": rv[0], "y": rv[1], "z": rv[2]}

    pbar = bar(len(info["data"]), "top harmonics")
    for idx, (t_ex, c_ex, umat) in info["data"].items():
        top, rank, irrep, mul, comp = idx

        pbar.set_postfix_str(f"top={top}, l={rank}, irrep={irrep}, n={mul}, comp={comp}")
        pbar.update(1)

        cs = {}
        for m in range(1, rank + 1):
            olm = Olm(rank, m, rv=rv)
            s = f"C{m}"
            cs[s] = sp.sqrt(2) * (-1) ** m * sp.re(olm)
            s = f"S{m}"
            cs[s] = sp.sqrt(2) * (-1) ** m * sp.im(olm)
        cs["C0"] = Olm(rank, 0, rv=rv)

        olm = [Olm(rank, m, rv=rv) for m in range(rank, -rank - 1, -1)]
        u_ex_e = sum([oi * ui for oi, ui in zip(olm, umat)], sp.S(0))

        t_ex_e = t_ex.subs(cs)
        c_ex_e = c_ex.subs(rvs)

        diff = sp.simplify(t_ex_e - c_ex_e)
        if diff != 0:
            err.append({"chk": "top_harm", "tag": "tesseral != cartesian", "idx": idx, "diff": diff})

        diff = sp.simplify(t_ex_e - u_ex_e)
        if diff != 0:
            err.append({"chk": "top_harm", "tag": "tesseral != u-matrix", "idx": idx, "diff": diff})
    pbar.close()

    return err


# ==================================================
def check_polar_axial_compatibility():
    # check consistency between polar and axial irrep.

    sconv = {1: "", -1: "-"}
    cconv = {-1: "", 0: "(x)", 1: "(y)", 2: "(z)"}

    pbar = bar(len(compatibility_relation), "compatibility relation")
    err = []
    for tag in ["Oh", "D6h"]:
        dp = {}
        da = {}
        for pg, (top, dic) in compatibility_relation.items():
            if top != tag:
                continue
            pbar.set_postfix_str(f"group={pg}")
            pbar.update(1)
            irrep_p = [sconv[vp] + ip + cconv[cp] for ((ip, cp, vp), (ia, ca, va)) in dic.values()]
            irrep_a = [sconv[va] + ia + cconv[ca] for ((ip, cp, vp), (ia, ca, va)) in dic.values()]
            dp[pg] = irrep_p
            da[pg] = irrep_a

        ct = {}
        for pg1, tbl1 in dp.items():
            for pg2, tbl2 in dp.items():
                ct[(pg1, pg2)] = [(i, j) for i, j in zip(tbl1, tbl2)]

        ct2 = {}
        for pg1, tbl1 in da.items():
            for pg2, tbl2 in da.items():
                ct2[(pg1, pg2)] = [(i, j) for i, j in zip(tbl1, tbl2)]

        for p, v1, v2 in zip(ct.keys(), ct.values(), ct2.values()):
            for i in v2:
                if i not in v1:
                    err.append({"chk": "pa_comp", "axial": i, "polar_list": p})
    pbar.close()

    return err


# ================================================== main
if __name__ == "__main__":
    print("- read binary file.")
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)
    so_matrix = BinaryManager("symmetry_operation_matrix", topdir=BIN_DIR)
    amg = BinaryManager("atomic_multipole_group", topdir=BIN_DIR)

    print("- start check.")

    err = []
    for gt in ["PG", "SG"]:
        if gt == "PG":
            gset = info["id_set"][gt]["all"] + info["id_set"][gt]["complex"]
        else:
            gset = info["id_set"][gt]["all"]
        for no in gset:
            tag = info["tag"][no]
            # check group.
            err += check_symmetry_operation_sg(tag, group, no)
            err += check_so_mapping_pg(tag, group, no)
            err += check_find_wyckoff_site(tag, group, no)
            err += check_find_wyckoff_bond(tag, group, no)
            err += check_character(tag, group, harmonics, so_matrix, no)
            err += check_rep_matrix(tag, group, harmonics, so_matrix, no)
            # check atomic multipole.
            err += check_amp_orthonormality_pg(tag, amg, no)
            err += check_amp_completeness_pg(tag, amg, no)
            err += check_amp_harmonics_consistency_pg(tag, amg, harmonics, no)

    # overall info.
    err += check_top_harmonics()
    err += check_polar_axial_compatibility()

    print(f"- end check.")
    if err:
        print(f"{len(err)} errors are found.")
        with open("check_error.txt", "w", encoding="utf-8") as f:
            f.write(f"{len(err)} errors are found.\n")
            for i in err:
                f.write(str(i) + "\n")
    else:
        print("All check is ok.")
