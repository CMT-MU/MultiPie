"""
Create binary data (cluster SAMB).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp

from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager
from multipie.util.util_crystal import convert_to_cartesian_hexagonal

from others.tools.utils.util_samb import expand_component, gather_component, orthogonalize

# ==================================================
h_cluster_samb = """
* Cluster SAMB for all real PG and SG ("PG"/"SG").
- PG_id/SG_id (str): (dict) SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (Dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_s" (str): (dict) symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_a" (str): (dict) anti-symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "vector" (str): (dict) bond-cluster vector SAMB data.
    - bw_tag (str): (Dict) vector SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,ns,sympy), ndarray(dim,sympy), ndarray(dim,2s+1,sympy)) [SAMB], [cartesian ex.], and [cartesian vector ex.] for each component.
NOTE:
  - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
  - multipoles are sorted in order of [Gamma, l, k, Q/G, n, p].
"""


# ==================================================
def create_mapping_rc(so, op_mapping, wyckoff):
    # so number -> rc site number.
    so_rc_map = [op_mapping[tag[: tag.rfind(":")] if tag.count(":") else tag] for tag in so]

    # symmetry-operation set.
    so_idx_list = sorted(set([tuple(i) for i in sum([v["mapping"] for v in wyckoff.values()], [])]), key=len)

    # symmetry-operation -> rc site.
    rc_idx_list = [(list(lst), [so_rc_map[abs(i) - 1] for i in lst], [-1 if i < 0 else 1 for i in lst]) for lst in so_idx_list]

    return rc_idx_list


# ==================================================
def create_cluster_basis_list(rc_idx_list, rc_harm, sign):
    lst = {}
    for rc_idx in rc_idx_list:
        dic = {}
        for idx, harm in rc_harm.items():
            if sign:
                b = np.sum(harm[0][:, rc_idx[1]] * rc_idx[2], axis=1)
            else:
                b = np.sum(harm[0][:, rc_idx[1]], axis=1)
            dic[idx] = b
        lst[str(rc_idx[0])] = dic
    return lst


# ==================================================
def create_cluster_list(so, op_mapping, wyckoff_site, wyckoff_bond, monopole, dipole):
    monopole_q = monopole.select(X="Q")  # polar only.
    dipole = dipole.sort("Gamma", "l", "k", "n")

    site_rc_idx_list = create_mapping_rc(so, op_mapping, wyckoff_site)
    bond_rc_idx_list = create_mapping_rc(so, op_mapping, wyckoff_bond)

    site_basis = create_cluster_basis_list(site_rc_idx_list, monopole_q, sign=False)
    bond_s_basis = create_cluster_basis_list(bond_rc_idx_list, monopole_q, sign=False)
    bond_a_basis = create_cluster_basis_list(bond_rc_idx_list, monopole_q, sign=True)
    vector_basis = create_cluster_basis_list(bond_rc_idx_list, dipole, sign=False)

    dic = {"site": site_basis, "bond_s": bond_s_basis, "bond_a": bond_a_basis, "vector": vector_basis}

    return dic


# ==================================================
def create_wyckoff_basis(basis, mapping, lst, vector):
    # create basis for given wyckoff.
    wp_basis = {}
    for idx, bs in basis.items():
        b = np.array([lst[str(ms)][idx] for ms in mapping])
        if not (b == 0).all():
            if vector:
                b = b.transpose(1, 0, 2)
            else:
                b = b.transpose(1, 0)
            wp_basis[idx] = (b, *bs[1:])

    ex_basis = expand_component(wp_basis, vector)

    n = len(mapping)
    if vector:
        n *= 3
    o_basis = orthogonalize(ex_basis, n)

    wp_basis = gather_component(o_basis, vector)

    return wp_basis


# ==================================================
def create_anti_sym_basis(bond_wyckoff, bond_a_samb, vec_samb, hexagonal, harmonics_idx, tol=1e-8):
    n = len(bond_wyckoff["conventional"])
    dic_vec = expand_component(vec_samb, vector=True)
    dic_bond_a = expand_component(bond_a_samb, vector=False)

    vec = bond_wyckoff["reference"][:n, 0:3]
    if hexagonal:
        vec = convert_to_cartesian_hexagonal(vec)
    vec = vec.astype(float)

    # vector rep. of bond_a.
    dic = {idx: np.asarray([i * j for i, j in zip(vec, val[0].astype(float))]) for idx, val in dic_bond_a.items()}

    # replace index to that of vector_samb with first independent finite overlap.
    bond_a0 = {}
    rep_idx = set()
    for idx, ba in dic.items():
        replaced = False
        ba = ba.reshape(-1)
        for idxv, val in dic_vec.items():
            if idxv in rep_idx:
                continue
            v = val[0].astype(float).reshape(-1)
            ip = ba @ v
            if idxv not in rep_idx and abs(ip) > tol:
                ex = val[1]
                ex_str = str(sp.expand(ex)).replace(" ", "")
                idx_r = harmonics_idx[(idxv[0][0], ex_str)]
                basis = dic_bond_a[idx][0]
                bond_a0[idx_r] = bond_a0.get(idx_r, []) + [(basis, ex)]
                replaced = True
                rep_idx.add(idxv)
                break
        if not replaced:  # no overlap (should not be reached).
            raise Exception(f"no overlap between bond_a and vector SAMB for {idx}.")

    # renumber.
    bond_a = {}
    for ((X, l, Gamma, mul, p, s, k, x), comp), vals in bond_a0.items():
        if len(vals) == 1:
            bond_a[((X, l, Gamma, mul, -1, s, k, x), comp)] = vals[0]
        else:
            for no, va in enumerate(vals):
                bond_a[((X, l, Gamma, mul, no + 1, s, k, x), comp)] = va

    if len(dic_bond_a.keys()) != len(bond_a.keys()):
        raise Exception("# of basis is different.")

    bond_a_basis = gather_component(bond_a, vector=False)

    return bond_a_basis


# ==================================================
def create_basis_gp(wyckoff_site, wyckoff_bond, monopole, dipole, lst, hexagonal, harmonics_idx):
    monopole_q = monopole.select(X="Q")
    dipole = dipole.sort("Gamma", "l", "k", "n")

    site = {}
    for s_wp in wyckoff_site.keys():
        print("  site =", s_wp)
        site[s_wp] = create_wyckoff_basis(monopole_q, wyckoff_site[s_wp]["mapping"], lst["site"], vector=False)

    bond_s = {}
    bond_a = {}
    vector = {}
    for b_wp in wyckoff_bond.keys():
        bond_wyckoff_mp = wyckoff_bond[b_wp]["mapping"]
        print("  bond_s =", b_wp)
        bond_s[b_wp] = create_wyckoff_basis(monopole_q, bond_wyckoff_mp, lst["bond_s"], vector=False)
        print("  bond_a =", b_wp)
        bond_a_pre = create_wyckoff_basis(monopole_q, bond_wyckoff_mp, lst["bond_a"], vector=False)
        print("  vector =", b_wp)
        vector[b_wp] = create_wyckoff_basis(dipole, bond_wyckoff_mp, lst["vector"], vector=True)
        # replace index to that of vector_samb with first independent finite overlap.
        bond_a[b_wp] = create_anti_sym_basis(wyckoff_bond[b_wp], bond_a_pre, vector[b_wp], hexagonal, harmonics_idx)

    dic = {"site": site, "bond_s": bond_s, "bond_a": bond_a, "vector": vector}

    return dic


# ==================================================
@timer
def create_cluster_samb():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)
    rc = BinaryManager("root_cluster", topdir=BIN_DIR)
    rc_harmonics = BinaryManager("harmonics_root_cluster", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)

    samb = BinaryManager(verbose=True, topdir=BIN_DIR)
    samb.add_comment(h_cluster_samb)

    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        gp = group[no]
        print("creating", group_tag, flush=True)
        so = gp["symmetry_operation"]["tag"]
        hexagonal = gp["info"].hexagonal_g
        parent = "D6h" if hexagonal else "Oh"
        op_mapping = rc["op_mapping"][parent]
        wyckoff_site = gp["wyckoff"]["site"]
        wyckoff_bond = gp["wyckoff"]["bond"]
        harmonics_idx = harmonics[no]["index"]

        monopole = rc_harmonics[no]["monopole"]
        dipole = rc_harmonics[no]["dipole"]

        lst = create_cluster_list(so, op_mapping, wyckoff_site, wyckoff_bond, monopole, dipole)

        samb[no] = create_basis_gp(wyckoff_site, wyckoff_bond, monopole, dipole, lst, hexagonal, harmonics_idx)

    for no in info["id_set"]["SG"]["all"]:
        group_tag = info["tag"][no]
        gp = group[no]
        print("creating", group_tag, flush=True)
        so = gp["symmetry_operation"]["tag"]
        hexagonal = gp["info"].hexagonal_g
        parent = "D6h" if hexagonal else "Oh"
        op_mapping = rc["op_mapping"][parent]
        wyckoff_site = gp["wyckoff"]["site"]
        wyckoff_bond = gp["wyckoff"]["bond"]
        pg_id = gp["info"].PG
        harmonics_idx = harmonics[pg_id]["index"]

        monopole = rc_harmonics[pg_id]["monopole"]
        dipole = rc_harmonics[pg_id]["dipole"]

        lst = create_cluster_list(so, op_mapping, wyckoff_site, wyckoff_bond, monopole, dipole)

        samb[no] = create_basis_gp(wyckoff_site, wyckoff_bond, monopole, dipole, lst, hexagonal, harmonics_idx)

    samb.save_binary("cluster_samb")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_cluster_samb()
