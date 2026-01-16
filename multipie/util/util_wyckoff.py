"""
For Wyckoff related.
"""

import numpy as np
from itertools import product

from multipie.util.util import str_to_sympy, replace
from multipie.util.util_crystal import shift_site, shift_bond, convert_to_primitive, TOL_SAME_SITE, DIGIT


# ==================================================
def find_vector(target, vec, tol=TOL_SAME_SITE):
    """
    Find vector (lowest index only).

    Args:
        target (array-like): vector.
        vec (array-like): list of vectors.
        tol (float, optional): tolerance for same vector.

    Returns:
        - (int) -- index (when not found, return None).
    """
    target = np.asarray(target, dtype=float)
    vec = np.asarray(vec, dtype=float)

    for no, v in enumerate(vec):
        if np.allclose(v, target, atol=tol):
            return no
    return None


# ==================================================
def create_mapping(all_vec, u_vec, tol=TOL_SAME_SITE):
    """
    Create mapping for vector list.

    Args:
        all_vec (ndarray): vector list.
        u_vec (ndarray): unique vector list.
        tol (float, optional): tolerance for same vector.

    Returns:
        - (list) -- mapping.
    """
    mapping = [[] for _ in range(len(u_vec))]
    for no, v in enumerate(all_vec):
        idx = find_vector(v, u_vec, tol)
        if idx is None:
            v[0:3] = -v[0:3]
            idx = find_vector(v, u_vec, tol)
            assert idx is not None
            mapping[idx].append(-(no + 1))
        else:
            mapping[idx].append(no + 1)

    return mapping


# ==================================================
def remove_duplicate_vector(vec, directional=True):
    """
    Remove duplicate vector.

    Args:
        vec (array-like): list of vectors.
        directional (bool, optional): directional vector ?

    Returns:
        - (ndarray) -- vectors.

    Note:
        - each vector is assumed to be site(3d) or bond(6d, vector+center).
    """
    vec = np.asarray(vec, dtype=float)

    vecs = []
    if directional:
        for v in vec:
            if not any(np.allclose(v, uniq_vec, atol=TOL_SAME_SITE) for uniq_vec in vecs):
                vecs.append(v)
    else:
        for v in vec:
            nv = v.copy()
            nv[0:3] = -v[0:3]
            if not any(np.allclose(v, uniq_vec, atol=TOL_SAME_SITE) for uniq_vec in vecs) and not any(
                np.allclose(nv, uniq_vec, atol=TOL_SAME_SITE) for uniq_vec in vecs
            ):
                vecs.append(v)
    return np.asarray(vecs)


# ==================================================
def find_xyz(pos, site, var=["x", "y", "z"]):
    """
    Find match in pos for given site, and obtain [x,y,z].

    Args:
        pos (ndarray): set of vectors in terms of var.
        site (ndarray): site to match.
        var (list, optional): variable.

    Returns:
        - (dict) -- (x,y,z) value (return None when not found).
    """
    sc = np.zeros(3)
    sxyz = np.eye(3)

    site = np.asarray(site, dtype=float)
    pos = np.asarray(pos, dtype=object)

    sub_b = dict(zip(var, sc))
    b = np.asarray([p.subs(sub_b) for p in pos], dtype=float)
    sub_A = [dict(zip(var, si)) for si in sxyz]
    A = np.array([[p.subs(subs) for subs in sub_A] for p in pos - b], dtype=float)

    solution = np.linalg.lstsq(A, site - b, rcond=None)[0].tolist()
    solution = dict(zip(var, solution))

    sol_site = replace(pos, solution).astype(float)
    if not np.allclose(sol_site, site, atol=TOL_SAME_SITE):
        return None

    return solution


# ==================================================
def create_equivalent_site(so_p, s0_p, shift):
    """
    Create equivalent site for given site, s0.

    Args:
        so_p (ndarray): symmetry operation (primitive) (PG:3x3) or (SG:4x4).
        s0_p (ndarray): representative site (primitive) (float).
        shift (bool): shift site for SG ?

    Returns:
        - (ndarray) -- set of equivalent site (primitive).
    """
    so_p = so_p.astype(float)

    if shift:
        site_p = (so_p @ np.pad(s0_p, (0, 1), constant_values=1))[:, 0:3]
        site_p = shift_site(site_p)
    else:
        site_p = so_p[:, 0:3, 0:3] @ s0_p

    site_p = remove_duplicate_vector(site_p)

    return site_p


# ==================================================
def create_equivalent_bond(so, so_p, b0_p, shift):
    """
    Create equivalent bond for given bond, v0@s0.

    Args:
        so (ndarray): symmetry operation (conventional) (PG:3x3) or (SG:4x4).
        so_p (ndarray): symmetry operation (primitive) (PG:3x3) or (SG:4x4).
        b0_p (ndarray): representative bond (primitive), (vector+center) (float).
        shift (bool): shift center for SG ?

    Returns:
        - (ndarray) -- set of equivalent bond (nondirectional, primitive) (vector+center).
    """
    so = so.astype(float)
    so_p = so_p.astype(float)
    v0, s0_p = b0_p[0:3], b0_p[3:6]

    vector = so[:, 0:3, 0:3] @ v0

    if shift:
        center_p = (so_p @ np.pad(s0_p, (0, 1), constant_values=1))[:, 0:3]
        center_p = shift_site(center_p)
    else:
        center_p = so_p[:, 0:3, 0:3] @ s0_p

    bond_p = np.hstack((vector[:, None], center_p[:, None])).reshape(-1, 6)
    bond_p = remove_duplicate_vector(bond_p, directional=False)

    return bond_p


# ==================================================
def wyckoff_site_for_search(site_ex):
    """
    Wyckoff site for search with surrounding cell (no plus set).

    Args:
        site_ex (ndarray): set of Wyckoff site in terms of (x,y,z).

    Returns:
        - (ndarray) -- set of Wyckoff site for search.
    """
    shift = np.array(list(product([-1, 0, 1], repeat=3)), dtype=object)
    pos = np.concatenate([site_ex + i for i in shift])
    return pos


# ==================================================
def wyckoff_bond_for_search(bond_ex):
    """
    Wyckoff bond for search with surrounding cell (no plus set).

    Args:
        bond_ex (ndarray): set of Wyckoff bond in terms of (X,Y,Z,x,y,z).

    Returns:
        - (ndarray) -- set of Wyckoff bond vector for search.
        - (ndarray) -- set of Wyckoff bond center for search.
    """
    vec, ctr = bond_ex[:, 0:3], bond_ex[:, 3:6]

    ctr = wyckoff_site_for_search(ctr)
    n = len(ctr) // len(bond_ex)
    vec = np.tile(vec, (n, 1))

    return vec, ctr


# ==================================================
def _evaluate_wyckoff_site(wyckoff_site, xyz, pset=None):
    """
    Evaluate Wyckoff site for given site for x, y, z (with plus set).

    Args:
        wyckoff_site (ndarray): set of Wyckoff site in terms of (x,y,z).
        xyz (dict): x, y, z value dict.
        pset (ndarray, optional): plus set, which must be given for SG.

    Returns:
        - (ndarray) -- set of Wyckoff site (plus_set0, plus_set1, ...).
    """
    all_site = replace(wyckoff_site, xyz).astype(float)
    if pset is not None:
        pset = pset.astype(float)
        all_site = np.concatenate([all_site + i for i in pset])
        all_site = shift_site(all_site)

    return all_site


# ==================================================
def evaluate_wyckoff_site(group, s_wp, xyz):
    """
    Evaluate Wyckoff site for given site for x, y, z (with plus set).

    Args:
        group (dict): group dict.
        s_wp (str): Wyckoff site.
        xyz (dict): x, y, z value dict.

    Returns:
        - (ndarray) -- set of Wyckoff site (plus_set0, plus_set1, ...).
    """
    pset = group["symmetry_operation"].get("plus_set", None)
    site_wyckoff = group["wyckoff"]["site"][s_wp]
    site = _evaluate_wyckoff_site(site_wyckoff["conventional"], xyz, pset)

    return site


# ==================================================
def _evaluate_wyckoff_bond(wyckoff_bond, XYZxyz, pset=None):
    """
    Evaluate Wyckoff bond for given bond for X, Y, Z, x, y, z (with plus set).

    Args:
        wyckoff_bond (ndarray): set of Wyckoff bond in terms of (X,Y,Z)@(x,y,z).
        XYZxyz (dict): X, Y, Z, x, y, z value dict.
        pset (ndarray, optional): plus set, which must be given for SG.

    Returns:
        - (ndarray) -- set of Wyckoff bond (vector+center) (plus_set0, plus_set1, ...).
    """
    all_bond = replace(wyckoff_bond, XYZxyz).astype(float)
    if pset is not None:
        pset = pset.astype(float)
        vec, ctr = all_bond[:, 0:3], all_bond[:, 3:6]
        ctr = np.concatenate([ctr + i for i in pset])
        ctr = shift_site(ctr)
        vec = np.tile(vec, (len(pset), 1))
        all_bond = np.hstack((vec[:, None], ctr[:, None])).reshape(-1, 6)

    return all_bond


# ==================================================
def evaluate_wyckoff_bond(group, b_wp, XYZxyz):
    """
    Evaluate Wyckoff bond for given bond for X, Y, Z, x, y, z (with plus set).

    Args:
        group (dict): group dict.
        b_wp (str): Wyckoff bond.
        XYZxyz (dict): X, Y, Z, x, y, z value dict.

    Returns:
        - (ndarray) -- set of Wyckoff bond (vector+center) (plus_set0, plus_set1, ...).
    """
    pset = group["symmetry_operation"].get("plus_set", None)
    bond_wyckoff = group["wyckoff"]["bond"][b_wp]
    bond = _evaluate_wyckoff_bond(bond_wyckoff["conventional"], XYZxyz, pset)

    return bond


# ==================================================
def _find_wyckoff_site_pg(so, wyckoff_site, site):
    """
    Find Wyckoff site for point group.

    Args:
        so (ndarray): symmetry operation (3x3).
        wyckoff_site (dict): Wyckoff site dict.
        site (ndarray): site to find.

    Returns:
        - (str) -- Wyckoff position (return None when not found).
        - (dict) -- x, y, z value dict.
        - (ndarray) -- first Wyckoff site (return None when not found).
    """

    def solve():
        for wp in wp_list:
            wp_pos = wyckoff_site[wp]["conventional"]
            for p in wp_pos:
                sol = find_xyz(p, site)
                if sol is None:
                    continue
                return wp, sol
        return None, None

    site_list = create_equivalent_site(so, site, shift=False)
    wp_list = [wp for wp in wyckoff_site.keys() if int(wp[:-1]) == len(site_list)]
    s_wp, sol = solve()
    if s_wp is None:
        return None, None, None
    sol = {k: round(v, DIGIT) for k, v in sol.items()}

    wyckoff_site_wp = wyckoff_site[s_wp]["conventional"][0]
    s0 = replace(wyckoff_site_wp, sol).astype(float).round(DIGIT)

    return s_wp, sol, s0


# ==================================================
def _find_wyckoff_site_sg(so_p, wyckoff_site, site, lattice, pset):
    """
    Find Wyckoff site for space group.

    Args:
        so_p (ndarray): symmetry operation (4x4, primitive).
        wyckoff_site (dict): Wyckoff site dict.
        site (ndarray): site to find.
        lattice (str): lattice.
        pset (ndarray): plus set.

    Returns:
        - (str) -- Wyckoff position (return None when not found).
        - (dict) -- x, y, z value dict.
        - (ndarray) -- first Wyckoff site (return None when not found).
    """

    def solve():
        for wp in wp_list:
            wp_pos = wyckoff_site_for_search(wyckoff_site[wp]["primitive"])
            for p in wp_pos:
                sol = find_xyz(p, site_p)
                if sol is None:
                    continue
                return wp, sol
        return None, None

    site = shift_site(site).astype(float)
    site_p = convert_to_primitive(lattice, site, shift=True).astype(float)

    site_p_list = create_equivalent_site(so_p, site_p, shift=True)
    wp_list = [wp for wp in wyckoff_site.keys() if int(wp[:-1]) == len(pset) * len(site_p_list)]
    s_wp, sol = solve()  # (x,y,z) in sol are in conventional coordinate.
    if s_wp is None:
        return None, None, None
    sol = {k: round(v, DIGIT) for k, v in sol.items()}

    wyckoff_site_wp = wyckoff_site[s_wp]["conventional"][0]
    s0 = shift_site(replace(wyckoff_site_wp, sol)).astype(float)

    return s_wp, sol, s0


# ==================================================
def _find_wyckoff_site_msg(so, wyckoff_site, site):
    """
    Find Wyckoff site for magnetic space group.

    Args:
        so (ndarray): symmetry operation (4x4).
        wyckoff_site (dict): Wyckoff site dict.
        site (ndarray): site to find.

    Returns:
        - (str) -- Wyckoff position (return None when not found).
        - (dict) -- x, y, z value dict.
        - (ndarray) -- first Wyckoff site (return None when not found).
    """

    def solve():
        for wp in wp_list:
            wp_pos = wyckoff_site_for_search(wyckoff_site[wp]["conventional"])
            for p in wp_pos:
                sol = find_xyz(p, site)
                if sol is None:
                    continue
                return wp, sol
        return None, None

    site = shift_site(site).astype(float)

    site_list = create_equivalent_site(so, site, shift=True)
    wp_list = [wp for wp in wyckoff_site.keys() if int(wp[:-1]) == len(site_list)]
    s_wp, sol = solve()  # (x,y,z) in sol are in conventional coordinate.
    if s_wp is None:
        return None, None, None
    sol = {k: round(v, DIGIT) for k, v in sol.items()}

    wyckoff_site_wp = wyckoff_site[s_wp]["conventional"][0]
    s0 = shift_site(replace(wyckoff_site_wp, sol)).astype(float)

    return s_wp, sol, s0


# ==================================================
def find_wyckoff_site(group, site, msg=False):
    """
    Find Wyckoff site.

    Args:
        group (dict): group dict.
        site (ndarray or str): site to find.
        msg (bool, optional): MSG ?

    Returns:
        - (str) -- Wyckoff position (return None when not found).
        - (ndarray) -- set of Wyckoff site (plus_set0, plus_set1, ...).

    Note:
        - site string is "[x,y,z]".
    """
    if type(site) == str:
        site = str_to_sympy(site)
    wyckoff_site = group["wyckoff"]["site"]

    if group["info"].lattice != "0":
        if msg:
            so = group["symmetry_operation"]["fractional"]
            s_wp, sol, s0 = _find_wyckoff_site_msg(so, wyckoff_site, site)
            if s_wp is None:
                return None, None
            sites = _evaluate_wyckoff_site(wyckoff_site[s_wp]["conventional"], sol)
            return s_wp, sites
        else:
            so_p = group["symmetry_operation"]["fractional_primitive"]
            pset = group["symmetry_operation"]["plus_set"]
            lattice = group["info"].lattice
            s_wp, sol, s0 = _find_wyckoff_site_sg(so_p, wyckoff_site, site, lattice, pset)
            if s_wp is None:
                return None, None
            sites = _evaluate_wyckoff_site(wyckoff_site[s_wp]["conventional"], sol, pset)
            return s_wp, sites
    else:
        so = group["symmetry_operation"]["fractional"]
        s_wp, sol, s0 = _find_wyckoff_site_pg(so, wyckoff_site, site)
        if s_wp is None:
            return None, None
        sites = _evaluate_wyckoff_site(wyckoff_site[s_wp]["conventional"], sol)
        return s_wp, sites


# ==================================================
def _find_wyckoff_bond_pg(so, wyckoff_site, wyckoff_bond, bond):
    """
    Find Wyckoff bond for point group.

    Args:
        so (ndarray): symmetry operation (3x3).
        wyckoff_site (dict): Wyckoff site dict.
        wyckoff_bond (dict): Wyckoff bond dict.
        bond (ndarray): bond (vector+center) to find.

    Returns:
        - (str) -- Wyckoff bond position (return None when not found).
        - (dict) -- X, Y, Z, x, y, z value dict.
        - (ndarray) -- first Wyckoff bond (vector+center) (return None when not found).
    """

    def solve():
        for wp in wp_list:
            bond = wyckoff_bond[wp]["conventional"]
            vec, ctr = bond[:, 0:3], bond[:, 3:6]
            for v, c in zip(vec, ctr):
                v_sol = find_xyz(v, vector, ["X", "Y", "Z"])
                if v_sol is None:
                    continue
                c_sol = find_xyz(c, center)
                if c_sol is None:
                    continue
                return wp, v_sol | c_sol
        return None, None

    vector, center = bond[0:3], bond[3:6]
    c_wp = _find_wyckoff_site_pg(so, wyckoff_site, center)[0]
    if c_wp is None:
        return None, None, None

    bond_list = create_equivalent_bond(so, so, bond, shift=False)
    wp_list = [wp for wp in wyckoff_site[c_wp]["bond"] if int(wp.split("@")[0][:-1]) == len(bond_list)]
    b_wp, sol = solve()
    if b_wp is None:
        return None, None, None, None
    sol = {k: round(v, DIGIT) for k, v in sol.items()}

    wyckoff_bond_wp = wyckoff_bond[b_wp]["conventional"][0]
    b0 = replace(wyckoff_bond_wp, sol).astype(float).round(DIGIT)

    return b_wp, sol, b0


# ==================================================
def _find_wyckoff_bond_sg(so, so_p, wyckoff_site, wyckoff_bond, bond, lattice, pset):
    """
    Find Wyckoff bond for space group.

    Args:
        so (ndarray): symmetry operation (4x4, conventional).
        so_p (ndarray): symmetry operation (4x4, primitive).
        wyckoff_site (dict): Wyckoff site dict.
        wyckoff_bond (dict): Wyckoff bond dict.
        bond (ndarray): bond (vector+center) to find.
        lattice (str): lattice.
        pset (ndarray): plus set.

    Returns:
        - (str) -- Wyckoff bond position (return None when not found).
        - (dict) -- X, Y, Z, x, y, z value dict.
        - (ndarray) -- first Wyckoff bond (vector+center) (return None when not found).
    """

    def solve():
        for wp in wp_list:
            vec, ctr = wyckoff_bond_for_search(wyckoff_bond[wp]["primitive"])
            for v, c in zip(vec, ctr):
                v_sol = find_xyz(v, vector, ["X", "Y", "Z"])
                if v_sol is None:
                    continue
                c_sol = find_xyz(c, center_p)
                if c_sol is None:
                    continue
                return wp, v_sol | c_sol
        return None, None

    vector, center = bond[0:3], bond[3:6]
    vector = vector.astype(float)
    center = shift_site(center).astype(float)
    center_p = convert_to_primitive(lattice, center, shift=True).astype(float)
    bond_p = np.concatenate([vector, center_p])

    c_wp = _find_wyckoff_site_sg(so_p, wyckoff_site, center, lattice, pset)[0]
    if c_wp is None:
        return None, None, None

    bond_p_list = create_equivalent_bond(so, so_p, bond_p, shift=True)
    wp_list = [wp for wp in wyckoff_site[c_wp]["bond"] if int(wp.split("@")[0][:-1]) == len(pset) * len(bond_p_list)]
    b_wp, sol = solve()  # (x,y,z) in sol are in conventional coordinate.
    if b_wp is None:
        return None, None, None
    sol = {k: round(v, DIGIT) for k, v in sol.items()}

    wyckoff_bond_wp = wyckoff_bond[b_wp]["conventional"][0]
    b0 = shift_bond(replace(wyckoff_bond_wp, sol))

    return b_wp, sol, b0


# ==================================================
def convert_to_bond(bond):
    """
    Convert to bond from str.

    Args:
        bond (str): bond.

    Returns:
        - (ndarray) -- bond (vector,center).

    Note:
        - bond string is "[tail];[head]/[vector]@[center]/[start]:[vector]".
    """
    if bond.count(";") > 0:
        t, h = bond.split(";")
        t = str_to_sympy(t)
        h = str_to_sympy(h)
        v = h - t
        c = (t + h) / 2
    elif bond.count("@") > 0:
        v, c = bond.split("@")
        v = str_to_sympy(v)
        c = str_to_sympy(c)
    elif bond.count(":") > 0:
        s, v = bond.split(":")
        s = str_to_sympy(s)
        v = str_to_sympy(v)
        c = s + v / 2

    b = np.concatenate([v, c])

    return b


# ==================================================
def find_wyckoff_bond(group, bond):
    """
    Find Wyckoff bond.

    Args:
        group (dict): group dict.
        bond (ndarray or str): bond (vector+center) to find.

    Returns:
        - (str) -- Wyckoff bond position (return None when not found).
        - (ndarray) -- set of Wyckoff bond (vector+center) (plus_set0, plus_set1, ...).

    Note:
        - bond string is "[tail];[head]/[vector]@[center]/[start]:[vector]".
    """
    if type(bond) == str:
        bond = convert_to_bond(bond)
    wyckoff_site = group["wyckoff"]["site"]
    wyckoff_bond = group["wyckoff"]["bond"]

    if group["info"].lattice != "0":
        so = group["symmetry_operation"]["fractional"]
        so_p = group["symmetry_operation"]["fractional_primitive"]
        pset = group["symmetry_operation"]["plus_set"]
        lattice = group["info"].lattice
        b_wp, sol, b0 = _find_wyckoff_bond_sg(so, so_p, wyckoff_site, wyckoff_bond, bond, lattice, pset)
        if b_wp is None:
            return None, None
        bonds = _evaluate_wyckoff_bond(wyckoff_bond[b_wp]["conventional"], sol, pset)
        return b_wp, bonds
    else:
        so = group["symmetry_operation"]["fractional"]
        b_wp, sol, b0 = _find_wyckoff_bond_pg(so, wyckoff_site, wyckoff_bond, bond)
        if b_wp is None:
            return None, None
        bonds = _evaluate_wyckoff_bond(wyckoff_bond[b_wp]["conventional"], sol)
        return b_wp, bonds


# ==================================================
def create_cell_site(group, site):
    """
    Create sites in unit cell (sorted and with plus set).

    Args:
        group (dict): group dict.
        site (ndarray or str): site.

    Returns:
        - (ndarray) -- sites in unit cell.
        - (list) -- SO mapping (index from 1).

    Note:
        - site string is "[x,y,z]".
    """
    if type(site) == str:
        site = str_to_sympy(site)
    if group["info"].lattice != "0":
        site = shift_site(site)

    so = group["symmetry_operation"]["fractional"].astype(float)

    if group["info"].lattice != "0":
        sites = (so @ np.pad(site, (0, 1), constant_values=1))[:, 0:3]
        if "plus_set" in group["symmetry_operation"].keys():
            pset = group["symmetry_operation"]["plus_set"]
            sites = np.concatenate([sites + i for i in pset.astype(float)])
        sites = shift_site(sites)
    else:
        sites = so[:, 0:3, 0:3] @ site

    u_sites = remove_duplicate_vector(sites)
    # rotate sites such that top site equals to given site.
    idx = np.where(np.linalg.norm(u_sites.astype(float) - site.astype(float), axis=1) < TOL_SAME_SITE)[0][0]
    u_sites = np.concatenate([u_sites[idx:], u_sites[:idx]], axis=0)

    mapping = create_mapping(sites, u_sites)
    mapping = [sorted([(j - 1) % len(so) + 1 for j in i]) for i in mapping]

    return u_sites, mapping


# ==================================================
def create_cell_bond(group, bond):
    """
    Create bonds in unit cell (sorted, nondirectional and with plus set).

    Args:
        group (dict): group dict.
        bond (ndarray or str): bond.

    Returns:
        - (ndarray) -- bonds in unit cell.
        - (list) -- SO mapping (index from 1).

    Note:
        - bond string is "[tail];[head]/[vector]@[center]/[start]:[vector]".
        - negative value in mapping represents reversed bond.
    """
    if type(bond) == str:
        bond = convert_to_bond(bond)
    v0, s0 = bond[0:3], bond[3:6]
    if group["info"].lattice != "0":
        s0 = shift_site(s0)
        bond = np.concatenate([v0, s0])

    so = group["symmetry_operation"]["fractional"].astype(float)

    vectors = so[:, 0:3, 0:3] @ v0

    if group["info"].lattice != "0":
        centers = (so @ np.pad(s0, (0, 1), constant_values=1))[:, 0:3]
        if "plus_set" in group["symmetry_operation"].keys():
            pset = group["symmetry_operation"]["plus_set"]
            centers = np.concatenate([centers + i for i in pset.astype(float)])
            vectors = np.tile(vectors, (len(pset), 1))
        centers = shift_site(centers)
    else:
        centers = so[:, 0:3, 0:3] @ s0

    bonds = np.hstack((vectors[:, None], centers[:, None])).reshape(-1, 6)
    u_bonds = remove_duplicate_vector(bonds, directional=False)
    # rotate bonds such that top bond equals to given bond.
    idx = np.where(np.linalg.norm(u_bonds.astype(float) - bond.astype(float), axis=1) < TOL_SAME_SITE)[0][0]
    u_bonds = np.concatenate([u_bonds[idx:], u_bonds[:idx]], axis=0)

    mapping = create_mapping(bonds, u_bonds)
    mapping = [sorted([(j - 1) % len(so) + 1 for j in i]) for i in mapping]

    return u_bonds, mapping
