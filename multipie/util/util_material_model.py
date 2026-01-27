"""
Utility for MaterialModel calss.
"""

import numpy as np

from multipie import RepSiteType, CellSiteType, BondInfoType, RepBondType, CellBondType, BraketInfoType
from multipie.util.util import progress_bar_step, progress_bar_done
from multipie.util.util_crystal import site_distance, shift_site, convert_to_primitive, TOL_SAME_SITE
from multipie.util.util_wyckoff import find_vector
from multipie.core.default_model import _site_property, _bond_property


TOL = 100.0 * TOL_SAME_SITE
DIGIT = 6


# ==================================================
def get_bond(tail, head, n, no):
    """
    Get bond name.

    Args:
        tail (str): tail atom.
        head (str): head atom.
        n (int): neighbor (from 1).
        no (int): multiplicity (from 1).

    Returns:
        - (str) -- bond name.
    """
    return f"{tail};{head}_{n:03d}_{no}"


# ==================================================
def get_tail_head(site_bond):
    """
    Get tail and head atoms.

    Args:
        site_bond (str): site or bond tag.

    Returns:
        - (str) -- tail atom.
        - (str) -- head atom.
    """
    if ";" not in site_bond:
        tail, head = site_bond, site_bond
    else:
        tail, head = site_bond.split(";")
        head = head.split("_")[0]

    return tail, head


# ==================================================
def get_neighbor_info(site_bond):
    """
    Get neighbor info.

    Args:
        site_bond (str): site or bond tag.

    Returns:
        - (int) -- neighbor (0 for site).
        - (int) -- multiplicity (-1 for site).
    """
    if ";" not in site_bond:
        return (0, -1)
    else:
        return tuple(map(int, site_bond.split(";")[1].split("_")[1:]))


# ==================================================
def unique_vector_index(vectors, so, tol=TOL):
    """
    Get unique vector indices (remove SO related vecrors).

    Args:
        vectors (ndarray): vectors.
        so (ndarray): set of symmetry operations (3x3) except for identity.
        tol (float, optional): tolerance to check same vector.

    Returns:
        - (ndarray) -- unique vector indices.
    """
    check_idx = np.arange(len(vectors))
    idx = []

    while len(check_idx) > 0:
        i = check_idx[0]
        v = vectors[i]
        idx.append(i)

        vd = so @ v

        tail_idx = check_idx[1:]
        tail_vecs = vectors[tail_idx]

        diff = tail_vecs[:, None, :] - vd[None, :, :]
        is_match = np.all(np.isclose(diff, 0, atol=tol), axis=2)
        mask = ~np.any(is_match, axis=1)

        check_idx = tail_idx[mask]

    return idx


# ==================================================
def create_equivalent_bond(so, pset, bond):
    """
    Create unique bonds.

    Args:
        so (ndarray): symmetry operation (conventional) (PG:3x3) or (SG:4x4).
        pset (ndarray): plus set.
        bond (ndarray): source bond.

    Returns:
        - (ndarray) -- unique bonds.
    """
    v0, c0 = bond[0:3], bond[3:6]
    vector = so[:, 0:3, 0:3] @ v0

    if pset is not None:
        center = (so @ np.pad(c0, (0, 1), constant_values=1))[:, 0:3]
        center = np.concatenate([center + i for i in pset.astype(float)])
        center = shift_site(center)
        vector = np.tile(vector, (len(pset), 1))
    else:
        center = so[:, 0:3, 0:3] @ c0

    bonds = np.hstack((vector[:, None], center[:, None])).reshape(-1, 6)

    return bonds


# ==================================================
def get_unique_bond(bond, tol=TOL):
    """
    Get unique bond.

    Args:
        bond (ndarray): set of bonds.
        tol (float, optional): tolerance.

    Returns:
        - (ndarray) -- unique bonds.
    """

    def regular_direction(v):
        d = np.array(v[:3])
        d[np.abs(d) < tol] = 0
        sign = next((1 if x > 0 else -1 for x in d if abs(x) >= tol), 1)
        return np.concatenate([sign * d, np.where(np.abs(v[3:]) < tol, 0, v[3:])])

    key = lambda v: tuple(int(round(x / tol)) for x in v)
    return np.asarray(sorted([regular_direction(v).tolist() for v in bond], key=key))


# ==================================================
def unique_bond_index(bond, tol=TOL):
    """
    Unique bond index.

    Args:
        bond (ndarray): bond set.
        tol (float, optional): tolerance.

    Returns:
        - (list) -- unique index.
    """
    n = len(bond)
    visited = np.zeros(n, bool)
    unique_indices = []

    for i in range(n):
        if not visited[i]:
            unique_indices.append(i)
            visited |= np.linalg.norm(bond - bond[i], axis=1) <= tol

    return unique_indices


# ==================================================
def get_basis_type(site_data, spinful):
    """
    Get atomic basis type.

    Args:
        site_data (dict): site data.
        spinful (bool): spinful basis ?

    Returns:
        - (str) -- atomic basis type, "jml/lgs/lg".
    """
    lst = []
    for pos, orb in site_data.values():
        tp = type(orb)
        if tp == str:
            lst.append(orb.count("(") > 0)
        elif tp == list:
            lst += [i.count("(") > 0 for i in orb]
        elif tp == dict:
            raise Exception("not implemented.")
        else:
            raise Exception(f"invalid orbital format, {orb}.")

    if all(lst):
        basis_type = "jml"
    elif not any(lst):
        if spinful:
            basis_type = "lgs"
        else:
            basis_type = "lg"
    else:
        raise Exception("jml, lgs, and lg formats coexist.")

    return basis_type


# ==================================================
def parse_orbital(orbital, basis_type, basis_info):
    """
    Parse atomic orbital information. ### input check needs to be implemented!!!

    Args:
        orbital (str or list): orbital information.
        basis_type (str): atomic orbital type.
        basis_info (dict): basis info.

    Returns:
        - ([[str]]) -- atomic basis in order of rank, s,p,d,f.

    Notes:
        - orbital format (str), each orbital (px, dxz, ...) or orbital set (p, (3/2,p), ...).
        - orbital format (list), set of the above format.
        - returned list is in full expression in fixed order compatible with atomic matrix element.
        - coexistence of spinless and spinfull, or different basis types (jml, lgs, lg) is not allowed.
    """
    str_rank = {"s": 0, "p": 1, "d": 2, "f": 3}

    def check_block(s):
        return len(s) > 0 and s[0] in str_rank.keys()

    def regularize(orb, basis_type, basis_info):
        orb = orb.replace(" ", "").lower()
        if basis_type == "jml":
            v = orb.replace("(", "").replace(")", "").split(",")
            if len(v[1]) == 1 and check_block(v[1]):
                rank = str_rank[v[1][0]]
                j = v[0]
                basis = [i for i in basis_info["jml"][rank] if i.split(",")[0][1:] == j]
            else:
                raise Exception(f"invalid orbital format, {orbital}.")
        elif basis_type == "lgs":
            if check_block(orb):
                rank = str_rank[orb[0]]
                if len(orb) == 1:
                    basis = basis_info["lgs"][rank]
                else:
                    basis = [f"({orb},u)", f"({orb},d)"]
            else:
                raise Exception(f"invalid orbital, {orb}.")
        elif basis_type == "lg":
            if check_block(orb):
                rank = str_rank[orb[0]]
                if len(orb) == 1:
                    basis = basis_info["lg"][rank]
                else:
                    basis = [orb]
            else:
                raise Exception(f"invalid orbital, {orb}.")
        else:
            raise Exception(f"invalid orbital format, {orbital}.")
        return rank, basis

    tp = type(orbital)
    basis_set = [[], [], [], []]  # s, p, d, f.
    if tp == str:
        rank, basis = regularize(orbital, basis_type, basis_info)
        basis_set[rank] += basis
    elif tp == list:
        for i in orbital:
            rank, basis = regularize(i, basis_type, basis_info)
            basis_set[rank] += basis
    elif tp == dict:
        raise Exception(f"to be available.")
    else:
        raise Exception(f"invalid orbital format, {orbital}.")

    basis_set = [sorted(block, key=lambda x: basis_info[basis_type][rank].index(x)) for rank, block in enumerate(basis_set)]

    return basis_set


# ==================================================
def convert_orbital_index(basis_set, basis_info_type):
    """
    Convert to orbital index.

    Args:
        basis_set ([[str]]): atomic basis.
        basis_info_type (dict): basis info. of given basis_type.

    Returns:
        - ([[int]]) -- orbital index corresponding to basis set.
    """
    bs_idx = [[basis_info_type[rank].index(t) for t in bsr if t in bsr] for rank, bsr in enumerate(basis_set)]

    return bs_idx


# ==================================================
def parse_neighbor(neighbor, tail_rank, head_rank):
    """
    Parse neighbor.

    Args:
        neighbor (int or [int] or tuple): neighbor info., (int or [int], [tail rank], [head rank]).
        tail_rank ([int] or [str]): tail rank.
        head_rank ([int] or [str]): head rank.

    Returns:
        - (tuple) -- ([neighbor], [tail rank], [head rank]).
    """
    rank = {"s": 0, "p": 1, "d": 2, "f": 3, 0: 0, 1: 1, 2: 2, 3: 3}

    if type(neighbor) == int:  # max neighbor.
        tail_rank = sorted(list(set(tail_rank)))
        head_rank = sorted(list(set(head_rank)))
        neighbor = (list(range(1, neighbor + 1)), tail_rank, head_rank)
    elif type(neighbor) == list:
        tail_rank = sorted(list(set(tail_rank)))
        head_rank = sorted(list(set(head_rank)))
        neighbor = (neighbor, tail_rank, head_rank)
    elif type(neighbor) == tuple:
        max_neighbor, tail, head = neighbor
        tail1 = set([rank[i] for i in tail])
        head1 = set([rank[i] for i in head])
        if not tail1.issubset(tail_rank):
            raise Exception(f"{tail} are not in {tail_rank}.")
        if not head1.issubset(head_rank):
            raise Exception(f"{head} are not in {head_rank}.")
        tail_rank = sorted(list(tail1))
        head_rank = sorted(list(head1))
        neighbor = (list(range(1, max_neighbor + 1)), tail_rank, head_rank)
    else:
        raise Exception(f"unknown format for neighbor {neighbor}")

    return neighbor


# ==================================================
def parse_samb_select(select, irreps):
    """
    Parse SAMB select dict.

    Args:
        select (dict): select condition dict.
        irreps (list): list of irreps.

    Returns:
        - (dict) -- regularized SAMB select.

    Note:
        - dict keys are "X", "l", "Gamma", "s".
        - null list represents all possible values.
        - "IR" in "Gamma" indicates the identity representation.
    """
    if select is None:
        select = {}

    for key in select.keys():
        if key not in ["X", "l", "Gamma", "s"]:
            raise ValueError(f"unknown key, {key}.")

    if "Gamma" in select.keys() and select["Gamma"] == "IR":
        select["Gamma"] = [irreps[0]]

    select = {k: [v] if type(v) != list else v for k, v in select.items()}

    if "X" not in select.keys() or len(select["X"]) == 0:
        select["X"] = ["Q", "G", "M", "T"]
    if "l" not in select.keys() or len(select["l"]) == 0:
        select["l"] = list(range(12))
    if "Gamma" not in select.keys() or len(select["Gamma"]) == 0:
        select["Gamma"] = irreps
    if "s" not in select.keys() or len(select["s"]) == 0:
        select["s"] = [0, 1]

    return select


# ==================================================
def parse_combined_select(select, irreps, default_samb_select, site_rep, bond_rep):
    """
    Parse select dict for combined SAMB.

    Args:
        select (dict): select condition dict.
        irreps (list): list of irreps.
        default_samb_select (dict): default SAMB select.
        site_rep (dict): site representative dict.
        bond_rep (dict): bond representative dict.

    Returns:
        - (dict) -- regularized SAMB select.
        - (dict) -- regularized other select.

    Note:
        - dict keys are "site" and "bond" in addition to SAMB_select.
        - null list represents all possible values.
    """
    if select is None:
        select = {}

    for key in select.keys():
        if key not in ["site", "bond", "X", "l", "Gamma", "s"]:
            raise ValueError(f"unkwon key, {key}.")

    # split samb select and others.
    samb_select = {key: val for key, val in select.items() if key in ["X", "l", "Gamma", "s"]}
    select = {key: val for key, val in select.items() if key in ["site", "bond"]}

    # regularize samb select.
    regularized_samb_select = parse_samb_select(samb_select, irreps)
    regularized_samb_select = {k: [x for x in v if x in default_samb_select[k]] for k, v in regularized_samb_select.items()}

    keys = list(regularized_samb_select.keys())
    for k, v in default_samb_select.items():
        if k not in keys:
            regularized_samb_select[k] = v

    # regularized combined select.
    select = {k: v if isinstance(v, list) else [v] for k, v in select.items()}

    regularized_select = {}
    for k, v in select.items():
        if k == "site":
            regularized_select["site"] = []
            for i in v:
                if isinstance(i, str):  # site.
                    regularized_select["site"].append((i, None))
                else:
                    a, b = i  # site, rank/[rank]
                    if not isinstance(b, list):
                        b = [b]
                    regularized_select["site"].append((a, b))
        elif (
            k == "bond"
        ):  # name, neighbor, [neighbor], (name, neighbor), (name, [neighbor]), (name, rank), (name, rank, neighbor), (name, rank, [neighbor]).
            regularized_select["bond"] = []
            for i in v:
                if isinstance(i, str):  # tail;head.
                    regularized_select["bond"].append((i, None, None))
                elif isinstance(i, int):  # neighbor.
                    regularized_select["bond"].append((None, None, [i]))
                elif isinstance(i, list) and all(isinstance(x, int) for x in i):  # neighbor.
                    regularized_select["bond"].append((None, None, i))
                elif len(i) == 2:
                    a, b = i
                    if isinstance(b, list):  # tail;head, neighbor.
                        regularized_select["bond"].append((a, None, b))
                    elif isinstance(b, str):  # tail;head, t_rank;h_rank.
                        regularized_select["bond"].append((a, b, None))
                    else:
                        regularized_select["bond"].append((a, None, [b]))  # tail;head, neighbor.
                elif len(i) == 3:
                    a, b, c = i
                    if isinstance(c, list):
                        regularized_select["bond"].append((a, b, c))  # tail;head, t_rank;h_rank, neighbor.
                    else:
                        regularized_select["bond"].append((a, b, [c]))  # tail;head, t_rank;h_rank, neighbor.

    # final filter.
    final_select = {}

    # site.
    default_site = [(k, [no for no, i in enumerate(v.orbital) if len(i) > 0]) for k, v in site_rep.items()]
    if "site" not in select.keys():
        final_select["site"] = default_site
    else:
        site = []
        for name, orb in regularized_select["site"]:
            if orb is None:
                d = [(s, o) for s, o in default_site if s == name]
            else:
                d = [(s, sorted(list(set(o) & set(orb)))) for s, o in default_site if s == name]
            site += [i for i in d if len(i[1]) > 0]
        site = list({tuple(tuple(x) if isinstance(x, list) else x for x in t) for t in site})
        site = [tuple(list(x) if isinstance(x, tuple) and all(isinstance(i, int) for i in x) else x for x in t) for t in site]
        final_select["site"] = sorted(site)
    site_list = [i[0] for i in final_select["site"]]

    # bond.
    default_bond = [
        (v.tail, v.head, v.neighbor, v.t_rank, v.h_rank) for v in bond_rep.values() if v.tail in site_list and v.head in site_list
    ]
    if "bond" not in select.keys():
        final_select["bond"] = default_bond
    else:
        bond = []
        for name, rank, neighbor in regularized_select["bond"]:
            if rank is None and neighbor is None:  # name only.
                tail, head = name.split(";")
                d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default_bond if t == tail and h == head]
                d += [(h, t, n, hr, tr) for t, h, n, tr, hr in default_bond if h == tail and t == head]
            elif name is None and rank is None:  # neighbor only.
                d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default_bond if n in neighbor]
            elif rank is None:  # name, neighbor.
                tail, head = name.split(";")
                d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default_bond if t == tail and h == head and n in neighbor]
                d += [(h, t, n, hr, tr) for t, h, n, tr, hr in default_bond if h == tail and t == head and n in neighbor]
            elif neighbor is None:  # name, rank.
                tail, head = name.split(";")
                t_rank, h_rank = rank.split(";")
                t_rank, h_rank = {int(t_rank)}, {int(h_rank)}
                d = [
                    (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                    for t, h, n, tr, hr in default_bond
                    if t == tail and h == head
                ]
                d += [
                    (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                    for t, h, n, tr, hr in default_bond
                    if h == tail and t == head
                ]
            else:
                tail, head = name.split(";")
                t_rank, h_rank = rank.split(";")
                t_rank, h_rank = set(int(t_rank)), set(int(h_rank))
                d = [
                    (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                    for t, h, n, tr, hr in default_bond
                    if t == tail and h == head and n in neighbor
                ]
                d += [
                    (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                    for t, h, n, tr, hr in default_bond
                    if h == tail and t == head and n in neighbor
                ]
            bond += d
        bond = list({tuple(tuple(x) if isinstance(x, list) else x for x in t) for t in bond})
        bond = [tuple(list(x) if isinstance(x, tuple) and all(isinstance(i, int) for i in x) else x for x in t) for t in bond]
        final_select["bond"] = sorted(bond)

    return regularized_samb_select, final_select


# ==================================================
def create_site_grid(site_dict, igrid=None):
    """
    Create site grid.

    Args:
        site_dict (dict): site dict.
        igrid (ndarray, optional): integer grid.

    Returns:
        - (dict) -- site grid (sorted), Dict[name, Dict[(#sublattice,#plus_set,i1,i2,i3), position]].
    Note:
        - if igrid is None, [[0,0,0]] is used.
    """
    if igrid is None:
        igrid = np.array([[0, 0, 0]], dtype=int)
    igrid_list = [i.tolist() for i in igrid]

    cell_site = site_dict["cell"]

    # add each grid point.
    site_grid = {}
    for name, cell_site_name in cell_site.items():
        site_grid_each = {(c.sublattice, c.plus_set, *i): c.position + i for c in cell_site_name for i in igrid_list}
        site_grid[name] = dict(sorted(site_grid_each.items()))

    return site_grid


# ==================================================
def create_site_so(group, site_dict):
    """
    Create symmetry operations for first Wyckoff site.

    Args:
        group (dict): group dict.
        site_dict (dict): site dict.

    Returns:
        - (dict) -- symmetry operations except identity, Dict[name, SOs].
    """
    cell_site = site_dict["cell"]
    so = group.symmetry_operation["fractional"][:, 0:3, 0:3].astype(float)

    site_so = {name: so[np.array(cell_site[name][0].mapping) - 1][1:] for name in cell_site.keys()}

    return site_so


# ==================================================
def create_wyckoff_dict(rep_site, rep_bond):
    """
    Create site_bond to wyckoff dict.

    Args:
        rep_site (dict): representative site dict.
        rep_bond (dict): representative bond dict.

    Returns:
        - (dict) -- site_bond to wyckoff dict, Dict[site_bond, wyckoff].
    """
    wyckoff_dict = {}
    for site, lst in rep_site.items():
        wyckoff_dict[site] = lst.wyckoff
    for bond, lst in rep_bond.items():
        wyckoff_dict[bond] = lst.wyckoff

    return wyckoff_dict


# ==================================================
def create_braket_dict(rep_site, rep_bond, basis_info_type):
    """
    Create site_bond to braket dict.

    Args:
        rep_site (dict): representative site dict.
        rep_bond (dict): representative bond dict.
        basis_info_type (str): basis info type.

    Returns:
        - (dict) -- site_bond to braket dict, Dict[site_bond, [BraketInfoType]].
    """
    braket_dict = {}

    # site.
    for name, lst in rep_site.items():
        rank = [no for no, o in enumerate(lst.orbital) if len(o) > 0]
        for tr in rank:
            for hr in rank:
                if hr > tr:  # skip if head_rank > tail_rank.
                    continue
                tidx = tuple(convert_orbital_index(lst.orbital, basis_info_type)[tr])
                hidx = tuple(convert_orbital_index(lst.orbital, basis_info_type)[hr])
                braket_dict[name] = braket_dict.get(name, []) + [BraketInfoType(hr, hidx, tr, tidx)]

    # bond.
    for name, lst in rep_bond.items():
        tail = lst.tail
        head = lst.head
        tail_rank = lst.t_rank
        head_rank = lst.h_rank
        for tr in tail_rank:
            for hr in head_rank:
                if head == tail and (hr > tr):  # skip if head_rank > tail_rank among same atom sites.
                    continue

                tidx = tuple(convert_orbital_index(rep_site[tail].orbital, basis_info_type)[tr])
                hidx = tuple(convert_orbital_index(rep_site[head].orbital, basis_info_type)[hr])

                braket_dict[name] = braket_dict.get(name, []) + [BraketInfoType(hr, hidx, tr, tidx)]

    return braket_dict


# ==================================================
def create_full_matrix_info(site_dict):
    """
    Create full matrix info.

    Args:
        site_dict (dict): site dict.

    Returns:
        - (dict) -- full matrix info, "ket": [(name, sublattice, rank, orbital)], "index": Dict[(atom,sublattice,rank), (top_idx, size)].
    """
    ket = []
    ket_dict = {}
    start_idx = 0
    for name, lst in site_dict["representative"].items():
        for c in site_dict["cell"][name]:
            if c.plus_set != 1:
                continue
            for rank, orbitals in enumerate(lst.orbital):
                num_orb = len(orbitals)
                for o in orbitals:
                    ket.append([name, c.sublattice, rank, o])

                if num_orb > 0:
                    ket_dict[(name, c.sublattice, rank)] = (start_idx, num_orb)
                    start_idx += num_orb

    dic = {"ket": ket, "index": ket_dict}

    return dic


# ==================================================
def parse_representative_site(group, site_data, basis_type, basis_info):
    """
    Parse representative site.

    Args:
        group (dict): group dict.
        site_data (dict): site data.
        basis_type (str): atomic basis type, "jml/lgs/lg".
        basis_info (dict): basis info dict.

    Returns:
        - (dict) -- site dict.
        - (list) -- atomic orbital rank list, [(tail,head,tail_rank,head_rank)].

    Notes:
        - "representative": representative site, Dict[name, RepSiteType].
        - "cell": cell site, Dict[name, [CellSiteType] ].
    """
    lattice = group.info.lattice
    ps = group.symmetry_operation.get("plus_set", None)
    npset = 1 if ps is None else len(ps)
    site_data = dict(sorted(site_data.items()))

    rep_site = {}
    cell_site = {}
    for c_no, (name, (pos, orb)) in enumerate(site_data.items()):
        pos = str(pos)
        wp, sites = group.find_wyckoff_site(pos)
        sites_primitive = convert_to_primitive(lattice, sites, shift=True)
        wyckoff_site = group.wyckoff["site"][wp]
        sym = wyckoff_site["symmetry"]
        pos = sites[0].tolist()
        orb = parse_orbital(orb, basis_type, basis_info)

        mapping = wyckoff_site["mapping"]
        n_sub = len(mapping)
        mapping = mapping * npset
        sublattice = [no % n_sub + 1 for no in range(len(sites))]
        pset = [no // n_sub + 1 for no in range(len(sites))]

        rep_site[name] = RepSiteType(c_no + 1, wp, sym, pos, orb)
        cell_site[name] = [
            CellSiteType(i + 1, s, sp, m, sl, ps)
            for i, (s, sp, m, sl, ps) in enumerate(zip(sites, sites_primitive, mapping, sublattice, pset))
        ]

    dic = {"representative": rep_site, "cell": cell_site}

    return dic


# ==================================================
def create_representative_bond(group, G, so, tail, heads, max_neighbor):
    """
    Create representative bonds.

    Args:
        group (dict): group dict.
        G (ndarray): metric tensor (3x3).
        so (ndarray): symmetry operations (nx3x3).
        tail (ndarray): tail position.
        heads (ndarray): head positions over symmetry related with grid.
        max_neighbor (int): max. neighbor.

    Returns:
        - (list) -- representative bonds for each neighbor.
    """
    all_bond = site_distance(tail, heads, G)
    if group.is_point_group:
        all_bond = list(all_bond.values())
    else:
        all_bond = list(all_bond.values())[:max_neighbor]
    rep_bond = []
    for lst in all_bond:
        vectors = lst - tail
        centers = 0.5 * (lst + tail)
        bonds = np.hstack((vectors[:, None], centers[:, None])).reshape(-1, 6)
        idx = unique_vector_index(vectors, so)
        bonds = bonds[idx]
        rep_bond.append(bonds)

    return rep_bond


# ==================================================
def remove_equivalent_representative_bond(so, pset, bonds):
    """
    Remove equivalnet representative bond.

    Args:
        so (ndarray): symmetry operation (conventional, fractional, 3x3 or 4x4).
        pset (ndarray): plus set.
        bonds (ndarray): set of representative bonds.

    Returns:
        - (ndarray) unique representative bonds.
    """
    unique_bonds = []
    for bn in bonds:
        if len(bn) == 1:
            unique_bonds.append(bn)
        else:
            bn_all = np.asarray([get_unique_bond(create_equivalent_bond(so, pset, i)).reshape(-1) for i in bn])
            idx = unique_bond_index(bn_all)
            bn = bn[idx]
            unique_bonds.append(bn)
    return unique_bonds


# ==================================================
def parse_representative_bond(group, G, site_grid, site_so, site_dict, bond_data, max_neighbor, verbose):
    """
    Parse representative bond.

    Args:
        group (dict): group dict.
        G (ndarray): metric tensor (3x3).
        site_grid (dict): site grid dict.
        site_so (dict): site SO dict.
        site_dict (dict): site dict.
        bond_data (list): bond data.
        max_neighbor (int): max. neighbor.
        verbose (bool): verbose progress bar ?

    Returns:
        - (dict) -- bond dict.
        - (list) -- atomic orbital rank list, [(tail,head,tail_rank,head_rank)].

    Notes:
        - "representative": representative bond, Dict[name, RepBondType].
        - "cell": cell bond, Dict[name, [CellBondType] ].
        - "neighbor": Dict[name, Dict[#neighbor, [rep_bond_tag] ]].
        - "info": [BondInfoType].
    """
    if verbose:
        progress = progress_bar_step(label="Analyzing ...")

    lattice = group.info.lattice
    ps = group.symmetry_operation.get("plus_set", None)
    npset = 1 if ps is None else len(ps)
    so_all = group.symmetry_operation["fractional"].astype(float)
    bond_data = sorted(bond_data)

    rep_bond = {}
    cell_bond = {}
    info_bond = []
    c_no = 0
    for tail_tag, head_tag, neighbor in bond_data:
        if tail_tag not in site_dict["representative"].keys():
            raise Exception(f"{tail_tag} is not found in sites.")
        if head_tag not in site_dict["representative"].keys():
            raise Exception(f"{head_tag} is not found in sites.")

        # swap if head_tag > tail_tag.
        if head_tag > tail_tag:
            tail_tag, head_tag = head_tag, tail_tag

        tail_rank = [no for no, orb in enumerate(site_dict["representative"][tail_tag].orbital) if len(orb) > 0]
        head_rank = [no for no, orb in enumerate(site_dict["representative"][head_tag].orbital) if len(orb) > 0]
        neighbor, tail_rank, head_rank = parse_neighbor(neighbor, tail_rank, head_rank)

        max_n = min(max(neighbor), max_neighbor)
        info_bond.append(BondInfoType(tail_tag, head_tag, neighbor, tail_rank, head_rank))

        # site idx.
        tail_info = site_dict["cell"][tail_tag]
        tail_pos = [i.position for i in tail_info]
        head_info = site_dict["cell"][head_tag]
        head_pos = [i.position for i in head_info]

        # create rep. bond.
        tail = site_grid[tail_tag][(1, 1, 0, 0, 0)]
        heads = np.asarray(list(site_grid[head_tag].values()))
        so = site_so[tail_tag]
        rep_bond_each = create_representative_bond(group, G, so, tail, heads, max_n)
        rep_bond_each = remove_equivalent_representative_bond(so_all, ps, rep_bond_each)

        for n, bonds in enumerate(rep_bond_each):
            if n + 1 not in neighbor:
                continue
            if verbose:
                next(progress)
            for bno, b in enumerate(bonds):
                name = get_bond(tail_tag, head_tag, n + 1, bno + 1)
                b_wp, all_bond = group.find_wyckoff_bond(b)
                wyckoff_bond = group.wyckoff["bond"][b_wp]
                vector, center = all_bond[:, 0:3], all_bond[:, 3:6]
                vector_p = convert_to_primitive(lattice, vector, shift=False)
                center_p = convert_to_primitive(lattice, center, shift=False)
                v0, c0 = vector[0], center[0]
                dist = float(np.sqrt(v0 @ G @ v0))
                v0, c0 = v0.tolist(), c0.tolist()
                mapping = wyckoff_bond["mapping"]
                d = str(mapping).count("-") == 0
                tail_p = center - 0.5 * vector
                head_p = center + 0.5 * vector
                if ps is not None:
                    tail_p = shift_site(tail_p, TOL)
                    head_p = shift_site(head_p, TOL)
                n_sub = len(mapping)
                mapping = mapping * npset
                sublattice = [no % n_sub + 1 for no in range(len(all_bond))]
                pset = [no // n_sub + 1 for no in range(len(all_bond))]

                tail_idx = [find_vector(i, tail_pos, TOL) + 1 for i in tail_p]
                tail_idx = [(tail_info[i - 1].sublattice, tail_info[i - 1].plus_set) for i in tail_idx]
                head_idx = [find_vector(i, head_pos, TOL) + 1 for i in head_p]
                head_idx = [(head_info[i - 1].sublattice, head_info[i - 1].plus_set) for i in head_idx]

                rep_bond[name] = RepBondType(c_no + 1, tail_tag, head_tag, n + 1, b_wp, d, v0, c0, dist, tail_rank, head_rank)
                cell_bond[name] = []
                for k, (v, vp, c, cp, m, sl, pi, ti, hi) in enumerate(
                    zip(vector, vector_p, center, center_p, mapping, sublattice, pset, tail_idx, head_idx)
                ):
                    s_tail = tail_pos[ti[0] - 1]
                    s_tail = convert_to_primitive(lattice, s_tail, shift=False)
                    s_head = head_pos[hi[0] - 1]
                    s_head = convert_to_primitive(lattice, s_head, shift=False)
                    n1, n2, n3 = -(vp - (s_head - s_tail))
                    n1, n2, n3 = round(n1), round(n2), round(n3)
                    cell_bond[name].append(CellBondType(k + 1, v, vp, c, cp, m, sl, pi, ti, hi, (n1, n2, n3)))
                c_no += 1

    dic = {"representative": rep_bond, "cell": cell_bond, "info": info_bond}

    if verbose:
        progress_bar_done(label="Analyzing ...")

    return dic


# ==================================================
def write_site_dict(site_dict):
    """
    Write site dict.

    Args:
        site_dict (dict): site dict.
    """
    rep_site = site_dict["representative"]
    cell_site = site_dict["cell"]

    print("--- representative site ---")
    for name, c in rep_site.items():
        print(f"{name}: #{c.no}, wyckoff = {c.wyckoff}, symmetry = {c.symmetry}, 1st = {c.position}, orbital = {c.orbital}")

    print("--- cell site ---")
    for name, cell_site_name in cell_site.items():
        print(f"tag = {name}")
        for i in cell_site_name:
            print(
                f"#{i.no}: position = {i.position.tolist()}, mapping = {i.mapping}, sublattice = {i.sublattice}, plus set = {i.plus_set}"
            )


# ==================================================
def write_site_grid(site_grid):
    """
    Write site grid.

    Args:
        site_grid (dict): site grid.
    """
    print("--- site grid ---")
    for name, val in site_grid.items():
        print(f"- {name} -")
        for idx, v in val.items():
            print("(sublattice, plus set, grid) =", idx, v.tolist())


# ==================================================
def write_bond_dict(bond_dict):
    """
    Write bond dict.

    Args:
        bond_dict (dict): bond dict.
    """
    rep_bond = bond_dict["representative"]
    cell_bond = bond_dict["cell"]
    info = bond_dict["info"]

    print("--- info ---")
    for i in info:
        print(f"{i.head}-{i.tail}: neighbor={i.neighbor}, head rank={i.h_rank}, tail rank={i.t_rank}")

    print("--- representative bond ---")
    for name, i in rep_bond.items():
        pos = str(i.vector) + "@" + str(i.center)
        print(
            f"{name}: #{i.no}, {i.neighbor}th, directional = {i.directional}, wyckoff = {i.wyckoff}, 1st = {pos}, distance = {i.distance}"
        )

    print("--- cell bond ---")
    for name, cell_bond_name in cell_bond.items():
        print(f"tag = {name}")
        for i in cell_bond_name:
            print(
                f"#{i.no}: bond = {i.vector.tolist()}@{i.center.tolist()}, mapping = {i.mapping}, sublattice = {i.sublattice}, plus set = {i.plus_set}, head(sublattice,plus_set) = {i.h_idx}, tail(sublattice,plus_set) = {i.t_idx}"
            )


# ==================================================
def qtdraw_site(qtdraw, site_dict, scale, mode, radius, show_rep_site):
    """
    Draw site by QtDraw.

    Args:
        qtdraw (QtDraw): QtDraw widget or application.
        site_dict (dict): site dict.
        scale (float): scale.
        mode (str): draw mode, "standard/detail".
        radius (float): base radius.
        show_rep_site (bool): show representative site?
    """
    rep_site = site_dict["representative"]
    cell_site = site_dict["cell"]

    for c_no, (name, cell_site_name) in enumerate(cell_site.items()):
        prop_no = min(c_no, len(_site_property) - 1)
        color, size, opacity = _site_property[prop_no]

        wp = rep_site[name].wyckoff
        sym = rep_site[name].symmetry
        for c in cell_site_name:
            c_name = f"{name}({c.plus_set})"

            if mode == "standard":
                label = f"#{c_no+1}({wp},{c.sublattice})"
            else:
                label = f"#{c_no+1}:{wp},{c.sublattice}) [{sym}]"

            if show_rep_site and c.sublattice == 1 and c.plus_set == 1:
                qtdraw.add_site(
                    position=c.position,
                    name=c_name + "*",
                    label=label,
                    color="gold",
                    size=size * scale * 1.2 * radius,
                    opacity=0.2,
                )

            qtdraw.add_site(
                position=c.position,
                name=c_name,
                label=label,
                color=color,
                size=size * scale * radius,
                opacity=opacity,
            )


# ==================================================
def qtdraw_bond(qtdraw, bond_dict, max_neighbor, scale, mode, width, show_rep_bond):
    """
    Draw bond by QtDraw.

    Args:
        qtdraw (QtDraw): QtDraw widget or application.
        bond_dict (dict): bond dict.
        max_neighbor (int): max. neighbor to draw.
        scale (float): scale.
        mode (str): draw mode, "standard/detail".
        width (float): base width.
        show_rep_bond (bool): show representative bond?
    """
    rep_bond = bond_dict["representative"]
    cell_bond = bond_dict["cell"]

    for c_no, (name, cell_bond_name) in enumerate(cell_bond.items()):
        n = rep_bond[name].neighbor
        wp = rep_bond[name].wyckoff
        directional = rep_bond[name].directional
        if n > max_neighbor:
            continue
        prop_no = min(c_no, len(_bond_property) - 1)
        ((color, color2), width1, opacity) = _bond_property[prop_no]
        if not directional:
            color2 = color
        for c in cell_bond_name:
            c_name = f"{name}({c.plus_set})"

            if mode == "standard":
                label = f"#{c_no+1}({wp},{c.sublattice})"
            else:
                label = f"#{c_no+1}({wp},{c.sublattice}) [{c.t_idx[0]}({c.t_idx[1]});{c.h_idx[0]}({c.h_idx[1]})]"

            if show_rep_bond:
                opacity1 = 0.5
                vcolor = "red" if c.sublattice == 1 and c.plus_set == 1 else "black"
                qtdraw.add_vector(
                    position=c.center,
                    direction=c.vector,
                    length=-0.22,
                    width=0.5 * width1 * width * scale,
                    cartesian=False,
                    name=c_name + "*",
                    label=label,
                    color=vcolor,
                )
            else:
                opacity1 = opacity

            qtdraw.add_bond(
                position=c.center,
                direction=c.vector,
                width=width1 * width * scale,
                cartesian=False,
                color=color,
                color2=color2,
                name=c_name,
                label=label,
                opacity=opacity1,
            )


# ==================================================
def create_qtdraw(qtdraw, group, name, cell_info, site_dict, bond_dict, prop):
    """
    Create QtDraw.

    Args:
        qtdraw (PyVistaWidget): PyVistaWidget.
        group (dict): group dict.
        name (str): model name.
        cell_info (dict): cell info. dict.
        site_dict (dict): site dict.
        bond_dict (dict): bond dict.
        prop (dict): property dict.
    """
    # setting.
    crystal = group.info.crystal
    cell = cell_info["cell"]
    max_neighbor = prop["max_neighbor"]
    cell_mode = prop["cell_mode"]
    mode = prop["mode"]
    scale = prop["scale"]
    radius = prop["site_radius"]
    width = prop["bond_width"]
    show_rep_site = prop["rep_site"]
    show_rep_bond = prop["rep_bond"]
    view = prop["view"]

    if cell_mode is None:
        cell_mode = "off" if group.is_point_group else "single"
    if scale is None:  # middle scale of [a,b,c].
        scale = list(sorted([cell["a"], cell["b"], cell["c"]]))[1]

    qtdraw.clear_data()
    qtdraw.set_model(name)
    qtdraw.set_crystal(crystal)
    qtdraw.set_unit_cell(cell)
    qtdraw.set_clip(False)
    qtdraw.set_cell(cell_mode)

    qtdraw_site(qtdraw, site_dict, scale, mode, radius, show_rep_site)
    qtdraw_bond(qtdraw, bond_dict, max_neighbor, scale, mode, width, show_rep_bond)

    qtdraw.set_view(view)
    qtdraw.mp_set_group(group=str(group))


# ==================================================
def create_atomic_samb_qtdraw(qtdraw, mm, name):
    """
    Create atomic SAMB QtDraw file.

    Args:
        qtdraw (PyVistaWidget): PyVistaWidget.
        mm (MaterialModel): material model.
        name (str): model name.
    """
    qtdraw.clear_data()
    qtdraw.set_model(name)
    qtdraw.set_cell("off")
    for xn in mm["atomic_id"]:
        mm.plot_atomic_samb(qtdraw, atomic_id=xn, label=False)


# ==================================================
def create_cluster_samb_qtdraw(qtdraw, mm, site_bond, name):
    """
    Create cluster SAMB QtDraw file.

    Args:
        qtdraw (PyVistaWidget): PyVistaWidget.
        mm (MaterialModel): material model.
        site_bond (str or list): (list of) site or bond.
        name (str): model name.
    """
    qtdraw.clear_data()
    qtdraw.set_model(name)
    qtdraw.set_crystal(mm["crystal"])
    if mm.group.is_point_group:
        qtdraw.set_cell("off")
    else:
        qtdraw.set_cell("single")

    for sb in site_bond:
        wp0 = mm["wyckoff"][sb]
        lst = [yn for yn, (wp, idx, comp) in mm["cluster_id"].items() if wp == wp0]
        for yn in lst:
            mm.plot_cluster_samb(qtdraw, sb, cluster_id=yn, label=False)
