"""
This file provides utility functions for symmetry adapted model.
"""
import sympy as sp
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

from gcoreutils.nsarray import NSArray

from multipie.tag.tag_list import TagList
from multipie.multipole.util.structure_samb_util import orthogonalize_fk, decompose_fk


# number of cpu cores
_cpu_num = multiprocessing.cpu_count()


# ==================================================
def create_atomic_samb_set(pg, atomic_braket, spinful, is_phonon, parallel=True):
    """
    create atomic multipole basis set for each atomic subspaces.

    Args:
        pg (PointGroup): point group.
        atomic_braket (dict): { "M_#": (bra_list, ket_list) }.
        spinful (bool): spinful ?
        is_phonon (bool): phonon system ?
        parallel (bool, optional): use parallel code ?

    Returns:
        tuple: information of atomic multipoles, (atomic_info, atomic_data).

    Note:
        atomic_info = { "M_#": ["amp_#"] }
        atomic_data = { "amp_#" : (TagMultipole, NSArray(matrix)) }
    """
    M_num = len(atomic_braket)
    n_jobs = max(abs(min(M_num + 1, _cpu_num - 2)), 1) if parallel else 1

    def proc(i, M_i):
        bra_list, ket_list = atomic_braket[M_i]
        atomic_samb = pg.atomic_samb(bra_list, ket_list, spinful)
        if is_phonon:
            atomic_samb = {tag: m for tag, m in atomic_samb.items() if tag.head == "Q"}
        return i, M_i, atomic_samb

    M_i_list = list(atomic_braket.keys())
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc)(i, M_i) for i, M_i in enumerate(M_i_list))
    res.sort(key=lambda x: x[0])

    atomic_info = {}
    atomic_data = {}
    i = 1
    for _, M_i, atomic_samb in res:
        lst = []
        for tag, m in atomic_samb.items():
            lst.append(f"amp_{i:03d}")
            atomic_data[f"amp_{i:03d}"] = (tag, m)
            i += 1
        atomic_info[M_i] = lst

    return atomic_info, atomic_data


# ==================================================
def create_site_cluster_samb_set(g, rep_site_dict, parallel=True):
    """
    create site-cluster multipole basis set for each S_#.

    Args:
        g (PointGroup/SpaceGroup): point/space group.
        rep_site_dict (dict): { cluster_tag: position }.
        parallel (bool, optional): use parallel code ?

    Returns:
        dict: { "S_#": (TagMultipole, NSArray(vector)) }
    """
    S_num = len(rep_site_dict)
    n_jobs = max(abs(min(S_num + 1, _cpu_num - 2)), 1) if parallel else 1

    def proc(i, S_i):
        rep_site = rep_site_dict[S_i]
        site_cluster_samb = g.site_cluster_samb(rep_site)[0]
        return i, S_i, site_cluster_samb

    S_i_list = list(rep_site_dict.keys())
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc)(i, S_i) for i, S_i in enumerate(S_i_list))
    res.sort(key=lambda x: x[0])

    scm_set = {S_i: site_cluster_samb for _, S_i, site_cluster_samb in res}

    return scm_set


# ==================================================
def create_bond_cluster_samb_set(g, rep_bond_dict, parallel=True):
    """
    create bond-cluster multipole basis set for each B_#.

    Args:
        g (PointGroup/SpaceGroup): point/space group.
        rep_bond_dict (dict): { cluster_tag: vector@center }.
        parallel (bool, optional): use parallel code ?

    Returns:
        dict: { "B_#": (TagMultipole, NSArray(vector)) }
    """
    B_num = len(rep_bond_dict)
    n_jobs = max(abs(min(B_num + 1, _cpu_num - 2)), 1) if parallel else 1

    def proc(i, B_i):
        rep_bond = rep_bond_dict[B_i]
        bond_cluster_samb = g.bond_cluster_samb(rep_bond)[0]
        return i, B_i, bond_cluster_samb

    B_i_list = list(rep_bond_dict.keys())
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc)(i, B_i) for i, B_i in enumerate(B_i_list))
    res.sort(key=lambda x: x[0])

    bcm_set = {B_i: bond_cluster_samb for _, B_i, bond_cluster_samb in res}

    return bcm_set


# ==================================================
def create_cluster_samb_set(g, rep_site_dict, rep_bond_dict, parallel=True):
    """
    create site-cluster multipole basis set for each S_#.

    Args:
        g (PointGroup/SpaceGroup): point/space group.
        rep_site_dict (dict): { cluster_tag: position }.
        rep_bond_dict (dict): { cluster_tag: vector@center }.
        parallel (bool, optional): use parallel code ?

    Returns:
        tuple: information of site/bond-cluster multipoles,
                (site_cluster_info, site_cluster_data, bond_cluster_info, bond_cluster_data).

    Note:
        site_cluster_info = { "S_#": ["smp_#"] }
        site_cluster_data = { "smp_#" : (TagMultipole, NSArray(vector)) ] }
        bond_cluster_info = { "B_#": ["bmp_#"] }
        bond_cluster_data = { "bmp_#" : (TagMultipole, NSArray(vector)) ] }
    """
    scm_set = create_site_cluster_samb_set(g, rep_site_dict, parallel)
    bcm_set = create_bond_cluster_samb_set(g, rep_bond_dict, parallel)

    site_cluster_info = {}
    site_cluster_data = {}
    i = 1
    for S_i, site_cluster_samb in scm_set.items():
        lst = []
        for tag, v in site_cluster_samb.items():
            lst.append(f"smp_{i:03d}")
            site_cluster_data[f"smp_{i:03d}"] = (tag, v)
            i += 1
        site_cluster_info[S_i] = lst

    bond_cluster_info = {}
    bond_cluster_data = {}
    for B_i, bond_cluster_samb in bcm_set.items():
        lst = []
        for tag, v in bond_cluster_samb.items():
            lst.append(f"bmp_{i:03d}")
            bond_cluster_data[f"bmp_{i:03d}"] = (tag, v)
            i += 1
        bond_cluster_info[B_i] = lst

    return site_cluster_info, site_cluster_data, bond_cluster_info, bond_cluster_data


# ==================================================
def uniform_samb(sbc_samb, braket_indexes, dim):
    """
    create uniform multipole basis set.

    Args:
        sbc_samb (list): site/bond-cluster multipole basis set, [(TagMultipole, NSArray(vector))].
        braket_indexes (list): [(bra_site_no, ket_site_no)]
        dim (int): dimension of matrix.

    Returns:
        dict: {TagMultipole: NSArray(matrix)}.
    """
    u_samb = {}
    for tag, v in sbc_samb:
        U = sp.Matrix.zeros(dim, dim)
        for vi, (bra_idx, ket_idx) in zip(v, braket_indexes):
            U[bra_idx, ket_idx] += vi

        if not U.is_diagonal():
            U = (U + U.adjoint()) / sp.sqrt(2)

        U = NSArray(str(U.tolist()), style="matrix", fmt="sympy").simplify()
        if not np.all(U == 0):
            u_samb[tag] = U

    return u_samb


# ==================================================
def create_uniform_samb_set(cluster_samb_set, braket_indexes_dict, dim, parallel=True):
    """
    create uniform multipole basis set.

    Args:
        cluster_samb_set (dict): { "S_#"/"B_#": [(TagMultipole, NSArray)] }.
        braket_indexes_dict (dict): { cluster_tag: [(bra_site_no, ket_site_no)] }.
        dim (int): dimension of matrix.
        parallel (bool, optional): use parallel code ?

    Returns:
        tuple: information of uniform multipoles, uniform_info, uniform_data.

    Note:
        uniform_info = { "S_#"/"B_#": ["ump_#"] }
        uniform_data = { "ump_#" : (TagMultipole, NSArray(matrix)) }
    """
    SB_num = len(cluster_samb_set)
    n_jobs = max(abs(min(SB_num + 1, _cpu_num - 2)), 1) if parallel else 1

    def proc(i, SB_i, braket_indexes):
        sbc_samb = cluster_samb_set[SB_i]
        is_diagonal = all([bra_idx == ket_idx for (bra_idx, ket_idx) in braket_indexes])
        head_list = ["Q"] if is_diagonal else ["Q", "T"]
        sbc_samb = [(tag, v) for tag, v in sbc_samb if tag.head in head_list]
        if len(sbc_samb) == 0:
            return i, SB_i, {}

        u_samb = uniform_samb(sbc_samb, braket_indexes, dim)
        return i, SB_i, u_samb

    res = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(proc)(i, SB_i, braket_indexes) for i, (SB_i, braket_indexes) in enumerate(braket_indexes_dict.items())
    )
    res.sort(key=lambda x: x[0])

    dic = {SB_i: u_samb for _, SB_i, u_samb in res}

    # orthogonalization
    def proc_ortho(i, head):
        tags = [(SB_i, tag) for SB_i, u_samb in dic.items() for tag in u_samb.keys() if tag.head == head]
        mats = [dic[SB_i][tag] for (SB_i, tag) in tags]
        if len(mats) == 0:
            return i, {}

        mats = NSArray(mats, style="matrix", fmt="sympy", real=False)
        mats, idx = NSArray.orthogonalize(mats)
        tags, mats = [tags[j] for j in idx], [mats[j] for j in idx]
        d = {}
        for (SB_i, tag), U in zip(tags, mats):
            if sp.Matrix(U).is_diagonal():
                tag = tag.replace(m_type="s")
            else:
                tag = tag.replace(m_type="u")
            d[(SB_i, tag)] = U

        um_orthogonalized = {
            (SB_i, tag): NSArray(mat, style="matrix", fmt="sympy", real=False) for (SB_i, tag), mat in d.items()
        }
        return i, um_orthogonalized

    head_list = ["Q", "T"]
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc_ortho)(i, head) for i, head in enumerate(head_list))
    res.sort(key=lambda x: x[0])

    uniform_info = {}
    uniform_data = {}
    i = 1
    for _, u_samb in res:
        for (SB_i, tag), v in u_samb.items():
            uniform_data[f"ump_{i:03d}"] = (tag, v)
            if SB_i in uniform_info:
                uniform_info[SB_i] += [f"ump_{i:03d}"]
            else:
                uniform_info[SB_i] = [f"ump_{i:03d}"]
            i += 1

    return uniform_info, uniform_data


# ==================================================
def structure_samb(bc_samb, bond_list, bond):
    """
    create structure multipole basis set.

    Args:
        bc_samb (list): bond-cluster multipole basis set, [(TagMultipole, NSArray(vector))].
        bond_list (list) : ["bond_#"].
        bond (dict, optional): { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }.

    Returns:
        dict: {TagMultipole: Symbol}.
    """
    n_sgn_list = []
    for bond_n in bond_list:
        vec = bond[bond_n][3]
        if "bond" in vec:
            n = int(vec.split("_")[1])
        else:
            n = int(bond_n.split("_")[1])

        if vec[0] == "-":
            n_sgn_list.append((n, +1))
        else:
            n_sgn_list.append((n, -1))

    fk_list = []
    tags = []
    for tag, v in bc_samb:
        fk = sp.S(0)
        for vn, (n, sgn) in zip(v, n_sgn_list):
            cn = sp.Symbol(f"c{n:03d}", real=True)
            sn = sp.Symbol(f"s{n:03d}", real=True)
            v = vn * (cn + sgn * sp.I * sn)
            fk += v + sp.conjugate(v)
        tags.append(tag)
        fk_list.append(sp.expand(fk))

    # orthogonalization
    fk_list, idx = orthogonalize_fk(fk_list)
    fk_list = [fk_list[i] for i in idx]
    tags = [tags[i] for i in idx]

    k_samb = {tag.replace(m_type="k"): fk for tag, fk in zip(tags, fk_list)}

    return k_samb


# ==================================================
def create_structure_samb_set(bc_samb_set, cluster_bond, bond, parallel=True):
    """
    create structure multipole basis set.

    Args:
        bc_samb_set (dict): { "B_#": [(TagMultipole, NSArray)] }.
        cluster_bond (dict, optional) : { cluster_tag: bond_list }.
        bond (dict, optional): { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }.
        parallel (bool, optional): use parallel code ?

    Returns:
        tuple: information of structure multipoles, structure_info, structure_data.

    Note:
        structure_info = { "B_#": ["kmp_#"] }
        structure_data = { "kmp_#" : (TagMultipole, Symbol) }
    """
    B_num = len(bc_samb_set)
    n_jobs = max(abs(min(B_num + 1, _cpu_num - 2)), 1) if parallel else 1

    def proc(i, B_i):
        bc_samb = bc_samb_set[B_i]
        bond_list = cluster_bond[B_i]
        k_samb = structure_samb(bc_samb, bond_list, bond)
        return (i, B_i, k_samb)

    B_i_list = list(cluster_bond.keys())
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc)(i, B_i) for i, B_i in enumerate(B_i_list))
    res.sort(key=lambda x: x[0])

    structure_info = {}
    structure_data = {}
    i = 1
    for _, B_i, k_samb in res:
        lst = []
        for tag, fk in k_samb.items():
            lst.append(f"kmp_{i:03d}")
            structure_data[f"kmp_{i:03d}"] = (tag, sp.expand(fk))
            i += 1
        structure_info[B_i] = lst

    return structure_info, structure_data


# ==================================================
def _check_complete_relation(x_tag_dict, y_tag_dict, z_info, z_data):
    """
    check complete relation of combined multipole basis set.

    Args:
        x_tag_dict (dict): multipole/harmonics tag dict, {M_#: TagList}.
        y_tag_dict (dict): multipole/harmonics tag dict, {S_#/B_#: TagList}.
        z_info (dict): { ("S_#"/"B_#", "M_#"): ["z_#"] }.
        z_data (dict): { "z_#" : [(coefficient, "amp_#", "smp_#"/"bmp_#"/"ump_#")] }.
    """
    for M_i, SB_i in z_info.keys():
        tag1_list = x_tag_dict[M_i]
        tag2_list = y_tag_dict[SB_i]
        Z_lst = [z_data[z_i] for z_i in z_info[(M_i, SB_i)]]
        n, n1, n2 = len(Z_lst), len(tag1_list), len(tag2_list)
        if n != n1 * n2:
            s = f"(M_i, SB_i) = {(M_i, SB_i)} \n"
            s += f"tag1_list = {[str(tag1) for tag1 in tag1_list]} \n"
            s += f"tag2_list = {[str(tag2) for tag2 in tag2_list]} \n"
            s += "tag2_list \n"
            for tag, lst in Z_lst:
                s += f"{tag} = {[(c, str(tag1), str(tag2)) for c, tag1, tag2 in lst]} \n"

            s += f"# of X = n1 = {n1} \n"
            s += f"# of Y = n2 = {n2} \n"
            s += f"# of Z = n = {n} \n"
            s += f"n1*n2 = n = {n} \n"
            s += f"n1*n2 = {n1*n2} must equal to n = {n}."
            raise Exception(s)


# ==================================================
def create_z_samb_set(
    g, x_tag_dict, y_tag_dict, M_SB_list, atomic_braket, toroidal_priority=False, parallel=True, **kwargs
):
    """
    create combined multipole basis set.

    Args:
        g (PointGroup/SpaceGroup): point/space group.
        x_tag_dict (dict): atomic multipole tag dict, {(M_#, TagMultipole): amp_i}.
        y_tag_dict (dict): site/bond-cluster tag dict, {(S_#/B_#, TagMultipole): smp_i/bmp_i/ump_i}.
        M_SB_list (list): [("M_#", "S_#"/"B_#")].
        atomic_braket (dict): { matrix_tag : (bra_list, ket_list) }.
        toroidal_priority (bool): create toroidal multipoles (G,T) in priority? else prioritize conventional multipoles (Q,M).
        parallel (bool, optional): use parallel code.
        kwargs (dict, optional): select conditions for multipoles,
                                 keywords in ["head", "rank", "irrep", "mul", "comp", "s", "k"].
    Returns:
        tuple: information of combined multipoles, z_info, z_data.

    Note:
        z_info = { ("S_#"/"B_#", "M_#"): ["z_#"] }
        z_data = { "z_#" : [(coefficient, "amp_#", "smp_#"/"bmp_#"/"ump_#")] }
    """
    E_am_dict = {}
    for M_i, SB_i in M_SB_list:
        o1, o2 = atomic_braket[M_i]
        if SB_i[0] == "B" and o1 != o2:
            E_am_dict[(M_i, SB_i)] = True
        else:
            E_am_dict[(M_i, SB_i)] = False

    M_SB_num = len(M_SB_list)
    n_jobs = max(abs(min(M_SB_num + 1, _cpu_num - 2)), 1) if parallel else 1

    irrep_list = kwargs.pop("irrep", None)
    if type(irrep_list) == str:
        irrep_list = [irrep_list]
    elif irrep_list is None:
        irrep_list = g.character.irrep_list if g.tag.pg == "" else g.pg.character.irrep_list
        irrep_list = [str(irrep) for irrep in irrep_list]

    irrep_M_SB_list = [(irrep, M_i, SB_i) for irrep in irrep_list for (M_i, SB_i) in M_SB_list]

    def proc(i, irrep, M_i, SB_i):
        x_tag_list = TagList([tag for (M_i_, tag) in x_tag_dict.keys() if M_i_ == M_i])
        if E_am_dict[(M_i, SB_i)]:
            x_tag_list = TagList([tag for tag in x_tag_list if tag.head in ("Q", "G")])
        y_tag_list = TagList([tag for (SB_i_, tag) in y_tag_dict.keys() if SB_i_ == SB_i])
        z_samb = g.z_samb(x_tag_list, y_tag_list, toroidal_priority, irrep=irrep, **kwargs)
        return i, irrep, M_i, SB_i, z_samb

    res = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(proc)(i, irrep, M_i, SB_i) for i, (irrep, M_i, SB_i) in enumerate(irrep_M_SB_list)
    )
    res.sort(key=lambda x: x[0])

    z_info = {}
    z_data = {}
    i = 1
    for _, irrep, M_i, SB_i, z_samb in res:
        lst = []
        for tags, c_xtag_ytag_list in z_samb.items():
            tag, _, _ = tags
            lst.append(f"z_{i:03d}")
            c_amp_cmp_list = [
                (c, x_tag_dict[(M_i, xtag)], y_tag_dict[(SB_i, ytag)]) for c, xtag, ytag in c_xtag_ytag_list
            ]
            z_data[f"z_{i:03d}"] = (tag, c_amp_cmp_list)
            i += 1
        if lst != []:
            z_info[(irrep, M_i, SB_i)] = lst

    # check basis (# of X)×(# of Y) = # of Z.
    if len(kwargs) == 0:
        _check_complete_relation(x_tag_dict, y_tag_dict, z_info, z_data)

    return z_info, z_data


# ==================================================
def _fourier_transform_bond_cluster_samb(bc_samb, u_samb_set, k_samb_set, braket_indexes, bond_list, bond, dim):
    """
    fourier series of bond-cluster multipoles.

    Yb_{a} => Σ_{bc} C_{abc} U_{b} fk_{c}

    Args:
        bc_samb (list): bond-cluster multipole basis set, [("bmp_#", NSArray(vector))].
        u_samb_set (list): uniform multipole basis set, [("ump_#", NSArray(matrix))].
        k_samb_set (list): structure multipole basis set, [("kmp_#", Symbol)].
        braket_indexes (list): [(bra_site_no, ket_site_no)].
        bond_list (list) : ["bond_#"].
        bond (dict, optional): { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }.
        dim (int): dimension of site-space matrix.

    Returns:
        dict: fourier series of bond-cluster multipoles, { "bmp_#" : [(coefficient, "ump_#", "kmp_#")] }.
    """

    def _check_fourier_transform(coeff_list):
        v = sp.Matrix(coeff_list)
        v = sp.adjoint(v).dot(v)
        norm = sp.sqrt(sp.simplify(v))

        if norm != 1:
            raise Exception(f"invalid fourier expansion coefficient = {v}, norm = {norm}.")

    n_sgn_list = []
    for bond_n in bond_list:
        vec = bond[bond_n][3]
        if "bond" in vec:
            n = int(vec.split("_")[1])
        else:
            n = int(bond_n.split("_")[1])

        if vec[0] == "-":
            n_sgn_list.append((n, +1))
        else:
            n_sgn_list.append((n, -1))

    Mfk_set = {}
    for bmp_i, v in bc_samb:
        Mk = sp.Matrix.zeros(dim, dim)
        for vi, (bra_idx, ket_idx), (n, sgn) in zip(v, braket_indexes, n_sgn_list):
            cn = sp.Symbol(f"c{n:03d}", real=True)
            sn = sp.Symbol(f"s{n:03d}", real=True)
            Mk[bra_idx, ket_idx] += vi * (cn + sgn * sp.I * sn)

        Mk = (Mk + Mk.adjoint()) / sp.sqrt(2)
        Mfk_set[bmp_i] = sp.simplify(Mk)

    bck_samb = {}
    for bmp_i, Mk in Mfk_set.items():
        coeff_ump_smp_list = []
        for ump_i, m in u_samb_set:
            fk = (sp.Matrix(m).adjoint() * Mk).trace()
            fk = sp.simplify(fk)
            coeffs = decompose_fk(fk, k_samb_set)
            for kmp_i, c in coeffs.items():
                coeff_ump_smp_list.append((c, ump_i, kmp_i))

        # check normalization
        coeff_list = [c for c, _, _ in coeff_ump_smp_list]
        _check_fourier_transform(coeff_list)

        bck_samb[bmp_i] = coeff_ump_smp_list

    return bck_samb


# ==================================================
def fourier_transform_bond_cluster_samb(bc_samb_set, u_samb_set, k_samb_set, cluster_bond, bond, dim, parallel=True):
    """
    fourier series of bond-cluster multipoles.

    Yb_{a} => Σ_{bc} C_{abc} U_{b} fk_{c}

    Args:
        bc_samb_set (dict): { "B_#": [("bmp_#", NSArray(vector))] }.
        u_samb_set (list): uniform multipole basis set, [("ump_#", NSArray(matrix))].
        k_samb_set (list): structure multipole basis set, [("kmp_#", Symbol)].
        cluster_bond (dict) : { cluster_tag: bond_list }.
        bond (dict): { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }.
        dim (int): dimension of matrix.
        parallel (bool, optional): use parallel code ?

    Returns:
        dict: { "bmp_#" : [(coefficient, "ump_#", "kmp_#")] }.
    """
    B_num = len(bc_samb_set)
    n_jobs = max(abs(min(B_num + 1, _cpu_num - 2)), 1) if parallel else 1

    braket_indexes_dict = {
        B_i: [bond[bond_tag][2] for bond_tag in bond_list] for B_i, bond_list in cluster_bond.items()
    }

    def proc(i, B_i):
        bc_samb = bc_samb_set[B_i]
        braket_indexes = braket_indexes_dict[B_i]
        bond_list = cluster_bond[B_i]
        bck_samb = _fourier_transform_bond_cluster_samb(
            bc_samb, u_samb_set, k_samb_set, braket_indexes, bond_list, bond, dim
        )
        return i, bck_samb

    B_i_list = list(cluster_bond.keys())
    res = Parallel(n_jobs=n_jobs, verbose=0)(delayed(proc)(i, B_i) for i, B_i in enumerate(B_i_list))
    res.sort(key=lambda x: x[0])

    bond_cluster_k_data = {}
    for _, bck_samb in res:
        bond_cluster_k_data |= bck_samb

    return bond_cluster_k_data


# ==================================================
def create_zk_samb_set(z_data, bc_samb_set, u_samb_set, k_samb_set, cluster_bond, bond, dim, parallel=True):
    """
    create k-space representation of combined multipole basis set.

    Zk = Σ_{abc} C_{abc} X_{a} U_{b} fk_{c}

    Args:
        z_data (dict): { "z_#": (TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#")]) }.
        bc_samb_set (dict): { "B_#": [("bmp_#", NSArray(vector))] }.
        u_samb_set (list): uniform multipole basis set, [("ump_#", NSArray(matrix))].
        k_samb_set (list): structure multipole basis set, [("kmp_#", Symbol)].
        cluster_bond (dict) : { cluster_tag: bond_list }.
        bond (dict): { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }.
        dim (int): dimension of matrix.
        parallel (bool, optional): use parallel code ?

    Returns:
        dict: { "z_#": (TagMultipole, [ (coeff, "amp_#", "ump_#", "kmp_#") ]) }.
    """
    bond_cluster_k_data = fourier_transform_bond_cluster_samb(
        bc_samb_set, u_samb_set, k_samb_set, cluster_bond, bond, dim, parallel
    )

    Zk_data = {}
    for z_i, (tag, lst) in z_data.items():
        lst_k = []
        for c1, amp_i, cmp_i in lst:
            if cmp_i[0] == "b":
                for c2, ump_i, kmp_i in bond_cluster_k_data[cmp_i]:
                    c = sp.simplify(sp.sympify(c1) * sp.sympify(c2))
                    lst_k.append((c, amp_i, ump_i, kmp_i))
            else:
                ump_i = cmp_i.replace("s", "u")
                lst_k.append((c1, amp_i, ump_i))

        Zk_data[z_i] = (tag, lst_k)

    return Zk_data


# ==================================================
def redefine_index(info, data, z_data, molecule=False, init_idx=1):
    """
    Redefine the serial number of multipoles.

    Args:
        info (dict): { "S_#"/"B_#"/"M_#": ["x_#"] }, x = amp/smp/bmp/ump/kmp.
        data (dict): { "x_#" : definition (NSArray/list) }.
        z_data (dict): { "z_#" : (TagMultipole, [(coefficient, "amp_#", "cmp_#")]) }/{ "z_#": (TagMultipole, [ (coeff, "amp_#", "ump_#", "kmp_#") ]) }.
        molecule (bool, optional): molecule ?
        init_idx (int): initial number of index.

    Returns:
        tuple: redefined information of multipoles, (info, data, z_data).
    """
    x = list(info.values())[0][0].split("_")[0]
    if molecule and x in ("smp", "bmp"):
        sb = "s" if x == "smp" else "b"
        x_i_list = [v[2].replace("u", sb) for lst in z_data.values() for v in lst[1] if v[2].split("_")[0] == "ump"]
    else:
        x_i_list = [x_i for lst in z_data.values() for v in lst[1] for x_i in v[1:] if x_i.split("_")[0] == x]
    x_i_list = sorted(list(set(x_i_list)))

    old_new_dict = {}
    new_info = {}
    new_data = {}
    i = init_idx
    for k, lst in info.items():
        new_lst = []
        for x_i_old in lst:
            if x_i_old in x_i_list:
                x_i_new = f"{x}_{i:03d}"
                old_new_dict[x_i_old] = x_i_new
                new_lst.append(x_i_new)
                new_data[x_i_new] = data[x_i_old]
                i += 1
        if len(new_lst) > 0:
            new_info[k] = new_lst

    new_z_data = {}
    for k, lst in z_data.items():
        tag = lst[0]
        new_lst = []
        for v in lst[1]:
            c = v[0]
            new_v = []
            for x_i_old in v[1:]:
                if x_i_old in old_new_dict:
                    x_i_new = old_new_dict[x_i_old]
                else:
                    x_i_new = x_i_old
                new_v.append(x_i_new)
            new_lst.append((c, *new_v))
        new_z_data[k] = (tag, new_lst)

    return new_info, new_data, new_z_data
