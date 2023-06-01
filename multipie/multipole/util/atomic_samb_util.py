"""
This file provides utility functions for calculation of atomic multipole basis set.
"""
import numpy as np

from gcoreutils.nsarray import NSArray

from multipie.multipole.util.atomic_orbital_util import rank, basis_type
from multipie.multipole.util.spin_orbital_basis import _standard_basis
from multipie.multipole.util.multipole_util import matrix_sum
from multipie.tag.tag_list import TagList
from multipie.tag.tag_multipole import TagMultipole
from multipie.const import __def_dict__


# ==================================================
def _extract_subspace(am_data, bra_list, ket_list, spinful, crystal):
    """
    extract subspace for given bra_list and ket_list.
    32✖️32 (spinful) or 16✖️16 (spinless) matrix Xlmsk => Xlmsk[bra_list, ket_list].

    Args:
        am_data (dict): data of atomic multipoles, { TagMultipole: Matrix (32✖️32/16✖️16)}.
        bra_list ([str]): orbital list of bra.
        ket_list ([str]): orbital list of ket.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        dict: atomic multipoles (subspace).
    """
    bra_rank_lst = [rank(o, spinful, crystal) for o in bra_list]
    ket_rank_lst = [rank(o, spinful, crystal) for o in ket_list]
    max_rank = max(bra_rank_lst) + max(ket_rank_lst)
    b_type = basis_type(bra_list, crystal)

    orbitals_dataset = _standard_basis[spinful][b_type]
    bra_idx_list = [orbitals_dataset.index(o) for o in bra_list]
    ket_idx_list = [orbitals_dataset.index(o) for o in ket_list]

    am_extracted = {}
    for tag, Xlm in am_data.items():
        if tag.rank > (max_rank + 1):
            continue

        Xlm_extracted = Xlm[bra_idx_list, :][:, ket_idx_list]

        if not np.all(Xlm_extracted == 0):
            am_extracted[TagMultipole(str(tag))] = Xlm_extracted

    return am_extracted


# ==================================================
def _unitary_transform(am_lm, hs):
    """
    unitary transform the atomic multipoles.
    XlΓnγ = Σ_{lm} U_{lm, lΓnγ} Xlm

    Args:
        am_lm (dict): atomic multipoles, { TagMultipole: Xlm }.
        hs (HarmonicsPGRSet): a set of point-group harmonics (real version).

    Returns:
        dict: transformed atomic multipoles (XlΓnγ).
    """
    shape = list(am_lm.values())[0].shape
    zm = NSArray.zeros(shape=shape, style="matrix", fmt="sympy")

    lXsk_list = [(t.l, t.head, t.s, t.k) for t in am_lm.keys()]
    lXsk_list = sorted(set(lXsk_list), key=lXsk_list.index)

    am_transformed = {}
    for l, X, s, k in lXsk_list:
        Xlm_list = []
        for m in reversed(range(-l, l + 1)):
            tag = TagMultipole.create_spherical(head=X, rank=l, mul=0, comp=m, s=s, k=k)
            Xlm_list.append(am_lm[tag]) if tag in am_lm else Xlm_list.append(zm)

        U = {h.tag: h.u_matrix() for h in hs.select(rank=l) if h.tag.i_type == __def_dict__["head_i"][X]}
        am_lg = {}
        for utag, u in U.items():
            Xlg = matrix_sum([Xlm_list[j] * u[j] for j in range(2 * l + 1)])
            if not np.all(Xlg == 0):
                tag = utag.replace(head=X, s=s, k=k, m_type="a")
                am_lg[tag] = Xlg.simplify()

        am_transformed.update(am_lg)

    return am_transformed


# ==================================================
def _orthogonalize(am_set, bra_list, ket_list, crystal):
    """
    orthogonalize atomic multipoles.

    Args:
        am_set (dict): dict of atomic multipoles in each subspace (bra_list, ket_list), { TagAtomicMultipole: NSArray }.
        bra_list ([str]): orbital list of bra.
        ket_list ([str]): orbital list of ket.
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        dict: orthogonalized multipoles.
    """

    def _ortho(**kwargs):
        t_type = kwargs.pop("t_type", None)
        tags = TagList(am_set.keys()).select(**kwargs)
        if t_type is None:
            am_set_ = {tag: m for tag, m in am_set.items() if tag in tags}
        else:
            am_set_ = {tag: m for tag, m in am_set.items() if tag in tags and tag.t_type == t_type}

        tags, mats = list(am_set_.keys()), list(am_set_.values())
        mats = NSArray(mats, style="matrix", fmt="sympy", real=False)

        if len(tags) > 0:
            mats, idx = NSArray.orthogonalize(mats)
            tags, mats = [tags[j] for j in idx], [mats[j] for j in idx]

        am_orthogonalized = {tag: NSArray(mat, style="matrix", fmt="sympy", real=False) for tag, mat in zip(tags, mats)}

        return am_orthogonalized

    b_type = basis_type(bra_list, crystal)
    if b_type == "jm":
        tags = [(tag.t_type, tag.irrep) for tag in am_set.keys()]
        tags = sorted(set(tags), key=tags.index)
        kwargs_lst = [{"t_type": t_type, "irrep": irrep} for t_type, irrep in tags]
    else:
        tags = [(tag.t_type, tag.s, tag.irrep) for tag in am_set.keys()]
        tags = sorted(set(tags), key=tags.index)
        kwargs_lst = [{"t_type": t_type, "s": s, "irrep": irrep} for t_type, s, irrep in tags]

    sub = [_ortho(**kwargs) for kwargs in kwargs_lst]
    am_set = {}
    for d in sub:
        am_set |= d

    return am_set


# ==================================================
def create_atomic_samb(bra_list, ket_list, spinful, crystal, bam, hs=None, ortho=True):
    """
    create atomic multipole basis set.

    Args:
        bra_list ([str]): orbital list of bra.
        ket_list ([str]): orbital list of ket.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.
        bam (BaseAtomicMultipoleDataset): atomic-multipole dataset for all four basis.
        hs (HarmonicsPG, optional): a set of point-group harmonics.
        parallel (bool, optional): use parallel code ?
        verbose (bool, optional): verbose parallel info ?
        ortho (bool, optional): orthogonalize ?

    Returns:
        dict: atomic SAMB, {TagMultipole: NSArray(matrix)}.
    """
    # if not spinful: 32✖️32 => 16✖️16
    U_rows = U_cols = [i for i in range(32) if i % 2 == 0]
    if not spinful:
        am_data = {tag: m[U_rows, :][:, U_cols] for tag, m in bam.items() if tag.s == 0}
    else:
        am_data = bam

    # extract subspace
    am_lm_set = _extract_subspace(am_data, bra_list, ket_list, spinful, crystal)

    # unitary transform
    if hs is None:
        am_set = am_lm_set
    else:
        am_set = _unitary_transform(am_lm_set, hs)

    # orthogonalization
    if ortho:
        am_set = _orthogonalize(am_set, bra_list, ket_list, crystal)

    return am_set
