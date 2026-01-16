"""
Utility for SAMB.
"""

import numpy as np
import sympy as sp

from multipie import PGMultipoleType
from multipie.util.util_dict import Dict
from multipie.util.util_gram_schmidt import gram_schmidt


# ==================================================
def orthogonalize_samb(samb, diagonal_block):
    """
    Orthogonalize SAMB.

    Args:
        samb (Dict): SAMB dict, Dict[idx, (expression, SAMB)].
        diagonal_block (bool): diagonal block ?

    Returns:
        - (Dict) -- orthogonalized SAMB.
    """
    if len(samb) == 0:
        return samb

    no_idx = [(k, ex) for k, (m, ex) in samb.items()]
    no_basis = np.asarray([m for m, ex in samb.values()])

    sh = list(no_basis.shape)
    no_basis = no_basis.reshape(len(no_idx), -1)

    df = 1 if diagonal_block else 1 / sp.sqrt(2)
    o_basis, nonzero = gram_schmidt(no_basis)
    o_basis *= sp.sqrt(sh[1]) * df  # multiply dimension of irrep. and off-diagonal factor.

    sh[0] = len(nonzero)

    o_basis = o_basis.reshape(tuple(sh))
    o_idx = [no_idx[i] for i in nonzero]

    samb = Dict(samb.key_type, {k: (v, ex) for (k, ex), v in zip(o_idx, o_basis)})

    # remove zero matrices
    samb = Dict(samb.key_type, {tag: (mat, ex) for tag, (mat, ex) in samb.items() if not np.all(mat == 0)})

    return samb


# ==================================================
def orthogonalize_multipole(samb, diagonal_block):
    """
    Orthogonalize mutlipole SAMB.

    Args:
        samb (Dict): SAMB dict, Dict[idx, (expression, SAMB)].
        diagonal_block (bool): diagonal block ?

    Returns:
        - (Dict) -- orthogonalized SAMB.
    """
    irreps = set([i.Gamma for i in samb.named_keys()])

    samb_ortho = Dict(PGMultipoleType)
    for X in [["Q", "G"], ["M", "T"]]:
        for Gamma in irreps:
            samb1 = samb.select(X=X, Gamma=Gamma)
            samb1 = samb1.sort(("X", X), "s", "k", "l", "n", "p")
            samb_ortho.update(orthogonalize_samb(samb1, diagonal_block))

    samb_ortho = samb_ortho.sort(("X", ["Q", "G", "M", "T"]), "s", "k", "Gamma", "l", "n", "p")

    return samb_ortho
