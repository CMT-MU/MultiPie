"""
For harmonics.
"""

import numpy as np

from multipie.core.group import Group


# ==================================================
def harmonics_decomposition(basis_pg, pg, rank, head="Q"):
    """
    Harmonics decomposition in terms of basis PG.

    Args:
        basis_pg (str): basis PG.
        pg (str): PG expressed by basis PG.
        rank (int): rank.
        head (str, optional): type, Q/G ?

    Returns:
        list: decomposition info., [(harm_tag, [(coeff,basis)].
    """
    pg = Group(pg)
    basis_pg = Group(basis_pg)

    basis_hset = basis_pg.harmonics.select(l=rank, X=head)
    hset = pg.harmonics.select(l=rank, X=head)

    decomp = []
    for h, v in hset.items():
        for comp, U in enumerate(v[1].T):
            tag = Group.tag_multipole(h, comp, latex=True)

            c = []
            for oh, ov in basis_hset.items():
                for ocomp, basis_U in enumerate(ov[1].T):
                    basis_tag = Group.tag_multipole(oh, ocomp, latex=True)
                    ci = np.vdot(basis_U, U)
                    if ci != 0:
                        c.append(((ci, basis_tag)))
            decomp.append((tag, c))

    return decomp
