"""
For SAMB creation.
"""

import numpy as np

from multipie import PGMultipoleType
from multipie.util.util_dict import Dict
from multipie.util.util_gram_schmidt import gram_schmidt, gram_schmidt_float


# ==================================================
def expand_component(basis, vector):
    ex_basis = {}
    for idx, v in basis.items():
        if len(v[0]) == 1:
            if vector:
                ex_basis[(idx, -1)] = (v[0][0], v[1][0], v[2][0])
            else:
                ex_basis[(idx, -1)] = (v[0][0], v[1][0])
        else:
            if vector:
                for no, (b, ex, ex1) in enumerate(zip(v[0], v[1], v[2])):
                    ex_basis[(idx, no)] = (b, ex, ex1)
            else:
                for no, (b, ex) in enumerate(zip(v[0], v[1])):
                    ex_basis[(idx, no)] = (b, ex)

    return ex_basis


# ==================================================
def gather_component(ex_basis, vector):
    dic = {}
    for (idx, no), v in ex_basis.items():
        dic[idx] = dic.get(idx, []) + [v]
    nc = 3 if vector else 2
    dic1 = {}
    for idx, v in dic.items():
        dic1[idx] = tuple([np.asarray([i[c] for i in v]) for c in range(nc)])

    basis = Dict(PGMultipoleType)
    for idx, v in dic1.items():
        basis[idx] = v

    return basis


# ==================================================
def orthogonalize(samb, n):  # real cluster samb only.
    no_idx = [(idx, v[1:]) for idx, v in samb.items()]
    no_basis = np.asarray([v[0] for v in samb.values()])
    sh = list(no_basis.shape)
    no_basis = no_basis.reshape(len(no_idx), -1)

    # first orthogonalize numerically, then independent bases are orthogonalized symbolically.
    _, nonzero = gram_schmidt_float(no_basis, n)
    o_idx = [no_idx[i] for i in nonzero]
    no_basis = [no_basis[i] for i in nonzero]
    o_basis, _ = gram_schmidt(no_basis, n)

    sh[0] = len(nonzero)
    o_basis = o_basis.reshape(tuple(sh))
    samb_ortho = {idx: (b, *ex) for (idx, ex), b in zip(o_idx, o_basis)}

    return samb_ortho
