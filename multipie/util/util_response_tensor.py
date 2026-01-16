"""
Response tensors up to 4th rank.
"""

import re
import numpy as np
import sympy as sp
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy import LeviCivita
from itertools import product


# ==================================================
def delta(i1, i2):
    """
    Kronecker delta.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2. 1-3.

    Returns:
        - (sympy) -- delta(i1,i2).
    """
    return KroneckerDelta(i1, i2)


# ==================================================
def epsilon(i1, i2, i3):
    """
    Levi-Civita epsilon.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2. 1-3.
        i3 (int): index 3, 1-3.

    Returns:
        - (sympy) -- epsilon(i1,i2,i3).
    """
    return LeviCivita(i1, i2, i3)


# ==================================================
def d1234(i1, i2, i3, i4):
    return delta(i1, i3) * delta(i2, i4) + delta(i1, i4) * delta(i2, i3)


# ==================================================
def d12345(i1, i2, i3, i4, i5):
    return (
        delta(i1, i3) * epsilon(i2, i4, i5)
        + delta(i2, i3) * epsilon(i1, i4, i5)
        + delta(i1, i4) * epsilon(i2, i3, i5)
        + delta(i2, i4) * epsilon(i1, i3, i5)
    )


# ==================================================
def d123456(i1, i2, i3, i4, i5, i6):
    return sp.S(delta(i1, i2) * d1234(i3, i4, i5, i6) + delta(i3, i4) * d1234(i1, i2, i5, i6)) / 2


# ==================================================
def d123456p(i1, i2, i3, i4, i5, i6):
    return (
        sp.S(
            delta(i1, i3) * d1234(i2, i4, i5, i6)
            + delta(i2, i3) * d1234(i1, i4, i5, i6)
            + delta(i1, i4) * d1234(i2, i3, i5, i6)
            + delta(i2, i4) * d1234(i1, i3, i5, i6)
        )
        / 2
    )


# ==================================================
def d123456pp(i1, i2, i3, i4, i5, i6):
    return sp.S(delta(i1, i2) * d1234(i3, i4, i5, i6) - delta(i3, i4) * d1234(i1, i2, i5, i6)) / 2


# ==================================================
def d1234567(i1, i2, i3, i4, i5, i6, i7):
    return (
        sp.S(
            epsilon(i1, i3, i7) * d1234(i2, i4, i5, i6)
            + epsilon(i2, i3, i7) * d1234(i1, i4, i5, i6)
            + epsilon(i2, i4, i7) * d1234(i1, i3, i5, i6)
            + epsilon(i1, i4, i7) * d1234(i2, i3, i5, i6)
        )
        / 2
    )


# ==================================================
def e1234(i1, i2, i3, i4):
    return delta(i1, i3) * delta(i2, i4) - delta(i1, i4) * delta(i2, i3)


# ==================================================
def e51234(i1, i2, i3, i4, i5):
    return (
        sp.S(
            delta(i1, i5) * epsilon(i2, i3, i4)
            - delta(i2, i5) * epsilon(i1, i3, i4)
            - delta(i3, i5) * epsilon(i1, i2, i4)
            + delta(i4, i5) * epsilon(i1, i2, i3)
        )
        / 2
    )


# ==================================================
def e561234(i1, i2, i3, i4, i5, i6):
    return epsilon(i1, i2, i5) * epsilon(i3, i4, i6) + epsilon(i1, i2, i6) * epsilon(i3, i4, i5)


# ==================================================
def g12534(i1, i2, i3, i4, i5):
    return delta(i1, i5) * epsilon(i2, i3, i4) + delta(i2, i5) * epsilon(i1, i3, i4)


# ==================================================
def g125634(i1, i2, i3, i4, i5, i6):
    return (
        sp.S(
            delta(i3, i5) * d1234(i1, i2, i4, i6)
            + delta(i3, i6) * d1234(i1, i2, i4, i5)
            - delta(i4, i5) * d1234(i1, i2, i3, i6)
            - delta(i4, i6) * d1234(i1, i2, i3, i5)
        )
        / 2
    )


# ==================================================
def mp_string(lst, head, sort=True):
    """
    Multipole string.

    Args:
        head (str): head of multipole.
        sort (bool, optional): sort subscript ?

    Returns:
        - (str) -- multipole string.
    """
    _dv = {1: "x", 2: "y", 3: "z"}

    lst = [_dv[i] for i in lst]
    if sort:
        s = head + "_{" + "".join(sorted("".join(lst))) + "}"
    else:
        s = head + "_{" + "".join("".join(lst)) + "}"
    return s


# ==================================================
def S0():
    """
    Rank 0 tensor component.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    return [(1, "Q^{(1)}")]


# ==================================================
def S1(i1):
    """
    Rank 1 tensor component.

    Args:
        i1 (int): index 1, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    return [(1, mp_string([i1], "Q^{(1)}"))]


# ==================================================
def S12(i1, i2):
    """
    Symmetric part of rank 2 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if i1 == i2:
        lst = [(1, mp_string([i1, i2], "Q^{(1)}")), (1, "Q^{(1)}")]
    else:
        lst = [(1, mp_string([i1, i2], "Q^{(1)}"))]

    return lst


# ==================================================
def A12(i1, i2):
    """
    Anti-symmetric part of rank 2 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    lst = []
    for i3 in range(1, 4):
        e = epsilon(i1, i2, i3)
        if e != 0:
            lst.append((e, mp_string([i3], "G^{(1)}")))

    return lst


# ==================================================
def S123(i1, i2, i3):
    """
    Symmetric part of rnak 3 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if i1 == i2:
        l1 = [(1, mp_string([i3], "Q^{(1)}"))]
    else:
        l1 = []

    l2 = []
    for i4 in range(1, 4):
        d = d1234(i1, i2, i3, i4)
        if d != 0:
            l2.append((d, mp_string([i4], "Q^{(2)}")))

    l3 = []
    for i4, i5 in product(range(1, 4), range(1, 4)):
        g = g12534(i1, i2, i3, i4, i5)
        if g != 0:
            l3.append((g, mp_string([i4, i5], "G^{(1)}")))

    l4 = [(1, mp_string([i1, i2, i3], "Q^{(1)}"))]

    return l1 + l2 + l3 + l4


# ==================================================
def A123(i1, i2, i3):
    """
    Anti-symmetric part of rank 3 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    e = epsilon(i1, i2, i3)
    if e != 0:
        l1 = [(e, "G^{(1)}")]
    else:
        l1 = []

    l2 = []
    for i4 in range(1, 4):
        e = e1234(i1, i2, i3, i4)
        if e != 0:
            l2.append((e, mp_string([i4], "Q^{(3)}")))

    l3 = []
    for i4 in range(1, 4):
        e = epsilon(i1, i2, i4)
        if e != 0:
            l3.append((e, mp_string([i3, i4], "G^{(2)}")))

    return l1 + l2 + l3


# ==================================================
def S1234(i1, i2, i3, i4):
    """
    SSS part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if i1 == i2 and i3 == i4:
        l1 = [(1, "Q^{(1)}")]
    else:
        l1 = []

    d = d1234(i1, i2, i3, i4)
    if d != 0:
        l2 = [(d, "Q^{(2)}")]
    else:
        l2 = []

    l3 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l3.append((d, mp_string([i5, i6], "Q^{(1)}")))

    l4 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456p(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l4.append((d, mp_string([i5, i6], "Q^{(2)}")))

    l5 = [(1, mp_string([i1, i2, i3, i4], "Q^{(1)}"))]

    return l1 + l2 + l3 + l4 + l5


# ==================================================
def Sb1234(i1, i2, i3, i4):
    """
    SSA part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    l1 = []
    for i5 in range(1, 4):
        d = d12345(i1, i2, i3, i4, i5)
        if d != 0:
            l1.append((d, mp_string([i5], "G^{(1)}")))

    l2 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456pp(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l2.append((d, mp_string([i5, i6], "Q^{(3)}")))

    l3 = []
    for i5, i6, i7 in product(range(1, 4), range(1, 4), range(1, 4)):
        d = d1234567(i1, i2, i3, i4, i5, i6, i7)
        if d != 0:
            l3.append((d, mp_string([i5, i6, i7], "G^{(1)}")))

    return l1 + l2 + l3


# ==================================================
def A1234(i1, i2, i3, i4):
    """
    AAS part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    e = e1234(i1, i2, i3, i4)
    if e != 0:
        l1 = [(e, "Q^{(3)}")]
    else:
        l1 = []

    l2 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        e = e561234(i1, i2, i3, i4, i5, i6)
        if e != 0:
            l2.append((e, mp_string([i5, i6], "Q^{(4)}")))

    return l1 + l2


# ==================================================
def Ab1234(i1, i2, i3, i4):
    """
    AAA part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    l1 = []
    for i5 in range(1, 4):
        e = e51234(i1, i2, i3, i4, i5)
        if e != 0:
            l1.append((e, mp_string([i5], "G^{(2)}")))

    return l1


# ==================================================
def M1234(i1, i2, i3, i4):
    """
    SA part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    l1 = []
    for i5 in range(1, 4):
        de = delta(i1, i2) * epsilon(i3, i4, i5)
        if de != 0:
            l1.append((de, mp_string([i5], "G^{(3)}")))

    l2 = []
    for i5 in range(1, 4):
        g = g12534(i1, i2, i3, i4, i5)
        if g != 0:
            l2.append((g, mp_string([i5], "G^{(4)}")))

    l3 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        g = g125634(i1, i2, i3, i4, i5, i6)
        if g != 0:
            l3.append((g, mp_string([i5, i6], "Q^{(5)}")))

    l4 = []
    for i5 in range(1, 4):
        e = epsilon(i3, i4, i5)
        if e != 0:
            l4.append((e, mp_string([i1, i2, i5], "G^{(2)}")))

    return l1 + l2 + l3 + l4


# ==================================================
def Mb1234(i1, i2, i3, i4):
    """
    AS part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    l1 = []
    for i5 in range(1, 4):
        de = delta(i3, i4) * epsilon(i1, i2, i5)
        if de != 0:
            l1.append((de, mp_string([i5], "G^{(5)}")))

    l2 = []
    for i5 in range(1, 4):
        g = g12534(i3, i4, i1, i2, i5)
        if g != 0:
            l2.append((g, mp_string([i5], "G^{(6)}")))

    l3 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        g = g125634(i3, i4, i1, i2, i5, i6)
        if g != 0:
            l3.append((g, mp_string([i5, i6], "Q^{(6)}")))

    l4 = []
    for i5 in range(1, 4):
        e = epsilon(i1, i2, i5)
        if e != 0:
            l4.append((e, mp_string([i3, i4, i5], "G^{(3)}")))

    return l1 + l2 + l3 + l4


# ==================================================
def T1234(i1, i2, i3, i4):
    """
    S part of rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if i1 == i2 and i3 == i4:
        d1 = 1
    else:
        d1 = 0
    d2 = d1234(i1, i2, i3, i4)
    if d1 + d2 != 0:
        l1 = [(d1 + d2, "Q^{(1)}")]
    else:
        l1 = []

    l2 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l2.append((d, mp_string([i5, i6], "Q^{(1)}")))

    l3 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456p(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l3.append((d, mp_string([i5, i6], "Q^{(2)}")))

    l4 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456pp(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l4.append((d, mp_string([i5, i6], "Q^{(3)}")))

    l5 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        g = g125634(i1, i2, i3, i4, i5, i6)
        if g != 0:
            l5.append((g, mp_string([i5, i6], "Q^{(5)}")))

    l6 = [(1, mp_string([i1, i2, i3, i4], "Q^{(1)}"))]

    return l1 + l2 + l3 + l4 + l5 + l6


# ==================================================
def P0():
    """
    Rank 0 tensor component.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    return S0()


# ==================================================
def P1(i1):
    """
    Rank 1 tensor component.

    Args:
        i1 (int): index 1, 1-3.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    return S1(i1)


# ==================================================
def P2(i1, i2, opt=None):
    """
    Rank 2 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        opt (str, optional): part, s/a.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if opt is None:
        return S12(i1, i2) + A12(i1, i2)
    elif opt == "s":
        return S12(i1, i2)
    elif opt == "a":
        return A12(i1, i2)
    else:
        raise ValueError(f"invalid option, {opt}")


# ==================================================
def P3(i1, i2, i3, opt=None):
    """
    Rank 3 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        opt (str, optional): part, s/a.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if opt is None:
        return S123(i1, i2, i3) + A123(i1, i2, i3)
    elif opt == "s":
        return S123(i1, i2, i3)
    elif opt == "a":
        return A123(i1, i2, i3)
    else:
        raise ValueError(f"invalid option, {opt}")


# ==================================================
def P4(i1, i2, i3, i4, opt=None):  # sss/ssa/aas/aaa/sa/as/s/a/t.
    """
    Rank 4 tensor component.

    Args:
        i1 (int): index 1, 1-3.
        i2 (int): index 2, 1-3.
        i3 (int): index 3, 1-3.
        i4 (int): index 4, 1-3.
        opt (str, optional): part, sss/ssa/aas/aaa/sa/as/s/a/t.

    Returns:
        - (list) -- expression list, (coeff, multipole).
    """
    if opt is None:
        return (
            S1234(i1, i2, i3, i4)
            + Sb1234(i1, i2, i3, i4)
            + A1234(i1, i2, i3, i4)
            + Ab1234(i1, i2, i3, i4)
            + M1234(i1, i2, i3, i4)
            + Mb1234(i1, i2, i3, i4)
        )
    elif opt == "sss":
        return S1234(i1, i2, i3, i4)
    elif opt == "ssa":
        return Sb1234(i1, i2, i3, i4)
    elif opt == "aas":
        return A1234(i1, i2, i3, i4)
    elif opt == "aaa":
        return Ab1234(i1, i2, i3, i4)
    elif opt == "sa":
        return M1234(i1, i2, i3, i4)
    elif opt == "as":
        return Mb1234(i1, i2, i3, i4)
    elif opt == "s":
        return S1234(i1, i2, i3, i4) + Sb1234(i1, i2, i3, i4) + M1234(i1, i2, i3, i4)
    elif opt == "a":
        return A1234(i1, i2, i3, i4) + Ab1234(i1, i2, i3, i4) + Mb1234(i1, i2, i3, i4)
    elif opt == "t":
        return T1234(i1, i2, i3, i4)
    else:
        raise ValueError(f"invalid option, {opt}")


# ==================================================
def create_active_dict(active, cartesian_mp):
    d = {}
    for cc, lst in cartesian_mp.items():
        for X in ["Q", "G", "T", "M"]:
            ex = []
            for c, v in lst:
                Xv = X + v
                if Xv in active:
                    ex.append((c, X + "_{" + v + "}"))
            if ex:
                d[X + cc] = ex
    return d


# ==================================================
def split_symbol(s):
    m = re.fullmatch(r"([A-Za-z]+)(?:\^\{\(([^)]*)\)\})?(?:_\{([^}]*)\})?", s)
    if m:
        base = m.group(1)
        sup = f"({m.group(2)})" if m.group(2) else ""
        sub = m.group(3) if m.group(3) else "s"
        return base, sup, sub
    else:
        return s, "", "s"


# ==================================================
def get_response_tensor_mp(rt, active_dict, axial_tensor, magnetic_tensor):
    to_axial_dict = {"Q": "G", "G": "Q", "T": "M", "M": "T"}
    to_magnetic_dict = {"Q": "T", "G": "M", "T": "Q", "M": "G"}

    lst = [(c, *split_symbol(v)) for c, v in rt]
    ex_lst = []
    for c, X, ss, comp in lst:
        if axial_tensor:
            X = to_axial_dict[X]
        if magnetic_tensor:
            X = to_magnetic_dict[X]
        if X + comp in active_dict.keys():
            mp = active_dict[X + comp]
            for mc, m in mp:
                ex_lst.append((c * mc, ss, m))

    s = sp.S(0)
    for c, ss, m in ex_lst:
        s += c * sp.Symbol(m + "^{" + ss + "}")

    return s


# ==================================================
def simplify_tensor(M):
    """
    Simplify tensor.

    Args:
        M (ndarray): tensor, [[(symbol, expression)]].

    Returns:
        - (ndarray): simplified tensor.
        - (dict): expression dict.
    """
    M = [[(sp.Symbol(c), m) for c, m in lst] for lst in M]
    Xs = [x for lst in M for _, m in lst for x in m.atoms(sp.Symbol)]
    Xs = list(set(Xs))
    Cs = list(reversed([c for lst in M for c, m in lst if m != sp.S(0)]))
    args = Xs + Cs
    eqs = [sp.Eq(c, m) for lst in M for c, m in lst]
    tensor_map = sp.solve(eqs, args)

    Ms = [[None for _ in lst] for lst in M]
    d = {}
    for i, lst in enumerate(M):
        for j, (c, m) in enumerate(lst):
            if c in tensor_map:
                c_ = tensor_map[c]
                if c_ != sp.S(0):
                    if all([i in Cs for i in c_.atoms(sp.Symbol)]):
                        c = c_
            if m != 0:
                Ms[i][j] = c
                sym = c.free_symbols
                if len(sym) > 1:
                    for s in sym:
                        if s not in d.keys():
                            raise Exception(f"fail to solve.")
                else:
                    if c.could_extract_minus_sign():
                        c = -c
                        m = -m
                    d[c] = m
            else:
                Ms[i][j] = sp.S(0)

    Ms = np.asarray(Ms)

    return Ms, d
