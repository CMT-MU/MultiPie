"""
This file provides utility functions for calculation of structure multipole basis set.
"""
import sympy as sp
import numpy as np

from gcoreutils.nsarray import NSArray


# ==================================================
def cs_list(s, real=False):
    """
    list of c### and s### contained in the given equation.

    Args:
        s (str/Symbol): string object to be converted to symbol.
        real (bool): real symbol ? If real=True, the object is recognized as a real number.

    Returns:
        list: [c###/s###]
    """
    s = sp.expand(sp.sympify(str(s)))
    cs_lst = []
    n = 1
    while s != 0:
        cn, sn = sp.Symbol(f"c{n:03d}"), sp.Symbol(f"s{n:03d}")
        c = sp.diff(s, cn)
        s = sp.limit(s, cn, 0)
        if c != 0:
            cs_lst.append(cn)
        c = sp.diff(s, sn)
        s = sp.limit(s, sn, 0)
        if c != 0:
            cs_lst.append(sn)
        n += 1

    if real:
        cs_lst = [sp.Symbol(str(csn), real=True) for csn in cs_lst]

    return cs_lst


# ==================================================
def to_symbol(s, real=False):
    """
    convert string to symbol.

    Args:
        s (str): string object to be converted to symbol.
        real (bool): real symbol ? If real=True, the object is recognized as a real number.

    Returns:
        Symbol: symbol.
    """
    cs_lst = cs_list(s)
    s = sp.expand(sp.sympify(str(s)))
    d = {sp.Symbol(str(csn), real=real): s.coeff(csn) for csn in cs_lst}
    return sp.expand(sum([c * csn for csn, c in d.items()]))


# ==================================================
def inner_product(f1, f2):
    """
    inner product of two structure factors.

    Args:
        f1 (str/Symbol): F1.
        f2 (str/Symbol): F2.

    Returns:
        NSArray: inner product of F1 and F2.

    Notes:
        - F1 and F2 are polynomials in terms of c### and s### (### are given by lst) only without const. terms.
    """
    f1 = str(to_symbol(f1))
    f2 = str(to_symbol(f2))

    ip_str = {}
    for csi in cs_list(f1):
        for csj in cs_list(f2):
            if csi != csj:
                ip_str[str(csi * csj)] = "(0)"
            else:
                ip_str[str(csi * csj)] = "(1/2)"

    prod = str((NSArray(f1) * NSArray(f2)).expand())
    for k, v in ip_str.items():
        prod = prod.replace(k, v)

    prod = sp.expand(sp.sympify(prod))

    return prod


# ==================================================
def orthogonalize_fk(v, nmax=None):
    """
    orthogonalize list of function of k by Gram-Schmidt orthogonalization method.

    Args:
        v (NSArray): list of function of k to be orthogonalized.
        nmax (int, optional): max. number of nonzero basis.

    Returns:
        NSArray: orthogonalized list of function of k.
    """
    ev = []

    cnt = 0
    for v0 in v:
        s = sp.S(0)
        for evi in ev:
            s += inner_product(evi, v0) * evi

        e = v0 - s
        e = sp.expand(e)
        d = sp.sqrt(inner_product(e, e))
        if d != 0:
            e = e / d
            e = sp.expand(e)
            cnt += 1
        else:
            e = sp.S(0)
        ev.append(e)
        if nmax is not None and cnt == nmax:
            break

    idx = np.array([i for i, b in enumerate(ev) if not b == sp.S(0)])

    return ev, idx


# ==================================================
def decompose_fk(fk, k_samb_set):
    """
    decompose k function into linear combination of structure multipoles.

    Args:
        fk (str/Symbol): function of k.
        k_samb_set (list): [ ("kmp_#", "formfactor") ].

    Returns:
        dict: dictionary of expansion coefficient {"smp_#": coeff}.
    """
    fk = sp.expand(to_symbol(fk))
    if fk == sp.S(0):
        return {}

    sm_set = {kmp_i: sp.expand(fk_) for kmp_i, fk_ in k_samb_set}

    coeffs = {}
    for kmp_i, fk_ in sm_set.items():
        c = inner_product(fk, fk_)
        if c != 0:
            coeffs[kmp_i] = c

    return coeffs


def test():
    c1, s1, c2, s2 = sp.symbols("c001 s001 c002 s002")
    eq = sp.expand(sp.sqrt(2) * (c1 - sp.I * s1) + sp.sqrt(3) * (c2 + sp.I * s2))
    print(f"eq = {eq}")

    cs_lst = cs_list(eq)
    print(cs_lst)

    eq_ = sp.conjugate(to_symbol(str(eq), real=True))
    print(f"eq_ = {eq_}")

    f1 = sp.expand(sp.sqrt(2) * (c1 - s1))
    f2 = sp.expand(sp.sqrt(3) * (c1 + s1))
    print(f"inner_product(f1, f1) = {inner_product(f1, f1)}")
    print(f"inner_product(f1, f2) = {inner_product(f1, f2)}")
    print(f"inner_product(f2, f2) = {inner_product(f2, f2)}")

    v = [f1, f2]
    v = orthogonalize_fk(v)[0]
    print(f"orthogonalize_fk = {v}")

    k_samb_set = [("kmp_1", v[0]), ("kmp_2", v[1])]
    print(f"decompose_fk = {decompose_fk(c1, k_samb_set)}")


# test()
