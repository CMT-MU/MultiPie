"""
calculate response tensors up to 4th rank.

Notes:
    - index number starts from 1.
"""
from itertools import product
import sympy as sp
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy import LeviCivita
from multipie.data.data_compatibility_relation import _data_compatibility_relation_cubic, _data_compatibility_relation_hexagonal
from multipie.data.data_cartesian_to_ch_harmoncis import (
    _data_cartesian_to_cubic_harmonics,
    _data_cartesian_to_hexagonal_harmonics,
)

_dv = {1: "x", 2: "y", 3: "z"}
_dq = {(1, 1): "1", (2, 2): "2", (3, 3): "3", (2, 3): "4", (3, 1): "5", (1, 2): "6", (3, 2): "-4", (1, 3): "-5", (2, 1): "-6"}
_dq1 = {(1, 1): "1", (2, 2): "2", (3, 3): "3", (2, 3): "4", (3, 1): "5", (1, 2): "6", (3, 2): "4", (1, 3): "5", (2, 1): "6"}
_mul = {"a": "1", "b": "2", "c": "3", "d": "4", "e": "5", "f": "6"}


# ==================================================
def delta(i1, i2):
    """
    wrapper to ``KroneckerDelta()``

    Args:
        i1 (Number, Symbol): The first index of the KroneckerDelta function.
        i2 (Number, Symbol): The second index of the KroneckerDelta function.

    Returns:
        Number, Symbol: 1(i1==i2)/0(i1!=i2)
    """
    return KroneckerDelta(i1, i2)


# ==================================================
def epsilon(i1, i2, i3):
    """
    wrapper to ``LeviCivita()``.

    Args:
        i1 (Number, Symbol): The first index of the Levi-Civita function.
        i2 (Number, Symbol): The second index of the Levi-Civita function.
        i3 (Number, Symbol): The third index of the Levi-Civita function.

    Returns:
        Number, Symbol: 1((i1-i2)*(i2-i3)*(i3-i1)/2==1)/0(otherwise)
    """
    return LeviCivita(i1, i2, i3)


# ==================================================
def d1234(i1, i2, i3, i4):
    return delta(i1, i3) * delta(i2, i4) + delta(i1, i4) * delta(i2, i3)


# ==================================================
def d12345(i1, i2, i3, i4, i5):
    return (
        epsilon(i1, i3, i5) * delta(i2, i4)
        + epsilon(i2, i3, i5) * delta(i1, i4)
        + epsilon(i1, i4, i5) * delta(i2, i3)
        + epsilon(i2, i4, i5) * delta(i1, i3)
    )


# ==================================================
def d123456(i1, i2, i3, i4, i5, i6):
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
            epsilon(i1, i3, i4) * delta(i2, i5)
            + epsilon(i1, i2, i4) * delta(i3, i5)
            - epsilon(i1, i2, i3) * delta(i4, i5)
            - epsilon(i2, i3, i4) * delta(i1, i5)
        )
        / 2
    )


# ==================================================
def e561234(i1, i2, i3, i4, i5, i6):
    return epsilon(i1, i2, i5) * epsilon(i3, i4, i6) + epsilon(i1, i2, i6) * epsilon(i3, i4, i5)


# ==================================================
def g12534(i1, i2, i3, i4, i5):
    return epsilon(i1, i3, i4) * delta(i2, i5) + epsilon(i2, i3, i4) * delta(i1, i5)


# ==================================================
def g125634(i1, i2, i3, i4, i5, i6):
    return (
        sp.S(
            d1234(i1, i2, i4, i5) * delta(i3, i6)
            + d1234(i1, i2, i4, i6) * delta(i3, i5)
            - d1234(i1, i2, i3, i5) * delta(i4, i6)
            - d1234(i1, i2, i3, i6) * delta(i4, i5)
        )
        / 2
    )


# ==================================================
def mp_string(lst, head, sup="", sort=True):
    """
    string of multipole

    Args:
        head (str): head of the multipole Q/G.
        sup (str): superscript.
        sort (bool, optional): sort ?

    Returns:
        str: string of multipole
    """
    lst = [_dv[i] for i in lst]
    if sort:
        s = head + "_{" + "".join(sorted("".join(lst))) + "}"
    else:
        s = head + "_{" + ("".join(lst)) + "}"
    if sup != "":
        s += "^{" + sup + "}"
    return s


# ==================================================
def mp_symbol(lst, head, sup="", sort=True):
    """
    symbol of multipole

    Args:
        head (str): head of the multipole Q/G.
        sup (str): superscript.
        sort (bool, optional): sort ?

    Returns:
        Symbol: symbol of multipole
    """
    return sp.Symbol(mp_string(lst, head, sup, sort))


# ==================================================
def S0():
    """
    component of the rank 0 tensor

    Returns:
        str: "Q_{0}"
    """
    return "Q_{0}"


# ==================================================
def S1(i1):
    """
    component of the rank 1 tensor

    Args:
        i1 (Number, Symbol): The index of the tensor.

    Returns:
        str: "Q_{i}" (i=x/y/z)
    """
    return mp_string([i1], "Q")


# ==================================================
def S12(i1, i2):
    """
    component of the symmetric part of the rank 2 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    if i1 == i2:
        lst = [(1, mp_string([i1, i2], "Q")), (1, "Q_{0}")]
    else:
        lst = [(1, mp_string([i1, i2], "Q"))]

    return lst


# ==================================================
def A12(i1, i2):
    """
    component of the asymmetric part of the rank 2 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    lst = []
    for i3 in range(1, 4):
        e = epsilon(i1, i2, i3)
        if e != 0:
            lst.append((e, mp_string([i3], "G")))

    return lst


# ==================================================
def S123(i1, i2, i3):
    """
    component of the symmetric part of the rank 3 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = [(1, mp_string([i1, i2, i3], "Q"))]
    l2 = []
    for i4, i5 in product(range(1, 4), range(1, 4)):
        g = g12534(i1, i2, i3, i4, i5)
        if g != 0:
            l2.append((g, mp_string([i4, i5], "G", "a")))
    l3 = []
    for i4 in range(1, 4):
        d = d1234(i1, i2, i3, i4)
        if d != 0:
            l3.append((d, mp_string([i4], "Q", "a")))

    if i1 == i2:
        l4 = [(1, mp_string([i3], "Q", "b"))]
    else:
        l4 = []

    return l1 + l2 + l3 + l4


# ==================================================
def A123(i1, i2, i3):
    """
    component of the asymmetric part of the rank 3 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = []
    for i4 in range(1, 4):
        e = epsilon(i1, i2, i4)
        if e != 0:
            l1.append((e, mp_string([i3, i4], "G", "b")))
    l2 = []
    for i4 in range(1, 4):
        e = e1234(i1, i2, i3, i4)
        if e != 0:
            l2.append((e, mp_string([i4], "Q", "c")))

    e = epsilon(i1, i2, i3)
    if e != 0:
        l3 = [(e, "G_{0}")]
    else:
        l3 = []

    return l1 + l2 + l3


# ==================================================
def S1234(i1, i2, i3, i4):
    """
    component of the (symmetric,symmetric,symmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = [(1, mp_string([i1, i2, i3, i4], "Q"))]
    l2 = []
    if i1 == i2:
        l2.append((1, mp_string([i3, i4], "Q", "a")))
    if i3 == i4:
        l2.append((1, mp_string([i1, i2], "Q", "a")))

    l3 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        d = d123456(i1, i2, i3, i4, i5, i6)
        if d != 0:
            l3.append((d, mp_string([i5, i6], "Q", "b")))

    if i1 == i2 and i3 == i4:
        l4 = [(1, "Q_{0}^{a}")]
    else:
        l4 = []

    d = d1234(i1, i2, i3, i4)
    if d != 0:
        l5 = [(d, "Q_{0}^{b}")]
    else:
        l5 = []

    return l1 + l2 + l3 + l4 + l5


# ==================================================
def Sb1234(i1, i2, i3, i4):
    """
    component of the (symmetric,symmetric,asymmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = []
    for i5, i6, i7 in product(range(1, 4), range(1, 4), range(1, 4)):
        d = d1234567(i1, i2, i3, i4, i5, i6, i7)
        if d != 0:
            l1.append((d, mp_string([i5, i6, i7], "G", "a")))
    l2 = []
    if i1 == i2:
        l2.append((1, mp_string([i3, i4], "Q", "c")))
    if i3 == i4:
        l2.append((-1, mp_string([i1, i2], "Q", "c")))
    l3 = []
    for i5 in range(1, 4):
        d = d12345(i1, i2, i3, i4, i5)
        if d != 0:
            l3.append((d, mp_string([i5], "G", "a")))

    return l1 + l2 + l3


# ==================================================
def A1234(i1, i2, i3, i4):
    """
    component of the (asymmetric,asymmetric,symmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        e = e561234(i1, i2, i3, i4, i5, i6)
        if e != 0:
            l1.append((e, mp_string([i5, i6], "Q", "f")))

    e = e1234(i1, i2, i3, i4)
    if e != 0:
        l2 = [(e, "Q_{0}^{c}")]
    else:
        l2 = []

    return l1 + l2


# ==================================================
def Ab1234(i1, i2, i3, i4):
    """
    component of the (asymmetric,asymmetric,asymmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    lst = []
    for i5 in range(1, 4):
        e = e51234(i1, i2, i3, i4, i5)
        if e != 0:
            lst.append((e, mp_string([i5], "G", "f")))

    return lst


# ==================================================
def M1234(i1, i2, i3, i4):
    """
    component of the (symmetric,asymmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = []
    for i5 in range(1, 4):
        e = epsilon(i3, i4, i5)
        if e != 0:
            l1.append((e, mp_string([i1, i2, i5], "G", "b")))
    l2 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        g = g125634(i1, i2, i3, i4, i5, i6)
        if g != 0:
            l2.append((g, mp_string([i5, i6], "Q", "d")))
    l3 = []
    for i5 in range(1, 4):
        g = g12534(i1, i2, i3, i4, i5)
        if g != 0:
            l3.append((g, mp_string([i5], "G", "b")))
    l4 = []
    for i5 in range(1, 4):
        de = delta(i1, i2) * epsilon(i3, i4, i5)
        if de != 0:
            l4.append((de, mp_string([i5], "G", "c")))

    return l1 + l2 + l3 + l4


# ==================================================
def Mb1234(i1, i2, i3, i4):
    """
    component of the (asymmetric,symmetric) part of the rank 4 tensor

    Args:
        i1 (Number, Symbol): The first index of the tensor.
        i2 (Number, Symbol): The second index of the tensor.
        i3 (Number, Symbol): The third index of the tensor.
        i4 (Number, Symbol): The fourth index of the tensor.

    Returns:
        list: [(coefficient, multipole)]
    """
    l1 = []
    for i5 in range(1, 4):
        e = epsilon(i1, i2, i5)
        if e != 0:
            l1.append((e, mp_string([i3, i4, i5], "G", "c")))
    l2 = []
    for i5, i6 in product(range(1, 4), range(1, 4)):
        g = g125634(i3, i4, i1, i2, i5, i6)
        if g != 0:
            l2.append((g, mp_string([i5, i6], "Q", "e")))
    l3 = []
    for i5 in range(1, 4):
        g = g12534(i3, i4, i1, i2, i5)
        if g != 0:
            l3.append((g, mp_string([i5], "G", "d")))
    l4 = []
    for i5 in range(1, 4):
        de = delta(i3, i4) * epsilon(i1, i2, i5)
        if de != 0:
            l4.append((de, mp_string([i5], "G", "e")))

    return l1 + l2 + l3 + l4


# ==================================================
def simplify_tensor(M):
    """
    simplify tensor symbol.

    ex) when Sxx = Syy, Syy is replaced by Sxx.

    Args:
        M ([list]):

    Returns:
        [list]: simplified tensor.
    """
    Xs = [x for lst in M for _, m in lst for x in m.atoms(sp.Symbol)]
    Xs = list(set(Xs))
    Cs = list(reversed([c for lst in M for c, m in lst if m != sp.S(0)]))
    args = Xs + Cs
    eqs = [sp.Eq(c, m) for lst in M for c, m in lst]
    tensor_map = sp.solve(eqs, args)

    Ms = [[None for _ in lst] for lst in M]
    for i, lst in enumerate(M):
        for j, (c, m) in enumerate(lst):
            if c in tensor_map:
                c_ = tensor_map[c]
                if c_ != sp.S(0):
                    if all([i in Cs for i in c_.atoms(sp.Symbol)]):
                        c = c_
            Ms[i][j] = (c, m)

    return Ms


# ==================================================
def create_tensor(pg_tag, tag):
    """
    create response tensors

    Args:
        pg_tag (TagGroup): tag of point group
        tag (TagResponseTensor): tag of response tensor.

    Returns:
        [list]: [[(tensor symbol, component)]]
    """
    rank, comp, i_type, t_type = tag.rank, tag.comp, tag.i_type, tag.t_type

    is_axial = i_type == "axial"
    is_magnetic_tensor = t_type == "magnetic"

    if pg_tag.is_hexagonal_subgroup():
        cr = _data_compatibility_relation_hexagonal
        mp_map = _data_cartesian_to_hexagonal_harmonics
    else:
        cr = _data_compatibility_relation_cubic
        mp_map = _data_cartesian_to_cubic_harmonics

    A_irrep = cr[str(pg_tag)][0]

    def is_A_irrep(idx):  # Is given irrep full symmetric (A) irrep ?
        return cr[str(pg_tag)][idx] == A_irrep

    def replace_head_axial(mp):
        if is_axial:
            mp = mp.translate(str.maketrans({"Q": "G", "G": "Q"}))
        return mp

    def replace_head_magnetic(mp):
        if is_magnetic_tensor:
            mp = mp.replace("Q", "T").replace("G", "M")
        return mp

    # create response tensors with full symmetric component
    # rank 0
    if rank == 0:
        mp = replace_head_axial(S0())
        c, mp, idx = mp_map[mp][0]
        x = sp.S(0)
        if is_A_irrep(idx):
            mp = replace_head_magnetic(mp)
            x += c * sp.Symbol(mp)
        M = [[(sp.Symbol(tag.latex()), x)]]
    # rank 1
    elif rank == 1:
        M = [[("", sp.S(0))] * 3] * 1
        for i1 in range(1, 4):
            x = sp.S(0)
            mp = replace_head_axial(S1(i1))
            for c, mp, idx in mp_map[mp]:
                if is_A_irrep(idx):
                    mp = replace_head_magnetic(mp)
                    x += c * sp.Symbol(mp)
            sub = _dv[i1]
            c = sp.Symbol(tag.latex() + "_{" + sub + "}")
            M[0][i1 - 1] = (c, x)
    # rank 2
    elif rank == 2:
        M = [[("", sp.S(0))] * 3 for _ in range(3)]
        if comp == "s":
            tf = S12
        else:
            tf = A12
        for i1, i2 in product(range(1, 4), range(1, 4)):
            x = sp.S(0)
            for c1, m in tf(i1, i2):
                m = replace_head_axial(m)
                for c2, mp, idx in mp_map[m]:
                    if is_A_irrep(idx):
                        mp = replace_head_magnetic(mp)
                        x += c1 * c2 * sp.Symbol(mp)
            sub = _dv[i1] + _dv[i2]
            c = sp.Symbol(tag.latex() + "_{" + sub + "}")
            M[i1 - 1][i2 - 1] = (c, x)
    # rank 3
    elif rank == 3:
        if comp == "s":
            tf = S123
            lst = [(i1, i2, i3) for i1, i2, i3 in product(range(1, 4), range(1, 4), range(1, 4)) if i1 <= i2]
            M = [[("", sp.S(0))] * 3 for _ in range(6)]
        else:
            tf = A123
            lst = [(i1, i2, i3) for i1, i2, i3 in product(range(1, 4), range(1, 4), range(1, 4)) if i1 != i2 and i1 <= i2]
            M = [[("", sp.S(0))] * 3 for _ in range(3)]
        for i1, i2, i3 in lst:
            x = sp.S(0)
            mul = ""
            for c1, m in tf(i1, i2, i3):
                if "^" in m:
                    m, mul = m.split("^")
                    mul = mul[1:-1]
                else:
                    mul = ""
                m = replace_head_axial(m)
                for c2, mp, idx in mp_map[m]:
                    if is_A_irrep(idx):
                        mp = replace_head_magnetic(mp)
                        if mul == "":
                            x += c1 * c2 * sp.Symbol(mp)
                        else:
                            x += c1 * c2 * sp.Symbol(mp + "[" + _mul[mul] + "]")

            if comp == "s":
                sub = _dq1[(i1, i2)] + _dv[i3]
                i = int(_dq1[(i1, i2)]) - 1
            else:
                sub = _dq1[(i1, i2)] + _dv[i3]
                i = int(_dq1[(i1, i2)]) - 4
                if int(_dq[(i1, i2)]) < 0:
                    x = -x
            j = i3 - 1
            c = sp.Symbol(tag.latex() + "_{" + sub + "}")
            M[i][j] = (c, x)
    # rank 4
    elif rank == 4:
        if comp == "sss":
            tf = S1234
        elif comp == "ssa":
            tf = Sb1234
        elif comp == "aas":
            tf = A1234
        elif comp == "aaa":
            tf = Ab1234
        elif comp == "sa":
            tf = M1234
        elif comp == "as":
            tf = Mb1234

        if comp in ("sss", "ssa"):
            lst = [
                (i1, i2, i3, i4)
                for i1, i2, i3, i4 in product(range(1, 4), range(1, 4), range(1, 4), range(1, 4))
                if i1 <= i2 and i3 <= i4
            ]
            M = [[("", sp.S(0))] * 6 for _ in range(6)]
        elif comp == "sa":
            lst = [
                (i1, i2, i3, i4)
                for i1, i2, i3, i4 in product(range(1, 4), range(1, 4), range(1, 4), range(1, 4))
                if i1 <= i2 and i3 <= i4 and i3 != i4
            ]
            M = [[("", sp.S(0))] * 3 for _ in range(6)]
        elif comp == "as":
            lst = [
                (i1, i2, i3, i4)
                for i1, i2, i3, i4 in product(range(1, 4), range(1, 4), range(1, 4), range(1, 4))
                if i1 <= i2 and i3 <= i4 and i1 != i2
            ]
            M = [[("", sp.S(0))] * 6 for _ in range(3)]
        else:
            lst = [
                (i1, i2, i3, i4)
                for i1, i2, i3, i4 in product(range(1, 4), range(1, 4), range(1, 4), range(1, 4))
                if i1 != i2 and i1 <= i2 and i3 != i4 and i3 <= i4
            ]
            M = [[("", sp.S(0))] * 3 for _ in range(3)]

        for i1, i2, i3, i4 in lst:
            x = sp.S(0)
            mul = ""
            for c1, m in tf(i1, i2, i3, i4):
                if "^" in m:
                    m, mul = m.split("^")
                    mul = mul[1:-1]
                else:
                    mul = ""
                m = replace_head_axial(m)
                for c2, mp, idx in mp_map[m]:
                    if is_A_irrep(idx):
                        mp = replace_head_magnetic(mp)
                        if mul == "":
                            x += c1 * c2 * sp.Symbol(mp)
                        else:
                            x += c1 * c2 * sp.Symbol(mp + "[" + _mul[mul] + "]")

            if comp in ("sss", "ssa"):
                i = int(_dq1[(i1, i2)]) - 1
                j = int(_dq1[(i3, i4)]) - 1
                sub = _dq1[(i1, i2)] + _dq1[(i3, i4)]
            elif comp == "sa":
                i = int(_dq1[(i1, i2)]) - 1
                j = int(_dq1[(i3, i4)]) - 4
                if int(_dq[(i3, i4)]) < 0:
                    x = -x
                sub = _dq1[(i1, i2)] + _dv[int(_dq1[(i3, i4)]) - 3]
            elif comp == "as":
                i = int(_dq1[(i1, i2)]) - 4
                j = int(_dq1[(i3, i4)]) - 1
                if int(_dq[(i1, i2)]) < 0:
                    x = -x
                sub = _dv[int(_dq1[(i1, i2)]) - 3] + _dq1[(i3, i4)]
            else:
                i = int(_dq1[(i1, i2)]) - 4
                j = int(_dq1[(i3, i4)]) - 4
                if int(_dq[(i1, i2)]) * int(_dq[(i3, i4)]) < 0:
                    x = -x
                sub = _dv[int(_dq1[(i1, i2)]) - 3] + _dv[int(_dq1[(i3, i4)]) - 3]

            c = sp.Symbol(tag.latex() + "_{" + sub + "}")
            M[i][j] = (c, x)

    M = [[(c, sp.simplify(m)) for c, m in lst] for lst in M]
    return M
