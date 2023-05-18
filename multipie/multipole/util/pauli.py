"""
This file provides Pauli matices.
"""
import sympy
from sympy.physics.wigner import wigner_3j


# ==================================================
"""
Pauli matrices
"""


def sm0():
    return sympy.Matrix([[1, 0], [0, 1]])


def smx():
    return sympy.Matrix([[0, 1], [1, 0]])


def smy():
    return sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])


def smz():
    return sympy.Matrix([[1, 0], [0, -1]])


def ss0():
    return sympy.symbols(r"\sigma_{0}", commutative=False)


def ssx():
    return sympy.symbols(r"\sigma_{x}", commutative=False)


def ssy():
    return sympy.symbols(r"\sigma_{y}", commutative=False)


def ssz():
    return sympy.symbols(r"\sigma_{z}", commutative=False)


# ==================================================
def pauli_matrix(i=4, matrix=True):
    """
    Pauli matrices (σ0, σx, σy, σz)

    Args:
        i (int, optional): = 0(0), 1(x), 2(y), 3(z), other(all)
        matrix (bool, optional): Matrix ?

    Returns:
        Matrix or sympy or dict: matrix or symbol or dict { str: Matrix }
    """
    if matrix:
        sm = [sm0(), smx(), smy(), smz()]
        if i in (0, 1, 2, 3):
            return sm[i]
        else:
            return {"0": sm[0], "x": sm[1], "y": sm[2], "z": sm[3]}
    else:
        ss = [ss0(), ssx(), ssy(), ssz()]
        if i in (0, 1, 2, 3):
            return ss[i]
        else:
            return {"0": ss[0], "x": ss[1], "y": ss[2], "z": ss[3]}


# ==================================================
def pauli_matrix_n(s, n, s1, s2):
    """
    Pauli matrices (σ00, σ11, σ10, σ1-1)

    Args:
        s (int): 0/1
        n (int): 0/+1/-1
        s1 (sympy): ±1/2
        s2 (sympy): ±1/2

    Returns:
       sympy: matrix element
    """
    if s == 0 and n == 0:
        return 1 if s1 == s2 else 0
    elif s == 1:
        return sympy.Rational((-1) ** (s1 - 1 / 2)) * sympy.sqrt(6) * wigner_3j(1 / 2, 1 / 2, 1, -s1, s2, n)
    else:
        return 0
