"""
equivalent operator.
"""
import sympy
from itertools import permutations
from gcoreutils.nsarray import NSArray


# ==================================================
def _symmetrized_polynomial(f):
    """
    convert polynomial to symmetrized sum.

    Args:
        f (sympy): polynomial.

    Returns:
        [(sympy, sympy, list)]: (coefficient, symmetrized factor, [variables])
    """
    tc = f.as_coefficients_dict()
    info = []
    for t, c in tc.items():
        tm = t.as_terms()
        s = ""
        fa = 1
        fs = 0
        for b, n in zip(tm[1], tm[0][0][1][1]):
            s += str(b) * n
            fa *= sympy.factorial(n)
            fs += n
        info.append((c, fa / sympy.factorial(fs), list(set(permutations(s)))))
    return info


# ==================================================
def _angular_momentum_matrix(j):
    """
    angular-momentum operator in desceding order of Jz.

    Args:
        j (str): magnitude of angular momentum, J (0, 1/2, 1, ...).

    Returns:
        - Matrix: Jx.
        - Matrix: Jy.
        - Matrix: Jz.
    """
    j = sympy.sympify(j)
    d = 2 * j + 1

    jz = sympy.zeros(d)
    for i in range(d):
        jz[i, i] = j - i

    jp = sympy.zeros(d)
    for i in range(1, d):
        m = j - i
        jp[i - 1, i] = sympy.sqrt((j - m) * (j + m + 1))
    jm = jp.T

    jx = (jp + jm) / 2
    jy = (jp - jm) / 2 / sympy.I

    return jx, jy, jz


# ==================================================
def equivalent_operator_from_poly(poly, j):
    """
    equivalent operator in descending order in Jz.

    Args:
        poly (sympy/str): (x,y,z) polynomial.
        j (str): magnitude of angular momentum, J (0, 1/2, 1, ...).

    Returns:
        NSArray: equivalent-operator matrix.
    """
    if type(poly) == str:
        poly = sympy.sympify(poly)
    poly = sympy.expand(poly)
    jx, jy, jz = _angular_momentum_matrix(j)
    jm = {"x": jx, "y": jy, "z": jz}
    sf = _symmetrized_polynomial(poly)

    m = sympy.zeros(jx.rows)
    for c1, c2, sv in sf:
        sf0 = sympy.zeros(jx.rows)
        for p in sv:
            sf1 = sympy.eye(jx.rows)
            for i in p:
                sf1 = sf1 * jm[i]
            sf0 += sf1
        m += c1 * c2 * sf0
    return NSArray(m.tolist(), "matrix")
