"""
This file contains utility functions to create matrix elements of atomic multipoles.
"""
import sympy as sp
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j
from joblib import Parallel, delayed

from multipie.data.data_tag_harmonics_alias import _data_alias_oh, _data_alias_d6h
from multipie.multipole.util.pauli import pauli_matrix_n
from multipie.multipole.util.atomic_orbital_util import rank
from multipie.multipole.util.spin_orbital_basis import _standard_spinless_basis, _standard_spinful_basis
from multipie.tag.tag_multipole import TagMultipole
from multipie.group.point_group import PointGroup


_lm_basis = [sp.sympify(o.replace("U", "1/2").replace("D", "-1/2")) for o in _standard_spinful_basis["lm"]]
_jm_basis = [sp.sympify(o) for o in _standard_spinful_basis["jm"]]


# ==================================================
def sort_matrix(M, idx_row_list=None, idx_col_list=None):
    """
    sort row and col of given matrix.

    Args:
        M (sympy): Matrix.
        idx_row_list (list, optional): list of row indexes after sorting.
        idx_col_list (list, optional): list of column indexes after sorting.

    Returns:
        sympy: sorted Matrix.
    """
    N_row = sp.S(M.rows)
    N_col = sp.S(M.cols)

    if idx_row_list is None:
        idx_row_list = [i for i in range(N_row)]
    if idx_col_list is None:
        idx_col_list = [i for i in range(N_col)]

    if len(idx_row_list) != N_row:
        raise Exception(f"invalid len(idx_row_list) = {len(idx_row_list)}")
    if len(idx_col_list) != N_col:
        raise Exception(f"invalid len(idx_col_list) = {len(idx_col_list)}")

    M_sorted = sp.zeros(N_row, N_col)
    for i, idx_row in enumerate(idx_row_list):
        for j, idx_col in enumerate(idx_col_list):
            M_sorted[i, j] = M[idx_row, idx_col]

    return M_sorted


# ==================================================
def Rl_l1l2(l, l1, l2):
    """
    Proportionality coefficients of magnetic (electric) toroidal multipoles.
    T(G) = Rl_l1l2 * Q(M)

    Args:
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: coefficient
    """
    return -sp.S(l1 * (l1 + 1) - l2 * (l2 + 1)) / (l + 1) / (l + 2)


# ==================================================
def Ql_l1l2(l, l1, l2):
    """
    reduced matrix element of Q.

    Args:
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: reduced matrix element.
    """
    return (-1) ** l1 * sp.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * wigner_3j(l1, l2, l, 0, 0, 0)


# ==================================================
def Ml_l1l2(l, l1, l2):
    """
    reduced matrix element of M.

    Args:
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: reduced matrix element.
    """
    return (
        (-1) ** l1
        * sp.sqrt(l * (2 * l + 1) * (2 * l - 1) * (2 * l1 + 1) * l2 * (l2 + 1))
        * sp.S(2 * (2 * l2 + 1))
        / (l + 1)
        * wigner_3j(l1, l2, l - 1, 0, 0, 0)
        * wigner_6j(l - 1, l, 1, l2, l2, l1)
    )


# ==================================================
def Tl_l1l2(l, l1, l2):
    """
    reduced matrix element of T.

    Args:
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: reduced matrix element.
    """
    return 0 if l == 0 else sp.I * Rl_l1l2(l, l1, l2) * Ql_l1l2(l, l1, l2)


# ==================================================
def Gl_l1l2(l, l1, l2):
    """
    reduced matrix element of G.

    Args:
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: reduced matrix element.
    """
    return 0 if l < 2 else sp.I * Rl_l1l2(l, l1, l2) * Ml_l1l2(l, l1, l2)


# ==================================================
def Xl_l1l2(X, l, l1, l2):
    """
    reduced matrix element of X (Q/M/T/G).

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        sympy: reduced matrix element.
    """
    qmtg_map = {"Q": Ql_l1l2, "M": Ml_l1l2, "T": Tl_l1l2, "G": Gl_l1l2}
    return qmtg_map[X](l, l1, l2)


# ==================================================
def is_zero(X, l, l1, l2):
    """
    matrix element is zero ?

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank of multipole.
        l1 (int): bra <l1|.
        l2 (int): ket |l2>.

    Returns:
        bool: zero ?
    """
    return (
        # odd parity multipoles are zero if (l1 + l2) % 2 == 0
        ((l1 + l2) % 2 == 0 and ((X in ("Q", "T") and l % 2 == 1) or (X in ("M", "G") and l % 2 == 0)))
        # even parity multipoles are zero if (l1 + l2) % 2 == 1
        or ((l1 + l2) % 2 == 1 and ((X in ("Q", "T") and l % 2 == 0) or (X in ("M", "G") and l % 2 == 1)))
    )


# ==================================================
def Xlm_lm(X, l, m, b1, b2):
    """
    matrix element of spinless atomic multipoles (LM basis).
    <b1|X_{lm}|b2>

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank.
        m (int): component.
        b1 (list): bra <l1, m1|.
        b2 (list): ket |l2, m2>.

    Returns:
        sympy: matrix element.
    """
    (l1, m1) = b1
    (l2, m2) = b2

    if is_zero(X, l, l1, l2) or (l < 0 or l1 < 0 or l2 < 0):
        return 0
    else:
        return (-1) ** (l1 + m1) * wigner_3j(l1, l2, l, -m1, m2, m) * Xl_l1l2(X, l, l1, l2)


# ==================================================
def Pl_sk_J1J2L1L2(l, s, k, J1, J2, l1, l2):
    """
    proportionality coefficient of spinful atomic multipole.
    P_{l}(s,k;J1,J2;L1,L2)

    Args:
        l (int) : rank of multipole.
        s (int) : 0 (spinless) or 1 (spinful).
        k (int) : -s <= k <= s.
        J1 (sympy) : J1.
        J2 (sympy) : J2.
        l1 (int) : l1.
        l2 (int) : l2.

    Returns:
        sympy: coefficient.
    """
    return (
        sp.I ** (s + k)
        * (-1) ** s
        * sp.sqrt((2 * l + 1) * (2 * J1 + 1) * (2 * J2 + 1) * sp.factorial(1 - s) * sp.factorial(2 + s))
        * wigner_9j(l1, J1, 1 / 2, l2, J2, 1 / 2, l + k, l, s)
    )


# ==================================================
def qmtg_map_spinful(X, s, k):
    """
    relation between spinful multipoles head (X) and spinless multipole (Q/M/T/G).

    Args:
        X (str): "Q"/"M"/"T"/"G".
        s (int) : 0 (spinless) or 1(spinful).
        k (int) : -s <= k <= s.

    Returns:
        str: "Q"/"M"/"T"/"G".
    """
    qmtg_map = {
        (0, 0): {"Q": "Q", "M": "M", "T": "T", "G": "G"},
        (1, 1): {"Q": "M", "M": "Q", "T": "G", "G": "T"},
        (1, 0): {"Q": "T", "M": "G", "T": "Q", "G": "M"},
        (1, -1): {"Q": "M", "M": "Q", "T": "G", "G": "T"},
    }
    return qmtg_map[(s, k)][X]


# ==================================================
def atomic_multipole_matrix_lm_basis(X, l, m, s, k, b1, b2):
    """
    matrix element of spinless/spinful multipoles (LMS basis).
    <b1|X_{lm}^{(s)}(k)|b2>

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank.
        m (int): component.
        s (int) : 0 (spinless) or 1(spinful).
        k (int) : -s <= k <= s.
        b1 (list): bra <l1, m1, s1|.
        b2 (list): ket |l2, m2, s2>.

    Returns:
        sympy: matrix element.
    """
    (l1, m1, s1) = b1
    (l2, m2, s2) = b2

    # spinless
    if s == 0 and k == 0:
        return Xlm_lm(X, l, m, (l1, m1), (l2, m2)) if s1 == s2 else 0
    # spinful
    else:
        X = qmtg_map_spinful(X, s, k)
        xlmsk = 0
        for n in range(-1, 2):
            xlm = Xlm_lm(X, l + k, m - n, (l1, m1), (l2, m2))
            if xlm == 0:
                continue

            wig3j = wigner_3j(l + k, l, s, m - n, -m, n)
            if wig3j == 0:
                continue

            xlmsk += wig3j * xlm * pauli_matrix_n(s, n, s1, s2)

        return xlmsk * sp.I ** (s + k) * (-1) ** (l + m) * sp.sqrt(2 * l + 1)


# ==================================================
def atomic_multipole_matrix_jm_basis(X, l, m, s, k, b1, b2):
    """
    matrix element of spinless/spinful multipoles (JM basis)
    <b1|X_{lm}^{(s)}(k)|b2>

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank.
        m (int): component.
        s (int) : 0 (spinless) or 1(spinful).
        k (int) : -s <= k <= s.
        b1 (list): bra <J1, M1, l1|.
        b2 (list): ket |J2, M2, l2>.

    Returns:
        Matrix: matrix element.
    """
    (J1, M1, l1) = b1
    (J2, M2, l2) = b2

    wig3j = wigner_3j(J1, J2, l, -M1, M2, m)
    if wig3j == 0:
        return 0

    P = Pl_sk_J1J2L1L2(l, s, k, J1, J2, l1, l2)
    if P == 0:
        return 0

    X = qmtg_map_spinful(X, s, k)
    return (-1) ** (J1 + M1) * wig3j * P * Xl_l1l2(X, l + k, l1, l2)


# ==================================================
def atomic_multipole_matrix(X, l, m, s, k, b1, b2, b_type="lm"):
    """
    matrix element of spinless/spinful multipoles (LMS/JM basis).
    <b1|X_{lm}^{(s)}(k)|b2>

    Args:
        X (str): "Q"/"M"/"T"/"G".
        l (int): rank.
        m (int): component.
        s (int) : 0 (spinless) or 1(spinful).
        k (int) : -s <= k <= s.
        b1 (list): bra <l1, m1, s1|.
        b2 (list): ket |l2, m2, s2>.
        b_type (str, optional): type of basis, ``lm/jm``.

    Returns:
        sympy: matrix element.
    """
    if b_type == "lm":
        return atomic_multipole_matrix_lm_basis(X, l, m, s, k, b1, b2)
    else:  # "jm"
        return atomic_multipole_matrix_jm_basis(X, l, m, s, k, b1, b2)


# ==================================================
def create_atomic_multipole(b_type="lm"):
    """
    create atomic multipoles X_{lm}^{(s)}(k) (LMS/JM basis).

    Args:
        b_type (str, optional): type of basis, ``lm/jm``.

    Returns:
        dict: atomic multipoles {MultipoleTag: Matrix}.
    """
    if b_type == "lm":
        basis = _lm_basis
    elif b_type == "jm":
        basis = _jm_basis
    else:
        raise Exception("b_type: invalid basis type. choose 'lm'/'jm'.")

    def proc(xlmsk):
        (X, l, m, s, k) = xlmsk
        tag = TagMultipole.create_spherical(head=X, rank=l, mul=0, comp=m, s=s, k=k)
        am = sp.Matrix.zeros(32)
        for i, b1 in enumerate(basis):
            for j, b2 in enumerate(basis):
                b1 = sp.sympify(b1)
                b2 = sp.sympify(b2)
                v = atomic_multipole_matrix(X, l, m, s, k, b1, b2, b_type)
                am[i, j] = sp.simplify(v)

        return tag, am

    max_l = 7
    max_s = 1
    Xlmsk_list = [
        (X, l, m, s, k)
        for X in ["Q", "M", "T", "G"]
        for l in range(max_l + 1)
        for m in reversed(range(-l, l + 1))
        for s in range(max_s + 1)
        for k in reversed(range(-s, s + 1))
    ]

    sub = Parallel(n_jobs=-1, verbose=10)([delayed(proc)(xlmsk) for xlmsk in Xlmsk_list])

    Xlmsk = {tag: am for tag, am in sub if not am.is_zero_matrix}

    return Xlmsk


# ==================================================
def create_atomic_multipole_harmonics_basis(Xlmsk_lm, harmonics_type):
    """
    create atomic multipoles X_{lm}^{(s)}(k) (harmonics basis).

    Args:
        Xlmsk_lm (dict): atomic multipoles (lm basis).
        harmonics_type (str): type of harmonics (cubic/hexagonal).

    Returns:
        dict: atomic multipoles {MultipoleTag: Matrix} (harmonics basis).
    """
    idx_row_list = [2 * i for i in range(16)] + [2 * i + 1 for i in range(16)]
    idx_col_list = idx_row_list
    Xlmsk_lm = {tag: sort_matrix(M, idx_row_list, idx_col_list) for tag, M in Xlmsk_lm.items()}

    if harmonics_type == "cubic":
        pg = "Oh"
        spdf_harmonics_dict = _data_alias_oh
    elif harmonics_type == "hexagonal":
        pg = "D6h"
        spdf_harmonics_dict = _data_alias_d6h
    else:
        raise Exception("harmonics_type: invalid harmonics type. choose 'cubic'/'hexagonal'.")

    spdf_basis = _standard_spinless_basis[harmonics_type]
    hs = PointGroup(pg).harmonics

    U_mat = sp.Matrix.zeros(16, 1)
    for i, o in enumerate(spdf_basis):
        l = rank(o, spinful=False, crystal=harmonics_type)
        h_tag = spdf_harmonics_dict[o]
        u = hs[h_tag].u_matrix().tolist()
        if l == 0:
            m = u + [0] * 15
        elif l == 1:
            m = [0] + u + [0] * 12
        elif l == 2:
            m = [0] * 4 + u + [0] * 7
        elif l == 3:
            m = [0] * 9 + u
        U_mat = U_mat.col_insert(i + 2, sp.Matrix(m))

    U_mat.col_del(0)

    U = sp.Matrix.zeros(32)
    U[:16, :16] = U[16:, 16:] = U_mat
    d = {tag: U.adjoint() * Xlmsk * U for tag, Xlmsk in Xlmsk_lm.items()}
    d = {tag: sp.simplify(m) for tag, m in d.items()}

    idx_row_list = []
    for i in range(16):
        idx_row_list += [i, i + 16]

    idx_col_list = idx_row_list
    Xlmsk_harmonics = {tag: sort_matrix(M, idx_row_list, idx_col_list) for tag, M in d.items()}

    return Xlmsk_harmonics
