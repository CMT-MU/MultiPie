"""
For creating atomic multipole data.
"""

import sympy as sp
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j
from joblib import Parallel, delayed

orbital = ["s", "p", "d", "f"]
spin = ["1/2", "-1/2"]

# ==================================================
""" standard spinful basis """
_standard_spinful_basis = {
    "jml": [
        f"({L + sp.S(s)/2},{sp.S(m) / 2},{L})"
        for L in range(len(orbital))
        for s in [-1, 1]
        for m in reversed(range(-(2 * L + s), 2 * L + s + 2, 2))
    ],
    "lms": [f"({l},{m},{s})" for l in range(len(orbital)) for m in reversed(range(-l, l + 1)) for s in spin],
}


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
        return sp.Rational((-1) ** (s1 - 1 / 2)) * sp.sqrt(6) * wigner_3j(1 / 2, 1 / 2, 1, -s1, s2, n)
    else:
        return 0


# ==================================================
_lm_basis = [sp.sympify(o.replace("U", "1/2").replace("D", "-1/2")) for o in _standard_spinful_basis["lms"]]
_jm_basis = [sp.sympify(o) for o in _standard_spinful_basis["jml"]]


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
def atomic_multipole_matrix(X, l, m, s, k, b1, b2, b_type="lms"):
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
        b_type (str, optional): type of basis, ``lms/jml``.

    Returns:
        sympy: matrix element.
    """
    if b_type == "lms":
        return atomic_multipole_matrix_lm_basis(X, l, m, s, k, b1, b2)
    else:  # "jml"
        return atomic_multipole_matrix_jm_basis(X, l, m, s, k, b1, b2)


# ==================================================
def create_atomic_multipole_matrix(b_type="lms"):
    """
    create atomic multipoles X_{lm}^{(s)}(k) (LMS/JM basis).

    Args:
        b_type (str, optional): type of basis, ``lms/jml``.

    Returns:
        dict: atomic multipoles, m=(l,l-1,...,-l), {(head,l,s,k): [[[sympy]]]}.
    """
    verbose = 0  # 0, 1, 10.

    if b_type == "lms":
        basis = _lm_basis
    elif b_type == "jml":
        basis = _jm_basis
    else:
        raise Exception("b_type: invalid basis type. choose 'lms'/'jml'.")

    def proc(xlmsk):
        (X, l, m, s, k) = xlmsk
        am = sp.Matrix.zeros(32)
        for i, b1 in enumerate(basis):
            for j, b2 in enumerate(basis):
                b1 = sp.sympify(b1)
                b2 = sp.sympify(b2)
                v = atomic_multipole_matrix(X, l, m, s, k, b1, b2, b_type)
                am[i, j] = sp.simplify(v)

        return xlmsk, am

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

    sub = Parallel(n_jobs=-1, verbose=verbose)([delayed(proc)(xlmsk) for xlmsk in Xlmsk_list])

    Xlmsk = {tag: am for tag, am in sub if not am.is_zero_matrix}

    Xlsk = {}
    for (X, l, m, s, k), am in Xlmsk.items():
        Xlsk[(X, l, s, k)] = Xlsk.get((X, l, s, k), []) + [am.tolist()]

    return Xlsk
