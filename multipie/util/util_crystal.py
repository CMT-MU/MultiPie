"""
For crystal (for sympy element).
"""

import numpy as np
import sympy as sp

from multipie.util.util import str_to_sympy

TOL_SAME_SITE = 1e-8
CHOP = 1e-8
DIGIT = 8
DISTANCE_DIGIT = 5

# Crystallographic definition.
# ==================================================
# transformation matrix from conventional to primitive cell (4x4).
# or from cartesian to fractional coordinate for hexagonal system (h).
# x' = P^-1 x, W' = P^-1WP.
_P_dict = {
    "A": "[[1,0,0,0],[0,1/2,-1/2,0],[0,1/2,1/2,0],[0,0,0,1]]",
    "B": "[[1/2,0,1/2,0],[0,1,0,0],[-1/2,0,1/2,0],[0,0,0,1]]",
    "C": "[[1/2,-1/2,0,0],[1/2,1/2,0,0],[0,0,1,0],[0,0,0,1]]",
    "P": "[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]",
    "I": "[[-1/2,1/2,1/2,0],[1/2,-1/2,1/2,0],[1/2,1/2,-1/2,0],[0,0,0,1]]",
    "F": "[[0,1/2,1/2,0],[1/2,0,1/2,0],[1/2,1/2,0,0],[0,0,0,1]]",
    "R": "[[2/3,-1/3,-1/3,0],[1/3,1/3,-2/3,0],[1/3,1/3,1/3,0],[0,0,0,1]]",
    "a": "[[0,1,0,0],[0,0,1,0],[1,0,0,0],[0,0,0,1]]",  # a-axis.
    "b": "[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]",  # b-axis (default setting).
    "c": "[[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]]",  # c-axis.
    "h": "[[1,-1/2,0,0],[0,sqrt(3)/2,0,0],[0,0,1,0],[0,0,0,1]]",  # hexagonal.
    "0": "[[1]]",  # point group.
}
P_dict = {lattice: str_to_sympy(m) for lattice, m in _P_dict.items()}
Pi_dict = {lattice: np.array(sp.Matrix(m).inv()) for lattice, m in P_dict.items()}

# ==================================================
# partial translation n x (1x3).
_t_dict = {
    "A": "[[0,0,0],[0,1/2,1/2]]",
    "B": "[[0,0,0],[1/2,0,1/2]]",
    "C": "[[0,0,0], [1/2,1/2,0]]",
    "P": "[[0,0,0]]",  # primitive.
    "I": "[[0,0,0],[1/2,1/2,1/2]]",
    "F": "[[0,0,0],[0,1/2,1/2],[1/2,0,1/2],[1/2,1/2,0]]",
    "R": "[[0,0,0],[2/3,1/3,1/3],[1/3,2/3,2/3]]",
    "0": None,  # point group.
}
t_dict = {lattice: None if t is None else str_to_sympy(t) for lattice, t in _t_dict.items()}


# ==================================================
def shift_site(site, tol=TOL_SAME_SITE):
    """
    Shift site within home unit cell.

    Args:
        site (ndarray): (set of) site (sympy or float).
        tol (float, optional): tolerance for same site.

    Returns:
        - (ndarray) -- shifted site.

    Note:
        - When using sympy elements, ensure that all elements are sympy. Creating an array like np.array([0, sp.S(1)/2, 1]) can be dangerous, as all elements may be treated as integers in some cases. To ensure that all elements are treated as sympy objects, it is recommended to use str_to_sympy.
    """
    site = np.mod(site, 1)
    site = np.vectorize(lambda i: i.args[0] if isinstance(i, sp.Mod) else i)(site)
    if site.dtype != object:
        site[np.abs(site - 1.0) < tol] = 0.0
    return site


# ==================================================
def shift_bond(bond):
    """
    Shift bond within home unit cell.

    Args:
        bond (ndarray): (set of) bond (sympy or float) [vector+center].

    Returns:
        - (ndarray) -- shifted bond.
    """
    if bond.ndim == 1:
        return np.concatenate([bond[0:3], shift_site(bond[3:6])])
    else:
        return np.concatenate([bond[:, 0:3], shift_site(bond[:, 3:6])])


# ==================================================
def convert_to_fractional_hexagonal(vec_c):
    """
    Convert vector to fractional coordinate for hexagonal system.

    Args:
        vec_c (ndarray): vector in cartesian coordinate.

    Returns:
        - (ndarray) -- in fractional coordinate.
    """
    return convert_to_primitive("h", vec_c, shift=False)


# ==================================================
def convert_to_fractional_hexagonal_matrix(mat_c):
    """
    Convert matrix to fractional coordinate for hexagonal systems.

    Args:
        mat_c (ndarray): matrix in fractioanl coordinate.

    Returns:
        - (ndarray) -- in fractional coordinate.
    """
    return convert_to_primitive_matrix("h", mat_c)


# ==================================================
def convert_to_cartesian_hexagonal(vec_f):
    """
    Convert vector to cartesian coordinate for hexagonal system.

    Args:
        vec_f (ndarray): vector in fractional coordinate.

    Returns:
        - (ndarray) -- in cartesian coordinate.
    """
    return convert_to_conventional("h", vec_f, plus_set=False, shift=False)


# ==================================================
def convert_to_cartesian_hexagonal_matrix(mat_f):
    """
    Convert matrix to cartesian coordinate for hexagonal system.

    Args:
        mat_f (ndarray): matrix in fractioanl coordinate.

    Returns:
        - (ndarray) -- in cartesian coordinate.
    """
    return convert_to_conventional_matrix("h", mat_f)


# ==================================================
def convert_to_primitive(lattice, vec_cf, shift=True):
    """
    Convert vector to primitive cell in fractional coordinate.

    Args:
        lattice (str): crystal lattice, (A/B/C/P/I/F/R/0). [0: point group].
        vec_cf (ndarray): vector of conventional cell in fractional coordinate.
        shift (bool, optional): shift to home cell ?

    Returns:
        - (ndarray) -- in primitive cell.
    """
    if lattice == "0":
        return vec_cf

    if lattice == "P":
        if shift:
            vec_cf = shift_site(vec_cf)
        return vec_cf
    else:
        if vec_cf.ndim == 1:
            vec_cf = np.pad(vec_cf, (0, 1), constant_values=1)
        else:
            vec_cf = np.pad(vec_cf, ((0, 0), (0, 1)), constant_values=1)

        Pi = Pi_dict[lattice]
        vec_pf = vec_cf @ Pi.T

        if shift:
            vec_pf = shift_site(vec_pf)

        if vec_pf.ndim == 1:
            vec_pf = vec_pf[0:3]
        else:
            vec_pf = vec_pf[:, 0:3]

        return vec_pf


# ==================================================
def convert_to_primitive_matrix(lattice, mat_cf):
    """
    Convert matrix to primitive cell in fractional coordinate.

    Args:
        lattice (str): crystal lattice, (A/B/C/P/I/F/R/0). [0: point group].
        mat_cf (ndarray): matrix of conventional cell in fractional coordinate.

    Returns:
        - (ndarray) -- in primitive cell.
    """
    if lattice in ["P", "0"]:
        return mat_cf
    else:
        n = mat_cf.shape[-1]
        if n == 3:
            if mat_cf.ndim == 2:
                mat_cf = np.pad(mat_cf, ((0, 1), (0, 1)), constant_values=0)
                mat_cf[n, n] = 1
            else:
                mat_cf = np.pad(mat_cf, ((0, 0), (0, 1), (0, 1)), constant_values=0)
                mat_cf[:, n, n] = 1

        P = P_dict[lattice]
        Pi = Pi_dict[lattice]
        mat_pf = Pi @ mat_cf @ P

        if n == 3:
            if mat_pf.ndim == 2:
                mat_pf = mat_pf[0:3, 0:3]
            else:
                mat_pf = mat_pf[:, 0:3, 0:3]

        return mat_pf


# ==================================================
def convert_to_conventional(lattice, vec_pf, plus_set=False, shift=True):
    """
    Convert vector to conventional cell in fractional coordinate.

    Args:
        lattice (str): crystal lattice, (A/B/C/P/I/F/R/0). [0: point group].
        vec_pf (ndarray): vector of primitive cell in fractional coordinate.
        plus_set (bool, optional): add partial translations ?
        shift (bool, optional): shift to home cell ?

    Returns:
        - (ndarray) -- in conventioanl cell.

    Note:
        - for plus_set, [set(t0), set(t1), ...].
    """
    if lattice == "0":
        return vec_pf

    if lattice == "P":
        if shift:
            vec_pf = shift_site(vec_pf)
        return vec_pf
    else:
        if vec_pf.ndim == 1:
            vec_pf = np.pad(vec_pf, (0, 1), constant_values=1)
        else:
            vec_pf = np.pad(vec_pf, ((0, 0), (0, 1)), constant_values=1)

        P = P_dict[lattice]
        vec_cf = vec_pf @ P.T

        if vec_cf.ndim == 1:
            vec_cf = vec_cf[0:3]
        else:
            vec_cf = vec_cf[:, 0:3]

        if plus_set:
            t = t_dict[lattice]
            vec_cf = np.concatenate([vec_cf + i for i in t])

        if shift:
            vec_cf = shift_site(vec_cf)

        return vec_cf


# ==================================================
def convert_to_conventional_matrix(lattice, mat_pf):
    """
    Convert matrix to conventional cell in fractional coordinate.

    Args:
        lattice (str): crystal lattice, (A/B/C/P/I/F/R/0). [0: point group].
        mat_pf (ndarray): matrix of primitive cell in fractional coordinate.

    Returns:
        - (ndarray) -- in conventioanl cell.
    """
    if lattice in ["P", "0"]:
        return mat_pf
    else:
        n = mat_pf.shape[-1]
        if n == 3:
            if mat_pf.ndim == 2:
                mat_pf = np.pad(mat_pf, ((0, 1), (0, 1)), constant_values=0)
                mat_pf[n, n] = 1
            else:
                mat_pf = np.pad(mat_pf, ((0, 0), (0, 1), (0, 1)), constant_values=0)
                mat_pf[:, n, n] = 1

        P = P_dict[lattice]
        Pi = Pi_dict[lattice]
        mat_cf = P @ mat_pf @ Pi

        if n == 3:
            if mat_cf.ndim == 2:
                mat_cf = mat_cf[0:3, 0:3]
            else:
                mat_cf = mat_cf[:, 0:3, 0:3]

        return mat_cf


# ==================================================
def create_igrid(cell_range):
    """
    Create integer grid.

    Args:
        cell_range (tuple): range, (a1,b1,a2,b2,a3,b3).

    Returns:
        - (ndarray) -- integer grid (float,N,3).
    """
    offset = cell_range[::2]
    N = [i - j for i, j in zip(cell_range[1::2], offset)]

    i_vals = np.arange(N[0])
    j_vals = np.arange(N[1])
    k_vals = np.arange(N[2])

    ii, jj, kk = np.meshgrid(i_vals, j_vals, k_vals, indexing="ij")

    igrid = np.stack([ii, jj, kk], axis=-1)

    igrid += np.array(offset)
    igrid = igrid.reshape(-1, 3)

    return igrid


# ==================================================
def site_distance(tail, heads, G=None, digit=DISTANCE_DIGIT):
    """
    group of heads with the same distance from tail (in increasing order).

    Args:
        tail (ndarray): tail vector.
        heads (ndarray): head vectros.
        G (ndarray, optional): metric matrix (None = unit matrix).
        digit (int, optional): accuracy of digit.

    Returns:
        - (dict) -- grouped heads (sorted), Dict[ distance(float), [heads] ].
    """
    if G is None:
        G = np.eye(3)

    d = heads - tail
    r = np.sqrt(np.einsum("ij,ij->i", d @ G, d))
    mask = r > 10**-digit
    r, heads = r[mask], heads[mask]
    keys = np.round(r, digit)
    unique, inverse = np.unique(keys, return_inverse=True)
    dic = {float(k): heads[inverse == i] for i, k in enumerate(unique)}
    dic = dict(sorted(dic.items()))

    return dic


# ==================================================
def get_cell_info(crystal, cell):
    """
    Get unit cell information.

    Args:
        crystal (str): crystal, "triclinic/monoclinic/orthorhombic/trigonal/hexagonal/tetragonal/cubic".
        cell (dict): cell, key="a/b/c/alpha/beta/gamma".

    Returns:
        - (dict) -- cell(dict), volume(float), A(ndarray,4,4,float), G(ndarray,4,4,float).

    Note:
        - key in cell is omittable.
        - alpha, beta, gamma in unit of degree.
    """
    a = float(cell.get("a", 1.0))
    b = float(cell.get("b", 1.0))
    c = float(cell.get("c", 1.0))
    alpha = float(cell.get("alpha", 90.0))
    beta = float(cell.get("beta", 90.0))
    gamma = float(cell.get("gamma", 90.0))

    if crystal == "monoclinic":
        alpha = gamma = 90.0
    elif crystal == "orthorhombic":
        alpha = beta = gamma = 90.0
    elif crystal in ["trigonal", "hexagonal"]:
        alpha = beta = 90.0
        gamma = 120.0
        b = a
    elif crystal == "tetragonal":
        alpha = beta = gamma = 90.0
        b = a
    elif crystal == "cubic":
        alpha = beta = gamma = 90.0
        b = c = a
    elif crystal == "triclinic":
        a = b = c = 1.0
        alpha = beta = gamma = 90.0

    ca = np.cos(alpha * np.pi / 180)
    cb = np.cos(beta * np.pi / 180)
    cc = np.cos(gamma * np.pi / 180)
    sc = np.sin(gamma * np.pi / 180)
    s = 1.0 - ca * ca - cb * cb - cc * cc + 2.0 * ca * cb * cc
    s = max(CHOP, np.sqrt(s))

    a1 = np.array([a, 0, 0])
    a2 = np.array([b * cc, b * sc, 0])
    a3 = np.array([c * cb, c * (ca - cb * cc) / sc, c * s / sc])

    A = np.eye(4)
    A[0:3, 0] = a1
    A[0:3, 1] = a2
    A[0:3, 2] = a3

    cell = {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}
    volume = float(a * b * c * s)

    A = A.round(DIGIT)
    G = (A.T @ A).round(DIGIT)
    dic = {"cell": cell, "volume": volume, "A": A, "G": G}

    return dic
