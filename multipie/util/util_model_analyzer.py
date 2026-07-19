"""
Utility for ModelAnalyzer calss.
"""

import os
import numpy as np
import sympy as sp
import subprocess
from collections import defaultdict
from multipie import Group
from multipie.util.util import str_to_sympy


# ==================================================
def grid_path(pts, gpath, N1=100, A=None):
    """
    Create grid points along path.

    Args:
        pts (dict): definitions of points, e.g., {"A":[1,2,3],"B":[3,4,5]}.
        gpath (str): path, e.g., "A-B|C-D-E".
        N1 (int, optional): number of divisions.
        A (ndarray, optional): conversion matrix.

    Returns:
        - (ndarray) -- grid points.
        - (ndarray) -- linear positions.
        - (dict) -- dict[disconnected linear position, label].
    """
    if A is None:
        d = np.array(list(pts.values())[0]).shape[0]
        A = np.eye(d)

    gpath = gpath.replace(" ", "").replace("\t", "").replace("\n", "")
    glabel = gpath.split("-")
    gpath = gpath.split("|")
    gpath = [i.split("-") for i in gpath]
    gpath = [[(i1, i2) for i1, i2 in zip(i[:-1], i[1:])] for i in gpath]
    gpath = sum(gpath, [])

    grid = []
    glin = []
    gdis = [0]
    x = 0
    for s, e in gpath:
        s = np.array(pts[s])
        e = np.array(pts[e])
        dv = (e - s) / N1
        d = np.linalg.norm(dv @ A.T)
        for j in range(N1):
            glin.append(x)
            grid.append(s + j * dv)
            x += d
        grid.append(s + N1 * dv)
        glin.append(x)
        gdis.append(x)

    grid = np.asarray(grid)
    glin = np.asarray(glin)
    gdis = dict(zip(gdis, glabel))

    return grid, glin, gdis


# ==================================================
def grid_all(N1, N2, N3):
    """
    Create N1 x N2 x N3 grid points.

    Args:
        N1 (int): number of divisions for a1.
        N2 (int): number of divisions for a2.
        N3 (int): number of divisions for a3.

    Returns:
        - (ndarray) -- grid points.
    """
    kpoints = np.array([[i / float(N1), j / float(N2), k / float(N3)] for i in range(N1) for j in range(N2) for k in range(N3)])

    return kpoints


# ==================================================
def fourier_r_to_k(OR, atom, kv, s=True):
    """
    Fourier transformation from R to k.

    Args:
        OR (dict): matrix dict, dict[(n1,n2,n3,m,n): value].
        atom (ndarray): atomic postions.
        kv (ndarray): list of k vector.
        s (bool, optional): include phase of atomic position ?

    Returns:
        - (ndarray) -- O(k) = sum_s Omn exp( 2pi i k*Rs ), Rs = (n1,n2,n3)-s(R[m]-R[n]).
    """
    eRa = np.exp(1j * (kv @ atom.T))

    Nk = len(kv)
    d = len(atom)
    S = np.zeros((Nk, d, d), dtype=complex)
    for (R, m, n), value in OR.items():
        eR = np.exp(1j * (2 * np.pi * kv @ np.asarray(R)))
        if s:
            S[:, m, n] += value * eR * np.conj(eRa[:, m]) * eRa[:, n]
        else:
            S[:, m, n] += value * eR

    return S


# ==================================================
def output_linear_dispersion_eig(filename, k, e, o=None, **kwargs):
    """
    Output band dispersion along high-symmetry lines.
    (only eigen values)

    Args:
        filename (str): file name.
        k (ndarray): k points along high-symmetry lines.
        e (ndarray): eigen values.
        o (list, optional): expectation value of any operator for each band, [[num_k, num_band]].
        kwargs (dict, optional): key words for generate_band_gnuplot.
    """
    kmax = np.max(k)
    emax = np.max(e)
    emin = np.min(e)
    num_wann = e.shape[1]

    fs = open(filename, "w")
    fs.write("# k Energy [eV] \n")
    fs.write(f"# Emax = {str(emax)}\n")
    fs.write(f"# Emin = {str(emin)}\n")
    fs.write(f"# num_wann = {str(num_wann)}\n")

    ef = kwargs.get("ef", None)

    if ef is None:
        ef = 0.0
        kwargs["ef"] = 0.0
        fs.write("# ef = ? (no shift) \n\n")
    else:
        fs.write(f"# shifted by fermi energy = {ef} [eV] \n\n")

    num_k, Nm = e.shape
    e = e.T

    for n in range(Nm):
        en = e[n] - ef
        s = ""
        for i in range(num_k):
            if o is None:
                s += "{k:0<20}   {e:<20} \n".format(k=k[i], e=en[i])
            else:
                s += "{k:0<20}   {e:<20}".format(k=k[i], e=en[i])
                for oi in o:
                    s += "   {o:<20}".format(o=oi[i, n])
                s += " \n"

        s += "\n"
        fs.write(s)

    fs.close()

    # generate gnuplot file
    generate_band_gnuplot_eig(filename, kmax, emax, emin, num_wann, **kwargs)


# ==================================================
def generate_band_gnuplot_eig(filename, kmax, emax, emin, num_wann, **kwargs):
    """
    Generate gnuplot file to plot band dispersion.
    (only eigen values)

    Args:
        filename (str): file name.
        kmax (float): maximum value in kpoints.
        emax (float): maximum value of eigen values.
        emin (float): minimum value of eigen values.
        num_wann (int): # of wannier functions.
        kwargs (dict, optional): key words for generate_band_gnuplot.
            - a (float): length of lattice vector.
            - ef (float): fermi energy.
            - k_dis_pos (dict): {disconnected linear position:label}.
            - lwidth (float): line width.
            - lc (str): line color.
            - ref_filename (str): file name of reference band data.
            - colormap (bool): with colormap of first expectation value of operator.
    """
    offset = (emax - emin) * 0.1

    a = kwargs.get("a", None)
    ef = kwargs.get("ef", 0.0)
    k_dis_pos = kwargs.get("k_dis_pos", None)
    ref_filename = kwargs.get("ref_filename", None)
    lwidth = kwargs.get("lwidth", 2)
    lc = kwargs.get("lc", "salmon")
    colormap = kwargs.get("colormap", False)

    fs = open("plot_band.gnu", "w")
    fs.write("unset key \n")
    fs.write("unset grid \n")
    fs.write(f"lwidth = {lwidth} \n")
    fs.write(f"set xrange [:{kmax}] \n")
    fs.write(f"set yrange [{emin-ef-offset}:{emax-ef+offset}] \n")
    fs.write("set tics font 'Times New Roman, 14' \n\n")
    fs.write("set size ratio 0.7 \n\n")

    fs.write('set palette defined ( -1.0 "royalblue", 0 "gray90", 1.0 "salmon")\n')
    fs.write("set cbrange [-1.0:1.0]\n\n")

    if k_dis_pos is not None:
        for pos, label in k_dis_pos.items():
            fs.write(f"set arrow from  {pos},  {emin-ef-offset} to {pos}, {emax-ef+offset} nohead \n")

        k_dis_pos = {pos: label.replace("G", "{/Symbol G}").replace("|", ":") for pos, label in k_dis_pos.items()}
        fs.write("set xtics (" + "".join([f'"{label}" {pos},' for pos, label in k_dis_pos.items()]) + ") \n\n")

    fs.write(f"ef = {ef} \n")

    if a is not None:
        fs.write(f"a = {a} \n")

    fs.write("set terminal postscript eps color enhanced \n\n")

    fn_eps = os.path.splitext(filename)[0] + ".eps"
    fs.write(f"set output '{fn_eps}' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        fs.write(f"'{ref_filename}' u ($1/a):($2-ef) w l lw lwidth lc 'dark-grey', ")

    if colormap:
        fs.write(f"'{filename}' u 1:2:3 w l lw lwidth lc palette, ")
    else:
        fs.write(f"'{filename}' u 1:2 w l lw lwidth lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 lc 'black'")

    fs.write(" \n\n")

    fs.write("set terminal pdf \n\n")

    fn_pdf = os.path.splitext(filename)[0] + ".pdf"
    fs.write(f"set output '{fn_pdf}' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        fs.write(f"'{ref_filename}' u ($1/a):($2-ef) w l lw lwidth lc 'dark-grey', ")

    if colormap:
        fs.write(f"'{filename}' u 1:2:3 w l lw lwidth lc palette, ")
    else:
        fs.write(f"'{filename}' u 1:2 w l lw lwidth lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 lc 'black'")

    fs.close()

    subprocess.run(f"gnuplot plot_band.gnu", shell=True)


# ==================================================
def parse_orbital(orbital_str):
    """
    Parse orbital string.

    Args:
        orbital_str (str): orbital string, (orbital, spin).

    Returns:
        - (tuple or str) -- orbital, spin or orbital.
    """
    if orbital_str.count("("):
        return orbital_str.strip("()").split(",")
    else:
        return orbital_str, "u"


# ==================================================
def create_operator_index(ket, mode=None):
    """
    Create operator index.

    Args:
        ket (list): ket list.
        mode (str, optional): "orbital" or "spin".

    Returns:
        - (list) -- index pair list.

    Note:
        - atom, sublattice, rank is the same in pair.
        - spin is the same for "orbital" mode.
        - orbital is the same for "spin" mode.
    """
    results = []

    for i, (atom_i, sub_i, rank_i, idx_i, orb_i) in enumerate(ket):
        o_i, s_i = parse_orbital(orb_i)

        for j, (atom_j, sub_j, rank_j, idx_j, orb_j) in enumerate(ket):
            o_j, s_j = parse_orbital(orb_j)

            if atom_i != atom_j or sub_i != sub_j or rank_i != rank_j:  # same atom, sub, rank.
                continue

            if mode == "spin" and o_i != o_j:  # same orbital.
                continue

            if mode == "orbital" and s_i != s_j:  # same spin.
                continue

            results.append((i, j, rank_i, idx_i, idx_j))

    return results


# ==================================================
def create_local_operator(ket, op_name, op_dict, spinful):
    """
    Create local operator.

    Args:
        ket (list): ket set.
        op_name (str): operator name.
        op_dict (dict): operator dict.
        spinful (bool): lgs basis ?

    Returns:
        - (ndarray) -- matrix for given operator.
    """
    if op_name[0] == "S":
        mode = "spin"
    else:
        mode = "orbital"

    dim = len(ket)
    mat = np.zeros((dim, dim), dtype=complex)
    lst = create_operator_index(ket, mode)
    if mode == "spin":
        for m, n, rank, b_idx, k_idx in lst:
            if spinful:
                b_idx = b_idx // 2
                k_idx = k_idx // 2
            mat[m, n] += complex(op_dict["s"][op_name][b_idx, k_idx])
    else:
        for m, n, rank, b_idx, k_idx in lst:
            if spinful:
                b_idx = b_idx // 2
                k_idx = k_idx // 2
            mat[m, n] += complex(op_dict[rank][op_name][b_idx, k_idx])

    return mat


# ==================================================
def create_all_local_operator():
    """
    Create local operators.

    Returns:
        - (dict) -- local operators, dict[rank, dict[name, matrix]].

    Note:
        - atomic basis is tesseral harmonics, s, px, py, pz, dv, dxy, dxz, dyz, du, f2, f1, fbz, f3, f3x, f3y, faz.
        - operator names are Lx, Ly, Lz, Qu, Qv, Qyz, Qxz, Qxy, which are not normalized as defined.
        - it is used for basis, lg and lgs, not for jml.
    """
    op_key = [("M", 1, "T1g", -1, -1, 0, 0, "q"), ("Q", 2, "Eg", -1, -1, 0, 0, "q"), ("Q", 2, "T2g", -1, -1, 0, 0, "q")]
    op_factor = {
        1: ["sqrt(2)", "sqrt(6)/5", "sqrt(6)/5"],
        2: ["sqrt(10)", "sqrt(14)/7", "sqrt(14)/7"],
        3: ["14/sqrt(7)", "sqrt(21)*2/15", "sqrt(21)*2/15"],
    }
    group = Group("Oh")

    lst = ["Lx", "Ly", "Lz", "Qu", "Qv", "Qyz", "Qxz", "Qxy"]
    d = {0: {name: str_to_sympy("[[0]]") for name in lst}}
    for rank in [1, 2, 3]:
        no = 0
        am = group.atomic_samb("lg", (rank, rank))
        dd = {}
        for l, f in zip(op_key, op_factor[rank]):
            f = str_to_sympy(f)
            v, ex = am[l]
            for e, m in zip(ex, v):
                name = lst[no]
                no = no + 1
                dd[name] = m * f
        d[rank] = dd

    spin = {"Sx": "[[0,1/2],[1/2,0]]", "Sy": "[[0,-I/2],[I/2,0]]", "Sz": "[[1/2,0],[0,-1/2]]"}
    d["s"] = {name: str_to_sympy(s) for name, s in spin.items()}

    return d


# ==================================================
def create_k_multipole(cluster_samb, cluster_vector):
    """
    Create momentum multipole.

    Args:
        cluster_samb (dict): cluster SAMB.
        cluster_vector (dict): cluster vector.

    Returns:
        - (dict) -- momentum multipole in terms of k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
        - (dict) -- cluster vector, dict[site/bond name, dict[kb, expression] ].
    """
    k_multipole = {}
    for k, v in cluster_samb.items():
        if "@" in k:  # bond.
            d = len(next(iter(v.values()))[0][0])
            kb = np.array(sp.symbols(" ".join([f"kb_{i+1}" for i in range(d)])), dtype=object)
            c = list(map(sp.cos, kb))
            s = list(map(sp.sin, kb))
            d_wp = {}
            for idx, (samb, sym) in v.items():
                if idx[0] == "Q":
                    d_wp[idx] = (sp.sqrt(2) * samb @ c, sym)
                else:
                    d_wp[idx] = (sp.sqrt(2) * sp.I * samb @ s, sym)
            k_multipole[k] = d_wp
        else:  # site.
            k_multipole[k] = v

    kb_dic = {}
    kv = np.array(sp.symbols("k_1 k_2 k_3"), dtype=object)
    for sb, lst in cluster_vector.items():
        if ";" in sb:
            d = len(lst[0])
            kb = np.array(sp.symbols(" ".join([f"kb_{i+1}" for i in range(d)])), dtype=object)
            kb_dic[sb] = {i: j @ kv for i, j in zip(kb, lst)}
        else:
            kb_dic[sb] = {}

    return k_multipole, kb_dic


# ==================================================
def create_k_matrix(matrix, cluster_dict, vector_dict):
    """
    Create momentum matrix.

    Args:
        matrix (dict): O(R) matrix dict.
        cluster_dict (dict): cluster dict.
        vector_dict (dict): vector dict.

    Returns:
        - (dict) -- momentum matrix, dict[tag, dict[(m,n), value] ].

    Notes:
        - only tight-binding gauge is supported.
    """
    k_matrix = {}
    for tag in matrix.keys():
        cluster = cluster_dict[tag]
        OR = matrix[tag]
        if ";" in cluster:
            vec = vector_dict[cluster]
            mat = defaultdict(lambda: sp.S(0))
            kb = np.array(sp.symbols(" ".join([f"kb_{i+1}" for i in range(len(vec))])), dtype=object)
            for (n1, n2, n3, m, n), (value, b_no) in OR.items():
                k = kb[b_no - 1] if b_no > 0 else -kb[-b_no - 1]
                mat[(m, n)] += value * sp.exp(sp.I * k)
            mat = {Rmn: v for Rmn, v in mat.items() if not v.is_zero}
            k_matrix[tag] = mat
        else:
            k_matrix[tag] = {(m, n): v for (n1, n2, n3, m, n), (v, b) in OR.items()}

    return k_matrix
