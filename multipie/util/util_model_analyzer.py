"""
Utility for ModelAnalyzer class.
"""

import os
import numpy as np
import sympy as sp
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from multipie import Group
from multipie.util.util import str_to_sympy
from multipie.util.util_constant import M_ZERO, k_B_SI, elem_charge_SI


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
        atom (ndarray): atomic positions.
        kv (ndarray): list of k vector.
        s (bool, optional): include phase of atomic position ?

    Returns:
        - (ndarray) -- O(k) = sum_s Omn exp( 2pi i k*Rs ), Rs = (n1,n2,n3)-s(R[m]-R[n]).
    """
    eRa = np.exp(1j * (2 * np.pi * kv @ atom.T))

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
def fourier_k_to_r(Ok, atom, kv, irvec, s=True):
    """
    Fourier transformation from k to R (inverse of fourier_r_to_k).

    Args:
        Ok (ndarray): O(k) matrix, shape (Nk, d, d).
        atom (ndarray): atomic positions.
        kv (ndarray): list of k vector.
        irvec (iterable): list/array of R = (n1, n2, n3) lattice vectors to
            evaluate O_R at.
        s (bool, optional): include phase of atomic position ? (must match the
            value used to generate `Ok`).

    Returns:
        dict: dict[(n1,n2,n3,m,n): value], O_R = (1/Nk) sum_k O(k) exp(-2pi i k*Rs),
        Rs = (n1,n2,n3)-s(R[m]-R[n]), i.e. the inverse of `fourier_r_to_k`.
    """
    eRa = np.exp(1j * (2 * np.pi * kv @ atom.T))

    Nk = len(kv)
    d = len(atom)

    OR = {}
    for R in irvec:
        R = tuple(int(v) for v in R)
        eR = np.exp(-1j * (2 * np.pi * kv @ np.asarray(R, dtype=float)))

        if s:
            phase = eR[:, None, None] * eRa[:, :, None] * np.conj(eRa)[:, None, :]
        else:
            phase = eR[:, None, None] * np.ones((Nk, d, d), dtype=complex)

        value = np.sum(Ok * phase, axis=0) / Nk

        for m in range(d):
            for n in range(d):
                OR[(R, m, n)] = value[m, n]

    return OR


# ==================================================
def output_dispersion(filename, k, e, o=None):
    """
    Output band dispersion along high-symmetry lines.

    Args:
        filename (str): file name.
        k (ndarray): k points along high-symmetry lines.
        e (ndarray): eigen values.
        o (list, optional): expectation value of any operator for each band, [[num_k, num_band]].
    """
    emax = np.max(e)
    emin = np.min(e)
    ef = 0.0

    fs = open(filename, "w")
    fs.write("# k Energy [eV] \n")
    fs.write(f"# Emax = {emax}\n")
    fs.write(f"# Emin = {emin}\n")
    fs.write(f"# ef = {ef}\n\n")

    num_k, Nm = e.shape
    e = e.T

    kv = np.append(k, -1.0)  # add dummy for last one.
    for n in range(Nm):
        en = e[n] - ef
        s = ""
        for i in range(num_k):
            s += "{k:0<20} {e:<20}".format(k=kv[i], e=en[i])
            if o is None:
                s += "\n"
            else:
                for oi in o:
                    s += " {o:<20}".format(o=oi[i, n])
                s += "\n"
            if np.abs(kv[i] - kv[i + 1]) < 1e-5:
                s += "\n"

        s += "\n"
        fs.write(s)

    fs.close()


# ==================================================
def create_gnuplot_cmd(filename, k_dis_pos, kmax, emax, emin, colormap=False, lwidth=2, lc="salmon"):
    """
    Create gnuplot file.

    Args:
        filename (str): file name.
        k_dis_pos (dict): info. of high-symmetry point in linear k, dict[disconnected position, label].
        kmax (float): maximum value in kpoints.
        emax (float): maximum value of eigen values.
        emin (float): minimum value of eigen values.
        colormap (bool, optional): with colormap of first expectation value of operator.
        lwidth (int, optional): line width.
        lc (str, optional): line color.
    """
    offset = (emax - emin) * 0.1
    ef = 0.0

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

    fs.write("set terminal postscript eps color enhanced \n\n")

    fn_eps = os.path.splitext(filename)[0] + ".eps"
    fs.write(f"set output '{fn_eps}' \n\n")
    fs.write("plot ")

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

    if colormap:
        fs.write(f"'{filename}' u 1:2:3 w l lw lwidth lc palette, ")
    else:
        fs.write(f"'{filename}' u 1:2 w l lw lwidth lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 lc 'black'")

    fs.close()


# ==================================================
def plot_save_dispersion(filename, k_dis_pos, colormap=False, lwidth=1, lc="salmon"):
    """
    Generate plot window for band dispersion (matplotlib).

    Args:
        filename (str): file name.
        k_dis_pos (dict): info. of high-symmetry point in linear k, dict[disconnected position, label].
        colormap (bool, optional): with colormap of first expectation value of operator.
        lwidth (int, optional): line width.
        lc (str, optional): line color.
    """
    # read file.
    bands = read_text_data(filename)
    emax = max(a[:, 1].max() for a in bands)
    emin = min(a[:, 1].min() for a in bands)
    kmax = max(a[:, 0].max() for a in bands)
    ef = 0.0

    offset = (emax - emin) * 0.1

    fig, ax = plt.subplots(figsize=(6, 3.5))  # aspect ratio 0.7
    # ax = fig.add_subplot(111)
    fig.subplots_adjust(right=0.84)

    ax.set_xlim(0, kmax)
    ax.set_ylim(emin - ef - offset, emax - ef + offset)

    ax.axhline(0, color="black", lw=0.5)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True)
    ax.tick_params(axis="y", which="both", left=True, right=True, direction="in")

    # vertical lines.
    positions = list(k_dis_pos.keys())
    labels = [v.replace("G", r"$\Gamma$").replace("|", ":") for v in k_dis_pos.values()]
    for x in positions:
        ax.axvline(x, color="black", lw=0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)

    if colormap:
        all_segments = []
        all_colors = []
        for band in bands:
            x = band[:, 0]
            y = band[:, 1] - ef
            c = band[:, 2]
            points = np.column_stack((x, y)).reshape(-1, 1, 2)
            seg = np.concatenate([points[:-1], points[1:]], axis=1)
            all_segments.extend(seg)
            all_colors.extend(c[:-1])

        lcobj = LineCollection(all_segments, cmap="coolwarm", norm=Normalize(-1.0, 1.0), linewidth=lwidth)
        lcobj.set_array(np.asarray(all_colors))
        ax.add_collection(lcobj)
        cbar = fig.colorbar(lcobj, ax=ax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0], fraction=0.035, pad=0.02)
        cbar.ax.tick_params(which="both", length=0)
        cbar.set_ticklabels(["-1", "-0.5", "0", "0.5", "1"])
    else:
        for band in bands:
            ax.plot(band[:, 0], band[:, 1] - ef, color=lc, lw=lwidth)

    fn_pdf = os.path.splitext(filename)[0] + ".pdf"
    fn_eps = os.path.splitext(filename)[0] + ".eps"
    fig.savefig(fn_pdf)
    fig.savefig(fn_eps)

    plt.close("all")


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
            kb = np.array([sp.Symbol(f"p_{i+1}", real=True) for i in range(d)], dtype=object)
            c = list(map(sp.cos, kb))
            s = list(map(sp.sin, kb))
            d_wp = {}
            for idx, (samb, sym) in v.items():
                if idx[0] == "Q":
                    d_wp[idx] = (np.sqrt(2) * np.vectorize(sp.re)(samb).astype(float) @ c, sym)
                else:
                    d_wp[idx] = (sp.I * (np.sqrt(2) * np.vectorize(sp.im)(samb).astype(float)) @ s, sym)
            k_multipole[k] = d_wp
        else:  # site.
            k_multipole[k] = v

    kb_dic = {}
    kv = np.array([sp.Symbol(f"k_{i}", real=True) for i in range(1, 4)], dtype=object)
    for sb, lst in cluster_vector.items():
        if ";" in sb:
            d = len(lst[0])
            kb = np.array([sp.Symbol(f"p_{i+1}", real=True) for i in range(d)], dtype=object)
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
            kb = np.array([sp.Symbol(f"p_{i+1}", real=True) for i in range(len(vec))], dtype=object)
            for (n1, n2, n3, m, n), (value, b_no) in OR.items():
                k = kb[b_no - 1] if b_no > 0 else -kb[-b_no - 1]
                mat[(m, n)] += value * sp.exp(sp.I * k)
            mat = {Rmn: v for Rmn, v in mat.items() if not v.is_zero}
            k_matrix[tag] = mat
        else:
            k_matrix[tag] = {(m, n): v for (n1, n2, n3, m, n), (v, b) in OR.items()}

    return k_matrix


# ==================================================
def add_local_parameter(matrix_info, parameter):
    """
    Add local parameter for given nonlocal parameter.

    Args:
        matrix_info (dict): matrix info.
        parameter (dict): parameter dict without site cluster.

    Returns:
        - (dict) -- parameter dict with adding local one.
    """
    if not len(set([i[2] for i in matrix_info["index"].keys()])) == 1:  # only for the same ranks.
        return
    rank = next(iter(matrix_info["index"]))[2]
    dim = matrix_info["dimension"]

    d = 2 * rank + 1
    nn = dim // d

    local = {}
    local_z = {}
    for tag, mat in matrix_info["matrix"].items():
        v = np.full((dim, dim), sp.S(0))
        for (n1, n2, n3, m, n), (val, no) in mat.items():
            v[m, n] += val
        if ";" in matrix_info["cluster"][tag]:  # nonlocal SAMB with parameter key.
            if tag not in parameter.keys():
                continue
            vd = v.reshape(nn, d, nn, d).sum(axis=0).transpose(1, 0, 2)
            s = np.full((dim, dim), sp.S(0))
            for i in range(nn):
                s[i * d : (i + 1) * d, i * d : (i + 1) * d] = vd[i]
            local[tag] = s
        else:  # local SAMB.
            local_z[tag] = v
    coeff = {tag: np.array([np.einsum("ij,ji->", z, s) for z in local_z.values()]) for tag, s in local.items()}

    zt = np.full(len(local_z.keys()), sp.S(0))
    for zj, val in parameter.items():
        zt -= val * coeff[zj]

    for zk, val in zip(local_z.keys(), zt):
        if not val.is_zero:
            parameter[zk] = val

    parameter = dict(sorted(parameter.items()))

    return parameter


# ==================================================
def build_and_solve_hermitian(mat, z_symbols, labels=None):
    """
    Build linear equations from a Hermitian matrix and solve for z_j.

    Args:
        mat (ndarray): hermitian matrix, only upper triangle is required.
        z_symbols (list): zj variables.
        labels (list, optional): name for bra-ket.

    Returns:
        - (dict) -- dict from matrix index to name, dict[(int,int), name].
        - (list) -- equations used to solve.
        - (dict) -- solution dict, if failed, empty.
    """
    rows, cols = mat.shape
    expr_to_positions = {}

    for m in range(rows):
        for n in range(m, cols):
            expr = sp.sympify(mat[m, n])
            if expr == 0:
                continue
            expr_expanded = sp.expand(expr)
            coeffs = tuple(expr_expanded.coeff(z) for z in z_symbols)
            const = expr_expanded - sum(c * z for c, z in zip(coeffs, z_symbols))
            key = coeffs + (const,)
            expr_to_positions.setdefault(key, []).append((m, n, expr_expanded))

    g_syms = {}
    lin_eqs = []
    for key, positions in expr_to_positions.items():
        m0, n0, expr0 = positions[0]

        if labels is not None:
            label_str = ",".join(f"({labels[m]},{labels[n]})" for m, n, _ in positions)
        else:
            label_str = ",".join(f"({m},{n})" for m, n, _ in positions)

        g = sp.Symbol(f"g_{{{label_str}}}")
        for m, n, _ in positions:
            g_syms[(m, n)] = g
        lin_eqs.append(sp.Eq(g, expr0))

    sol = sp.linsolve(lin_eqs, z_symbols)

    if sol == sp.EmptySet:
        sol_dict = {}
    else:
        (values,) = sol  # take unique tuple in FiniteSet.
        sol_dict = {str(zj): val for zj, val in zip(z_symbols, values)}

    return g_syms, lin_eqs, sol_dict


# ==================================================
def convert_zj_atomic_var(matrix_info, combined_cluster, combined_id, IR):
    """
    Convert from zj to atomic variable.

    Args:
        matrix_info (dict): matrix info.
        combined_cluster (dict): combined cluster.
        combined_id (dict): combined id.
        IR (str): identity irrep.

    Returns:
        - (dict) -- zj to var for each cluster, dict[bond name, dict[zj, var] ].
    """
    dim = matrix_info["dimension"]
    ket = [i.replace("@", "_").replace("(", "").replace(")", "") for i in matrix_info["ket_site"].keys()]

    # classify zj for each cluster.
    cluster = {}
    for tag, name in combined_cluster.items():
        if ";" in name and combined_id[tag][2][2] == IR:
            cluster.setdefault(name, []).append(tag)

    # construct upper-triangle matrix.
    dic = {}
    for name, tags in cluster.items():
        var = []
        mat = np.full((dim, dim), sp.S(0))
        for zj in tags:
            zjv = sp.Symbol(zj, real=True)
            var.append(zjv)
            for (n1, n2, n3, m, n), (val, no) in matrix_info["matrix"][zj].items():
                if no in (1, -1) and m <= n:  # only upper triangle.
                    mat[m, n] += val * zjv
        dic[name] = (mat, var)

    # solve zj for each cluster.
    result = {}
    for name, (mat, var) in dic.items():
        g_syms, lin_eqs, sol = build_and_solve_hermitian(mat, var, ket)
        result[name] = sol

    return result


# ==================================================
def read_text_data(filename):
    """
    Read text data.

    Args:
        filename (str): file name.

    Returns:
        - (list) -- block of data (ndarray) separated by empty line, [block][x, y].
    """
    bands = []
    band = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            # skip comment.
            if line.startswith("#"):
                continue
            # each block ends with empty line.
            if line == "":
                if band:
                    bands.append(np.array(band))
                    band = []
                continue
            band.append([float(x) for x in line.split()])

    # last block.
    if band:
        bands.append(np.array(band))

    return bands


# ==================================================
def is_almost_zero(x):
    """check if x is numerically zero (within a tolerance)."""
    return np.abs(x) < M_ZERO * 100


# ==================================================
def kelvin_to_ev(T_kelvin):
    """convert temperature from Kelvin to eV (k_B * T)."""
    return T_kelvin * k_B_SI / elem_charge_SI


# ==================================================
def fermi_dirac(x, T=0.0, unit="Kelvin"):
    """
    Fermi-Dirac distribution function, f(x) = 1 / (1 + exp(x/T)).

    Args:
        x (ndarray): energy relative to the chemical potential (E - mu), in eV.
        T (float, optional): temperature (Kelvin or eV, see `unit`).
        unit (str, optional): unit of T, "Kelvin" or "eV".

    Returns:
        ndarray: Fermi-Dirac occupation, 0 <= f(x) <= 1.
    """
    if T == 0.0:
        return np.where(x < 0.0, 1.0, np.where(x > 0.0, 0.0, 0.5))

    T_eV = kelvin_to_ev(T) if unit == "Kelvin" else T
    return 0.5 * (1.0 - np.tanh(0.5 * x / T_eV))


# ==================================================
def fermi_dirac_deriv(x, T=0.01, unit="Kelvin"):
    """
    Minus the derivative of the Fermi-Dirac distribution w.r.t. x, -df/dx = f(x)*f(-x)/T.
    A positive, bell-shaped function peaked at x=0 (used e.g. for smearing/DOS broadening).

    Args:
        x (ndarray): energy relative to the chemical potential (E - mu), in eV.
        T (float, optional): temperature (Kelvin or eV, see `unit`). Must be nonzero.
        unit (str, optional): unit of T, "Kelvin" or "eV".

    Returns:
        ndarray: -df/dx.
    """
    T_eV = kelvin_to_ev(T) if unit == "Kelvin" else T
    return fermi_dirac(x, T, unit) * fermi_dirac(-x, T, unit) / T_eV


# ==================================================
def fermi_dirac_deriv2(x, T=0.01, unit="Kelvin"):
    """
    Second derivative of the Fermi-Dirac distribution w.r.t. x, -d^2f/dx^2
    (i.e. the derivative of `fermi_dirac_deriv`).

    Args:
        x (ndarray): energy relative to the chemical potential (E - mu), in eV.
        T (float, optional): temperature (Kelvin or eV, see `unit`). Must be nonzero.
        unit (str, optional): unit of T, "Kelvin" or "eV".

    Returns:
        ndarray: -d^2f/dx^2.
    """
    T_eV = kelvin_to_ev(T) if unit == "Kelvin" else T
    return (1 - 2 * fermi_dirac(-x, T, unit)) * fermi_dirac(x, T, unit) * fermi_dirac(-x, T, unit) / T_eV / T_eV
