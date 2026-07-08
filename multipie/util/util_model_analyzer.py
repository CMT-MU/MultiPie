"""
Utility for ModelAnalyzer calss.
"""

import numpy as np
import subprocess


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
    S = np.zeros((Nk, d, d), dtype=np.complex128)
    for (R, m, n), value in OR.items():
        eR = np.exp(1j * (2 * np.pi * kv @ np.asarray(R)))
        if s:
            S[:, m, n] += value * eR * np.conj(eRa[:, m]) * eRa[:, n]
        else:
            S[:, m, n] += value * eR

    return S


# ==================================================
def output_linear_dispersion_eig(outdir, filename, k, e, o=None, **kwargs):
    """
    Output band dispersion along high-symmetry lines.
    (only eigen values)

    Args:
        outdir (str): input and output files are found in this directory.
        filename (str): file name.
        k (ndarray): k points along high-symmetry lines.
        e (ndarray): eigen values.
        o (ndarray, optional): expectation value of any operator for each band, [num_band, num_k, component].
        kwargs (dict, optional): key words for generate_band_gnuplot.
    """
    kmax = np.max(k)
    emax = np.max(e)
    emin = np.min(e)
    num_wann = e.shape[1]

    filename = filename[:-4] if filename[-4:] == ".txt" else filename

    fs = open(outdir + "/" + filename + ".txt", "w")
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
                if o.ndim == 2:
                    s += "   {o:<20}".format(o=o[n, i])
                else:
                    for a in range(o.shape[2]):
                        s += "   {o:<20}".format(o=o[n, i, a])

                s += " \n"

        s += "\n"
        fs.write(s)

    fs.close()

    # generate gnuplot file
    generate_band_gnuplot_eig(outdir, filename, kmax, emax, emin, num_wann, **kwargs)


# ==================================================
def generate_band_gnuplot_eig(outdir, filename, kmax, emax, emin, num_wann, **kwargs):
    """
    Generate gnuplot file to plot band dispersion.
    (only eigen values)

    Args:
        outdir (str): input and output files are found in this directory.
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
    """
    offset = (emax - emin) * 0.1

    a = kwargs.get("a", None)
    ef = kwargs.get("ef", 0.0)
    k_dis_pos = kwargs.get("k_dis_pos", None)
    ref_filename = kwargs.get("ref_filename", None)
    lwidth = kwargs.get("lwidth", 3)
    lc = kwargs.get("lc", "salmon")

    fs = open(f"{outdir}/plot_band.gnu", "w")
    fs.write("unset key \n")
    fs.write("unset grid \n")
    fs.write(f"lwidth = {lwidth} \n")
    fs.write(f"set xrange [:{kmax}] \n")
    fs.write(f"set yrange [{emin-ef-offset}:{emax-ef+offset}] \n")
    fs.write("set tics font 'Times New Roman, 30' \n\n")
    fs.write("set size ratio 0.7 \n\n")

    if k_dis_pos is not None:
        for pos, label in k_dis_pos.items():
            fs.write(f"set arrow from  {pos},  {emin-ef-offset} to {pos}, {emax-ef+offset} nohead \n")

        k_dis_pos = {pos: "{/Symbol G}" if label == "G" else label for pos, label in k_dis_pos.items()}
        fs.write("set xtics (" + "".join([f"'{label}' {pos}," for pos, label in k_dis_pos.items()]) + ") \n\n")

    fs.write(f"ef = {ef} \n")

    if a is not None:
        fs.write(f"a = {a} \n")

    fs.write("set terminal postscript eps color enhanced \n\n")

    fs.write(f"set output '{filename}.eps' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        # fs.write(f"'{ref_filename}' u ($1/(2*pi)):($2-ef) w l lw lwidth lc 'dark-grey', ")
        fs.write(f"'{ref_filename}' u ($1/a):($2-ef) w l lw lwidth lc 'dark-grey', ")

    fs.write(f"'{filename}.txt' u 1:2 w l lw lwidth dt (3,1) lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 dt (2,1) lc 'black'")

    fs.write(" \n\n")

    fs.write("set terminal pdf \n\n")

    fs.write(f"set output '{filename}.pdf' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        # fs.write(f"'{ref_filename}' u ($1/(2*pi)):($2-ef) w l lw lwidth lc 'dark-grey', ")
        fs.write(f"'{ref_filename}' u ($1/a):($2-ef) w l lw lwidth lc 'dark-grey', ")

    fs.write(f"'{filename}.txt' u 1:2 w l lw lwidth dt (3,1) lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 dt (2,1) lc 'black'")

    fs.close()

    subprocess.run(f"cd {outdir} ; gnuplot plot_band.gnu", shell=True)


# ==================================================
def output_linear_dispersion(outdir, filename, k, e, u, **kwargs):
    """
    Output band dispersion along high-symmetry lines.

    Args:
        outdir (str): input and output files are found in this directory.
        filename (str): file name.
        k (ndarray): k points along high-symmetry lines.
        e (ndarray): eigen values.
        u (ndarray): eigen vectors.
        kwargs (dict, optional): key words for generate_band_gnuplot.
    """
    kmax = np.max(k)
    emax = np.max(e)
    emin = np.min(e)
    num_wann = e.shape[1]

    filename = filename[:-4] if filename[-4:] == ".txt" else filename

    fs = open(outdir + "/" + filename + ".txt", "w")
    fs.write("# n = band, j = orbital, E_n: energy, W_jn: weight\n")
    fs.write("# k E_0 W_00 W_10 ... W_M0\n")
    fs.write("# k E_1 W_01 W_11 ... W_M1\n")
    fs.write("# ...\n")
    fs.write("# k E_d W_0d W_1d ... W_dd\n")
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
    u = np.abs(u) ** 2

    for n in range(Nm):
        en = e[n] - ef
        s = ""
        for i in range(num_k):
            s += "{k:0<20}   {e:<20}".format(k=k[i], e=en[i])
            for j in range(Nm):
                s += " " + str(u[i, j, n])

            s += "\n"

        s += "\n"
        fs.write(s)

    fs.close()

    # generate gnuplot file
    generate_band_gnuplot(outdir, filename, kmax, emax, emin, num_wann, **kwargs)


# ==================================================
def generate_band_gnuplot(outdir, filename, kmax, emax, emin, num_wann, **kwargs):
    """
    Generate gnuplot file to plot band dispersion.

    Args:
        outdir (str): input and output files are found in this directory.
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
    """
    offset = (emax - emin) * 0.1

    a = kwargs.get("a", None)
    ef = kwargs.get("ef", 0.0)
    k_dis_pos = kwargs.get("k_dis_pos", None)
    ref_filename = kwargs.get("ref_filename", None)
    lwidth = kwargs.get("lwidth", 3)
    lc = kwargs.get("lc", "salmon")

    fs = open(f"{outdir}/plot_band_detail.gnu", "w")
    fs.write("unset key \n")
    fs.write("unset grid \n")
    fs.write(f"lwidth = {lwidth} \n")
    fs.write(f"set xrange [:{kmax}] \n")
    fs.write(f"set yrange [{emin-ef-offset}:{emax-ef+offset}] \n")
    fs.write("set tics font 'Times New Roman, 30' \n\n")
    fs.write("set size ratio 0.7 \n\n")

    if k_dis_pos is not None:
        for pos, label in k_dis_pos.items():
            fs.write(f"set arrow from  {pos},  {emin-ef-offset} to  {pos}, {emax-ef+offset} nohead \n")

        k_dis_pos = {pos: "{/Symbol G}" if label == "G" else label for pos, label in k_dis_pos.items()}
        fs.write("set xtics (" + "".join([f"'{label}' {pos}," for pos, label in k_dis_pos.items()]) + ") \n\n")

    fs.write(f"ef = {ef} \n")

    if a is not None:
        fs.write(f"a = {a} \n")

    fs.write("set terminal postscript eps color enhanced \n\n")

    fs.write(f"set output '{filename}.eps' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        fs.write(f"'{ref_filename}' u ($1/(2*pi)):($2-ef) w l lw lwidth lc 'dark-grey', ")

    fs.write(f"'{filename}.txt' u 1:2 w l lw lwidth dt (3,1) lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 dt (2,1) lc 'black'")

    fs.write(" \n\n")

    fs.write("set terminal pdf \n\n")

    fs.write(f"set output '{filename}.pdf' \n\n")
    fs.write("plot ")

    if ref_filename is not None:
        # fs.write(f"'{ref_filename}' u ($1/(2*pi)):($2-ef) w l lw lwidth lc 'dark-grey', ")
        fs.write(f"'{ref_filename}' u ($1/a):($2-ef) w l lw lwidth lc 'dark-grey', ")

    fs.write(f"'{filename}.txt' u 1:2 w l lw lwidth dt (3,1) lc '{lc}', ")

    fs.write(f"{0.0} lw 0.5 dt (2,1) lc 'black'")

    fs.close()

    subprocess.run(f"cd {outdir} ; gnuplot plot_band_detail.gnu", shell=True)
