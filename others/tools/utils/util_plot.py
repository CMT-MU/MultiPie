"""
For plot.
"""

import numpy as np
import sympy as sp

from multipie.util.util_binary import BinaryManager

from others.tools.utils.util_spherical_harmonics import create_harmonics_set, is_transform_real, create_internal_basis


# ==================================================
def create_spherical_point(N=8):
    """
    Create discrete points on sphere.

    Args:
        N (int, optional): the number of divisions of pi.

    Returns:
        - (ndarray) -- points, [[x(float), y(float), z(float)]].
    """
    # create spherical points.
    pts = []
    for th in range(1, N):
        t = np.pi * th / N
        for phi in range(2 * N):
            p = np.pi * phi / N
            pt = (t, p)
            pts.append(pt)
    pts = [(0, 0)] + pts + [(np.pi, 0)]

    # convert to cartesian coordinate.
    pos = []
    for th, phi in pts:
        x = np.sin(th) * np.cos(phi)
        y = np.sin(th) * np.sin(phi)
        z = np.cos(th)
        pos.append([x, y, z])
    pos = np.asarray(pos)

    return pos


# ==================================================
def plot_harmonics(qtdraw, l, s=0, k=0, point=None, v=None):
    """
    Plot harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        l (int): rank.
        s (int, optional): internal rank.
        k (int, optional): internal component.
        point (array-like, optional): spherical points.
        v (array-like, optional): weight vector.
    """
    if v is None:
        v = np.full(2 * l + 1, 0)
        v[l] = 1  # C0.
    else:
        v = np.array(v, ndmin=1)

    if not is_transform_real(v):
        raise Exception(f"v={v.tolist()} gives complex harmonics.")

    rv = sp.symbols("X Y Z", real=True)
    ex_form = (create_harmonics_set(l) @ v).simplify()

    X = "Q" if (s + k) % 2 == 0 else "G"
    tag = f"({X},{l},{s},{k},q),v={v.tolist()}: " + f"${sp.latex(ex_form)}$"

    ex = v @ create_harmonics_set(l, s, k, rv)
    ex = np.vectorize(sp.factor)(ex)

    _plot_harmonics(qtdraw, ex, rv, s, point, tag)


# ==================================================
def plot_multipole_harmonics(qtdraw, harmonics_pg, idx, comp=-1, point=None):
    X, l, Gamma, n, p, s, k, x = idx
    tag = f"{X}({l},{Gamma},{n},{p},{s},{k},{comp})[{x}]"

    if comp == -1:
        comp = 0

    if s == 0:
        ex = harmonics_pg["harmonics"][idx][0][comp]
    elif s < 4:
        ex = harmonics_pg[idx][1][comp]
    else:
        return

    info = BinaryManager("info")

    R = sp.symbols("X Y Z", real=True)
    rv = info["harmonics"]["variable"]
    sub = dict(zip(rv, R))
    ex = np.vectorize(lambda i: i.subs(sub))(ex)

    _plot_harmonics(qtdraw, ex, R, s, point, tag)


# ==================================================
def _plot_harmonics(qtdraw, ex, rv, s=0, point=None, tag=""):
    """
    Plot harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        ex (sympy): expression (internal)
        l (int): rank.
        s (int, optional): internal rank.
        k (int, optional): internal component.
        point (array-like, optional): spherical points.
        v (array-like, optional): weight vector.
    """
    if point is None:
        point = create_spherical_point()

    sv = sp.symbols("x y z", real=True)
    if s > 1:
        basis = create_internal_basis(s, sv, factor=False)
        ex = ex @ basis

    ex = np.asarray([np.vectorize(lambda i: sp.radsimp(i.subs({rv[0]: p[0], rv[1]: p[1], rv[2]: p[2]})))(ex) for p in point])
    if s < 2:
        ex = ex.astype(float)

    qtdraw.clear_data()
    if s == 0:
        _plot_scalar_harmonics(qtdraw, point, ex)
    elif s == 1:
        _plot_vector_harmonics(qtdraw, point, ex)
    else:
        _plot_multipolar_harmonics(qtdraw, point, ex)

    qtdraw.add_text2d(tag)
    qtdraw.set_view()


# ==================================================
def _plot_scalar_harmonics(qtdraw, point, tlm):
    """
    Plot scalar harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        point (ndarray): spherical points.
        tlm (ndarray): scalar harmonics at given points.
    """
    for s, v in zip(point, tlm):
        if v > 0:
            c = "salmon"
        elif v < 0:
            c = "aqua"
        else:
            c = "white"
            v = 0.01
        qtdraw.add_site(position=s, size=0.1 * abs(v), color=c)


# ==================================================
def _plot_vector_harmonics(qtdraw, point, tlm):
    """
    Plot vector harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        point (ndarray): spherical points.
        tlm (ndarray): vector harmonics at given points.
    """
    for s, v in zip(point, tlm):
        if np.linalg.norm(v) > 1e-6:
            qtdraw.add_vector(position=s, direction=v, length=-0.3)


# ==================================================
def _plot_multipolar_harmonics(qtdraw, point, tlm):
    """
    Plot multipolar harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        point (ndarray): spherical points.
        tlm (ndarray): multipolar harmonics at given points.
    """
    for s, q in zip(point, tlm):
        qtdraw.add_orbital(position=s, shape=str(q), size=0.1, color="coolwarm")


# ==================================================
def plot_site(qtdraw, site, rep=True):
    """
    Plot site cluster.

    Args:
        qtdraw (QtDraw): QtDraw object.
        site (ndarray): site cluster.
        rep (bool, optional): highlight for 0th site ?
    """
    for no, s in enumerate(site):
        color = "pink" if rep and no == 0 else "silver"
        qtdraw.add_site(position=s, color=color, size=0.07, label=f"s{no+1}", name="site_cluster")


# ==================================================
def plot_bond(qtdraw, bond):
    """
    Plot bond cluster.

    Args:
        qtdraw (QtDraw): QtDraw object.
        bond (ndarray): bond cluster.
    """
    for no, b in enumerate(bond):
        color = "white"
        v, c = b[0:3], b[3:6]
        qtdraw.add_bond(
            position=c,
            direction=v,
            color=color,
            color2=color,
            width=0.03,
            cartesian=False,
            opacity=0.7,
            label=f"b{no+1}",
            name="bond_cluster",
        )
        acolor = "red" if no == 0 else "black"
        qtdraw.add_vector(
            position=c,
            direction=v,
            length=-0.2,
            color=acolor,
            width=0.01,
            cartesian=False,
            label=f"b{no+1}",
            name="bond_cluster",
        )


# ==================================================
def plot_site_samb(qtdraw, site, samb):
    """
    Plot site-cluster SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        site (ndarray): site cluster.
        samb (ndarray): site-cluster samb.
    """
    if isinstance(samb, np.ndarray):
        samb = samb.astype(float)
    for s, v in zip(site, samb):
        if v > 0:
            c = "salmon"
        elif v < 0:
            c = "aqua"
        else:
            c = "silver"
            v = 0.5
        qtdraw.add_site(position=s, size=0.2 * abs(v), color=c, name="site_samb")


# ==================================================
def plot_bond_samb(qtdraw, bond, samb, sym=True):
    """
    Plot bond-cluster SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        bond (ndarray): bond cluster.
        samb (ndarray): bond-cluster samb.
        sym (bool, optional): symmetric bond ?
    """
    if isinstance(samb, np.ndarray):
        samb = samb.astype(float)
    if sym:
        for b, h in zip(bond, samb):
            v, c = b[0:3], b[3:6]
            if abs(h) < 1e-6:
                qtdraw.add_bond(
                    position=c, direction=v, color="white", color2="white", width=0.03, cartesian=False, name="bond_s_samb"
                )
            else:
                cl = "aqua" if h < 0 else "salmon"
                width = 0.07 * abs(h)
                qtdraw.add_bond(position=c, direction=v, color=cl, color2=cl, width=width, cartesian=False, name="bond_s_samb")
    else:
        for b, h in zip(bond, samb):
            v, c = b[0:3], b[3:6]
            if abs(h) < 1e-6:
                qtdraw.add_bond(
                    position=c,
                    direction=0.7 * v,
                    color="white",
                    color2="white",
                    width=0.03,
                    cartesian=False,
                    name="bond_s_samb",
                )
            else:
                if h < 0:
                    v = -v
                cl = "salmon"
                width = 0.07 * abs(h)
                qtdraw.add_vector(
                    position=c, direction=v, length=-0.7, color=cl, width=width, cartesian=False, name="bond_a_samb"
                )


# ==================================================
def plot_vector_samb(qtdraw, bond, samb):
    """
    Plot bond-cluster vector SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        bond (ndarray): bond cluster.
        samb (ndarray): bond-cluster vector samb.
    """
    if isinstance(samb, np.ndarray):
        samb = samb.astype(float)
    for b, h in zip(bond, samb):
        s = b[3:6]
        qtdraw.add_site(position=s, size=0.07, color="white", opacity=0.4)
        if not (np.abs(h) < 1e-6).all():
            qtdraw.add_vector(position=s, direction=h, length=-0.8, width=0.02, name="bond_v_samb")
