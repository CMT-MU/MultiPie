"""
Utility for plot.
"""

import numpy as np


# ==================================================
def plot_site(qtdraw, name, site, rep=True):
    """
    Plot site cluster.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): cluster name.
        site (ndarray): site cluster.
        rep (bool, optional): highlight for 0th site ?
    """
    for no, s in enumerate(site):
        color = "pink" if rep and no == 0 else "silver"
        qtdraw.add_site(position=s, color=color, size=0.07, label=f"s{no+1}", name=name)


# ==================================================
def plot_bond(qtdraw, name, bond, rep=True):
    """
    Plot bond cluster.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): cluster name.
        bond (ndarray): bond cluster.
        rep (bool, optional): highlight for 0th bond ?
    """
    opt = 0.7
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
            opacity=opt,
            label=f"b{no+1}",
            name=name,
        )
        acolor = "red" if rep and no == 0 else "black"
        qtdraw.add_vector(
            position=c, direction=v, length=-0.2, color=acolor, width=0.01, cartesian=False, label=f"b{no+1}", name=name
        )


# ==================================================
def plot_site_samb(qtdraw, name, site, samb, label=None):
    """
    Plot site-cluster SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): SAMB name.
        site (ndarray): site cluster.
        samb (ndarray): site-cluster samb.
        label (str, optional): SAMB label.
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
        qtdraw.add_site(position=s, size=0.2 * abs(v), color=c, name=name, label=label)


# ==================================================
def plot_bond_samb(qtdraw, name, bond, samb, sym=True, label=None):
    """
    Plot bond-cluster SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): SAMB name.
        bond (ndarray): bond cluster.
        samb (ndarray): bond-cluster samb.
        sym (bool, optional): symmetric bond ?
        label (str, optional): SAMB label.
    """
    if isinstance(samb, np.ndarray):
        samb = samb.astype(float)
    if sym:
        for b, h in zip(bond, samb):
            v, c = b[0:3], b[3:6]
            if abs(h) < 1e-6:
                qtdraw.add_bond(
                    position=c, direction=v, color="white", color2="white", width=0.03, cartesian=False, name=name, label=label
                )
            else:
                cl = "aqua" if h < 0 else "salmon"
                width = 0.07 * abs(h)
                qtdraw.add_bond(
                    position=c, direction=v, color=cl, color2=cl, width=width, cartesian=False, name=name, label=label
                )
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
                    name=name,
                    label=label,
                )
            else:
                if h < 0:
                    v = -v
                cl = "salmon"
                width = 0.07 * abs(h)
                qtdraw.add_vector(
                    position=c, direction=v, length=-0.7, color=cl, width=width, cartesian=False, name=name, label=label
                )


# ==================================================
def plot_vector_samb(qtdraw, name, bond, samb, label=None):
    """
    Plot bond-cluster vector SAMB.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): SAMB name.
        bond (ndarray): bond cluster.
        samb (ndarray): bond-cluster vector samb.
        label (str, optional): SAMB label.
    """
    if isinstance(samb, np.ndarray):
        samb = samb.astype(float)
    for b, h in zip(bond, samb):
        s = b[3:6]
        qtdraw.add_site(position=s, size=0.07, color="white", opacity=0.4, label=label)
        if not (np.abs(h) < 1e-6).all():
            qtdraw.add_vector(position=s, direction=h, length=-0.8, width=0.02, name=name, label=label)


# ==================================================
def plot_harmonics(qtdraw, name, samb, point=[[0, 0, 0]], label=None):
    """
    Plot harmonics.

    Args:
        qtdraw (QtDraw): QtDraw object.
        name (str): SAMB name.
        samb (sympy): samb expression.
        point (ndarray, optional): center of harmonics.
        label (str, optional): SAMB label.
    """
    point = np.asarray(point, dtype=float)
    for s in point:
        qtdraw.add_orbital(position=s, shape=str(samb), size=0.1, color="coolwarm", name=name, label=label)
