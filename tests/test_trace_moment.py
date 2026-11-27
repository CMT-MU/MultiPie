# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Tests for the trace moments (util_trace_moment, ModelAnalyzer.trace_moment).

The SnTe model and the reference polynomials are Eqs. (131), (140), (146), (150) of R. Oiwa and H. Kusunose,
J. Phys. Soc. Jpn. 91, 014701 (2022).
"""

import os

import numpy as np
import pytest
import sympy as sp

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict
from multipie.util.util_model_analyzer import fourier_r_to_k, grid_all
from multipie.util.util_trace_moment import TraceMoment, order_tuples

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "docs", "src", "examples")

# ==================================================
# two-orbital (px, py) x spin model of SnTe monolayer, Eq. (131).
_s0 = np.eye(2)
_sz = np.diag([1, -1])
_Q0 = np.eye(2)
_Qxy = np.array([[0, 1], [1, 0]])
_Qv = np.diag([1, -1])
_Mz = np.array([[0, -1j], [1j, 0]])


def _bond(b, orb, spin=_s0):
    # b: dict[r, value] in the convention O(k) = sum_r O(r) exp(-ik.r) of the paper, R = -r in MultiPie.
    M = np.kron(orb, spin)
    return {(-r[0], -r[1], 0, m, n): v * M[m, n] for r, v in b.items() for m in range(4) for n in range(4) if M[m, n] != 0}


_c = {(0, 0): 1}
_Qbx = {(1, 0): 1, (-1, 0): 1}
_Qby = {(0, 1): 1, (0, -1): 1}
_Tyby = {(0, 1): 1j, (0, -1): -1j}
_Qbxy = {(1, 1): 1, (-1, 1): 1, (-1, -1): 1, (1, -1): 1}
_Qxybxy = {(1, 1): 1, (-1, 1): -1, (-1, -1): 1, (1, -1): -1}
_Tybxy = {(1, 1): 1j, (-1, 1): 1j, (-1, -1): -1j, (1, -1): -1j}

SNTE = {
    "e1": _bond(_c, _Q0),
    "e2": _bond(_c, _Qv),
    "tx1": _bond(_Qbx, _Q0),
    "tx2": _bond(_Qbx, _Qv),
    "ty1": _bond(_Qby, _Q0),
    "ty2": _bond(_Qby, _Qv),
    "ty3": _bond(_Tyby, _Mz),
    "t1p": _bond(_Qbxy, _Q0),
    "t2p": _bond(_Qbxy, _Qv),
    "t3p": _bond(_Qxybxy, _Qxy),
    "t4p": _bond(_Tybxy, _Mz),
    "lam": _bond(_c, _Mz, _sz),
    "mz": _bond(_c, _Mz),
}
SNTE_PARAMETER = dict(e1=0.0, e2=0.05, tx1=-0.2, tx2=-0.05, ty1=-0.3, ty2=-0.1, ty3=-0.05)
SNTE_PARAMETER |= dict(t1p=0.1, t2p=0.05, t3p=0.05, t4p=0.05, lam=0.1, mz=0.0)


# ==================================================
@pytest.fixture(scope="module")
def snte():
    op = {"vx": "v_x", "vy": "v_y", "vxy": "v_xy", "Mz": _bond(_c, _Mz)}
    return TraceMoment(SNTE, np.zeros((4, 3)), np.eye(3), SNTE_PARAMETER, op)


# ==================================================
def _naive(tm, ops, order, grid=(24, 24, 1)):
    """
    (1/N) sum_k Tr[O_1 h^i1 ...] by matrix products.
    """
    kv = grid_all(*grid)
    atom = tm._atom
    par = tm._parameter

    def to_k(form):
        return sum(fourier_r_to_k(M, atom, kv, False) * (1.0 if key is None else par[key]) for key, M in form.items())

    H = to_k(tm._h)
    X = np.broadcast_to(np.eye(tm._dim), H.shape).astype(complex)
    for o, i in zip(ops, order):
        X = X @ to_k(tm._operator[o])
        for _ in range(i):
            X = X @ H
    return np.mean(np.trace(X, axis1=1, axis2=2))


# ==================================================
def test_snte_value(snte):
    g = snte.value(("vx", "vy", "vy"), 2)
    for order in [(0, 0, 0), (1, 0, 0), (0, 1, 1), (2, 1, 2)]:
        assert np.isclose(g[order], _naive(snte, ("vx", "vy", "vy"), order), atol=1e-12)
    g = snte.value(("vxy", "Mz"), 3)
    assert np.isclose(g[3, 2], _naive(snte, ("vxy", "Mz"), (3, 2)), atol=1e-12)
    # Gamma^{ij}_{A;B}* = Gamma^{ji}_{B;A}, Eq. (40).
    gab = snte.value(("vx", "Mz"), 3)
    gba = snte.value(("Mz", "vx"), 3)
    assert np.allclose(np.conj(gab), gba.T, atol=1e-12)
    # exact grid: a denser grid gives the same value.
    assert np.allclose(snte.value(("vx", "vy"), 3), snte.value(("vx", "vy"), 3, grid=(30, 30, 1)), atol=1e-12)


# ==================================================
def test_snte_chebyshev(snte):
    W, Ec = 2.0, 0.1
    gp = snte.value(("vx", "vy"), 3, scale=(W, Ec))
    gc = snte.value(("vx", "vy"), 3, polynomial="chebyshev", scale=(W, Ec))
    T = [[1, 0, 0, 0], [0, 1, 0, 0], [-1, 0, 2, 0], [0, -3, 0, 4]]  # T_i(x) = sum_l c_l x^l.
    assert np.allclose(gc, np.einsum("il,jm,lm->ij", T, T, gp), atol=1e-12)


# ==================================================
def test_snte_polynomial(snte):
    S = {k: sp.Symbol(k, real=True) for k in SNTE}
    e2, tx1, tx2, ty1, ty2, ty3 = S["e2"], S["tx1"], S["tx2"], S["ty1"], S["ty2"], S["ty3"]
    t1, t2, t3, t4, mz = S["t1p"], S["t2p"], S["t3p"], S["t4p"], S["mz"]
    var = [k for k in SNTE if k != "mz"]

    # nonlinear Hall effect, Eq. (140) (Tr[v_x h v_y v_y]).
    r = snte.lowest(("vx", "vy", "vy"), n_term=2, max_order=4, part="im", variable=var)
    ref = -128 * t3 * (t2 * ty3 * ty1 + ty2 * ty3 * t1 + 2 * ty2 * t4 * ty1 + 4 * t2 * t4 * t1)
    assert [x["order"] for x in r] == [(0, 0, 1), (1, 0, 0)]
    assert sp.expand(r[1]["polynomial"] - ref) == 0
    assert sp.expand(r[0]["polynomial"] + ref) == 0
    assert set(r[1]["variable"]) == {"t1p", "t2p", "t3p", "t4p", "ty1", "ty2", "ty3"}
    assert r[1]["essential"] == ["t3p"]  # Table II (t_y3 or t'_4 is essential as a pair).
    r = snte.lowest(("vx", "vy", "vy"), n_term=1, max_order=4, part="im", variable=var, expand=False)
    assert r[0]["polynomial"] is None and r[0]["essential"] == ["t3p"]
    # Gamma^{010} = Tr[v_x v_y h v_y] is real.
    assert snte.polynomial(("vx", "vy", "vy"), (0, 1, 0), part="im", variable=var) == 0

    # orbital magneto-current effect, Eq. (146) (Tr[M_z h^2 v_y]).
    p = snte.polynomial(("Mz", "vy"), (2, 0), part="re", variable=var)
    ref = 16 * (ty3 * (e2 * ty2 + 2 * tx2 * t2) + 2 * t4 * (e2 * t2 + tx2 * ty2))
    assert sp.expand(p - ref) == 0

    # anomalous Hall effect under m_z, Eq. (150).
    r = snte.lowest(("vx", "vy"), n_term=1, max_order=4, part="im")
    assert r[0]["order"] == (0, 2)
    p = snte.polynomial(("vx", "vy"), (2, 0), part="im")
    assert sp.expand(p - 32 * mz * t3 * (tx1 * ty2 - ty1 * tx2)) == 0

    # the polynomial at the current parameter is the numerical value.
    p = snte.polynomial(("vx", "Mz"), (1, 2))
    v = complex(p.subs({S[k]: v for k, v in SNTE_PARAMETER.items()}))
    assert np.isclose(v, snte.value(("vx", "Mz"), 2)[1, 2], atol=1e-12)
    # fixed parameters.
    p = snte.polynomial(("vx", "Mz"), (1, 2), variable=["tx1", "ty3"])
    v = complex(p.subs({S["tx1"]: SNTE_PARAMETER["tx1"], S["ty3"]: SNTE_PARAMETER["ty3"]}))
    assert np.isclose(v, snte.value(("vx", "Mz"), 2)[1, 2], atol=1e-12)


# ==================================================
def test_order_tuples():
    assert order_tuples(2, 2) == [(0, 2), (1, 1), (2, 0)]
    assert len(order_tuples(3, 3)) == 10


# ==================================================
@pytest.fixture(scope="module")
def graphene(tmp_path_factory):
    model = read_dict(os.path.abspath(os.path.join(EXAMPLE_DIR, "graphene_in.py")))
    model["pdf"] = {"create": False}
    model["qtdraw"] = {"create": False}
    control = read_dict(os.path.abspath(os.path.join(EXAMPLE_DIR, "graphene_ctrl.py")))
    control["samb"]["samb_figure"] = False
    control["samb"]["k_multipole"] = False
    control["grid"] = (4, 4, 1)
    control["output"] = {}

    topdir = str(tmp_path_factory.mktemp("graphene_trace"))
    cwd = os.getcwd()
    try:
        create_model(model, topdir=topdir)
        os.chdir(topdir)
        ma = ModelAnalyzer(topdir)
        ma.analyze(control)
    finally:
        os.chdir(cwd)
    return topdir, ma


# ==================================================
def _write_hr(filename, OR, dim):
    """
    Write an operator in the format of seedname_hr.dat (degeneracy 1).
    """
    Rs = sorted({k[:3] for k in OR} | {tuple(-i for i in k[:3]) for k in OR})
    with open(filename, "w") as f:
        print("test", file=f)
        print(dim, file=f)
        print(len(Rs), file=f)
        print(" ".join(["1"] * len(Rs)), file=f)
        for R in Rs:
            for n in range(dim):
                for m in range(dim):
                    v = complex(OR.get((*R, m, n), 0))
                    print(f"{R[0]} {R[1]} {R[2]} {m + 1} {n + 1} {v.real:.15e} {v.imag:.15e}", file=f)


# ==================================================
def test_graphene(graphene):
    topdir, ma = graphene
    tm = ma.trace_moment({"vx": "v_x", "vy": "v_y"})
    assert tm.variable == list(ma.get_samb_matrix().keys())

    # velocity: d_k h(k) (Cartesian) by finite difference in the gauge with the positions, k_frac = A k_cart / 2pi.
    HR, atom, _ = ma.get_hk_data()
    A = np.asarray(ma["info"]["A"])
    k = np.array([[0.13, 0.27, 0.0]])
    dk = 1e-6
    kc = A @ np.array([dk, 0, 0]) / (2 * np.pi)
    dH = (fourier_r_to_k(HR, atom, k + kc, True) - fourier_r_to_k(HR, atom, k - kc, True)) / (2 * dk)
    vx = sum(tm._parameter[z] * fourier_r_to_k(M, atom, k, True) for z, M in tm._operator["vx"].items())
    assert np.allclose(vx, dH, atol=1e-6)

    g = tm.value(("vx", "vy"), 3)
    for order in [(0, 0), (1, 1), (3, 2)]:
        assert np.isclose(g[order], _naive(tm, ("vx", "vy"), order, grid=(20, 20, 1)), atol=1e-12)

    # operator from a file of seedname_hr.dat and from a matrix dict.
    vx_r = {(*R, m, n): v * tm._parameter["z2"] for (R, m, n), v in tm._operator["vx"]["z2"].items()}
    _write_hr(os.path.join(topdir, "vx_hr.dat"), vx_r, tm._dim)
    tm.add_operator("vx_file", "vx_hr.dat")
    tm.add_operator("vx_dict", vx_r)
    g1 = tm.value(("vx_file", "vy"), 2)
    g2 = tm.value(("vx_dict", "vy"), 2)
    assert np.allclose(g1, tm.value(("vx", "vy"), 2), atol=1e-12)
    assert np.allclose(g2, g1, atol=1e-12)

    # Drude weight: Tr[v_x h^0 v_x] is the lowest term; the on-site SAMB (z1) does not appear.
    r = tm.lowest(("vx", "vx"), n_term=1, part="re")
    assert r[0]["order"] == (0, 0)
    assert "z1" not in r[0]["variable"]
    v = r[0]["polynomial"].subs({sp.Symbol(z, real=True): tm._parameter[z] for z in r[0]["variable"]})
    assert np.isclose(float(v), tm.value(("vx", "vx"), 0, part="re")[0, 0], atol=1e-12)

    # exact coefficients (the cell vectors are kept in full precision, sqrt(3)/2 not rounded): Drude weight and the
    # anomalous Hall effect of the Haldane type (z4: the imaginary hopping of the second neighbors, A2g).
    z1, z2, z3, z4 = (sp.Symbol(f"z{i}", real=True) for i in range(1, 5))
    var = ["z1", "z2", "z3", "z4"]
    r = tm.lowest(("vx", "vx"), n_term=1, part="re", variable=var)
    assert sp.expand(r[0]["polynomial"] - (z2**2 + 3 * z3**2 + 3 * z4**2) / 6) == 0
    r = tm.lowest(("vx", "vy"), n_term=2, part="im", variable=var)
    assert [x["order"] for x in r] == [(0, 1), (1, 0)]
    assert sp.expand(r[0]["polynomial"] - z2**2 * z4 / 4) == 0
    assert sp.expand(r[1]["polynomial"] + z2**2 * z4 / 4) == 0
    assert r[0]["essential"] == ["z2", "z4"]
