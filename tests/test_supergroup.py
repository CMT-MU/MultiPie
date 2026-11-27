# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Supergroups other than the direct parent (util_supergroup), and S_kj = <Z_k, g Z_j>, F_j = sum_k (delta_kj - S_kj) <Z_k>
(MaterialModel.image_overlap, parent_forbidden_part).

Reference values: right-handed Te (P3_121, No.152) with the parents P-3m1 (No.164, x -> 0) and R-3m (No.166, x -> 1/3),
and No.198 (P2_13) with Pa-3 (No.205), with Slater-Koster models (p orbitals of Te, s and p of No.198).
"""

import contextlib
import io
import itertools
import logging

import numpy as np
import pytest

from multipie import Group
from multipie.core.material_model import MaterialModel
from multipie.util import util_supergroup as us

TE = {
    "group": "SG:152",
    "cell": {"a": 4.457, "c": 5.929},
    "site": {"Te": ("[0.274,0,1/3]", ["p"])},
    "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A1"]},
    "atomic_select": {"X": ["Q"]},
}
N198 = {
    "group": "SG:198",
    "cell": {"a": 4.5},
    "site": {"Co": ("[0.14,0.14,0.14]", ["s", "p"])},
    "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A"]},
}


def _build(base, bonds, name="M", **kw):
    m = {"model": name, "spinful": False, "qtdraw": {"create": False}, "pdf": {"create": False}, "coset": False}
    m.update(base)
    m["bond"] = bonds
    m.update(kw)
    mm = MaterialModel()
    with contextlib.redirect_stdout(io.StringIO()):
        mm.analyze(m)
    return mm


@pytest.fixture(scope="module")
def te164():
    return _build(TE, [("Te", "Te", [1, 2, 6, 10])], "te164")


@pytest.fixture(scope="module")
def te166():
    return _build(TE, [("Te", "Te", [1, 2])], "te166")


@pytest.fixture(scope="module")
def n198():
    return _build(N198, [("Co", "Co", [1, 3])], "n198")


# ==================================================
# Slater-Koster model and expectation values (label space {(a, b, n): <a,0|O|b,n>}).
class _TB:
    def __init__(self, mm):
        self.A = np.array(mm["unit_vector"], float).T
        kets = mm["full_matrix"]["ket"]
        self.site_of = [(k[0], k[1]) for k in kets]
        self.orb = [k[4] for k in kets]
        self.pos = {
            (n, v.sublattice): np.array(v.position, float) for n, c in mm["site"]["cell"].items() for v in c if v.plus_set == 1
        }
        d = mm.get_combined_samb_matrix(fmt="value", bond=False)
        self.Z = {z: {(a, b, (n1, n2, n3)): complex(v) for (n1, n2, n3, a, b), v in dd.items()} for z, dd in d.items()}
        self.nk = len(kets)

    def positions(self, formula, x0):
        out = {}
        for s, p in self.pos.items():
            for f in formula:
                d = p - f(x0)
                if np.allclose(d, np.rint(d), atol=1e-6):
                    out[s] = (f, np.rint(d))
                    break
        return lambda x: {s: f(x) + off for s, (f, off) in out.items()}

    def hamiltonian(self, pos, V, cutoff, onsite):
        ax = {"px": 0, "py": 1, "pz": 2}

        def sk(o1, o2, l, v):
            if o1 == "s" and o2 == "s":
                return v["sss"]
            if o1 == "s":
                return l[ax[o2]] * v["sps"]
            if o2 == "s":
                return -l[ax[o1]] * v["sps"]
            i, j = ax[o1], ax[o2]
            return l[i] * l[j] * (v["pps"] - v["ppp"]) + (v["ppp"] if i == j else 0)

        H = {(a, a, (0, 0, 0)): onsite[self.orb[a]] for a in range(self.nk)}
        for a, b in itertools.product(range(self.nk), repeat=2):
            ra, rb = pos[self.site_of[a]], pos[self.site_of[b]]
            for n in itertools.product(range(-2, 3), repeat=3):
                dv = self.A @ (rb + np.array(n) - ra)
                d = np.linalg.norm(dv)
                if 1e-6 < d <= cutoff:
                    H[(a, b, n)] = H.get((a, b, n), 0) + sk(self.orb[a], self.orb[b], dv / d, V(d))
        return H

    def _k(self, O, k):
        M = np.zeros((self.nk, self.nk), complex)
        for (a, b, n), v in O.items():
            M[a, b] += v * np.exp(2j * np.pi * np.dot(k, n))
        return M

    def expectation(self, H, nel, nkg, T=0.02):
        ks = [np.array(k) / nkg for k in itertools.product(range(nkg), repeat=3)]
        data = [np.linalg.eigh(self._k(H, k)) for k in ks]
        E = np.concatenate([e for e, _ in data])
        lo, hi = E.min() - 1, E.max() + 1
        f = lambda mu: np.sum(1 / (np.exp(np.clip((E - mu) / T, -60, 60)) + 1)) / len(ks)
        for _ in range(100):
            mu = (lo + hi) / 2
            lo, hi = (mu, hi) if f(mu) < nel else (lo, mu)
        rhos = []
        for e, u in data:
            occ = 1 / (np.exp(np.clip((e - mu) / T, -60, 60)) + 1)
            rhos.append((u * occ) @ u.conj().T)
        return lambda O: (sum(np.sum(r.T * self._k(O, k)) for k, r in zip(ks, rhos)) / len(ks)).real


F_TE = [lambda x: np.array([x, 0, 1 / 3]), lambda x: np.array([0, x, 2 / 3]), lambda x: np.array([-x, -x, 0])]
V_TE = lambda d: {"pps": 2.5 * np.exp(-2 * (d - 2.9)), "ppp": -0.6 * np.exp(-2 * (d - 2.9)), "sss": 0, "sps": 0}
F_198 = [
    lambda x: np.array([x, x, x]),
    lambda x: np.array([-x + 0.5, -x, x + 0.5]),
    lambda x: np.array([-x, x + 0.5, -x + 0.5]),
    lambda x: np.array([x + 0.5, -x + 0.5, -x]),
]
V_198 = lambda d: {k: v * np.exp(-2 * (d - 2.76)) for k, v in {"sss": -1.0, "sps": 1.2, "pps": 1.5, "ppp": -0.4}.items()}


def _te_state(tb, x):
    return tb.expectation(tb.hamiltonian(tb.positions(F_TE, 0.274)(x), V_TE, 3.7, {"px": 0, "py": 0, "pz": 0}), 6, 10)


def _label(gz):
    return {z: {(a, b, (n1, n2, n3)): complex(v) for (n1, n2, n3, a, b), v in d.items()} for z, d in gz.items()}


# ==================================================
def test_supergroups():
    r = Group("SG:152").supergroups(max_index=6)
    d = {x["group"]: x for x in r}
    assert d["SG:155"]["direct"] and (d["SG:155"]["index"], d["SG:155"]["index_parent"]) == (3, 1)
    assert d["SG:166"]["path"] == ["SG:152", "SG:155", "SG:166"] and (d["SG:166"]["index"], d["SG:166"]["index_parent"]) == (6, 2)
    assert d["SG:164"]["path"] == ["SG:152", "SG:150", "SG:164"] and d["SG:164"]["index"] == 6
    assert [x["group"] for x in Group("SG:152").supergroups(chain=True)] == ["SG:164", "SG:191"]
    assert {x["group"] for x in Group("SG:198").supergroups(max_index=2)} >= {"SG:205"}
    assert Group("SG:152").embedding("SG:143") == []  # a subgroup.


# ==================================================
def test_embedding():
    # 152 -> 166: R = 1, origin (2/3, 2/3, 0) (the atom 3 of x = 1/3); {3+|0,0,1/3} of H is {3+|tau} of R-3m.
    o = np.array([2 / 3, 2 / 3, 0])
    A = np.eye(4)
    A[:3, 3] = -o
    e = us.resolve_embedding("SG:152", "SG:166", A)
    assert np.allclose(e["R"], np.eye(3))
    W = np.eye(4)
    W[:3, :3] = [[0, -1, 0], [1, -1, 0], [0, 0, 1]]
    W[:3, 3] = [0, 0, 1 / 3]
    t = (A @ W @ np.linalg.inv(A))[:3, 3]
    assert np.allclose(t % 1, [2 / 3, 1 / 3, 1 / 3])
    assert len(Group("SG:152").embedding("SG:166")) == 2  # the other: x -> 2/3.
    # 152 -> 164 (the representative parent of the binaries): R = diag(1, 1, 3).
    assert np.allclose(Group("SG:152").embedding("SG:164")[0]["R"], np.diag([1, 1, 3]))


# ==================================================
def test_coset_decomposition():
    G = Group("SG:152")
    assert len(G.coset_decomposition("SG:164")) == 2
    cs = G.coset_decomposition("SG:164", modulo="self")
    assert len(cs) == 6 and "1:[0,0,1/3]" in [c["rep"] for c in cs]
    cs = G.coset_decomposition("SG:166", modulo="self")
    assert len(cs) == 6 and [c["kind"] for c in cs].count("translation") == 2 and [c["kind"] for c in cs].count("inversion") == 3
    cs = Group("SG:198").coset_decomposition("SG:205")
    assert [(c["rep"], c["kind"], c["order"]) for c in cs] == [("1:[0,0,0]", "identity", 1), ("-1:[0,0,0]", "inversion", 2)]


# ==================================================
def test_parent_limit():
    G = Group("SG:152")
    assert G.parent_limit("SG:164", "3a")["3a"] == [{"param": {"x": "0"}, "position": "[0,0,1/3]", "parent_wyckoff": "1a"}]
    lims = [G.parent_limit("SG:166", "3a", k)["3a"][0] for k in range(2)]
    assert [(d["param"]["x"], d["parent_wyckoff"]) for d in lims] == [("1/3", "3a"), ("2/3", "3b")]
    lims = Group("SG:198").parent_limit("SG:205", "4a")["4a"]
    assert [(d["param"]["x"], d["parent_wyckoff"]) for d in lims] == [("0", "4a"), ("1/2", "4b")]


# ==================================================
def test_map_site_bond(te166, te164):
    # g exchanges the atoms 1 and 2 and keeps the atom 3 for both parents.
    G = Group("SG:152")
    for p, op in [("SG:164", "-1:[0,0,0]"), ("SG:166", "-1:[1/3,1/3,0]")]:
        assert {i: j for i, (j, _) in G.map_site(op, "3a", p).items()} == {1: 2, 2: 1, 3: 3}
    # the nearest neighbors moved by g: the 1st neighbors (164); the 1st and the 2nd neighbors (166).
    nb = lambda mm, p: sorted(d["neighbor"] for d in mm.map_bond(parent=p)["Te;Te_001_1"])
    assert nb(te164, "SG:164") == [1, 1]
    assert nb(te166, "SG:166") == [1, 2]


# ==================================================
def test_parent_shells(te164, te166, n198):
    sh = te164.parent_shells(parent="SG:164")[("Te", "Te")]
    assert [2, 6, 10] in [d["neighbor"] for d in sh] and all(d["complete"] for d in sh)
    assert te166.parent_shells(parent="SG:166")[("Te", "Te")][0]["neighbor"] == [1, 2]
    assert n198.parent_shells(parent="SG:205")[("Co", "Co")][0]["neighbor"] == [1, 3]


# ==================================================
def test_image_overlap(te164, te166, n198):
    count = lambda r: {v: int(np.sum(np.isclose(np.real(r["eigenvalue"]), v))) for v in (1, 0, -1)}
    r = te164.image_overlap(parent="SG:164")
    assert len(r["samb"]) == 26 and r["complete"] and count(r) == {1: 10, 0: 16, -1: 0}
    r = te166.image_overlap(parent="SG:166")
    assert len(r["samb"]) == 14 and r["complete"] and count(r) == {1: 6, 0: 8, -1: 0}
    i = {z: k for k, z in enumerate(r["samb"])}
    S = r["S"]
    assert np.allclose(S[np.ix_([i["z5"], i["z9"]], [i["z5"], i["z9"]])], [[1 / 3, np.sqrt(2) / 3], [np.sqrt(2) / 3, 2 / 3]])
    col = {z: S[i[z], i["z7"]] for z in r["samb"] if abs(S[i[z], i["z7"]]) > 1e-8}
    assert set(col) == {"z7", "z10", "z12"} and np.allclose(
        [col["z7"], col["z10"], col["z12"]], [1 / 3, -0.3086, 0.3563], atol=1e-4
    )
    r = n198.image_overlap(parent="SG:205")
    assert len(r["samb"]) == 36 and r["complete"] and count(r)[0] == 0
    i = {z: k for k, z in enumerate(r["samb"])}
    assert np.isclose(r["S"][i["z21"], i["z5"]], 1)  # g Z_5 = Z_21 (s-s Q0 of the 3rd neighbors).
    # the diagonal is w_+ - w_- of "weight"; incomplete model (bonds [1, 2] for P-3m1).
    mm = _build(TE, [("Te", "Te", [1, 2])], coset=None)
    w = mm["coset"][0]["weight"]
    assert np.allclose(np.diag(mm.image_overlap()["S"]), [w[z][0] - w[z][1] for z in mm["combined_id"]])
    assert mm["coset"][0]["complete"] is False and not mm.image_overlap(parent="SG:164")["complete"]


# ==================================================
@pytest.mark.parametrize(
    "case, x, expected",
    [
        ("te164", 0.274, {"z3": 0.3378, "z4": -0.3317, "z7": -0.9093, "z8": 0.5268, "z9": -0.1333, "z12": -0.1396}),
        ("te166", 0.274, {"z5": -0.1565, "z7": -0.4669, "z9": 0.1107, "z12": 0.2322}),
        ("te166", 1 / 3, None),
        ("n198", 0.14, {"z2": -0.0362, "z5": 0.1737}),
    ],
)
def test_parent_forbidden_part(case, x, expected, request):
    mm = request.getfixturevalue(case)
    parent = {"te164": "SG:164", "te166": "SG:166", "n198": "SG:205"}[case]
    tb = _TB(mm)
    if case == "n198":
        pos = tb.positions(F_198, 0.14)(x)
        E = tb.expectation(tb.hamiltonian(pos, V_198, 4.6, {"s": -6.0, "px": 1.0, "py": 1.0, "pz": 1.0}), 4, 6)
    else:
        E = _te_state(tb, x)
    e = {z: E(tb.Z[z]) for z in tb.Z}
    F = mm.parent_forbidden_part(e, parent=parent)
    if expected is None:  # parent structure.
        assert all(abs(v) < 1e-8 for v in F.values())
    else:
        assert all(abs(F[z] - v) < 1e-4 for z, v in expected.items())
    # F_j = chi_j = <Z_j> - <g Z_j> with g Z_j of transform_samb.
    gz = _label(mm.transform_samb(parent=parent))
    assert max(abs(F[z] - (e[z] - E(gz[z]))) for z in tb.Z) < 1e-9


# ==================================================
def test_bond_order(te166):
    # F_5 (R-3m) = 2 sqrt(2)/3 (b1 - b2), b_n: bond order of the n-th neighbors.
    tb = _TB(te166)
    L = lambda a, b, n: np.linalg.norm(tb.A @ (tb.pos[tb.site_of[b]] + np.array(n) - tb.pos[tb.site_of[a]]))
    keys = {k for z in tb.Z for k in tb.Z[z]}
    d1 = min(L(*k) for k in tb.Z["z5"])
    d2 = min(L(*k) for k in tb.Z["z9"])
    O = [{k: 1.0 for k in keys if tb.orb[k[0]] == tb.orb[k[1]] and abs(L(*k) - d) < 1e-3} for d in (d1, d2)]
    nb = [len({(tb.site_of[a], tb.site_of[b], n) for (a, b, n) in o}) for o in O]
    for x in (0.274, 0.30, 0.36):
        E = _te_state(tb, x)
        e = {z: E(tb.Z[z]) for z in tb.Z}
        F = te166.parent_forbidden_part(e, parent="SG:166")
        b1, b2 = E(O[0]) / nb[0], E(O[1]) / nb[1]
        assert np.isclose(F["z5"] / (b1 - b2), 2 * np.sqrt(2) / 3, atol=1e-6)


# ==================================================
@pytest.mark.parametrize("parent", ["SG:164", "SG:166"])
def test_inversion_model(parent):
    # <Z'_k, g Z_j> = delta_jk eps_j (Te: all +1) with the atoms of the inversion model matched on the parent structure.
    mm = _build(TE, [("Te", "Te", [1, 2])], "TeR")
    mi = _build(TE, [("Te", "Te", [1, 2])], "TeR", create_inversion_model=True)
    rc, _ = mm._coset_space(None, parent, None)
    g = rc["g"]
    tb, ti = _TB(mm), _TB(mi)
    img = us.map_site([rc["lim_pos"][s] for s in tb.pos], g, np.eye(3))
    sites = list(tb.pos)
    # the image crystal g(model) is the inversion model translated by T.
    q = {s: g[:3, :3] @ tb.pos[s] + g[:3, 3] for s in sites}
    lab = None
    for s0 in ti.pos:
        T = q[sites[0]] - ti.pos[s0]
        cand = {}
        for s in sites:
            for s2, p2 in ti.pos.items():
                Ln = q[s] - T - p2
                if np.allclose(Ln, np.rint(Ln), atol=1e-6):
                    j, m = img[sites.index(s)]
                    cand[s2] = (sites[j], m - np.rint(Ln))
        if len(cand) == len(sites):
            lab = cand
            break
    assert lab is not None
    idx = {(tb.site_of[a], tb.orb[a]): a for a in range(tb.nk)}

    def tr(O):
        out = {}
        for (a, b, n), v in O.items():
            (sa, ma), (sb, mb) = lab[ti.site_of[a]], lab[ti.site_of[b]]
            key = (idx[(sa, ti.orb[a])], idx[(sb, ti.orb[b])], tuple(int(x) for x in mb + np.array(n) - ma))
            out[key] = out.get(key, 0) + v
        return out

    gz = _label(mm.transform_samb(parent=parent))
    zp = {z: tr(ti.Z[z]) for z in ti.Z}
    A = np.array([[sum(np.conj(v) * gz[j].get(k, 0) for k, v in zp[kk].items()).real for kk in ti.Z] for j in tb.Z])
    assert np.allclose(A, np.eye(len(A)), atol=1e-8)


# ==================================================
def test_coset_input(caplog):
    mm = _build(TE, [("Te", "Te", [1, 2])], coset={"parent": "SG:166"})
    (r,) = mm["coset"]
    assert r["parent"] == "SG:166" and r["kind"] == "inversion" and r["image_group"] == "SG:154"
    assert r["complete"] and r["embedding"]["path"] == ["SG:152", "SG:155", "SG:166"]
    assert np.allclose(sorted(r["S_eigenvalue"]), [0] * 8 + [1] * 6, atol=1e-8)
    assert all(v == [(z, 1.0)] for z, v in r["samb"].items())
    # incomplete: warning.
    mi = _build(TE, [("Te", "Te", [1, 2])])
    with caplog.at_level(logging.WARNING):
        mi.parent_forbidden_part({z: 0 for z in mi["combined_id"]}, parent="SG:164")
    assert "neighbor shells" in caplog.text


# ==================================================
def test_order_3():
    # D4h in Oh, g = 3+[111] (n = 3): F_j^(m) = n <Pi_m Z_j>, F_j = sum_{m>0} F_j^(m) = n <Z_j> - sum_t <g^t Z_j>.
    md = _build(
        {"group": "PG:15", "site": {"A": ("[0.30,0.30,0.25]", ["s"])}, "SAMB_select": {"X": ["Q", "G"], "Gamma": ["A1g"]}},
        [("A", "A", [1, 2])],
    )
    r = md.image_overlap(parent="PG:32")
    assert r["order"] == 3 and r["complete"] and len(r["projection"]) == 3
    zs = r["samb"]
    # a state invariant under H: rho = exp(-H) / Tr, H = sum_j c_j Z_j.
    d = md.get_combined_samb_matrix(fmt="value", bond=False)
    nk = len(md["full_matrix"]["ket"])

    def mat(O):
        M = np.zeros((nk, nk), complex)
        for (n1, n2, n3, a, b), v in O.items():
            M[a, b] += complex(v)
        return M

    rng = np.random.default_rng(3)
    H = sum(c * mat(d[z]) for c, z in zip(rng.normal(size=len(zs)), zs))
    w, u = np.linalg.eigh(H)
    rho = (u * np.exp(-w)) @ u.conj().T
    rho /= np.trace(rho)
    E = lambda O: np.sum(rho.T * mat(O)).real
    e = {z: E(d[z]) for z in zs}
    Fm = md.parent_forbidden_part(e, parent="PG:32", by_eigenvalue=True)
    F = md.parent_forbidden_part(e, parent="PG:32")
    g1 = md.transform_samb(parent="PG:32")
    g2 = md.transform_samb(op="3-[111]", parent="PG:32")  # g^2 = g^-1.
    for z in zs:
        assert np.isclose(sum(Fm[z]), 3 * e[z])
        assert np.isclose(F[z], 3 * e[z] - e[z] - E(g1[z]) - E(g2[z]), atol=1e-10)
        assert np.isclose(Fm[z][1], np.conj(Fm[z][2]))  # w^1 and w^2 = w^-1.
