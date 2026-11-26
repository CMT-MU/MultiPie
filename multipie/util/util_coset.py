# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Coset of a model: the SAMB of a model of H moved by a representative g of a coset g H in a parent P.

The data of the cosets are those of the group binaries (MultiPieGen): "coset" of the group file (the first coset in
the direct parent) and "parent", "coset" of info (the parents and the cosets in them).

- Reduction of g (in the coordinates of H).

  - "rotation" (det g > 0): the image system is written in the frame moved by g; its data are those of the model, and
    g Z_j = +Z'_j (the SAMB of the image system with the same label).
  - "improper" (det g < 0, H not chiral): g = g' m^-1 with an improper operation m of H and the rotation g' = g m:
    g Z_j = sum_k D(m^-1)_kj Z'_k in the frame moved by g' (D = 1 for the identity irrep).
  - "inversion" (det g < 0, H chiral): g = r I with the inversion I of the first coset in the direct parent and a
    rotation r: I maps H to H' = I H I^-1 (map M), the sites and bonds of each Wyckoff position W to those of W'
    (row i -> row perm[i], bond reversed if sign = -1), and the cluster SAMB by eps (label by label) or by U (within an
    irrep). The atomic SAMB of the orbitals of ranks (l1, l2) are multiplied by (-1)^(l1+l2). The frame is moved by r.

- Weight. On the parent structure (the limit of the Wyckoff parameters at which the structure is invariant under the
  direct parent, "parent_limit"), g is a symmetry, and the weight of the eigenvalue w^m (w = exp(2 pi i/n), n the
  order of g modulo H) of Z_j is w_m(j) = (1/n) sum_t w^(-m t) <Z_j, g^t Z_j>/<Z_j, Z_j>. Given only for the direct
  parent.
"""

import itertools

import numpy as np
import sympy as sp

from multipie.core.group import Group
from multipie.util.util_crystal import P_dict

TOL = 1e-6
DIG = 6
_x, _y, _z = sp.symbols("x y z", real=True)


# ==================================================
def _float(M):
    return np.array(sp.Matrix(M).evalf(), dtype=float)


# ==================================================
def _mat4(M):
    M = np.asarray(M, dtype=float)
    if M.shape == (3, 3):
        W = np.eye(4)
        W[:3, :3] = M
        return W
    return M


# ==================================================
def _lattice(group):
    """primitive vectors (columns) in the conventional coordinates, None for PG."""
    if group.is_point_group:
        return None
    return np.array(P_dict[group.info.lattice], dtype=float)[:3, :3]


# ==================================================
def _reduce(v, L):
    """v modulo the lattice L (None: no lattice)."""
    v = np.asarray(v, dtype=float)
    if L is None:
        return np.where(np.abs(v) < 1e-9, 0.0, v)
    f = np.linalg.solve(L, v)
    f = f - np.floor(f + TOL)
    f = np.where(np.abs(f) < 1e-9, 0.0, f)
    return L @ f


# ==================================================
def _okey(W, L):
    return (tuple(np.rint(W[:3, :3] * 1e6).astype(int).ravel()), tuple(np.rint(_reduce(W[:3, 3], L) * 1e5).astype(int)))


# ==================================================
def _ops(group):
    """operations (4x4 float, conventional) of a group, with the plus set."""
    ops = [_mat4(_float(W)) for W in group.symmetry_operation["fractional"]]
    if group.is_point_group:
        return ops
    out = []
    for c in group.symmetry_operation["plus_set"]:
        for W in ops:
            W = W.copy()
            W[:3, 3] += np.array(c, dtype=float)
            out.append(W)
    return out


# ==================================================
def _parse_tag(tag):
    """ "point part:[t1,t2,t3]" -> (point part, translation or None)."""
    if ":" not in tag:
        return tag, None
    p, t = tag.split(":")
    return p, np.array([float(sp.sympify(v)) for v in t.strip("[]").split(",")])


# ==================================================
def _op_matrix(group, parent, R, tag):
    """the operation tag of the parent in the coordinates of H (4x4 float); SG: translation of the tag (in H)."""
    p, t = _parse_tag(tag)
    Pg = Group(parent)
    for W, tg in zip(Pg.symmetry_operation["fractional"], Pg.symmetry_operation["tag"]):
        if tg.split(":")[0] == p:
            W = np.linalg.inv(R) @ _mat4(_float(W)) @ R
            if t is not None:
                W[:3, 3] = t
            return W
    raise ValueError(f"symmetry operation {tag} is not found in {parent}.")


# ==================================================
def group_id(group):
    """id of a PG or SG, e.g. "PG:32", "SG:152"."""
    return group.info.PG if group.group_type == "PG" else group.info.SG


# ==================================================
def parse_coset(group, spec):
    """
    Representatives of the cosets to analyze.

    Args:
        group (Group): group H of the model.
        spec (None, False, dict or list): "coset" of the model input.

    Returns:
        - (list) -- [(parent id, representative tag)], [] for False or a group without a parent.

    Note:
        - None: the representative of the first coset in the direct parent.
        - {"parent": id, "rep": tag}: the operation must be one of the cosets of info["coset"][H][parent] (not a
          translation only, nor an operation of H).
    """
    if spec is False:
        return []
    first = group.coset
    if first is None:
        return []
    if spec is None:
        return [(first["parent"], first["rep"])]
    specs = [spec] if isinstance(spec, dict) else list(spec)
    info = Group.global_info()
    gid = group_id(group)
    out = []
    for s in specs:
        if not isinstance(s, dict) or set(s) - {"parent", "rep"} or "rep" not in s:
            raise ValueError(f"coset must be None, False, {{'parent': id, 'rep': tag}} or a list of them, but {s!r} is given.")
        parent = s.get("parent", first["parent"])
        parent = group_id(Group(parent))
        cosets = info["coset"][gid].get(parent)
        if cosets is None:
            raise ValueError(f"{parent} is not a parent of {gid}: {info['parent'][gid]}.")
        if not any(s["rep"] in c for c in cosets):
            raise ValueError(f"{s['rep']} is not an operation of the cosets of {gid} in {parent}.")
        out.append((parent, s["rep"]))
    return out


# ==================================================
def _chain_R(group, parent):
    """x_parent = R x_H along the chain of parents (4x4 float)."""
    info = Group.global_info()
    R = np.eye(4)
    gid = group_id(group)
    for p in info["parent"][gid]:
        R = _mat4(_float(Group(gid).coset["R"])) @ R
        gid = p
        if p == parent:
            return R
    raise ValueError(f"{parent} is not a parent of {group_id(group)}.")


# ==================================================
def _order(g, Hkeys, L):
    n, G = 1, g.copy()
    while _okey(G, L) not in Hkeys:
        G = G @ g
        n += 1
        if n > 12:
            raise ValueError("order of the representative is not found.")
    return n


# ==================================================
def _samb_vectors(mm, wyckoff_map, transform):
    """
    Combined SAMB in the label space: Z_j = sum cg X(t1,c1) (x) Y(W; t2,c2).

    Returns:
        - (dict) -- {z: {(UniqueSAMBType with W replaced, t1, c1, t2, c2): coefficient}} of the image transform(Z_j).
    """
    out = {}
    for z, (_, u, idx, comp) in mm["combined_id"].items():
        st = u.samb_type
        lc = mm["combined_samb"][st][idx][0][comp]
        out[z] = transform(u, lc)
    return out


# ==================================================
def _value(a):
    """complex -> float (or complex) rounded to 10 digits."""
    a = complex(a)
    if abs(a.imag) < 1e-10:
        return round(a.real, 10) + 0.0
    return complex(round(a.real, 10), round(a.imag, 10))


# ==================================================
def samb_map(mm, coef):
    """
    {z_j: [(z_k, A_jk)]} from the image vectors in the label space.

    Args:
        mm (MaterialModel): model.
        coef (dict): {z: {(UniqueSAMBType, t1, c1, t2, c2): coefficient}}, image of Z_j in the labels of the image
            system (the same UniqueSAMBType as the model for the SAMB of the image system with the same label).
    """
    basis = {}
    for z, (_, u, idx, comp) in mm["combined_id"].items():
        lc = mm["combined_samb"][u.samb_type][idx][0][comp]
        basis[z] = {(u, t1, c1, t2, c2): complex(cg) for cg, t1, c1, t2, c2 in lc}
    res = {}
    for zj, v in coef.items():
        terms = []
        for zk, b in basis.items():
            a = sum(np.conj(b[k]) * c for k, c in v.items() if k in b)
            if abs(a) > 1e-8:
                terms.append((zk, _value(a)))
        res[zj] = terms
    return res


# ==================================================
def _identity_coef(mm):
    out = {}
    for z, (_, u, idx, comp) in mm["combined_id"].items():
        lc = mm["combined_samb"][u.samb_type][idx][0][comp]
        out[z] = {(u, t1, c1, t2, c2): complex(cg) for cg, t1, c1, t2, c2 in lc}
    return out


# ==================================================
def _inversion_coef(mm, inv):
    """image of the combined SAMB under the inversion, in the labels of the image system."""
    out = {}
    for z, (_, u, idx, comp) in mm["combined_id"].items():
        st = u.samb_type
        kind = "bond" if "@" in st.wyckoff else "site"
        t = inv[kind][st.wyckoff]
        parity = (-1) ** (st.bk_info.bh_rank + st.bk_info.kt_rank)
        lc = mm["combined_samb"][st][idx][0][comp]
        v = {}
        for cg, t1, c1, t2, c2 in lc:
            if t["diag"]:
                images = [((t2, c2), complex(sp.S(t["eps"][t2][c2])))]
            else:
                blk = next(b for b in t["U"].values() if (t2, c2) in b[0])
                rows, rows2, U = blk
                i = rows.index((t2, c2))
                images = [(rows2[j], complex(sp.S(U[i, j]))) for j in range(len(rows2)) if complex(sp.S(U[i, j])) != 0]
            for (k2, cc2), e in images:
                key = (u, t1, c1, k2, cc2)
                v[key] = v.get(key, 0) + complex(cg) * parity * e
        out[z] = v
    return out


# ==================================================
# orbitals and matrix elements
_poly = {}


def _orbital_poly(name):
    if name not in _poly:
        expr = Group.global_info()["harmonics"]["basis_function"][name][0]
        _poly[name] = sp.lambdify((_x, _y, _z), sp.sympify(expr).subs({"x": _x, "y": _y, "z": _z}), "numpy")
    return _poly[name]


def _split_orbital(o):
    if o.startswith("("):
        name, s = o.strip("()").split(",")
        return name, s.lower()
    return o, None


def _spin_matrix(C):
    """SU(2) matrix of the proper part of the orthogonal matrix C (basis u, d)."""
    Cp = C * np.linalg.det(C)
    w = np.sqrt(max(0.0, 1 + np.trace(Cp))) / 2
    if w > 1e-6:
        v = np.array([Cp[2, 1] - Cp[1, 2], Cp[0, 2] - Cp[2, 0], Cp[1, 0] - Cp[0, 1]]) / (4 * w)
    else:
        k = np.argmax(np.diag(Cp))
        v = np.zeros(3)
        v[k] = np.sqrt(max(0.0, (1 + Cp[k, k]) / 2))
        for j in range(3):
            if j != k:
                v[j] = (Cp[k, j] + Cp[j, k]) / (4 * v[k])
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], complex)
    return w * np.eye(2) - 1j * (v[0] * sx + v[1] * sy + v[2] * sz)


_pts = np.random.default_rng(1).normal(size=(64, 3))


def orbital_matrix(C, orbs):
    """D[b, a]: the orbital a rotated by the Cartesian orthogonal matrix C is the sum of D[b, a] (orbital b), (g f)(r) = f(C^T r)."""
    names = [_split_orbital(o) for o in orbs]
    onames = list(dict.fromkeys(n for n, _ in names))
    F = np.array([[complex(_orbital_poly(n)(*p)) for n in onames] for p in _pts])
    G = np.array([[complex(_orbital_poly(n)(*(C.T @ p))) for n in onames] for p in _pts])
    Dorb = np.linalg.lstsq(F, G, rcond=None)[0]
    if not np.allclose(F @ Dorb, G, atol=1e-8):
        raise ValueError("the orbitals of a site are not closed under the operation.")
    U = _spin_matrix(C) if names[0][1] is not None else None
    si = {"u": 0, "d": 1}
    D = np.zeros((len(orbs), len(orbs)), complex)
    for a, (na, sa) in enumerate(names):
        for b, (nb, sb) in enumerate(names):
            D[b, a] = Dorb[onames.index(nb), onames.index(na)] * (1 if U is None else U[si[sb], si[sa]])
    return D


# ==================================================
def model_data(mm):
    """kets, site positions (conventional), and the SAMB matrices with conventional lattice vectors."""
    L = _lattice(mm.group)
    Pm = np.eye(3) if L is None else L
    kets = mm["full_matrix"]["ket"]
    pos = {}
    for name, val in mm["site"]["cell"].items():
        for v in val:
            if v.plus_set == 1:
                pos[(name, v.sublattice)] = np.array(v.position, dtype=float)
    site_of = [(k[0], k[1]) for k in kets]
    orb_of = [k[4] for k in kets]
    Z = {}
    for z, d in mm.get_combined_samb_matrix(fmt="value", bond=False).items():
        Z[z] = {
            (m, n, tuple(np.round(Pm @ np.array([n1, n2, n3], float), DIG) + 0.0)): complex(v)
            for (n1, n2, n3, m, n), v in d.items()
        }
    return {"site_of": site_of, "orb_of": orb_of, "pos": pos, "Z": Z}


# ==================================================
class _Frame:
    """operations of H and of the parent modulo the lattice of H, in the conventional coordinates of H."""

    def __init__(self, group, parent, R):
        self.group = group
        self.R = R
        self.L = _lattice(group)
        self.Hops = _ops(group)
        Ri = np.linalg.inv(R)
        Pg = Group(parent)
        Pops = [Ri @ W @ R for W in _ops(Pg)]
        trans = [np.zeros(3)]
        if self.L is not None:
            LP = _lattice(Pg)
            gens = [Ri[:3, :3] @ LP[:, i] for i in range(3)]
            T = {self.tkey(np.zeros(3)): np.zeros(3)}
            frontier = [np.zeros(3)]
            while frontier:
                new = []
                for t in frontier:
                    for g in gens:
                        u = _reduce(t + g, self.L)
                        k = self.tkey(u)
                        if k not in T:
                            T[k] = u
                            new.append(u)
                frontier = new
            trans = list(T.values())
        els = {}
        for W in Pops:
            for u in trans:
                Wt = W.copy()
                Wt[:3, 3] = _reduce(W[:3, 3] + u, self.L)
                els.setdefault(_okey(Wt, self.L), Wt)
        self.G = list(els.values())

    def tkey(self, v):
        return tuple(np.round(_reduce(v, self.L), DIG) + 0.0)


# ==================================================
def _find_wyckoff(frame, p, w):
    """parameter values, operation h, lattice vector n with p = h rep(params) + n, and (A, b) of the representative."""
    rep = sp.Matrix(frame.group.wyckoff["site"][w]["conventional"][0])
    syms = sorted(rep.free_symbols, key=str)
    A = np.array([[float(sp.diff(r, s)) for s in syms] for r in rep], dtype=float).reshape(3, len(syms))
    b = np.array([float(r.subs({s: 0 for s in syms})) for r in rep], dtype=float)
    rng = (0, -1, 1, -2, 2) if frame.L is not None else (0,)
    for h in frame.Hops:
        for n in itertools.product(rng, repeat=3):
            rhs = p - h[:3, :3] @ b - h[:3, 3] - np.array(n, float)
            M = h[:3, :3] @ A
            x = np.linalg.lstsq(M, rhs, rcond=None)[0] if M.shape[1] else np.zeros(0)
            if np.allclose(M @ x, rhs, atol=TOL):
                return dict(zip([str(s) for s in syms], x)), h, np.array(n, float), (A, b)
    raise ValueError(f"site {p} is not on the Wyckoff position {w}.")


# ==================================================
def _limit_point(lim, values):
    """parameter values on the limit nearest to the model values (free parameters by least squares)."""
    exprs = {k: sp.sympify(v) for k, v in lim["param"].items()}
    free = sorted({s for e in exprs.values() for s in e.free_symbols}, key=str)
    names = list(values)

    def vec(fv):
        sub = dict(zip(free, fv))
        return np.array([float(exprs[k].subs(sub)) if k in exprs else values[k] for k in names])

    fv = []
    if free:
        v0 = vec(np.zeros(len(free)))
        J = np.array([vec(np.eye(len(free))[i]) - v0 for i in range(len(free))]).T
        fv = np.linalg.lstsq(J, np.array([values[k] for k in names]) - v0, rcond=None)[0]
    return dict(zip(names, vec(fv)))


# ==================================================
def parent_structure(mm, md, frame, limits):
    """
    Positions of every site in the parent structure (the nearest limit of each site).

    Returns:
        - (dict) -- {site: {"wyckoff", "param", "limit": index, "limit_param", "parent_wyckoff"}}.
        - (dict) -- {(site, sublattice): position}.
    """
    A_cart = np.array(mm["unit_vector"], dtype=float).T
    info, lim_pos = {}, {}
    for name, rep in mm["site"]["representative"].items():
        p = np.array(rep.position, dtype=float)
        vals, h, n, (A, b) = _find_wyckoff(frame, p, rep.wyckoff)
        cands = limits.get(rep.wyckoff, [])
        if vals:
            if not cands:
                raise ValueError(f"site '{name}' ({rep.wyckoff}): no parent structure.")
            pts = []
            for lim in cands:
                lv = _limit_point(lim, vals)
                q = h[:3, :3] @ (A @ np.array([lv[k] for k in vals]) + b) + h[:3, 3] + n
                d = min(
                    np.linalg.norm(A_cart @ (W[:3, :3] @ q + W[:3, 3] + np.array(m) - p))
                    for W in frame.Hops
                    for m in (itertools.product(range(-1, 2), repeat=3) if frame.L is not None else [(0, 0, 0)])
                )
                pts.append((d, lim, lv, q))
            k = int(np.argmin([t[0] for t in pts]))
            _, lim, lv, q0 = pts[k]
            info[name] = {
                "wyckoff": rep.wyckoff,
                "param": {s: float(v) for s, v in vals.items()},
                "limit": k,
                "limit_param": dict(lim["param"]),
                "parent_wyckoff": lim["parent_wyckoff"],
            }
        else:
            q0 = p
            info[name] = {"wyckoff": rep.wyckoff, "param": {}, "limit": None, "limit_param": {}, "parent_wyckoff": None}
        for (nm, sl), ps in md["pos"].items():
            if nm != name:
                continue
            for W in frame.Hops:
                d = ps - (W[:3, :3] @ p + W[:3, 3])
                if np.allclose(d, np.rint(d) if frame.L is not None else 0, atol=TOL):
                    q = W[:3, :3] @ q0 + W[:3, 3]
                    if frame.L is not None:
                        q = q + np.rint(ps - q)
                    lim_pos[(nm, sl)] = q
                    break
            else:
                raise ValueError(f"site {nm}({sl}) is not generated by the group.")
    keys = {frame.tkey(v) for v in lim_pos.values()}
    if len(keys) != len(lim_pos):
        raise ValueError("sites coincide in the parent structure.")
    for g in frame.G:
        for k, v in lim_pos.items():
            if frame.tkey(g[:3, :3] @ v + g[:3, 3]) not in keys:
                raise ValueError(f"parent structure is not invariant under the parent: site {k}.")
    return info, lim_pos


# ==================================================
class _KeySpace:
    """matrix elements on the parent structure: (bra position, ket position, bra orbital, ket orbital), modulo L_H."""

    def __init__(self, frame, C):
        self.frame = frame
        self.C = C
        self.Ci = np.linalg.inv(C)
        self.keys, self.idx = [], {}

    def set_sites(self, positions, orbitals):
        """sites of the structure (reduced) and their orbitals."""
        self.sites = [_reduce(p, self.frame.L) for p in positions]
        self.site_orbs = list(orbitals)

    def orbitals(self, r):
        """orbitals of the site at r (modulo the lattice)."""
        L = self.frame.L
        red = _reduce(np.asarray(r, dtype=float), L)
        for p, o in zip(self.sites, self.site_orbs):
            d = red - p
            dd = d if L is None else L @ (np.linalg.solve(L, d) - np.rint(np.linalg.solve(L, d)))
            if np.allclose(dd, 0, atol=1e-5):
                return o
        raise ValueError(f"position {r} is not a site of the structure.")

    def snap(self, r):
        L = self.frame.L
        r = np.asarray(r, dtype=float)
        red = _reduce(r, L)
        for p in self.sites:
            d = red - p
            dd = d if L is None else L @ (np.linalg.solve(L, d) - np.rint(np.linalg.solve(L, d)))
            if np.allclose(dd, 0, atol=1e-5):
                t = np.zeros(3) if L is None else L @ np.rint(np.linalg.solve(L, r - p))
                return p, t
        raise ValueError(f"position {r} is not a site of the parent structure.")

    def canon(self, rb, rk):
        pb, tb = self.snap(rb)
        pk, tk = self.snap(rk)
        return tuple(np.round(pb, DIG) + 0.0), tuple(np.round(pk + tk - tb, DIG) + 0.0)

    def add(self, k):
        if k not in self.idx:
            self.idx[k] = len(self.keys)
            self.keys.append(k)
        return self.idx[k]

    def vector(self, d, n):
        v = np.zeros(n, complex)
        for k, val in d.items():
            v[self.idx[k]] += val
        return v

    def close(self, G):
        pairs = {(k[0], k[1]) for k in self.keys}
        frontier = list(pairs)
        while frontier:
            new = []
            for rb, rk in frontier:
                for g in G:
                    q = self.canon(g[:3, :3] @ np.array(rb) + g[:3, 3], g[:3, :3] @ np.array(rk) + g[:3, 3])
                    if q not in pairs:
                        pairs.add(q)
                        new.append(q)
            frontier = new
        for rb, rk in sorted(pairs):
            for a in self.orbitals(rb):
                for b in self.orbitals(rk):
                    self.add((rb, rk, a, b))

    def action(self, g):
        """matrix of the action of g on the key space."""
        Cg = self.C @ g[:3, :3] @ self.Ci
        Dcache = {}
        n = len(self.keys)
        M = np.zeros((n, n), complex)
        for i, (rb, rk, a, b) in enumerate(self.keys):
            ob = self.orbitals(rb)
            ok = self.orbitals(rk)
            for lst in (ob, ok):
                if tuple(lst) not in Dcache:
                    Dcache[tuple(lst)] = orbital_matrix(Cg, lst)
            Db, Dk = Dcache[tuple(ob)], Dcache[tuple(ok)]
            nb, nk = self.canon(g[:3, :3] @ np.array(rb) + g[:3, 3], g[:3, :3] @ np.array(rk) + g[:3, 3])
            ia, ib = ob.index(a), ok.index(b)
            for a2, x in enumerate(Db[:, ia]):
                if abs(x) < 1e-12:
                    continue
                for b2, y in enumerate(Dk[:, ib]):
                    if abs(y) < 1e-12:
                        continue
                    M[self.idx[(nb, nk, ob[a2], ok[b2])], i] += x * np.conj(y)
        return M


# ==================================================
def weight(mm, frame, limits, g, n):
    """
    Weights of the eigenvalues of g on the parent structure.

    Returns:
        - (dict) -- parent structure of each site (parent_structure).
        - (dict) -- {z: [w_m for m = 0..n-1]}.
    """
    md = model_data(mm)
    info, lim_pos = parent_structure(mm, md, frame, limits)
    C = np.array(mm["unit_vector"], dtype=float).T
    ks = _KeySpace(frame, C)
    ks.set_sites(
        list(lim_pos.values()), [[o for o, s in zip(md["orb_of"], md["site_of"]) if s == name_sl] for name_sl in lim_pos]
    )
    Zp = {}
    for z, d in md["Z"].items():
        e = {}
        for (m, k, Rc), v in d.items():
            cb, ck = ks.canon(lim_pos[md["site_of"][m]], lim_pos[md["site_of"][k]] + np.array(Rc))
            key = (cb, ck, md["orb_of"][m], md["orb_of"][k])
            e[key] = e.get(key, 0) + v
            ks.add(key)
        Zp[z] = e
    ks.close([g])
    G = ks.action(g)
    N = len(ks.keys)
    zs = list(md["Z"])
    X = np.array([ks.vector(Zp[z], N) for z in zs]).T
    T = [X]
    for _ in range(1, n):
        T.append(G @ T[-1])
    nx2 = np.einsum("ij,ij->j", X.conj(), X).real
    ov = [np.einsum("ij,ij->j", X.conj(), Tt) / nx2 for Tt in T]
    w = np.exp(2j * np.pi / n)
    res = {
        z: [round(float(np.real(sum(w ** (-m * t) * ov[t][j] for t in range(n)) / n)), 10) + 0.0 for m in range(n)]
        for j, z in enumerate(zs)
    }
    return info, res


# ==================================================
def _site_bond_image(mm, inv):
    """images of the sites and bonds of the model (Wyckoff position, permutation, bond reversal)."""
    site = {}
    for name, rep in mm["site"]["representative"].items():
        t = inv["site"][rep.wyckoff]
        site[name] = {"wyckoff": rep.wyckoff, "image": t["image"][1], "perm": t["perm"]}
    bond = {}
    for name, rep in mm["bond"]["representative"].items():
        t = inv["bond"][rep.wyckoff]
        bond[name] = {"wyckoff": rep.wyckoff, "image": t["image"][1], "perm": t["perm"], "sign": t["sign"]}
    return site, bond


# ==================================================
def coset_analysis(mm, spec=None, with_weight=True):
    """
    Coset analysis of a model.

    Args:
        mm (MaterialModel): model (analyzed).
        spec (None, False, dict or list, optional): "coset" of the model input.
        with_weight (bool, optional): calculate the weights on the parent structure (direct parent only) ?

    Returns:
        - (list) -- result for each representative, dict with "parent", "rep", "rep_matrix", "order", "kind", "frame",
          "image_group", "site", "bond", "samb", "parent_structure", "weight", "weight_note" (the reason when the
          weights are not given, e.g. no parent structure of a site).
    """
    group = mm.group
    reps = parse_coset(group, spec)
    if not reps:
        return []
    first = group.coset
    L = _lattice(group)
    Hops = _ops(group)
    Hkeys = {_okey(W, L) for W in Hops}
    chiral = all(np.linalg.det(W[:3, :3]) > 0 for W in Hops)
    out = []
    for parent, tag in reps:
        R = _chain_R(group, parent)
        g = _op_matrix(group, parent, R, tag)
        n = _order(g, Hkeys, L)
        res = {"parent": parent, "rep": tag, "rep_matrix": g.tolist(), "order": n}
        if np.linalg.det(g[:3, :3]) > 0:
            res.update(kind="rotation", frame=g.tolist(), image_group=group_id(group), site=None, bond=None)
            res["samb"] = samb_map(mm, _identity_coef(mm))
        elif not chiral:
            m = next(W for W in Hops if np.linalg.det(W[:3, :3]) < 0)
            g2 = g @ m
            frame = _Frame(group, group_id(group), np.eye(4))
            md = model_data(mm)
            res.update(kind="improper", frame=g2.tolist(), image_group=group_id(group), site=None, bond=None)
            res["samb"] = _model_action(mm, md, frame, np.linalg.inv(m))
        else:
            inv = first["inversion"]
            I0 = _mat4(_float(first["rep_matrix"]))
            r = g @ np.linalg.inv(I0)
            site, bond = _site_bond_image(mm, inv)
            res.update(kind="inversion", frame=r.tolist(), image_group=inv["image"], map=_float(inv["map"]).tolist())
            res.update(site=site, bond=bond)
            res["samb"] = samb_map(mm, _inversion_coef(mm, inv))
        res["parent_structure"], res["weight"], res["weight_note"] = None, None, None
        if with_weight and parent == first["parent"]:
            frame = _Frame(group, parent, R)
            try:
                res["parent_structure"], res["weight"] = weight(mm, frame, first["parent_limit"], g, n)
            except ValueError as e:  # no parent structure for the sites of the model.
                res["weight_note"] = str(e)
        out.append(res)
    return out


# ==================================================
def _model_action(mm, md, frame, m):
    """{z_j: [(z_k, A_jk)]} of the operation m of H acting on the SAMB of the model."""
    C = np.array(mm["unit_vector"], dtype=float).T
    ks = _KeySpace(frame, C)
    ks.set_sites(
        list(md["pos"].values()), [[o for o, s in zip(md["orb_of"], md["site_of"]) if s == name_sl] for name_sl in md["pos"]]
    )
    Zp = {}
    for z, d in md["Z"].items():
        e = {}
        for (a, b, Rc), v in d.items():
            key = ks.canon(md["pos"][md["site_of"][a]], md["pos"][md["site_of"][b]] + np.array(Rc)) + (
                md["orb_of"][a],
                md["orb_of"][b],
            )
            e[key] = e.get(key, 0) + v
            ks.add(key)
        Zp[z] = e
    ks.close([m])
    G = ks.action(m)
    N = len(ks.keys)
    zs = list(md["Z"])
    X = np.array([ks.vector(Zp[z], N) for z in zs]).T
    A = X.conj().T @ (G @ X)
    res = {}
    for j, zj in enumerate(zs):
        res[zj] = [(zk, _value(A[k, j])) for k, zk in enumerate(zs) if abs(A[k, j]) > 1e-8]
    return res
