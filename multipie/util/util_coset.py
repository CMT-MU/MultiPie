# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Coset of a model: the SAMB of a model of H moved by a representative g of a coset g H in a supergroup P.

The data of the cosets are those of the group binaries (MultiPieGen): "coset" of the group file (the first coset in
the direct parent) and "parent", "coset" of info (the parents and the cosets in them); any other supergroup is reached
through the direct supergroups ("supergroup" of info, util_supergroup).

- Reduction of g (in the coordinates of H).

  - "rotation" (det g > 0): the image system is written in the frame moved by g; its data are those of the model, and
    g Z_j = +Z'_j (the SAMB of the image system with the same label).
  - "improper" (det g < 0, H not chiral): g = g' m^-1 with an improper operation m of H and the rotation g' = g m:
    g Z_j = sum_k D(m^-1)_kj Z'_k in the frame moved by g' (D = 1 for the identity irrep).
  - "inversion" (det g < 0, H chiral): g = r I with the inversion I of the first coset in the direct parent and a
    rotation r: I maps H to H' = I H I^-1 (map M), the sites and bonds of each Wyckoff position W to those of W'
    (row i -> row perm[i], bond reversed if sign = -1), and the cluster SAMB by eps (label by label) or by U (within an
    irrep). The atomic SAMB of the orbitals of ranks (l1, l2) are multiplied by (-1)^(l1+l2). The frame is moved by r.

- Image system of the inversion. Its representative (row 0) of each Wyckoff position is the image of row i0 with
  perm[i0] = 0 (image_representative); a model of H' created from its own input may have other representative bonds,
  and the inversion model (model input "create_inversion_model", inversion_model) has these representatives.

- Parent structure (resolve_coset). The limit of the Wyckoff parameters at which the structure is invariant under P:
  "parent_limit" of the binaries for the direct parent, util_supergroup.parent_limit otherwise; for a supergroup that
  is not a parent of the binaries, the embedding whose parent structure is the nearest to the model. On it, g is a
  symmetry and acts on the matrix elements of the model (image_space, the atoms matched by their parent positions).

- Weight. The weight of the eigenvalue w^m (w = exp(2 pi i/n), n the order of g modulo H) of Z_j is
  w_m(j) = (1/n) sum_t w^(-m t) <Z_j, g^t Z_j>/<Z_j, Z_j> (weight_from_space).

- S and F. S^(t)_kj = <Z_k, g^t Z_j> (overlap_matrix) and n Pi^(m) = sum_t w^(-m t) S^(t) (projected_matrix) give
  <g^t Z_j> and n <Pi_m Z_j> for a state invariant under H if every g^t Z_j lies on the bonds of the model
  ("complete", parent_shells); F_j = sum_{m>0} n <Pi_m Z_j> (MaterialModel.parent_forbidden_part).
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
        - (list) -- [(parent id, representative tag or None, embedding or None)], [] for False, or for None and a
          group without a direct parent.

    Note:
        - None: the representative of the first coset in the direct parent.
        - {"parent": id, "rep": tag}: for a parent of the binaries (info["parent"][H]), the operation must be one of
          the cosets of info["coset"][H][parent] (not a translation only, nor an operation of H; default: the first).
        - any supergroup P (Group.supergroups), or with "embedding" (index of Group.embedding): "rep" is an element
          of P in the coordinates of H (default: the first coset other than H modulo the lattice of H, the inversion
          if any); the embedding is the one whose parent structure is the nearest to the model if not given.
    """
    if spec is False:
        return []
    first = group.coset
    if spec is None:
        return [] if first is None else [(first["parent"], first["rep"], None)]
    specs = [spec] if isinstance(spec, dict) else list(spec)
    info = Group.global_info()
    gid = group_id(group)
    out = []
    for s in specs:
        if not isinstance(s, dict) or set(s) - {"parent", "rep", "embedding"} or ("rep" not in s and "parent" not in s):
            raise ValueError(
                f"coset must be None, False, {{'parent': id, 'rep': tag, 'embedding': index}} or a list of them, but {s!r} is given."
            )
        if "parent" not in s and first is None:
            raise ValueError(f"{gid} has no direct parent: give the parent in {s!r}.")
        parent = s.get("parent", first["parent"] if first else None)
        parent = group_id(Group(parent))
        cosets = info["coset"].get(gid, {}).get(parent)
        emb = s.get("embedding")
        if cosets is None or emb is not None:  # any supergroup (util_supergroup).
            from multipie.util.util_supergroup import embedding

            if not embedding(gid, parent):
                raise ValueError(f"{parent} is not a supergroup of {gid}.")
            out.append((parent, s.get("rep"), emb))
            continue
        rep = s.get("rep", cosets[0][0])
        if not any(rep in c for c in cosets):
            raise ValueError(f"{rep} is not an operation of the cosets of {gid} in {parent}.")
        out.append((parent, rep, None))
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

    Note:
        - g Z_j = sum_k A_kj Z'_k is solved in the label space by least squares: the label vectors of the SAMB are
          not normalized in general (e.g., the norm is 1/2 for a hybrid of the orbitals of different ranks, whose
          hermitian partner has no label of its own).
    """
    zs = list(mm["combined_id"])
    keys = {}
    basis = {}
    for z, (_, u, idx, comp) in mm["combined_id"].items():
        lc = mm["combined_samb"][u.samb_type][idx][0][comp]
        basis[z] = {(u, t1, c1, t2, c2): complex(cg) for cg, t1, c1, t2, c2 in lc}
        for k in basis[z]:
            keys.setdefault(k, len(keys))
    B = np.zeros((len(keys), len(zs)), dtype=complex)
    for j, z in enumerate(zs):
        for k, c in basis[z].items():
            B[keys[k], j] = c
    V = np.zeros((len(keys), len(zs)), dtype=complex)
    for j, z in enumerate(zs):
        for k, c in coef[z].items():
            if k not in keys:
                raise ValueError(f"the image of {z} is not in the span of the SAMB of the model (label {k}).")
            V[keys[k], j] = c
    A = np.linalg.lstsq(B, V, rcond=None)[0]
    if not np.allclose(B @ A, V, atol=1e-8):
        raise ValueError("the image of the SAMB is not in the span of the SAMB of the model.")
    res = {}
    for j, zj in enumerate(zs):
        res[zj] = [(zk, _value(A[k, j])) for k, zk in enumerate(zs) if abs(A[k, j]) > 1e-8]
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
def _round(v):
    """vector rounded to DIG digits (list of float)."""
    return [round(float(x), DIG) + 0.0 for x in v]


# ==================================================
def _mod1(v):
    """vector reduced to [0,1) (components within 1e-10 of an integer are snapped to it)."""
    v = np.asarray(v, dtype=float)
    v = np.where(np.abs(v - np.rint(v)) < 1e-10, np.rint(v), v)
    return np.mod(v, 1) + 0.0


# ==================================================
def image_representative(inv, site_dict, bond_dict, lattice=True):
    """
    Representative sites and bonds of the image system under the inversion (first coset of a chiral group).

    The eps/U of the cluster SAMB (MultiPieGen) map row i of each Wyckoff position W to row perm[i] of W' (bond
    reversed if sign = -1). The representative (row 0) of the image is therefore M (row i0) with perm[i0] = 0, which
    is M (representative) when perm[0] = 0.

    Args:
        inv (dict): "inversion" of the first coset (Group.coset).
        site_dict (dict): site dict of the model ("representative", "cell").
        bond_dict (dict): bond dict of the model ("representative", "cell").
        lattice (bool, optional): space group (positions reduced to [0,1)) ?

    Returns:
        - (dict) -- {site: position} of the image (fractional, in [0,1), ndarray).
        - (dict) -- {bond: (vector, center)} of the image (fractional, center in [0,1), ndarray).
    """
    M = _mat4(_float(inv["map"]))
    W, w = M[:3, :3], M[:3, 3]
    red = _mod1 if lattice else (lambda v: np.asarray(v, dtype=float) + 0.0)
    site = {}
    for name, rep in site_dict["representative"].items():
        i0 = inv["site"][rep.wyckoff]["perm"].index(0)
        p = np.asarray(site_dict["cell"][name][i0].position, dtype=float)
        site[name] = red(W @ p + w)
    bond = {}
    for name, rep in bond_dict["representative"].items():
        t = inv["bond"][rep.wyckoff]
        i0 = t["perm"].index(0)
        b = bond_dict["cell"][name][i0]
        v = t["sign"][i0] * (W @ np.asarray(b.vector, dtype=float))
        c = W @ np.asarray(b.center, dtype=float) + w
        bond[name] = (v, red(c))
    return site, bond


# ==================================================
def _num_str(x):
    """number as a short string (fraction with a denominator up to 24, else the decimal)."""
    f = sp.Rational(x).limit_denominator(24)
    if abs(float(f) - x) < 1e-10:
        return str(f)
    return repr(round(float(x), 12))


# ==================================================
def inversion_model(group, site_data, site_dict, bond_dict):
    """
    Input of the inversion model: the image of a model of a chiral group under the first coset (inversion).

    Args:
        group (Group): group of the model (its first coset is the inversion).
        site_data (dict): "site" of the model input.
        site_dict (dict): site dict of the model.
        bond_dict (dict): bond dict of the model.

    Returns:
        - (str) -- image group.
        - (dict) -- "site" of the image, {name: (position, orbital)}.
        - (dict) -- representative bonds of the image, {(tail, head): [bonds of the n-th neighbor]}
          (rep_bond_given of parse_representative_bond).
        - (dict) -- information ("group", "image_group", "rep", "map", "site", "bond") of the source and the image.
    """
    f = group.coset
    inv = f["inversion"]
    rs, rb = image_representative(inv, site_dict, bond_dict, not group.is_point_group)
    site = {
        name: ("[" + ",".join(_num_str(x) for x in rs[name]) + "]", site_data[name][1]) for name in site_dict["representative"]
    }

    rep_bond = {}
    for name, rep in bond_dict["representative"].items():
        lst = rep_bond.setdefault((rep.tail, rep.head), [])
        while len(lst) < rep.neighbor:
            lst.append([])
        v, c = rb[name]
        lst[rep.neighbor - 1].append(np.concatenate([v, c]))
    for info in bond_dict["info"]:  # the neighbors not selected.
        lst = rep_bond.setdefault((info.tail, info.head), [])
        while len(lst) < max(info.neighbor):
            lst.append([])

    info = {
        "group": group_id(group),
        "image_group": inv["image"],
        "rep": f["rep"],
        "map": _mat4(_float(inv["map"])).tolist(),
        "site": {
            name: {"source": _round(rep.position), "image": _round(rs[name])} for name, rep in site_dict["representative"].items()
        },
        "bond": {
            name: {"source": (_round(rep.vector), _round(rep.center)), "image": (_round(rb[name][0]), _round(rb[name][1]))}
            for name, rep in bond_dict["representative"].items()
        },
    }
    return inv["image"], site, rep_bond, info


# ==================================================
def check_inversion_model(info, site_dict, bond_dict, lattice=True):
    """
    Check that the representative sites and bonds of the inversion model are the images (info of inversion_model).
    """

    def same(a, b, lattice):
        d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
        return np.allclose(d - np.rint(d), 0, atol=TOL) if lattice else np.allclose(d, 0, atol=TOL)

    for name, rep in site_dict["representative"].items():
        if not same(rep.position, info["site"][name]["image"], lattice):
            raise RuntimeError(f"inversion model: the representative of site '{name}' is not the image.")
    for name, rep in bond_dict["representative"].items():
        v, c = info["bond"][name]["image"]
        if not (same(rep.vector, v, False) and same(rep.center, c, lattice)):
            raise RuntimeError(f"inversion model: the representative of bond '{name}' is not the image.")
    if set(bond_dict["representative"]) != set(info["bond"]):
        raise RuntimeError("inversion model: the bonds differ from those of the model.")


# ==================================================
def _site_bond_image(mm, inv):
    """images of the sites and bonds of the model (Wyckoff position, permutation, bond reversal, representative)."""
    rs, rb = image_representative(inv, mm["site"], mm["bond"], not mm.group.is_point_group)
    site = {}
    for name, rep in mm["site"]["representative"].items():
        t = inv["site"][rep.wyckoff]
        site[name] = {"wyckoff": rep.wyckoff, "image": t["image"][1], "perm": t["perm"], "position": _round(rs[name])}
    bond = {}
    for name, rep in mm["bond"]["representative"].items():
        t = inv["bond"][rep.wyckoff]
        bond[name] = {"wyckoff": rep.wyckoff, "image": t["image"][1], "perm": t["perm"], "sign": t["sign"]}
        bond[name]["vector"], bond[name]["center"] = _round(rb[name][0]), _round(rb[name][1])
    return site, bond


# ==================================================
def coset_analysis(mm, spec=None, with_weight=True):
    """
    Coset analysis of a model.

    Args:
        mm (MaterialModel): model (analyzed).
        spec (None, False, dict or list, optional): "coset" of the model input.
        with_weight (bool, optional): calculate the weights, S and the shells on the parent structure ?

    Returns:
        - (list) -- result for each representative, dict with "parent", "rep", "rep_matrix", "order", "kind", "frame",
          "image_group", "site", "bond", "samb", "embedding" (None for the parents of the binaries), "parent_structure",
          "weight", "weight_note" (the reason when the weights are not given, e.g. no parent structure of a site), "S"
          (S[k][j] = <Z_k, g Z_j>), "S_eigenvalue", "complete" (the model contains the whole neighbor shells of the
          parent structure), "shell" (parent_shells).
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
    for parent, tag, emb in reps:
        rc = None
        if emb is None and tag is not None and parent in Group.global_info()["parent"].get(group_id(group), []):
            R = _chain_R(group, parent)
            g = _op_matrix(group, parent, R, tag)
        else:
            rc = resolve_coset(mm, parent, tag, emb, require_structure=False)
            g, tag = rc["g"], rc["rep"]
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
        res["embedding"] = None
        if rc is not None and rc["embedding"] is not None:
            e = rc["embedding"]
            res["embedding"] = {"A": np.asarray(e["A"]).tolist(), "origin": _round(e["origin"]), "path": e["path"]}
        res["parent_structure"], res["weight"], res["weight_note"] = None, None, None
        res["S"], res["S_eigenvalue"], res["complete"], res["shell"] = None, None, None, None
        if with_weight:
            try:
                if rc is None:
                    rc = resolve_coset(mm, parent, tag)
                if rc["lim_pos"] is None:
                    raise ValueError(rc["note"])
                sp_ = image_space(mm, rc)
                res["parent_structure"] = rc["parent_structure"]
                res["weight"] = weight_from_space(sp_, n)
                S, _ = overlap_matrix(sp_)
                _, complete = projected_matrix(sp_, n)
                res["S"] = np.asarray(S).tolist()
                ev = np.linalg.eigvalsh((S + np.conj(S.T)) / 2) if n == 2 else np.linalg.eigvals(S)
                res["S_eigenvalue"] = [_value(v) for v in ev]
                res["complete"] = complete
                res["shell"] = parent_shells(mm, rc)
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


# ==================================================
# any supergroup (util_supergroup): the image g Z_j on the parent structure, S, F, shells.
def _model_wyckoff(mm):
    return sorted({rep.wyckoff for rep in mm["site"]["representative"].values()})


# ==================================================
def _displacement(mm, md, lim_pos):
    """largest distance between the sites of the model and those of the parent structure (Cartesian)."""
    C = np.array(mm["unit_vector"], dtype=float).T
    return max(np.linalg.norm(C @ (lim_pos[k] - p)) for k, p in md["pos"].items())


# ==================================================
def resolve_coset(mm, parent=None, rep=None, embedding=None, require_structure=True):
    """
    The parent, the embedding, the representative g, and the parent structure used by image_space.

    Args:
        mm (MaterialModel): model (analyzed).
        parent (str, optional): supergroup P (default: the direct parent of the binaries).
        rep (str or array, optional): element g of P ("point operation of P:[translation in H]" or 4x4 in the
            coordinates of H; default: the first coset of the binaries for its parents, the first coset other than
            H modulo the lattice of H (the inversion if any, a translation only if no other) otherwise).
        embedding (None, int, dict or array, optional): embedding of H in P (util_supergroup.resolve_embedding);
            default: the parents of the binaries with their transformation, otherwise the embedding whose parent
            structure is the nearest to the model.
        require_structure (bool, optional): raise an error if no embedding gives a parent structure of the model
            (otherwise the first embedding, with "parent_structure" and "lim_pos" None and the reason in "note").

    Returns:
        - (dict) -- "parent", "rep", "g", "order", "frame" (_Frame), "limits", "embedding" (None for the binaries),
          "parent_structure", "lim_pos", "md".
    """
    from multipie.util import util_supergroup as us

    group = mm.group
    gid = group_id(group)
    first = group.coset
    info = Group.global_info()
    if parent is None:
        if first is None:
            raise ValueError(f"{gid} has no parent: give the parent.")
        parent = first["parent"]
    parent = group_id(Group(parent))
    L = _lattice(group)
    Hops = _ops(group)
    Hkeys = {_okey(W, L) for W in Hops}
    md = model_data(mm)
    if embedding is None and parent in info["parent"].get(gid, []):
        R = _chain_R(group, parent)
        if rep is None:
            rep = first["rep"] if parent == first["parent"] else info["coset"][gid][parent][0][0]
        g = (
            _op_matrix(group, parent, R, rep)
            if isinstance(rep, str)
            else us._mat4(np.array(rep, float)[:3, :3], np.array(rep, float)[:3, 3])
        )
        limits = first["parent_limit"] if parent == first["parent"] else us.parent_limit(gid, parent, _model_wyckoff(mm), R)[0]
        frame = _Frame(group, parent, R)
        pstruct, lim_pos = parent_structure(mm, md, frame, limits)
        emb = None
        note = None
    else:
        cands = us.embedding(gid, parent) if embedding is None else [us.resolve_embedding(gid, parent, embedding)]
        if not cands:
            raise ValueError(f"{parent} is not a supergroup of {gid}.")
        best, err = None, None
        for e in cands:
            limits = us.parent_limit(gid, parent, _model_wyckoff(mm), e)[0]
            frame = _Frame(group, parent, e["A"])
            try:
                pstruct, lim_pos = parent_structure(mm, md, frame, limits)
            except ValueError as ex:
                err = ex
                continue
            d = _displacement(mm, md, lim_pos)
            if best is None or d < best[0] - 1e-8:
                best = (d, e, limits, frame, pstruct, lim_pos)
        if best is None:
            if require_structure:
                raise ValueError(f"no parent structure of the model for {parent}: {err}")
            e = cands[0]
            best = (None, e, None, _Frame(group, parent, e["A"]), None, None)
            note = f"no parent structure of the model for {parent}: {err}"
        else:
            note = None
        _, emb, limits, frame, pstruct, lim_pos = best
        R = emb["A"]
        if rep is None:
            cs, _ = us.coset_decomposition(gid, parent, emb, modulo="self")
            nonid = [c for c in cs if c["kind"] not in ("identity", "translation")]
            inv = [c for c in nonid if c["kind"] == "inversion"]
            trans = [c for c in cs if c["kind"] == "translation"]  # a k-supergroup of the same point group.
            rep = (inv or nonid or trans)[0]["rep"]
        g = us.op_matrix(gid, parent, rep, emb)
    n = _order(g, Hkeys, L)
    rep_tag = rep if isinstance(rep, str) else us.vec_str(np.array(rep).ravel())
    return {
        "parent": parent,
        "rep": rep_tag,
        "g": g,
        "order": n,
        "frame": frame,
        "limits": limits,
        "embedding": emb,
        "parent_structure": pstruct,
        "lim_pos": lim_pos,
        "md": md,
        "note": note,
    }


# ==================================================
def image_space(mm, rc):
    """
    The SAMB of the model and their images g Z_j on the parent structure (keys of the matrix elements).

    Args:
        mm (MaterialModel): model.
        rc (dict): resolve_coset.

    Returns:
        - (dict) -- "ks" (_KeySpace), "G" (action of g), "X" (columns: Z_j), "zs", "n0" (number of keys of the
          model; the keys beyond are reached by g only).
    """
    md, lim_pos, frame, g = rc["md"], rc["lim_pos"], rc["frame"], rc["g"]
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
    n0 = len(ks.keys)
    ks.close([g])
    G = ks.action(g)
    N = len(ks.keys)
    zs = list(md["Z"])
    X = np.array([ks.vector(Zp[z], N) for z in zs]).T
    return {"ks": ks, "G": G, "X": X, "zs": zs, "n0": n0}


# ==================================================
def overlap_matrix(sp_, t=1):
    """
    S^(t)_kj = <Z_k, g^t Z_j> and whether every g^t Z_j lies on the bonds (pairs of sites) of the model.

    Args:
        sp_ (dict): image_space.
        t (int, optional): power of g.

    Returns:
        - (ndarray) -- S^(t) (real if the imaginary part vanishes).
        - (bool) -- complete.
    """
    X, G = sp_["X"], sp_["G"]
    GX = X
    for _ in range(t):
        GX = G @ GX
    nrm = np.sqrt(np.einsum("ij,ij->j", X.conj(), X).real)
    S = (X.conj().T @ GX) / np.outer(nrm, nrm)
    if np.allclose(S.imag, 0, atol=1e-10):
        S = S.real
    S = np.where(np.abs(S) < 1e-12, 0, S) + 0
    # every g^t Z_j on the pairs of sites (bonds) of the model.
    keys = sp_["ks"].keys
    pairs = {(k[0], k[1]) for k in keys[: sp_["n0"]]}
    out = [i for i, k in enumerate(keys) if (k[0], k[1]) not in pairs]
    complete = bool(np.allclose(GX[out], 0, atol=1e-8)) if out else True
    return S, complete


# ==================================================
def projected_matrix(sp_, n):
    """
    Matrices of the projections on the eigenvalues w^m of g (order n, w = exp(2 pi i/n)) in the SAMB of the model.

    Returns:
        - (list) -- [n Pi^(m)_kj = sum_t w^(-m t) S^(t)_kj for m = 0..n-1], so that n <Pi_m Z_j> = sum_k (n Pi^(m))_kj <Z_k>
          for a state invariant under H.
        - (bool) -- complete (for every power of g).
    """
    Ss, complete = [], True
    for t in range(n):
        S, c = overlap_matrix(sp_, t)
        Ss.append(S)
        complete = complete and c
    w = np.exp(2j * np.pi / n)
    out = []
    for m in range(n):
        M = sum(w ** (-m * t) * Ss[t] for t in range(n))
        if np.allclose(np.imag(M), 0, atol=1e-10):
            M = np.real(M)
        out.append(np.where(np.abs(M) < 1e-12, 0, M) + 0)
    return out, complete


# ==================================================
def weight_from_space(sp_, n):
    """weights of the eigenvalues w^m of g (order n), {z: [w_m]}."""
    X, G = sp_["X"], sp_["G"]
    T = [X]
    for _ in range(1, n):
        T.append(G @ T[-1])
    nx2 = np.einsum("ij,ij->j", X.conj(), X).real
    ov = [np.einsum("ij,ij->j", X.conj(), Tt) / nx2 for Tt in T]
    w = np.exp(2j * np.pi / n)
    return {
        z: [round(float(np.real(sum(w ** (-m * t) * ov[t][j] for t in range(n)) / n)), 10) + 0.0 for m in range(n)]
        for j, z in enumerate(sp_["zs"])
    }


# ==================================================
def transformed_matrix(mm, rc, sp_, z=None):
    """
    g Z_j in the basis of the model (the atoms matched on the parent structure).

    Returns:
        - (dict) -- {z: {(n1, n2, n3, a, b): value}} (as get_combined_samb_matrix(fmt="value")), a, b: kets, (n1, n2, n3):
          primitive lattice vector of the ket.
    """
    md, lim_pos, ks = rc["md"], rc["lim_pos"], sp_["ks"]
    L = ks.frame.L
    Pm = np.eye(3) if L is None else L
    red = {k: ks.snap(p)[0] for k, p in lim_pos.items()}
    ket_index = {(s, o): i for i, (s, o) in enumerate(zip(md["site_of"], md["orb_of"]))}

    def site_at(r):
        for k, p in red.items():
            if np.allclose(p, r, atol=1e-5):
                return k
        raise ValueError(f"position {r} is not a site of the parent structure.")

    zs = sp_["zs"] if z is None else ([z] if isinstance(z, str) else list(z))
    GX = sp_["G"] @ sp_["X"]
    out = {}
    for zz in zs:
        j = sp_["zs"].index(zz)
        d = {}
        for i, (rb, rk, a, b) in enumerate(ks.keys):
            v = GX[i, j]
            if abs(v) < 1e-12:
                continue
            m = site_at(np.array(rb))
            pk, tk = ks.snap(np.array(rk))
            k = site_at(pk)
            Rc = np.array(rk) - lim_pos[k] + (lim_pos[m] - np.array(rb))
            n = np.rint(np.linalg.solve(Pm, Rc)).astype(int)
            key = (int(n[0]), int(n[1]), int(n[2]), ket_index[(m, a)], ket_index[(k, b)])
            d[key] = d.get(key, 0) + (v.real if abs(v.imag) < 1e-12 else v)
        out[zz] = d
    return out


# ==================================================
def parent_shells(mm, rc, cells=3):
    """
    The neighbor shells of the parent structure, and the neighbors of the structure of the model in each.

    Returns:
        - (dict) -- {(tail, head): [{"shell" (from 1), "distance" (parent structure), "neighbor" ([n] of the model
          structure), "model" ([n] in the bonds of the model), "complete" (no neighbor or all neighbors in the model)}]}
          for the pairs of the bonds of the model, up to the last shell with a neighbor of the model.
    """
    from multipie.util.util_material_model import DISTANCE_DIGIT

    md, lim_pos = rc["md"], rc["lim_pos"]
    C = np.array(mm["unit_vector"], dtype=float).T
    L = _lattice(mm.group)
    Pm = np.eye(3) if L is None else L
    cell = list(itertools.product(range(-cells, cells + 1), repeat=3)) if L is not None else [(0, 0, 0)]
    out = {}
    for b in mm["bond"]["info"]:
        tail, head, nb = b.tail, b.head, list(b.neighbor)
        t = (tail, 1)
        rows = []
        for k, p in md["pos"].items():
            if k[0] != head:
                continue
            for n in cell:
                R = Pm @ np.array(n, dtype=float)
                d = round(float(np.linalg.norm(C @ (p + R - md["pos"][t]))), DISTANCE_DIGIT)
                dp = round(float(np.linalg.norm(C @ (lim_pos[k] + R - lim_pos[t]))), DISTANCE_DIGIT)
                if d > 0:
                    rows.append((d, dp))
        dist = sorted({d for d, _ in rows})
        nidx = {d: i + 1 for i, d in enumerate(dist)}
        pdist = sorted({dp for _, dp in rows if dp > 0})
        shells = []
        for s, dp in enumerate(pdist):
            ns = sorted({nidx[d] for d, x in rows if x == dp})
            inm = [x for x in ns if x in nb]
            shells.append({"shell": s + 1, "distance": dp, "neighbor": ns, "model": inm, "complete": not inm or inm == ns})
        last = max((i for i, d in enumerate(shells) if d["model"]), default=-1)
        out[(tail, head)] = shells[: last + 1]
    return out


# ==================================================
def map_bond(mm, rc, cells=2):
    """
    The bonds of the model moved by g on the parent structure.

    Returns:
        - (dict) -- {bond name: [{"bond": (tail atom, head atom, n), "image": (tail atom, head atom, n), "neighbor" (of
          the image in the structure of the model), "distance"}]} for the bonds from the first tail atom (atoms:
          (site, sublattice), n: conventional lattice vector).
    """
    from multipie.util.util_material_model import DISTANCE_DIGIT
    from multipie.util.util_supergroup import map_site

    md, lim_pos, g = rc["md"], rc["lim_pos"], rc["g"]
    C = np.array(mm["unit_vector"], dtype=float).T
    L = _lattice(mm.group)
    keys = list(lim_pos)
    img = map_site([lim_pos[k] for k in keys], g, L)
    out = {}
    shells = {}
    for name, rep in mm["bond"]["representative"].items():
        tail, head = rep.tail, rep.head
        t = (tail, 1)
        if (tail, head) not in shells:
            Pm = np.eye(3) if L is None else L
            cell = list(itertools.product(range(-3, 4), repeat=3)) if L is not None else [(0, 0, 0)]
            ds = set()
            for k, p in md["pos"].items():
                if k[0] == head:
                    for n in cell:
                        d = round(float(np.linalg.norm(C @ (p + Pm @ np.array(n, float) - md["pos"][t]))), DISTANCE_DIGIT)
                        if d > 0:
                            ds.add(d)
            shells[(tail, head)] = {d: i + 1 for i, d in enumerate(sorted(ds))}
        nidx = shells[(tail, head)]
        res = []
        Pm = np.eye(3) if L is None else L
        cell = list(itertools.product(range(-cells, cells + 1), repeat=3)) if L is not None else [(0, 0, 0)]
        for k, p in md["pos"].items():
            if k[0] != head:
                continue
            for n in cell:
                R = Pm @ np.array(n, dtype=float)
                d = round(float(np.linalg.norm(C @ (p + R - md["pos"][t]))), DISTANCE_DIGIT)
                if d == 0 or nidx.get(d) != rep.neighbor:
                    continue
                it, st = img[keys.index(t)]
                ih, sh = img[keys.index(k)]
                ti, hi = keys[it], keys[ih]
                Rimg = (
                    np.round(
                        g[:3, :3] @ (lim_pos[k] + R)
                        + g[:3, 3]
                        - (g[:3, :3] @ lim_pos[t] + g[:3, 3])
                        - (lim_pos[hi] - lim_pos[ti]),
                        6,
                    )
                    + 0.0
                )
                di = round(float(np.linalg.norm(C @ (md["pos"][hi] + Rimg - md["pos"][ti]))), DISTANCE_DIGIT)
                res.append(
                    {
                        "bond": (t, k, tuple(R.tolist())),
                        "image": (ti, hi, tuple(Rimg.tolist())),
                        "neighbor": nidx.get(di),
                        "distance": di,
                    }
                )
        out[name] = res
    return out
