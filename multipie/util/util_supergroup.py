# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Supergroups of a group H (PG or SG), not only the direct parents of the group binaries.

- Direct supergroups: Group.global_info()["supergroup"] of the binaries (MultiPieGen, supergroup.pkl): the supergroups
  P in which H is a maximal subgroup (t- or k-subgroup, [T_P : T_H] <= 4, not isomorphic), with an embedding A
  (x_P = A x_H, conventional coordinates), given up to the normalizers of P and of H.
- Supergroups (supergroups): the chains of direct supergroups, A = A_k ... A_1, the index multiplied, with
  [T_P : T_H] <= 4 (MAX_LATTICE_INDEX; the lattices of P finer and finer give infinitely many supergroups); one entry
  for each supergroup up to the normalizer of H, with its chain ("path").
- Embeddings (embedding): A n for every element n of the affine normalizer of H (n H n^-1 = H, within the entries
  searched; translations along a continuous direction of a polar group are not used), two embeddings being the same
  when they give the same set of the operations of P in the coordinates of H (the parent group seen from H).
  Different embeddings give different parent structures for the same atoms of H.
- Index: [P : H] modulo the lattice of H ("index") and modulo the lattice of P ("index_parent").
- Cosets (coset_decomposition): the left cosets g H of P (modulo the lattice of P or of H); the representative is the
  pure translation or the inversion if the coset has one, else the proper rotation of the highest order, with the
  smallest translation.
- Parent structure (parent_limit): the values of the Wyckoff parameters at which the orbit of H is invariant under P:
  c p(x) lies on the orbit for the representative c of every right coset H c (exact linear equations in x); the
  orbits whose atoms coincide are excluded.
"""

import functools
import itertools
from fractions import Fraction

import numpy as np
import sympy as sp

from multipie.core.group import Group
from multipie.util.util_crystal import P_dict

TOL = 1e-6
GRID = 24
MAX_LATTICE_INDEX = 4  # upper limit of [T_P : T_H] of the chains.
_HEX = np.array([[1, -0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1]])


# ==================================================
def gid(group):
    """id of a PG or SG, e.g. "PG:32", "SG:152"."""
    group = group if isinstance(group, Group) else Group(group)
    return group.info.PG if group.group_type == "PG" else group.info.SG


# ==================================================
def _float(M):
    return np.array(sp.Matrix(M).evalf(), dtype=float)


# ==================================================
def _mat4(W, t=None):
    A = np.eye(4)
    A[:3, :3] = W
    if t is not None:
        A[:3, 3] = t
    return A


# ==================================================
def _lattice(group):
    """primitive vectors (columns) in the conventional coordinates, None for PG."""
    if group.is_point_group:
        return None
    return np.array(P_dict[group.info.lattice], dtype=float)[:3, :3]


# ==================================================
def _frac(x, den=GRID * 4):
    return Fraction(x).limit_denominator(den)


# ==================================================
def _str(x):
    f = _frac(x)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


# ==================================================
def vec_str(v):
    return "[" + ",".join(_str(x + 0.0) for x in v) + "]"


# ==================================================
@functools.lru_cache(maxsize=None)
def _group_data(g):
    """point operations, translations, tags, and the operations (with the plus set) of a group (conventional)."""
    group = Group(g)
    so = group.symmetry_operation
    ops = [_float(W) for W in so["fractional"]]
    if group.is_point_group:
        ops = [_mat4(W) for W in ops]
    tags = list(so["tag"])
    L = _lattice(group)
    plus = [np.zeros(3)] if group.is_point_group else [np.array(c, dtype=float) for c in _float(so["plus_set"])]
    full = []
    for c in plus:
        for W in ops:
            W = W.copy()
            W[:3, 3] = W[:3, 3] + c
            full.append(W)
    point = {}
    for W, tag in zip(ops, tags):
        point.setdefault(_wkey(W[:3, :3]), (W[:3, :3], W[:3, 3], tag.split(":")[0]))
    return {"group": group, "L": L, "ops": ops, "full": full, "tags": tags, "plus": plus, "point": point}


# ==================================================
def _wkey(W):
    return tuple(np.rint(np.asarray(W, dtype=float) * 1e6).astype(np.int64).ravel())


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
def okey(W, L):
    """key of an operation (4x4) modulo the lattice L."""
    return _wkey(W[:3, :3]), tuple(np.rint(_reduce(W[:3, 3], L) * 1e5).astype(np.int64))


# ==================================================
def _generators(mats):
    """a small generating set of a finite matrix group (list of 3x3)."""
    keys = {_wkey(W): W for W in mats}
    gens, span = [], {_wkey(np.eye(3)): np.eye(3)}
    for W in sorted(keys.values(), key=lambda W: (-_order(W), _wkey(W))):
        if _wkey(W) in span:
            continue
        gens.append(W)
        frontier = list(span.values())
        while frontier:
            new = []
            for A in frontier:
                for G in gens:
                    B = A @ G
                    k = _wkey(B)
                    if k not in span:
                        span[k] = B
                        new.append(B)
            frontier = new
        if len(span) == len(keys):
            break
    return gens


# ==================================================
def _order(W):
    n, A = 1, np.asarray(W, dtype=float)
    while not np.allclose(A, np.eye(3)):
        A = A @ W
        n += 1
    return n


# ==================================================
def _nullspace(K):
    """rational basis (rows, RREF) of the null space of the integer matrix K."""
    M = sp.Matrix(np.rint(K).astype(int))
    ns = M.nullspace()
    if not ns:
        return None, []
    B = sp.Matrix.hstack(*ns).T.rref()[0]
    piv = [next(j for j in range(B.cols) if B[i, j] != 0) for i in range(B.rows)]
    return B, piv


# ==================================================
_SIGNED_PERM = [
    np.array(P, dtype=float) * np.array(sg, dtype=float)[:, None]
    for P in itertools.permutations(np.eye(3).tolist())
    for sg in itertools.product([1, -1], repeat=3)
]


# ==================================================
def _point_candidates(Hd, Pd):
    """linear parts N (primitive bases of H and P) mapping the point group of H into that of P (right-handed)."""
    LH = np.eye(3) if Hd["L"] is None else Hd["L"]
    LP = np.eye(3) if Pd["L"] is None else Pd["L"]
    Hp = [np.linalg.inv(LH) @ W @ LH for W, _, _ in Hd["point"].values()]
    Pp = [np.linalg.inv(LP) @ V @ LP for V, _, _ in Pd["point"].values()]
    Pkeys = {_wkey(V) for V in Pp}
    gens = _generators(Hp)
    sig = lambda W: (round(np.linalg.det(W)), round(np.trace(W)), _order(W))
    cands = [[V for V in Pp if sig(V) == sig(W)] for W in gens]
    I3 = np.eye(3)
    sols = {}

    def rec(i, K):
        if K is not None and np.linalg.matrix_rank(K, tol=1e-8) >= 9:
            return
        if i == len(gens):
            B, piv = _nullspace(K)
            if B is not None:
                sols.setdefault(tuple(np.array(B.evalf(), dtype=float).ravel().round(8)), (B, piv))
            return
        for V in cands[i]:
            Ki = np.kron(I3, gens[i].T) - np.kron(V, I3)  # vec(N W - V N) = Ki vec(N) (row major).
            rec(i + 1, Ki if K is None else np.vstack([K, Ki]))

    pg = Hd["L"] is None
    Ns = []
    if not gens:
        Ns = list(_SIGNED_PERM)
    else:
        rec(0, None)
        for B, piv in sols.values():
            d = len(piv)
            Bf = np.array(B.evalf(), dtype=float)
            if d > 5:
                Ns += [N for N in _SIGNED_PERM if np.allclose(np.linalg.lstsq(Bf.T, N.ravel(), rcond=None)[0] @ Bf, N.ravel())]
                continue
            vals = [-1, -0.5, 0, 0.5, 1] if pg else (list(range(-3, 4)) if d <= 3 else list(range(-2, 3)))
            for c in itertools.product(vals, repeat=d):
                Ns.append((np.array(c, dtype=float) @ Bf).reshape(3, 3))
    out = {}
    for N in Ns:
        if np.linalg.det(N) < 1e-6:
            continue
        if not pg and not np.allclose(N, np.rint(N), atol=1e-8):
            continue
        Ni = np.linalg.inv(N)
        if not all(_wkey(N @ W @ Ni) in Pkeys for W in Hp):
            continue
        if pg:
            CH = _HEX if Hd["group"].is_hexagonal_subgroup else np.eye(3)
            CP = _HEX if Pd["group"].is_hexagonal_subgroup else np.eye(3)
            O = CP @ N @ np.linalg.inv(CH)
            if not np.allclose(O @ O.T, np.eye(3), atol=1e-6):
                continue
        out.setdefault(_wkey(N), N)
    return list(out.values())


# ==================================================
def _translation_part(Pd, V):
    """translation of P for the point operation V (conventional), None if V is not in P."""
    t = Pd["point"].get(_wkey(V))
    return None if t is None else t[1]


# ==================================================
_grid_cache = {}


def _grid():
    if GRID not in _grid_cache:
        r = np.arange(GRID) / GRID
        _grid_cache[GRID] = np.array(list(itertools.product(r, r, r)))
    return _grid_cache[GRID]


# ==================================================
def _shifts(Hd, Pd, R):
    """shifts s (conventional of P) with x_P = R x_H + s mapping H into P (on the grid of 1/24, primitive of P)."""
    LP = Pd["L"]
    LPi = np.linalg.inv(LP)
    S = _grid()
    ok = np.ones(len(S), dtype=bool)
    Ri = np.linalg.inv(R)
    for W in Hd["full"]:
        V = R @ W[:3, :3] @ Ri
        v = _translation_part(Pd, V)
        if v is None:
            return []
        a = LPi @ (R @ W[:3, 3] - v)
        M = LPi @ (np.eye(3) - V) @ LP
        res = a[None, :] + S[ok] @ M.T
        good = np.all(np.abs(res - np.rint(res)) < 1e-6, axis=1)
        idx = np.nonzero(ok)[0]
        ok[idx[~good]] = False
        if not ok.any():
            return []
    return [LP @ s for s in S[ok]]


# ==================================================
def _parent_ops_in_h(Pd, A):
    """operations of P (with the plus set) in the coordinates of H, and the lattice of P in the coordinates of H."""
    Ai = np.linalg.inv(A)
    ops = [Ai @ W @ A for W in Pd["full"]]
    L = None if Pd["L"] is None else Ai[:3, :3] @ Pd["L"]
    return ops, L


# ==================================================
def _class_key(Hd, Pd, A):
    """the operations of P in the coordinates of H modulo the lattice of P, and that lattice modulo the lattice of H."""
    ops, L = _parent_ops_in_h(Pd, A)
    lat = None
    if L is not None:
        LH = Hd["L"]
        gens = [np.linalg.solve(LH, L[:, i]) for i in range(3)]
        pts = {tuple(np.zeros(3))}
        frontier = [np.zeros(3)]
        while frontier:
            new = []
            for t in frontier:
                for g in gens:
                    u = t + g
                    u = u - np.floor(u + TOL)
                    k = tuple(np.round(u, 5) + 0.0)
                    if k not in pts:
                        pts.add(k)
                        new.append(u)
            frontier = new
        lat = frozenset(pts)
        shifts = [np.array(t) for t in pts]
        keys = set()
        for W in ops:
            c = np.linalg.solve(LH, W[:3, 3])
            cand = []
            for t in shifts:
                u = c + t
                u = u - np.floor(u + TOL)
                cand.append(tuple(np.rint(u * 1e5).astype(np.int64)))
            keys.add((_wkey(W[:3, :3]), min(cand)))
        return lat, frozenset(keys)
    return lat, frozenset(okey(W, L) for W in ops)


# ==================================================
def _cost(A):
    """ordering of equivalent embeddings: the simplest R, then the smallest shift."""
    R, s = A[:3, :3], A[:3, 3]
    return (
        int(np.sum(np.diag(R) < 0)),
        float(np.sum(np.abs(R - np.diag(np.diag(R))))),
        float(np.sum(np.abs(R))),
        float(np.sum(np.minimum(s % 1, 1 - s % 1))),
        tuple(np.round(s % 1, 6)),
    )


# ==================================================
def _direct(h):
    """direct supergroups of the binaries: [(parent, A (4x4 float), index, index_parent)]."""
    data = Group.global_info().get("supergroup")
    if data is None:
        raise ValueError("the group binaries have no supergroup data (MultiPieGen supergroup.pkl).")
    return [(d["parent"], _float(d["A"]), int(d["index"]), int(d["index_parent"])) for d in data.get(h, [])]


# ==================================================
@functools.lru_cache(maxsize=None)
def _normalizer(h):
    """elements n of the affine normalizer of H (4x4, coordinates of H), modulo the lattice of H."""
    Hd = _group_data(h)
    LH = np.eye(3) if Hd["L"] is None else Hd["L"]
    out = {}
    if Hd["L"] is None:
        for M in _point_candidates(Hd, Hd):
            out.setdefault(_wkey(M), _mat4(M))
        return list(out.values())
    # continuous translations (polar directions): not used.
    K = np.vstack([np.eye(3) - W for W, _, _ in Hd["point"].values()])
    _, sv, Vt = np.linalg.svd(K)
    C = Vt[int(np.sum(sv > 1e-8)) :]
    for M in _point_candidates(Hd, Hd):
        if round(abs(np.linalg.det(M))) != 1:
            continue
        R = LH @ M @ np.linalg.inv(LH)
        R = np.where(np.abs(R - np.rint(R * 6) / 6) < 1e-9, np.rint(R * 6) / 6, R) + 0.0
        for t in _shifts(Hd, Hd, R):
            if len(C) and np.linalg.norm(C @ t) > 1e-6:
                continue
            n = _mat4(R, _reduce(t, LH))
            out.setdefault(okey(n, LH), n)
    return list(out.values())


# ==================================================
@functools.lru_cache(maxsize=None)
def _reach(h):
    """ids of the groups reached from H by the direct supergroups."""
    out, frontier = set(), [h]
    while frontier:
        new = []
        for q in frontier:
            for p, *_ in _direct(q):
                if p not in out:
                    out.add(p)
                    new.append(p)
        frontier = new
    return frozenset(out)


# ==================================================
@functools.lru_cache(maxsize=None)
def _chains(h, max_index, target=None):
    """
    supergroups along the chains of direct supergroups: [(P, A, index, index_parent, path)], one for each supergroup up
    to the normalizer of H (only the chains to the target if given).
    """
    Hd = _group_data(h)
    nodes, seen = [], set()
    frontier = [(h, np.eye(4), 1, 1, [h])]
    while frontier:
        new = []
        for q, A, n, n_p, path in frontier:
            for p, Ad, nd, ndp in _direct(q):
                if target is not None and p != target and target not in _reach(p):
                    continue
                n2, np2 = n * nd, n_p * ndp
                if (max_index is not None and n2 > max_index) or n2 // np2 > MAX_LATTICE_INDEX:
                    continue
                A2 = Ad @ A
                Pd = _group_data(p)
                A2[:3, 3] = _reduce(A2[:3, 3], Pd["L"])
                if (p, _class_key(Hd, Pd, A2)) in seen:
                    continue
                for m in _normalizer(h):
                    seen.add((p, _class_key(Hd, Pd, A2 @ m)))
                node = (p, A2, n2, np2, path + [p])
                nodes.append(node)
                new.append(node)
        frontier = new
    if target is not None:
        nodes = [v for v in nodes if v[0] == target]
    return sorted(nodes, key=lambda v: (v[2], int(v[0].split(":")[1]), _cost(v[1])))


# ==================================================
@functools.lru_cache(maxsize=None)
def _embeddings(h, p):
    Hd, Pd = _group_data(h), _group_data(p)
    classes = {}
    for q, A, n, n_p, path in _chains(h, None, p):
        for m in sorted(_normalizer(h), key=_cost):
            B = A @ m
            B[:3, 3] = _reduce(B[:3, 3], Pd["L"])
            k = _class_key(Hd, Pd, B)
            if k not in classes:
                classes[k] = (B, n, n_p, tuple(path))
    return list(classes.values())


# ==================================================
def embedding(group, parent):
    """
    Embeddings of the group H in a supergroup P.

    Args:
        group (Group or str): group H.
        parent (Group or str): supergroup P.

    Returns:
        - (list) -- [{"parent", "A" (x_P = A x_H, 4x4), "R", "shift", "origin" (the origin of P in the coordinates
          of H), "index", "index_parent", "path"}]; [] if P is not a supergroup.
    """
    h, p = gid(group), gid(parent)
    L = _lattice(Group(h))
    out = []
    for A, n, n_p, path in _embeddings(h, p):
        R, s = A[:3, :3], A[:3, 3]
        out.append(
            {
                "parent": p,
                "A": A,
                "R": R,
                "shift": s,
                "origin": _reduce(-np.linalg.solve(R, s), L),
                "index": n,
                "index_parent": n_p,
                "path": list(path),
            }
        )
    return out


# ==================================================
def resolve_embedding(group, parent, embedding=None):
    """
    The embedding of H in P used by the other functions.

    Args:
        group (Group or str): group H.
        parent (Group or str): supergroup P.
        embedding (None, int, dict or array, optional): None (the first one), the index in embedding(), a dict of
            embedding(), or A (4x4, x_P = A x_H).

    Returns:
        - (dict) -- embedding (see embedding()).
    """
    if isinstance(embedding, dict):
        return embedding
    if embedding is not None and not isinstance(embedding, (int, np.integer)):
        A = np.array(embedding, dtype=float)
        if A.shape != (4, 4):
            raise ValueError(f"embedding must be None, an index, a dict or a 4x4 matrix, but {embedding!r} is given.")
        h, p = gid(group), gid(parent)
        Hd, Pd = _group_data(h), _group_data(p)
        k = _class_key(Hd, Pd, A)
        for e in _emb_list(h, p):
            if _class_key(Hd, Pd, e["A"]) == k:
                return e
        raise ValueError(f"the matrix {A.tolist()} is not an embedding of {h} in {p}.")
    es = _emb_list(gid(group), gid(parent))
    if not es:
        raise ValueError(f"{gid(parent)} is not a supergroup of {gid(group)}.")
    k = 0 if embedding is None else int(embedding)
    if not 0 <= k < len(es):
        raise ValueError(f"embedding {k} is out of range (0-{len(es) - 1}).")
    return es[k]


def _emb_list(h, p):
    return embedding(h, p)


# ==================================================
def _translations(Hd, Pd, A):
    """translations of the lattice of P modulo that of H (coordinates of H)."""
    if Hd["L"] is None:
        return [np.zeros(3)]
    LH = Hd["L"]
    L = np.linalg.inv(A)[:3, :3] @ Pd["L"]
    gens = [L[:, i] for i in range(3)]
    T = {_tkey(np.zeros(3), LH): np.zeros(3)}
    frontier = [np.zeros(3)]
    while frontier:
        new = []
        for t in frontier:
            for g in gens:
                u = _reduce(t + g, LH)
                k = _tkey(u, LH)
                if k not in T:
                    T[k] = u
                    new.append(u)
        frontier = new
    return list(T.values())


def _tkey(v, L):
    return tuple(np.rint(_reduce(v, L) * 1e5).astype(np.int64))


# ==================================================
def parent_elements(group, parent, embedding=None):
    """
    Elements of P modulo the lattice of H, in the coordinates of H.

    Returns:
        - (list) -- [4x4 matrix].
        - (dict) -- embedding.
    """
    h, p = gid(group), gid(parent)
    e = resolve_embedding(h, p, embedding)
    Hd, Pd = _group_data(h), _group_data(p)
    ops, _ = _parent_ops_in_h(Pd, e["A"])
    trans = _translations(Hd, Pd, e["A"])
    out = {}
    for W in ops:
        for t in trans:
            M = W.copy()
            M[:3, 3] = _reduce(W[:3, 3] + t, Hd["L"])
            out.setdefault(okey(M, Hd["L"]), M)
    return list(out.values()), e


# ==================================================
def op_tag(group, parent, e, W):
    """tag of an element W (coordinates of H) of P: point operation of P and the translation in H."""
    Pd = _group_data(gid(parent))
    R = e["R"]
    V = R @ W[:3, :3] @ np.linalg.inv(R)
    t = Pd["point"][_wkey(V)][2]
    if Pd["L"] is None:
        return t
    return t + ":" + vec_str(W[:3, 3])


# ==================================================
def op_matrix(group, parent, tag, embedding=None):
    """
    Element of P (coordinates of H) from its tag "point operation of P:[translation in H]" (or a 4x4 matrix).

    Returns:
        - (ndarray) -- 4x4 matrix.
    """
    h, p = gid(group), gid(parent)
    e = resolve_embedding(h, p, embedding)
    Hd, Pd = _group_data(h), _group_data(p)
    if not isinstance(tag, str):
        W = np.array(tag, dtype=float)
        if W.shape == (3, 3):
            W = _mat4(W)
    else:
        if ":" in tag:
            pt, t = tag.split(":")
            t = np.array([float(sp.Rational(sp.sympify(v))) for v in t.strip("[]").split(",")])
        else:
            pt, t = tag, np.zeros(3)
        R = e["R"]
        Vs = [V for V, _, tg in Pd["point"].values() if tg == pt]
        if not Vs:
            raise ValueError(f"symmetry operation {tag} is not found in {p}.")
        W = _mat4(np.linalg.inv(R) @ Vs[0] @ R, t)
    els, _ = parent_elements(h, p, e)
    keys = {okey(M, Hd["L"]) for M in els}
    if okey(W, Hd["L"]) not in keys:
        raise ValueError(f"{tag} is not an element of {p} (embedding origin {vec_str(e['origin'])}).")
    return W


# ==================================================
def _h_keys(Hd, L):
    return {okey(W, L) for W in Hd["full"]}


# ==================================================
def _simplicity(W, L):
    """ordering of the elements for the choice of the representatives."""
    t = _reduce(W[:3, 3], L)
    det = np.linalg.det(W[:3, :3])
    Wp = W[:3, :3]
    kind = 0 if np.allclose(Wp, np.eye(3)) else (1 if np.allclose(Wp, -np.eye(3)) else (2 if det > 0 else 3))
    return (
        kind,
        -_order(Wp),
        float(np.sum(np.minimum(np.abs(t), 1 - np.abs(t)) ** 2)),
        _wkey(Wp),
        tuple(np.round(t, 6)),
    )


# ==================================================
def coset_decomposition(group, parent, embedding=None, modulo="parent"):
    """
    Left cosets g H of P.

    Args:
        group (Group or str): group H.
        parent (Group or str): supergroup P.
        embedding (None, int, dict or array, optional): embedding (resolve_embedding).
        modulo (str, optional): "parent" (modulo the lattice of P) or "self" (modulo the lattice of H, the cosets of
            the translations of P included).

    Returns:
        - (list) -- [{"rep" (tag), "rep_matrix" (coordinates of H), "ops" (tags), "order" (smallest n with g^n in H),
          "kind" ("identity", "rotation", "improper", "inversion", "translation")}], H first.
        - (dict) -- embedding.
    """
    if modulo not in ("parent", "self"):
        raise ValueError(f"modulo must be 'parent' or 'self', but {modulo!r} is given.")
    h, p = gid(group), gid(parent)
    Hd = _group_data(h)
    els, e = parent_elements(h, p, embedding)
    L = Hd["L"]
    if modulo == "parent" and L is not None:
        L = np.linalg.inv(e["A"])[:3, :3] @ _group_data(p)["L"]
    elem = {}
    for W in els:
        M = W.copy()
        M[:3, 3] = _reduce(W[:3, 3], L)
        elem.setdefault(okey(M, L), M)
    hk = _h_keys(Hd, L)
    Hel = [W for k, W in elem.items() if k in hk]
    chiral = all(np.linalg.det(W[:3, :3]) > 0 for W in Hd["full"])
    order = sorted(elem.values(), key=lambda W: (okey(W, L) not in hk, _simplicity(W, L)))
    done, out = set(), []
    for g in order:
        k = okey(g, L)
        if k in done:
            continue
        cos = []
        for hh in Hel:
            M = g @ hh
            M[:3, 3] = _reduce(M[:3, 3], L)
            kk = okey(M, L)
            if kk not in done:
                done.add(kk)
                cos.append(M)
        cos.sort(key=lambda W: _simplicity(W, L))
        rep = cos[0] if k not in hk else _mat4(np.eye(3))
        n, G = 1, rep.copy()
        while okey(G, L) not in hk:
            G = G @ rep
            n += 1
        if k in hk:
            kind = "identity"
        elif np.allclose(rep[:3, :3], np.eye(3)):
            kind = "translation"
        elif np.linalg.det(rep[:3, :3]) > 0:
            kind = "rotation"
        else:
            kind = "inversion" if chiral else "improper"
        out.append(
            {
                "rep": op_tag(h, p, e, rep),
                "rep_matrix": rep,
                "ops": [op_tag(h, p, e, W) for W in cos],
                "order": n,
                "kind": kind,
            }
        )
    return out, e


# ==================================================
def _right_coset_reps(Hd, els):
    """representatives c of the right cosets H c (modulo the lattice of H)."""
    L = Hd["L"]
    done, reps = set(), []
    for c in sorted(els, key=lambda W: _simplicity(W, L)):
        if okey(c, L) in done:
            continue
        reps.append(c)
        for hh in Hd["full"]:
            M = hh @ c
            done.add(okey(M, L))
    return reps


# ==================================================
def _wyckoff_affine(group, w):
    """symbols, and (A, b) of the representative of the Wyckoff site w (conventional)."""
    rep = sp.Matrix(group.wyckoff["site"][w]["conventional"][0])
    syms = sorted(rep.free_symbols, key=str)
    A = np.array([[float(sp.diff(r, s)) for s in syms] for r in rep], dtype=float).reshape(3, len(syms))
    b = np.array([float(r.subs({s: 0 for s in syms})) for r in rep], dtype=float)
    return syms, A, b


# ==================================================
def _orbit(Hd, q):
    L = Hd["L"]
    pts = {}
    for W in Hd["full"]:
        r = _reduce(W[:3, :3] @ q + W[:3, 3], L)
        pts.setdefault(tuple(np.rint(r * 1e5).astype(np.int64)), r)
    return list(pts.values())


# ==================================================
def _rref(rows, k):
    """exact RREF (tuple of rows, Fractions) of the augmented rows [E | f]; None if inconsistent."""
    rows = [list(r) for r in rows]
    out, col = [], 0
    m = len(rows)
    r0 = 0
    for col in range(k):
        piv = next((i for i in range(r0, m) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[r0], rows[piv] = rows[piv], rows[r0]
        pv = rows[r0][col]
        rows[r0] = [x / pv for x in rows[r0]]
        for i in range(m):
            if i != r0 and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r0])]
        r0 += 1
    for i in range(r0, m):
        if rows[i][k] != 0:
            return None
    return tuple(tuple(r) for r in rows[:r0])


def _rat(x):
    return sp.Rational(_frac(x, 10**4).numerator, _frac(x, 10**4).denominator)


# ==================================================
def parent_limit(group, parent, wyckoff=None, embedding=None):
    """
    Parent structures of the Wyckoff sites of H for a supergroup P.

    Args:
        group (Group or str): group H.
        parent (Group or str): supergroup P.
        wyckoff (str or list, optional): Wyckoff sites (default: all with free parameters).
        embedding (None, int, dict or array, optional): embedding (resolve_embedding).

    Returns:
        - (dict) -- {wyckoff: [{"param": {name: expression}, "position" (representative), "parent_wyckoff"}]}, the
          same format as Group.coset["parent_limit"] (parameters modulo 1).
        - (dict) -- embedding.
    """
    h, p = gid(group), gid(parent)
    G = Group(h)
    Hd = _group_data(h)
    els, e = parent_elements(h, p, embedding)
    reps = [c for c in _right_coset_reps(Hd, els)]
    L = Hd["L"]
    LH = np.eye(3) if L is None else L
    if wyckoff is None:
        ws = [w for w in G.wyckoff["site"] if sp.Matrix(G.wyckoff["site"][w]["conventional"][0]).free_symbols]
    else:
        ws = [wyckoff] if isinstance(wyckoff, str) else list(wyckoff)
    box = list(itertools.product(range(-2, 3), repeat=3)) if L is not None else [(0, 0, 0)]
    Pg = Group(p)
    res = {}
    for w in ws:
        syms, A, b = _wyckoff_affine(G, w)
        k = len(syms)
        rng = np.random.default_rng(7)
        generic = rng.uniform(0.05, 0.45, size=k)
        mult = len(_orbit(Hd, A @ generic + b))
        # candidate affine subspaces for each representative c: exact RREF (rows of [E | f]) of E x = f.
        systems = {(): ()}
        for c in reps:
            cand = {}
            for W in Hd["full"]:
                M = (c[:3, :3] - W[:3, :3]) @ A
                for n in box:
                    rhs = (W[:3, :3] - c[:3, :3]) @ b + W[:3, 3] - c[:3, 3] + LH @ np.array(n, dtype=float)
                    rows = [[_frac(M[i, j], 48) for j in range(k)] + [_frac(rhs[i], 48)] for i in range(3)]
                    red = _rref(rows, k)
                    if red is not None:
                        cand[red] = True
            new = {}
            for sys0 in systems:
                for sys1 in cand:
                    red = _rref([list(r) for r in sys0] + [list(r) for r in sys1], k)
                    if red is not None:
                        new[red] = True
            systems = new
            if not systems:
                break
        lims = {}
        for red in systems:
            if k == 0 or not red:
                continue  # no constraint: no parameter is fixed.
            piv = [next(j for j in range(k) if r[j] != 0) for r in red]
            free = [s for j, s in enumerate(syms) if j not in piv]
            expr = {s: s for s in free}
            for r, j in zip(red, piv):
                c0 = r[k] - (r[k] // 1)  # modulo 1.
                v = sp.Rational(c0.numerator, c0.denominator) - sum(
                    sp.Rational(r[jj].numerator, r[jj].denominator) * syms[jj] for jj in range(k) if jj not in piv
                )
                expr[syms[j]] = v
            key = tuple(str(expr[s]) for s in syms)
            if key in lims:
                continue
            gv = {s: sp.Rational(int(rng.integers(1, 97)), 97) for s in free}
            xv = sp.Matrix([expr[s] for s in syms]).subs(gv)
            q = np.array((sp.Matrix(A) * xv).evalf(), dtype=float).ravel() + b
            if len(_orbit(Hd, q)) != mult:
                continue  # atoms coincide.
            pos = sp.Matrix(A) * sp.Matrix([expr[s] for s in syms]) + sp.Matrix([_rat(v) for v in b])
            pos = [sp.nsimplify(v) for v in pos]
            qP = e["A"] @ np.append(q, 1.0)
            pw = Pg.find_wyckoff_site("[" + ",".join(f"{v:.12f}" for v in qP[:3]) + "]")[0]
            lims[key] = {
                "param": {str(s): str(expr[s]) for s in syms},
                "position": "[" + ",".join(str(v).replace(" ", "") for v in pos) + "]",
                "parent_wyckoff": pw,
            }
        res[w] = sorted(
            lims.values(), key=lambda d: [float(sp.sympify(v).subs({s: 0 for s in syms})) for v in d["param"].values()]
        )
    return res, e


# ==================================================
def map_site(positions, op, lattice=None):
    """
    Atoms of a structure mapped by an operation.

    Args:
        positions (list): positions of the atoms (coordinates of H).
        op (ndarray): operation (4x4, coordinates of H).
        lattice (ndarray, optional): primitive vectors of H (columns, conventional), None for PG.

    Returns:
        - (dict) -- {i: (j, n)}: op r_i = r_j + n (n: lattice vector, conventional).
    """
    pos = [np.asarray(p, dtype=float) for p in positions]
    out = {}
    for i, p in enumerate(pos):
        q = op[:3, :3] @ p + op[:3, 3]
        for j, r in enumerate(pos):
            d = q - r
            if lattice is None:
                if np.allclose(d, 0, atol=TOL):
                    out[i] = (j, np.zeros(3))
                    break
            else:
                c = np.linalg.solve(lattice, d)
                if np.allclose(c, np.rint(c), atol=TOL):
                    out[i] = (j, lattice @ np.rint(c) + 0.0)
                    break
        else:
            raise ValueError(f"the structure is not invariant under the operation: atom {i}.")
    return out


# ==================================================
def supergroups(group, max_index=None, chain=False):
    """
    Supergroups of a group H (not only the direct parents).

    Args:
        group (Group or str): group H (PG or SG).
        max_index (int, optional): upper limit of the index (modulo the lattice of H).
        chain (bool, optional): only the chain of the representative parents of the group binaries
            (Group.global_info()["parent"]).

    Returns:
        - (list) -- [{"group", "A", "R", "shift", "origin", "index", "index_parent", "path" (H, ..., P: a chain of
          direct supergroups), "direct"}], in order of the index; one entry for each supergroup up to the normalizer
          of H (embedding() gives all).
    """
    h = gid(group)
    L = _lattice(Group(h))
    out = []
    for p, A, n, n_p, path in _chains(h, max_index):
        R, s = A[:3, :3], A[:3, 3]
        out.append(
            {
                "group": p,
                "A": A,
                "R": R,
                "shift": s,
                "origin": _reduce(-np.linalg.solve(R, s), L),
                "index": n,
                "index_parent": n_p,
                "path": path,
                "direct": len(path) == 2,
            }
        )
    if chain:
        ch = Group.global_info()["parent"].get(h, [])
        return [next(d for d in out if d["group"] == p) for p in ch if any(d["group"] == p for d in out)]
    return out
