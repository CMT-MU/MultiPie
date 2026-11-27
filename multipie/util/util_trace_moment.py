# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Trace moments of a SAMB Hamiltonian, the model-dependent part of the response tensors.

For the Hamiltonian h = sum_j z_j Z_j of a model and operators O_1, ..., O_n,

    Gamma^{i_1 ... i_n}_{O_1; ...; O_n} = (1/N) sum_k Tr[O_1(k) h(k)^{i_1} O_2(k) h(k)^{i_2} ... O_n(k) h(k)^{i_n}]

(per unit cell, N: number of k points), Eqs. (31)-(33) of R. Oiwa and H. Kusunose, J. Phys. Soc. Jpn. 91, 014701 (2022).
Gamma is a polynomial in z_j: its lowest nonzero terms show the essential parameters of a response.

- Operators. A matrix dict {(n1,n2,n3,m,n): value} (as the SAMB matrices of MultiPie, the indices of the full matrix
  from 0, value may be (value, bond_no)), a file in the format of seedname_hr.dat of Wannier90 (indices from 1,
  divided by the Wigner-Seitz degeneracy), the name of a SAMB of the model (zj), or the velocity operator of the
  Hamiltonian "v_a1a2...", v_{a1...an}(k) = d_{a1} ... d_{an} h(k) (a = x, y, z, Cartesian; hbar = 1; the length unit
  of the cell). In real space, v_{a1...an}(R)_{mn} = (i x_{a1}) ... (i x_{an}) h(R)_{mn} with the bond vector
  x = (R + r_n - r_m) A, so that v is linear in z_j.
- Exact k sum. O(k) = sum_R O(R) exp(2 pi i k.R) (lattice gauge; a trace does not depend on the gauge). A product of
  operators whose ranges are r^(1), r^(2), ... in a_i reaches at most sum r^(.), and the uniform grid of N_i = sum r^(.)
  + 1 points gives the exact k sum (the R = 0 component of the product).
- Numerical values. With h(k) = U e U^dagger at each k, O~ = U^dagger O U and V_{i,b} = e_b^i (or T_i(e_b) in the
  Chebyshev form), Gamma^{i_1...i_n} = sum_{b_0...b_{n-1}} O~1_{b0 b1} V_{i1 b1} O~2_{b1 b2} V_{i2 b2} ... O~n_{b_{n-1} b0} V_{in b0}:
  one diagonalization gives all orders.
- Polynomials. The word O_1 h^{i_1} O_2 ... is multiplied as a dict {monomial: matrix(k)} with h = h_0 + sum_{j in V} z_j Z_j
  (V: the variables, h_0: the fixed part), and the trace is summed on the exact grid. A polynomial is identically
  zero if its value at random z vanishes; the variables on which Gamma does not depend are dropped before the expansion.
"""

import itertools
import re
from fractions import Fraction

import numpy as np
import sympy as sp

from multipie.util.util_model_analyzer import fourier_r_to_k, grid_all
from multipie.util.util_wannier import read_hr, apply_ws_degeneracy

_AXIS = {"x": 0, "y": 1, "z": 2}
_VELOCITY = re.compile(r"v_([xyz]+)")
_BUDGET = 64e6  # bytes of a work array.
_TOL = 1e-9  # relative tolerance of zero.


# ==================================================
def matrix_dict(O, dim):
    """
    Normalize an operator matrix dict.

    Args:
        O (dict): dict[(n1,n2,n3,m,n), value or (value, bond_no)] or dict[((n1,n2,n3),m,n), value].
        dim (int): dimension of the full matrix.

    Returns:
        - (dict) -- dict[((n1,n2,n3),m,n), complex], without zeros.
    """
    out = {}
    for key, v in O.items():
        if isinstance(v, tuple):
            v = v[0]
        if len(key) == 5:
            R, m, n = tuple(int(i) for i in key[:3]), int(key[3]), int(key[4])
        elif len(key) == 3 and len(key[0]) == 3:
            R, m, n = tuple(int(i) for i in key[0]), int(key[1]), int(key[2])
        else:
            raise ValueError(f"invalid key of operator matrix, {key}.")
        if not (0 <= m < dim and 0 <= n < dim):
            raise ValueError(f"index ({m},{n}) of operator matrix is out of the dimension {dim}.")
        v = complex(v)
        if v != 0:
            out[(R, m, n)] = out.get((R, m, n), 0) + v

    return out


# ==================================================
def read_operator_hr(filename, topdir=None):
    """
    Read an operator in the format of seedname_hr.dat.

    Args:
        filename (str): file name.
        topdir (str, optional): directory of the file.

    Returns:
        - (dict) -- dict[(n1,n2,n3,m,n), value] (indices from 0, divided by the Wigner-Seitz degeneracy).
    """
    hr, irvec, ndegen = read_hr(filename, topdir)
    return apply_ws_degeneracy(hr, irvec, ndegen)


# ==================================================
def derivative_matrix(OR, atom, A, axis):
    """
    k derivative of an operator in real space.

    Args:
        OR (dict): dict[((n1,n2,n3),m,n), value].
        atom (ndarray): positions of the kets (fractional).
        A (ndarray): primitive vectors, [a1, a2, a3].
        axis (str): Cartesian axes, e.g., "xy" for d_x d_y.

    Returns:
        - (dict) -- d_{a1} ... d_{an} O, dict[((n1,n2,n3),m,n), value], (i x_{a1}) ... (i x_{an}) O(R)_{mn}, x = (R + r_n - r_m) A.
    """
    idx = [_AXIS[a] for a in axis]
    out = {}
    for (R, m, n), v in OR.items():
        x = (np.asarray(R, dtype=float) + atom[n] - atom[m]) @ A
        f = np.prod([1j * x[i] for i in idx])
        if abs(f) > 1e-14:
            out[(R, m, n)] = complex(f * v)

    return out


# ==================================================
def matrix_range(OR):
    """
    Range of an operator.

    Args:
        OR (dict): dict[((n1,n2,n3),m,n), value].

    Returns:
        - (ndarray) -- maximum of abs(n_i) for i = 1, 2, 3.
    """
    if not OR:
        return np.zeros(3, dtype=int)
    return np.max(np.abs(np.asarray([R for R, _, _ in OR], dtype=int)), axis=0)


# ==================================================
def order_tuples(n, total):
    """
    Orders (i_1, ..., i_n) with sum i = total in lexicographic order.

    Args:
        n (int): number of operators.
        total (int): total order.

    Returns:
        - (list) -- list of tuples.
    """
    return [t for t in itertools.product(range(total + 1), repeat=n) if sum(t) == total]


# ==================================================
def _simplify(c, simplify):
    """
    Simplify a real coefficient.
    """
    if simplify:
        for n in (1, 2, 3, 5, 6, 7):  # c = (p/q) sqrt(n) with q <= 1000.
            x = c / np.sqrt(n)
            f = Fraction(x).limit_denominator(1000)
            if abs(float(f) - x) <= 1e-9 * max(1, abs(x)):
                return sp.Rational(f.numerator, f.denominator) * sp.sqrt(n)
    return sp.Float(c, 12)


# ==================================================
class TraceMoment:
    # ==================================================
    def __init__(self, samb, atom, A, parameter=None, operator=None, topdir=None):
        """
        Trace moments of a SAMB Hamiltonian.

        Args:
            samb (dict): SAMB matrices, dict[zj, dict[(n1,n2,n3,m,n), value or (value, bond_no)]].
            atom (array-like): positions of the kets (fractional, primitive).
            A (array-like): primitive vectors, [a1, a2, a3].
            parameter (dict, optional): values of zj (missing: 0).
            operator (dict, optional): operators, dict[name, spec] (see add_operator).
            topdir (str, optional): directory of the operator files.
        """
        self._atom = np.asarray(atom, dtype=float)
        self._A = np.asarray(A, dtype=float)
        self._dim = len(self._atom)
        self._topdir = topdir
        self._Z = {zj: matrix_dict(m, self._dim) for zj, m in samb.items()}
        self._names = list(self._Z.keys())
        self._parameter = {zj: 0.0 for zj in self._names}
        if parameter:
            self.set_parameter(parameter)
        self._h = {zj: Z for zj, Z in self._Z.items()}  # form of h, dict[zj, matrix].
        self._operator = {}
        if operator:
            for name, spec in operator.items():
                self.add_operator(name, spec)

    # ==================================================
    @property
    def variable(self):
        """
        Names of the SAMB, zj.

        Returns:
            - (list) -- list of zj.
        """
        return list(self._names)

    # ==================================================
    @property
    def operator(self):
        """
        Names of the operators.

        Returns:
            - (list) -- list of names.
        """
        return list(self._operator.keys())

    # ==================================================
    def set_parameter(self, parameter):
        """
        Set values of zj.

        Args:
            parameter (dict): dict[zj, value].
        """
        for zj, v in parameter.items():
            if zj not in self._parameter:
                raise KeyError(f"unknown parameter {zj}.")
            self._parameter[zj] = float(v)

    # ==================================================
    def add_operator(self, name, spec):
        """
        Add an operator.

        Args:
            name (str): name of the operator.
            spec (dict or str): matrix dict {(n1,n2,n3,m,n): value} (indices from 0), file name of the format of
                seedname_hr.dat (indices from 1), SAMB name zj, or the velocity "v_a1a2..." (d_{a1}...d_{an} h).
        """
        if isinstance(spec, dict):
            form = {None: matrix_dict(spec, self._dim)}
        elif isinstance(spec, str):
            m = _VELOCITY.fullmatch(spec)
            if m:
                form = {zj: derivative_matrix(Z, self._atom, self._A, m.group(1)) for zj, Z in self._Z.items()}
            elif spec in self._Z:
                form = {None: self._Z[spec]}
            else:
                form = {None: matrix_dict(read_operator_hr(spec, self._topdir), self._dim)}
        else:
            raise ValueError(f"invalid operator spec for {name}.")
        self._operator[name] = form

    # ==================================================
    def _check_ops(self, ops):
        """
        Check names of operators.
        """
        if isinstance(ops, str):
            ops = (ops,)
        ops = tuple(ops)
        if not ops:
            raise ValueError("no operator is given.")
        for o in ops:
            if o not in self._operator:
                raise KeyError(f"unknown operator {o}, add it by add_operator.")
        return ops

    # ==================================================
    @staticmethod
    def _reduce(form, value, variable):
        """
        Reduce a form to {None: fixed part, zj: matrix (zj in variable)}.
        """
        out = {}
        const = {}
        for key, M in form.items():
            if key is not None and key in variable:
                out[key] = M
                continue
            c = 1.0 if key is None else value[key]
            if c == 0:
                continue
            for Rmn, v in M.items():
                const[Rmn] = const.get(Rmn, 0) + c * v
        const = {Rmn: v for Rmn, v in const.items() if v != 0}
        if const:
            out[None] = const
        return out

    # ==================================================
    @staticmethod
    def _range(form):
        """
        Range of a form.
        """
        r = np.zeros(3, dtype=int)
        for M in form.values():
            r = np.maximum(r, matrix_range(M))
        return r

    # ==================================================
    def _grid(self, forms, n_h, h_form, grid):
        """
        k grid (exact if grid is None).
        """
        if grid is None:
            N = sum((self._range(f) for f in forms), np.zeros(3, dtype=int)) + n_h * self._range(h_form) + 1
        else:
            N = np.asarray(grid, dtype=int)
        return grid_all(*[int(i) for i in N])

    # ==================================================
    def _to_k(self, form, kv):
        """
        Form in k space, dict[key, ndarray(k,d,d)].
        """
        return {key: fourier_r_to_k(M, self._atom, kv, False) for key, M in form.items()}

    # ==================================================
    def _values(self, ops, order, value, total=None, grid=None, polynomial="power", scale=None):
        """
        Numerical values, ndarray[(order+1,)*n], and the bound sum_k d prod |O| rho^i (for the zero test).
        """
        n = len(ops)
        forms = [self._reduce(self._operator[o], value, ()) for o in ops]
        h_form = self._reduce(self._h, value, ())
        n_h = n * order if total is None else total
        kv = self._grid(forms, n_h, h_form, grid)
        Nk = len(kv)
        d = self._dim

        def h_k(k):
            return self._to_k(h_form, k).get(None, np.zeros((len(k), d, d), dtype=complex))

        if polynomial not in ("power", "chebyshev"):
            raise ValueError(f"unknown polynomial {polynomial}, power/chebyshev.")
        if polynomial == "chebyshev" and scale is None:  # Eqs. (19), (20) with delta = 0.01.
            c = max(1, int(_BUDGET / (16 * d * d)))
            e = np.concatenate([np.linalg.eigvalsh(h_k(kv[i : i + c])).ravel() for i in range(0, Nk, c)])
            emin, emax = e.min(), e.max()
            scale = ((emax - emin) / (2 * (1 - 0.01)) if emax > emin else 1.0, (emax + emin) / 2)
        W, Ec = (1.0, 0.0) if scale is None else (float(scale[0]), float(scale[1]))

        letters = "abcdefghijklmnopqrstuvw"[: 2 * n]
        b = letters[:n]
        idx = letters[n:]
        sub = []
        for a in range(n):
            sub.append(f"k{b[a]}{b[(a + 1) % n]}")
            sub.append(f"k{idx[a]}{b[(a + 1) % n]}")
        expr = ",".join(sub) + "->" + idx

        chunk = max(1, int(_BUDGET / (16 * d ** max(n, 2))))
        gamma = np.zeros((order + 1,) * n, dtype=complex)
        rho = 0.0
        onorm = np.zeros(n)
        path = None
        for i in range(0, Nk, chunk):
            k = kv[i : i + chunk]
            e, U = np.linalg.eigh(h_k(k))
            e = (e - Ec) / W
            rho = max(rho, np.abs(e).max())
            Ud = np.conj(np.swapaxes(U, 1, 2))
            O = []
            for a, f in enumerate(forms):
                Ok = self._to_k(f, k).get(None, np.zeros((len(k), d, d), dtype=complex))
                onorm[a] = max(onorm[a], np.sqrt(np.max(np.sum(np.abs(Ok) ** 2, axis=(1, 2)))))
                O.append(Ud @ Ok @ U)
            if polynomial == "chebyshev":
                V = np.empty((len(k), order + 1, d))
                V[:, 0] = 1.0
                if order > 0:
                    V[:, 1] = e
                for j in range(2, order + 1):
                    V[:, j] = 2 * e * V[:, j - 1] - V[:, j - 2]
            else:
                V = e[:, None, :] ** np.arange(order + 1)[None, :, None]
            args = []
            for a in range(n):
                args += [O[a], V]
            if path is None:
                path = np.einsum_path(expr, *args, optimize="greedy")[0]
            gamma += np.einsum(expr, *args, optimize=path)

        gamma /= Nk
        tot = np.zeros((order + 1,) * n, dtype=int)
        for a in range(n):
            shape = [1] * n
            shape[a] = order + 1
            tot = tot + np.arange(order + 1).reshape(shape)
        bound = d * np.prod(onorm) * np.maximum(rho, 1.0) ** tot

        return gamma, bound

    # ==================================================
    def value(self, ops, order, parameter=None, polynomial="power", scale=None, grid=None, part=None):
        """
        Numerical trace moments at the current parameter.

        Args:
            ops (tuple): names of the operators, (O_1, ..., O_n).
            order (int): maximum order of each power of h.
            parameter (dict, optional): values of zj overriding the current ones.
            polynomial (str, optional): "power" (h^i) or "chebyshev" (T_i(h~), h~ = (h - Ec)/W).
            scale (tuple, optional): (W, Ec) of h~ = (h - Ec)/W. [default: (1, 0) for power, the spectrum for chebyshev]
            grid (tuple, optional): k grid. [default: the exact grid]
            part (str, optional): "re"/"im"/None.

        Returns:
            - (ndarray) -- Gamma[i_1, ..., i_n] per unit cell.
        """
        ops = self._check_ops(ops)
        value = dict(self._parameter)
        if parameter:
            value.update({zj: float(v) for zj, v in parameter.items()})
        gamma, _ = self._values(ops, order, value, grid=grid, polynomial=polynomial, scale=scale)
        return self._part(gamma, part)

    # ==================================================
    @staticmethod
    def _part(x, part):
        if part is None:
            return x
        if part == "re":
            return x.real
        if part == "im":
            return x.imag
        raise ValueError(f"unknown part {part}, re/im/None.")

    # ==================================================
    def _split_variable(self, variable, parameter):
        """
        Variables and the values of the others.
        """
        value = dict(self._parameter)
        if parameter:
            value.update({zj: float(v) for zj, v in parameter.items()})
        if variable is None:
            variable = list(self._names)
        for zj in variable:
            if zj not in self._Z:
                raise KeyError(f"unknown variable {zj}.")
        return list(variable), value

    # ==================================================
    def polynomial(self, ops, order, part=None, variable=None, parameter=None, simplify=True, grid=None):
        """
        Trace moment as a polynomial in zj.

        Args:
            ops (tuple): names of the operators, (O_1, ..., O_n).
            order (tuple): orders (i_1, ..., i_n).
            part (str, optional): "re"/"im"/None.
            variable (list, optional): zj kept as symbols; the others are fixed to the current values. [default: all]
            parameter (dict, optional): values of the fixed zj overriding the current ones.
            simplify (bool, optional): simplify the coefficients (sqrt of small integers) ?
            grid (tuple, optional): k grid. [default: the exact grid]

        Returns:
            - (sympy) -- Gamma^{i_1...i_n} per unit cell.
        """
        ops = self._check_ops(ops)
        order = tuple(int(i) for i in order)
        if len(order) != len(ops):
            raise ValueError("the numbers of operators and orders differ.")
        variable, value = self._split_variable(variable, parameter)
        coef = self._coefficient(ops, order, variable, value, grid)
        return self._to_sympy(coef, variable, part, simplify)

    # ==================================================
    def _coefficient(self, ops, order, variable, value, grid=None):
        """
        Coefficients of the monomials, dict[monomial (sorted indices of variable), complex].
        """
        vset = set(variable)
        vidx = {zj: i for i, zj in enumerate(variable)}
        h_form = self._reduce(self._h, value, vset)
        op_forms = [self._reduce(self._operator[o], value, vset) for o in ops]
        kv = self._grid(op_forms, sum(order), h_form, grid)
        Nk = len(kv)
        d = self._dim

        # word O_1 h^i1 O_2 ... as forms (index of the distinct forms) and the keys of each form.
        forms = [h_form] + op_forms
        word = []
        for a in range(len(ops)):
            word.append(a + 1)
            word += [0] * order[a]
        keys = [list(f.keys()) for f in forms]

        # structure of the products: monomials after each factor, and the target of each pair (monomial, key).
        def mono(m, key):
            return m if key is None else tuple(sorted(m + (vidx[key],)))

        monos = [mono((), key) for key in keys[word[0]]]
        steps = []
        for w in word[1:]:
            new = {}
            target = np.array([new.setdefault(mono(m, key), len(new)) for m in monos for key in keys[w]], dtype=int)
            perm = np.argsort(target, kind="stable")
            starts = np.flatnonzero(np.r_[True, np.diff(target[perm]) != 0])  # every target appears.
            steps.append((w, target, len(new), perm, starts))
            monos = list(new.keys())
        size = max([len(keys[word[0]])] + [len(st[1]) for st in steps])
        chunk = max(1, int(_BUDGET / (16 * d * d * size)))

        coef = np.zeros(len(monos), dtype=complex)
        for i in range(0, Nk, chunk):
            k = kv[i : i + chunk]
            Fk = [np.asarray(list(self._to_k(f, k).values())) for f in forms]  # [key, k, d, d].
            P = Fk[word[0]]
            for n, (w, target, n_new, perm, starts) in enumerate(steps):
                F = Fk[w]
                if n == len(steps) - 1:  # last factor: trace of the product.
                    t = np.einsum("pkab,fkba->pf", P, F, optimize=True).ravel()
                    coef += np.bincount(target, weights=t.real, minlength=n_new)
                    coef += 1j * np.bincount(target, weights=t.imag, minlength=n_new)
                    break
                X = (P[:, None] @ F[None]).reshape(-1, *P.shape[1:])
                P = np.add.reduceat(X[perm], starts, axis=0)
            else:  # a word of one factor.
                coef += np.einsum("pkaa->p", P)

        coef = dict(zip(monos, coef))
        return {m: c / Nk for m, c in coef.items()}

    # ==================================================
    def _to_sympy(self, coef, variable, part, simplify):
        """
        Polynomial from the coefficients.
        """
        sym = [sp.Symbol(zj, real=True) for zj in variable]
        cmax = max((abs(c) for c in coef.values()), default=0.0)
        terms = []
        for m, c in coef.items():
            c = complex(self._part(np.asarray(c), part))
            if abs(c) <= _TOL * cmax or abs(c) == 0:
                continue
            mono = sp.Mul(*[sym[i] for i in m])
            re = _simplify(c.real, simplify) if abs(c.real) > _TOL * cmax else 0
            im = _simplify(c.imag, simplify) if abs(c.imag) > _TOL * cmax else 0
            terms.append((re + sp.I * im) * mono)
        return sp.expand(sp.Add(*terms))

    # ==================================================
    def lowest(
        self,
        ops,
        n_term=1,
        max_order=10,
        part=None,
        variable=None,
        parameter=None,
        simplify=True,
        factor=False,
        expand=True,
        seed=0,
    ):
        """
        Lowest nonzero trace moments as polynomials in zj.

        Args:
            ops (tuple): names of the operators, (O_1, ..., O_n).
            n_term (int, optional): number of nonzero terms from the lowest total order sum i.
            max_order (int, optional): maximum total order.
            part (str, optional): "re"/"im"/None.
            variable (list, optional): zj kept as symbols; the others are fixed to the current values. [default: all]
            parameter (dict, optional): values of the fixed zj overriding the current ones.
            simplify (bool, optional): simplify the coefficients ?
            expand (bool, optional): expand the polynomials ? if False, only "variable" and "essential" are given (fast for many zj).
            factor (bool, optional): factorize the polynomials fully (sympy.factor, slow for many terms) ? if False, only the common factors are taken out (sympy.factor_terms).
            seed (int, optional): seed of the random values for the zero test.

        Returns:
            - (list) -- list of dict with "order" (tuple), "variable" (zj appearing), "essential" (zj whose zero makes the term vanish), "polynomial" (sympy), "factor" (sympy, common factors taken out).

        Notes:
            - the terms are in the order of the total order sum i, and in lexicographic order of (i_1, ..., i_n) for the same total order.
            - a term is nonzero if its value at random zj (variables) is nonzero; "variable" and "essential" are also determined numerically at random zj.
            - the number of terms of a polynomial grows as (number of zj)^(degree): use variable (a subset of zj) or expand=False for a model with many zj.
        """
        ops = self._check_ops(ops)
        n = len(ops)
        variable, value = self._split_variable(variable, parameter)
        rng = np.random.default_rng(seed)

        def random_value():
            v = dict(value)
            for zj in variable:
                v[zj] = float(rng.choice([-1, 1]) * rng.uniform(0.5, 1.5))
            return v

        trial = [random_value(), random_value()]
        vals = [self._values(ops, max_order, v, total=max_order) for v in trial]

        found = []
        for t in range(max_order + 1):
            for order in order_tuples(n, t):
                nz = any(abs(self._part(g[order], part)) > _TOL * b[order] for g, b in vals)
                if nz:
                    found.append(order)
                if len(found) >= n_term:
                    break
            if len(found) >= n_term:
                break

        if not found:
            return []

        # variables on which the found terms depend, and those whose zero makes the term vanish (essential).
        base, (g0, b0) = trial[0], vals[0]
        tmax = max(sum(o) for o in found)
        depend = {o: [] for o in found}
        essential = {o: [] for o in found}
        for zj in variable:
            v = dict(base)
            v[zj] = base[zj] + float(rng.choice([-1, 1]) * rng.uniform(0.5, 1.5))
            g1, b1 = self._values(ops, tmax, v, total=tmax)
            v[zj] = 0.0
            g2, b2 = self._values(ops, tmax, v, total=tmax)
            for o in found:
                if abs(self._part(g1[o] - g0[o], part)) > _TOL * max(b0[o], b1[o]):
                    depend[o].append(zj)
                    if abs(self._part(g2[o], part)) <= _TOL * max(b0[o], b2[o]):
                        essential[o].append(zj)

        result = []
        for o in found:
            var = depend[o]
            r = {"order": o, "variable": var, "essential": essential[o], "polynomial": None, "factor": None}
            if expand:
                fixed = dict(value)
                for zj in variable:
                    if zj not in var:
                        fixed[zj] = 0.0
                coef = self._coefficient(ops, o, var, fixed)
                poly = self._to_sympy(coef, var, part, simplify)
                r["polynomial"] = poly
                r["factor"] = sp.factor(poly) if factor else sp.factor_terms(poly)
            result.append(r)

        return result
