# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create the list of the atomic multipole operators up to rank 2 (Group.atomic_multipole_operator), included in
docs/src/atomic_multipole_table.md and atomic_multipole_table_ja.md.

    python docs/src/examples/atomic_multipole/create_table.py

It writes table_en.md and table_ja.md next to this file. tests/test_atomic_multipole_operator.py checks that they are
up to date.
"""

import os

import sympy as sp

from multipie import Group
from multipie.util.util_atomic_multipole_operator import canonical_form, op_angular_momentum, op_pauli, op_position

MAX_RANK = 2
TERMS_PER_LINE = 4

# real components (tesseral): C_l0 = X_l0, C_lm = ((-1)^m X_lm + X_l,-m)/sqrt(2), S_lm = ((-1)^m X_lm - X_l,-m)/(sqrt(2) i),
# named as the orbitals, in the order of the rank.
COMPONENT = {
    0: [(0, "")],
    1: [(1, "x"), (-1, "y"), (0, "z")],
    2: [(0, "u"), (2, "v"), (-1, "yz"), (1, "xz"), (-2, "xy")],
}
SECTION = [(0, 0), (1, -1), (1, 0), (1, 1)]  # (s, k).
NAME = {
    "en": {
        "Q": "Electric (Q)",
        "G": "Electric toroidal (G)",
        "T": "Magnetic toroidal (T)",
        "M": "Magnetic (M)",
        (0, 0): "$s = 0$, $k = 0$",
        (1, -1): "$s = 1$, $k = -1$",
        (1, 0): "$s = 1$, $k = 0$",
        (1, 1): "$s = 1$, $k = 1$",
    },
    "ja": {
        "Q": "電気（Q）",
        "G": "電気トロイダル（G）",
        "T": "磁気トロイダル（T）",
        "M": "磁気（M）",
        (0, 0): "$s = 0$、$k = 0$",
        (1, -1): "$s = 1$、$k = -1$",
        (1, 0): "$s = 1$、$k = 0$",
        (1, 1): "$s = 1$、$k = 1$",
    },
}


# ==================================================
# vectors of operators, keeping the order of the products.
def _dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def _cross(a, b):
    return [a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3] for i in range(3)]


def _times(*factors):
    # product of scalars and vectors (in this order), as a vector (or a scalar if no vector is given).
    vec = [f for f in factors if isinstance(f, list)]
    if not vec:
        return sp.Mul(*factors)
    return [sp.Mul(*[f[i] if isinstance(f, list) else f for f in factors]) for i in range(3)]


def _add(*terms):
    return [sum(t[i] for t in terms) for i in range(3)]


r, l, sigma = list(op_position), list(op_angular_momentum), list(op_pauli)
rl_lr = _add(_cross(r, l), _times(-1, _cross(l, r)))  # r x l - l x r.
l2 = _dot(l, l)
rs, ls = _dot(r, sigma), _dot(l, sigma)

# compact forms of the multiplets of rank 0 and 1, in vectors (X, l, s, k): (LaTeX, expression[, vanishing]).
# create() checks that expression (+ vanishing) agrees with the components term by term. vanishing are the terms with
# r.l or l.r, which vanish as operators (l = -i r x nabla) and are not shown (tests/test_atomic_multipole_operator.py
# checks that their matrix elements vanish).
R, L, S = r"\boldsymbol{r}", r"\boldsymbol{l}", r"\boldsymbol{\sigma}"
COMPACT = {
    ("Q", 0, 0, 0): ("1", sp.S(1)),
    ("Q", 1, 0, 0): (R, r),
    ("T", 1, 0, 0): (rf"\frac{{1}}{{6}} ({R} \times {L} - {L} \times {R})", _times(sp.Rational(1, 6), rl_lr)),
    ("M", 1, 0, 0): (L, l),
    ("M", 1, 1, -1): (S, sigma),
    ("Q", 1, 1, 0): (
        rf"\frac{{\sqrt{{2}}}}{{12}} ({L} \times {R} - {R} \times {L}) \times {S}",
        _times(-sp.sqrt(2) / 12, _cross(rl_lr, sigma)),
    ),
    ("G", 1, 1, 0): (rf"-\frac{{\sqrt{{2}}}}{{2}} {L} \times {S}", _times(-sp.sqrt(2) / 2, _cross(l, sigma))),
    ("T", 1, 1, 0): (rf"-\frac{{\sqrt{{2}}}}{{2}} {R} \times {S}", _times(-sp.sqrt(2) / 2, _cross(r, sigma))),
    ("Q", 0, 1, 1): (rf"\frac{{\sqrt{{3}}}}{{3}} {L} \cdot {S}", _times(sp.sqrt(3) / 3, ls)),
    ("G", 0, 1, 1): (
        rf"\frac{{\sqrt{{3}}}}{{18}} ({R} \times {L} - {L} \times {R}) \cdot {S}",
        _times(sp.sqrt(3) / 18, _dot(rl_lr, sigma)),
    ),
    ("M", 0, 1, 1): (rf"\frac{{\sqrt{{3}}}}{{3}} {R} \cdot {S}", _times(sp.sqrt(3) / 3, rs)),
    ("Q", 1, 1, 1): (
        rf"\frac{{\sqrt{{10}}}}{{10}} [{R} ({L} \cdot {S}) + ({L} \cdot {S}) {R}"
        " \\\\\n"
        rf"&\quad + {L} ({R} \cdot {S}) + ({R} \cdot {S}) {L}]",
        _times(sp.sqrt(10) / 10, _add(_times(r, ls), _times(ls, r), _times(l, rs), _times(rs, l))),
        _times(-sp.sqrt(10) / 15, _dot(r, l) + _dot(l, r), sigma),
    ),
    ("G", 1, 1, 1): (
        rf"\frac{{\sqrt{{10}}}}{{40}} \{{({R} \cdot {S}) ({R} \times {L}) + {R} [({R} \times {L}) \cdot {S}]"
        " \\\\\n"
        rf"&\quad - ({L} \times {R}) ({R} \cdot {S}) - [({L} \times {R}) \cdot {S}] {R}\}}",
        _times(
            sp.sqrt(10) / 40,
            _add(
                _times(rs, _cross(r, l)),
                _times(r, _dot(_cross(r, l), sigma)),
                _times(-1, _cross(l, r), rs),
                _times(-1, _dot(_cross(l, r), sigma), r),
            ),
        ),
    ),
    ("T", 1, 1, 1): (
        rf"\frac{{\sqrt{{10}}\, i}}{{60}} [({R} \cdot {S}) {L} {L}^2 + {R} ({L} \cdot {S}) {L}^2"
        " \\\\\n"
        rf"&\quad - {L}^2 {L} ({R} \cdot {S}) - {L}^2 ({L} \cdot {S}) {R}]",
        _times(
            sp.sqrt(10) * sp.I / 60,
            _add(_times(rs, l, l2), _times(r, ls, l2), _times(-1, l2, l, rs), _times(-1, l2, ls, r)),
        ),
        _times(-sp.sqrt(10) * sp.I / 90, _dot(r, l) * l2 - l2 * _dot(l, r), sigma),
    ),
    ("M", 1, 1, 1): (
        rf"\frac{{\sqrt{{10}}}}{{10}} [3 {R} ({R} \cdot {S}) - r^2 {S}]",
        _times(sp.sqrt(10) / 10, _add(_times(3, r, rs), _times(-_dot(r, r), sigma))),
    ),
}


# ==================================================
def real_components(X, l, s, k):
    """
    Real components of X_l^(s)(k).

    Returns:
        - (list) -- [(name, expression)] in the order of COMPONENT[l].
    """
    f = Group.atomic_multipole_operator
    res = []
    for m, name in COMPONENT[l]:
        if m == 0:
            ex = f(X, l, 0, s, k)
        else:
            a, b = f(X, l, abs(m), s, k), f(X, l, -abs(m), s, k)
            sign = (-1) ** abs(m)
            if m > 0:
                ex = canonical_form((sign * a + b) / sp.sqrt(2))
            else:
                ex = canonical_form((sign * a - b) / (sp.sqrt(2) * sp.I))
        res.append((name, ex))
    return res


# ==================================================
def symbol(X, l, name, s, k):
    # X_0, X_x, X_u, ... (s = 0), X^(1)_0(k), X^(1)_x(k), ... (s = 1), and the vector of rank 1 (name = None).
    if name is None:
        head, sub = rf"\boldsymbol{{{X}}}", ""
    else:
        head, sub = X, f"_{{{name or l}}}"
    return f"{head}{sub}" if s == 0 else f"{head}^{{(1)}}{sub}({k})"


# ==================================================
def compact(X, l, s, k, components):
    """
    Row of the compact form (rank 0 and 1), checked against the components.
    """
    tex, ex, *vanishing = COMPACT[(X, l, s, k)]
    if vanishing:
        ex = _add(ex, vanishing[0])
    ex = ex if isinstance(ex, list) else [ex]
    for (c, a), b in zip(components, ex):
        if sp.expand(canonical_form(sp.expand(b)) - a) != 0:
            raise ValueError(f"compact form of {X}, l={l}, s={s}, k={k} disagrees ({c}).")
    return f"{symbol(X, l, None if l == 1 else '', s, k)} &= {tex}"


# ==================================================
def latex(ex):
    """
    LaTeX of an expression, as lines of TERMS_PER_LINE terms.
    """
    terms = sorted(sp.Add.make_args(ex), key=sp.default_sort_key)
    lines, line = [], ""
    for no, t in enumerate(terms):
        s = sp.latex(t)
        if no > 0:
            s = s if s.startswith("-") else "+ " + s
        line = (line + " " + s).strip()
        if (no + 1) % TERMS_PER_LINE == 0 and no + 1 < len(terms):
            lines.append(line)
            line = ""
    if line:
        lines.append(line)
    return lines


# ==================================================
def create(lang):
    """
    Markdown of the list for a language.
    """
    name = NAME[lang]
    out = []
    for s, k in SECTION:
        out.append(f"## {name[(s, k)]}\n")
        for X in ["Q", "G", "T", "M"]:
            blocks = []
            for l in range(MAX_RANK + 1):
                if l + k < 0:
                    continue
                components = real_components(X, l, s, k)
                if all(ex == 0 for _, ex in components):  # the vanishing ones are not listed.
                    continue
                if l <= 1:
                    rows = [compact(X, l, s, k, components)]
                else:
                    rows = [f"{symbol(X, l, c, s, k)} &= " + " \\\\\n&\\quad ".join(latex(ex)) for c, ex in components if ex != 0]
                if rows:
                    blocks.append("$$\n\\begin{aligned}\n" + " \\\\\n".join(rows) + "\n\\end{aligned}\n$$\n")
            if blocks:
                out += [f"### {name[X]}\n"] + blocks
    return "\n".join(out)


# ==================================================
def output_file(lang):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), f"table_{lang}.md")


# ==================================================
if __name__ == "__main__":
    for lang in ["en", "ja"]:
        with open(output_file(lang), "w", encoding="utf-8") as f:
            f.write(create(lang))
        print(f"write {output_file(lang)}")
