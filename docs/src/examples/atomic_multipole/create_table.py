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
from multipie.util.util_atomic_multipole_operator import canonical_form

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
        (0, 0): "Spinless, $s = 0$",
        (1, -1): "Spin-dependent, $s = 1$, $k = -1$",
        (1, 0): "Spin-dependent, $s = 1$, $k = 0$",
        (1, 1): "Spin-dependent, $s = 1$, $k = 1$",
    },
    "ja": {
        "Q": "電気（Q）",
        "G": "電気トロイダル（G）",
        "T": "磁気トロイダル（T）",
        "M": "磁気（M）",
        (0, 0): "スピンによらないもの、$s = 0$",
        (1, -1): "スピンに依存するもの、$s = 1$、$k = -1$",
        (1, 0): "スピンに依存するもの、$s = 1$、$k = 0$",
        (1, 1): "スピンに依存するもの、$s = 1$、$k = 1$",
    },
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
    sub = f"{l}" if name == "" else f"{l},{name}"
    return f"{X}_{{{sub}}}" if s == 0 else f"{X}^{{(1)}}_{{{sub}}}({k})"


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
                rows = [
                    f"{symbol(X, l, c, s, k)} &= " + " \\\\\n&\\quad ".join(latex(ex))
                    for c, ex in real_components(X, l, s, k)
                    if ex != 0  # the vanishing ones are not listed.
                ]
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
