# multipolar_harmonics_expression.py
"""
Expressions of multipolar spherical harmonics with internal degree of freedom (orbital part always polar).

    python multipolar_harmonics_expression.py
"""

import numpy as np
import sympy as sp

from multipie import Group

f = Group.multipolar_harmonics_expression

# r.Q (electric monopole) and Q x r (electric toroidal dipole, z component) with a polar vector Q (s = 1, x = q).
print("r.Q/sqrt(3)       :", f("Q", 0, 0, s=1, k=1, x="q"))
print("(Q x r)_z/sqrt(2) :", f("G", 1, 0, s=1, k=0, x="q"))

# the same with a spin (magnetic axial vector): x = g, and X = M (G -> M).
print("r.sigma (M0)      :", f("M", 0, 0, s=1, k=1, x="g"))

# the wrong parity of X and x gives 0, and s = 0 is the harmonics (with G_{0} for an axial one).
print("wrong parity      :", f("Q", 0, 0, s=1, k=1, x="g"))
print("s = 0, Q          :", f("Q", 2, 0))
print("s = 0, G          :", f("G", 1, 0, x="g"))

# a rank beyond the multipolar harmonics of Group (l <= 3): rank 4 with s = 1.
print("l = 4, s = 1      :", f("G", 4, 0, s=1, k=0, x="q"))

# symmetry adaptation to a point group: sum_m U_{m gamma} O_{lm}^{(s)}(k), with U of the harmonics (X, l, Gamma, n).
g = Group("C4v")
U = g.harmonics[("Q", 1, "A1", -1, -1, 0, 0, "q")][1]
ex = np.array([f("Q", 1, m, s=1, k=0, x="g") for m in range(1, -2, -1)])
print("Rashba (C4v)      :", sp.expand((U.T @ ex)[0]))
print("multipolar_harmonics:", g.multipolar_harmonics[("Q", 1, "A1", -1, -1, 1, 0, "g")][1][0])
