# atomic_multipole_operator.py
"""
Operator expressions of atomic multipoles (JPSJ 89, 104704 (2020)) in the position, l, and sigma.

    python atomic_multipole_operator.py
"""

from multipie import Group

f = Group.atomic_multipole_operator

# spinless: charge distribution, orbital angular momentum, and magnetic toroidal dipole.
print("Q_20 (quadrupole)         :", f("Q", 2, 0))
print("M_10 (orbital moment)     :", f("M", 1, 0))
print("T_10 (toroidal dipole)    :", f("T", 1, 0))

# spinful (s = 1): spin, spin-orbit coupling, magnetic monopole, and electric toroidal monopole.
print("M_10^(1)(-1) (spin)       :", f("M", 1, 0, s=1, k=-1))
print("Q_00^(1)(1) (spin-orbit)  :", f("Q", 0, 0, s=1, k=1))
print("M_00^(1)(1) (r.sigma)     :", f("M", 0, 0, s=1, k=1))
print("G_00^(1)(1) (t.sigma)     :", f("G", 0, 0, s=1, k=1))
