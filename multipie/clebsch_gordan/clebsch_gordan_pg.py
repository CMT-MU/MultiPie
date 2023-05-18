"""
This class manages Clebsch-Gordan coefficient in point group.
"""
import sympy as sp
from sympy.physics.quantum.cg import CG

from multipie.harmonics.harmonics_pg import HarmonicsPG
from multipie.tag.tag_group import TagGroup


# ==================================================
class ClebschGordanPG:
    """
    point-group Clebsch-Gordan coefficient.

    Attributes:
        tag (TagGroup): point-group tag.
        harmonics (HarmonicsPG): the class to manage harmonics.
    """

    # ==================================================
    def __init__(self, pg_tag):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str, optional): point-group tag.
        """
        pg_tag = TagGroup(str(pg_tag))
        self.tag = pg_tag
        """point-group tag."""
        self.harmonics = HarmonicsPG(pg_tag)
        """the class to manage harmonics."""

    # ==================================================
    def __str__(self):
        return str(self.tag)

    # ==================================================
    def __repr__(self):
        return repr(self.tag)

    # ==================================================
    def latex(self):
        return self.tag.latex()

    # ==================================================
    def _tag_info(self, tag):
        t = +1 if tag.t_type == "electric" else -1
        p = 0 if tag.i_type == "polar" else 1
        return t, p, tag.rank, tag.irrep, tag.mul, tag.comp

    # ==================================================
    def cg(self, tag1, tag2, tag):
        """
        Clebsch-Gordan (CG) coefficient for direct product of 2 irreps., < tag1,comp1;tag2,comp2|tag>.

        Args:
            tag1 (TagMultipole): tag1.
            tag2 (TagMultipole): tag2.
            tag (TagMultipole): tag.

        Returns:
            sympy; CG coefficient
        """
        t1, p1, l1, _, _, _ = self._tag_info(tag1)
        t2, p2, l2, _, _, _ = self._tag_info(tag2)
        t, p, l, _, _, _ = self._tag_info(tag)

        if t1 * t2 * t != 1:
            return sp.S(0)
        if (l1 + l2 - l + p1 + p2 - p) % 2 != 0:
            return sp.S(0)
        if l < abs(l1 - l2) or l > l1 + l2:
            return sp.S(0)

        for h in self.harmonics.select(rank=l1):
            if self._tag_info(h.tag)[1:] == self._tag_info(tag1)[1:]:
                u1 = h.u_matrix()

        for h in self.harmonics.select(rank=l2):
            if self._tag_info(h.tag)[1:] == self._tag_info(tag2)[1:]:
                u2 = h.u_matrix()

        for h in self.harmonics.select(rank=l):
            if self._tag_info(h.tag)[1:] == self._tag_info(tag)[1:]:
                u = h.u_matrix()

        s = sp.S(0)
        for i, m1 in enumerate(reversed(range(-l1, l1 + 1))):
            for j, m2 in enumerate(reversed(range(-l2, l2 + 1))):
                for k, m in enumerate(reversed(range(-l, l + 1))):
                    if m != m1 + m2:
                        continue
                    c = CG(l1, m1, l2, m2, l, m).doit()
                    if c != 0:
                        s += (-sp.I) ** (l1 + l2 - l) * c * sp.conjugate(u1[i] * u2[j]) * u[k]

        return sp.expand(s)
