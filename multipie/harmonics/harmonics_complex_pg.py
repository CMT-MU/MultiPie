"""
HarmonicsPG manages a set of point-group harmonics.
"""
from multipie.harmonics.harmonics import Harmonics
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_multipole import TagMultipole
from multipie.tag.tag_list import TagList
from multipie.data.data_harmonics import (
    _data_harmonics_polar,
    _data_harmonics_axial,
)


# ==================================================
class HarmonicsComplexPG(dict):  # dict of (multipole tag, harmonics), {TagMultipole: Harmonics}.
    """
    a set of point-group harmonics.

    Attributes:
        tag (TagGroup): point-group tag.
    """

    # ==================================================
    def __init__(self, pg_tag):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str): tag of point group.
        """
        pg_tag = TagGroup(str(pg_tag))
        self.tag = pg_tag
        """point-group tag."""

        hset = (_data_harmonics_polar, _data_harmonics_axial)

        d = {}
        for i in [0, 1]:
            for tag, hstr in hset[i][str(self.tag)].items():
                m_tag = TagMultipole(tag)
                d[m_tag] = Harmonics(tag, *hstr)

        self.update(sorted(d.items()))

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
    def group(self, rank, axial=False):
        """
        group for given rank and i_type.

        Args:
            rank (int): rank.
            axial (bool, optional): axial ?

        Returns:
            [(str,int)]: irreps. with given rank and i_type, (irrep., multiplicity).
        """
        head = "G" if axial else "Q"
        tags = self.key_list().select(rank=rank, head=head)
        g = list(sorted({(i.irrep, i.mul) for i in tags}))
        return g

    # ==================================================
    def select(self, **kwargs) -> list[Harmonics]:
        """
        select harmonics with given keywords.

        Args:
            kwargs (dict): select conditions for harmonics, acceptable keywords are common to TagMultipole.

        Returns:
            [Harmonics]: selected harmonics.
        """
        tags = self.key_list().select(**kwargs)
        hs = [self.get(i) for i in tags]
        hs = sorted(hs, key=lambda i: i.tag.rank)

        return hs

    # ==================================================
    def __getitem__(self, tag) -> Harmonics:
        return self.get(TagMultipole(str(tag)).to_harmonics())

    # ==================================================
    def key_list(self):
        return TagList(self.keys())
