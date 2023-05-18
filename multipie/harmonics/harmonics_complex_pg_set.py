"""
HarmonicsPGSet manages a set of point-group harmonics.
"""
from multipie.harmonics.harmonics_complex_pg import HarmonicsComplexPG
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_list import TagList


# ==================================================
class HarmonicsComplexPGSet(dict):  # dict of (group tag, hamonics set), {TagGroup: HarmonicsComplexPG}.
    """
    a set of point-group harmonics.

    Attributes:
        tag (str): class name tag.
    """

    # ==================================================
    def __init__(self):
        """
        initialize the class.
        """
        self.tag = __class__.__name__
        """class name tag."""
        for tag in TagGroup.create():
            self[tag] = HarmonicsComplexPG(tag)

    # ==================================================
    def __str__(self):
        return self.tag

    # ==================================================
    def __repr__(self):
        return self.tag

    # ==================================================
    def latex(self):
        return self.tag

    # ==================================================
    def __getitem__(self, tag):
        if type(tag) == str:
            return self.get(TagGroup(tag))
        else:
            return self.get(tag)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())
