"""
This class manages a set of point-group CG coefficients.
"""
from multipie.tag.tag_group import TagGroup
from multipie.clebsch_gordan.clebsch_gordan_pg import ClebschGordanPG


# ==================================================
class ClebschGordanPGSet(dict):  # {TagGroup: ClebschGordanPG}
    """
    a set of point-group CG coefficients.

    Attributes:
        tag (str): tag
    """

    # ==================================================
    def __init__(self):
        """
        initialize the class.
        """
        self.tag = __class__.__name__
        """class name tag."""
        for tag in TagGroup.create():
            self[tag] = ClebschGordanPG(tag)

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
