"""
SymmetryOperationGSet manages a set of point/space-group symmetry operation.
"""
from multipie.symmetry_operation.symmetry_operation_g import SymmetryOperationG
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_list import TagList


# ==================================================
class SymmetryOperationGSet(dict):  # dict of (group tag, symmetry operation set), {TagGroup: SymmetryOperationG}.
    """
    a set of point/space-group symmetry operations.

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
            self[tag] = SymmetryOperationG(str(tag))
        for tag in TagGroup.create(space_group=True):
            self[tag] = SymmetryOperationG(str(tag))

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
