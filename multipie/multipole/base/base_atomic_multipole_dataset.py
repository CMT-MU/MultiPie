"""
BaseAtomicMultipoleDataset manages atomic-multipole dataset for all four spinful basis.
"""
from multipie.multipole.base.base_atomic_multipole_set import BaseAtomicMultipoleSet


# ==================================================
class BaseAtomicMultipoleDataset(dict):  # dict of (basis-type, multipole set) { str: BaseAtomicMultipoleSet }.
    """
    atomic-multipole dataset for all four basis.

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

        b_type_lst = ["lm", "jm", "cubic", "hexagonal"]
        for b_type in b_type_lst:
            self[b_type] = BaseAtomicMultipoleSet(b_type)

    # ==================================================
    def __str__(self):
        return self.tag

    # ==================================================
    def __repr__(self):
        return self.tag

    # ==================================================
    def latex(self):
        return self.tag
