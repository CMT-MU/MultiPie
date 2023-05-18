"""
BaseAtomicMultipoleSet manages a set of atomic multipoles for lm, jm, cubic, and hexagonal basis.
"""
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_multipole import TagMultipole
from multipie.tag.tag_list import TagList
from multipie.data.data_atomic_multipoles import (
    _data_atomic_multipoles_lm,
    _data_atomic_multipoles_lm,
    _data_atomic_multipoles_cubic,
    _data_atomic_multipoles_hexagonal,
)


_data_atomic_multipoles_set = {
    "lm": _data_atomic_multipoles_lm,
    "jm": _data_atomic_multipoles_lm,
    "cubic": _data_atomic_multipoles_cubic,
    "hexagonal": _data_atomic_multipoles_hexagonal,
}


# ==================================================
class BaseAtomicMultipoleSet(dict):  # dict of (multipole tag, matrix) { TagMultipole: Matrix }.
    """
    a set of atomic multipoles for lm, jm, cubic, and hexagonal basis.

    Attributes:
        tag (str): tag
    """

    # ==================================================
    def __init__(self, b_type):
        """
        initialize the class.

        Args:
            b_type (str): basis type, (lm/jm/cubic/hexagonal).
        """
        self.tag = b_type
        """basis type, (lm/jm/cubic/hexagonal)."""

        for mtag, m in _data_atomic_multipoles_set[b_type].items():
            self[TagMultipole(mtag)] = NSArray(m, style="matrix", fmt="sympy", real=False)

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
            return self.get(TagMultipole(tag))
        else:
            return self.get(tag)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())
