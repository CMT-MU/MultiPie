"""
SymmetryOperation manages symmetry operation.
"""
from gcoreutils.nsarray import NSArray
from multipie.symmetry_operation.util.symmetry_operation_util import (
    to_cartesian,
    to_reduced,
    to_axial,
    reflection_matrix,
    rotation_matrix,
    rotoinversion_matrix,
)
from multipie.tag.tag_symmetry_operation import TagSymmetryOperation


# ==================================================
class SymmetryOperation:
    """
    symmetry operation.

    Attributes:
        tag (TagSymmetryOperation): symmetry-operation tag.
        m_polar (NSArray): symmetry-operation matices for polar vector, (3x3:point group) or (4x4:space group).
        m_axial (NSArray): symmetry-operation matices for axial vector, (3x3:point group) or (4x4:space group).
    """

    # ==================================================
    def __init__(self, so_tag, crystal):
        """
        initialize the class.

        Args:
            so_tag (TagSymmetryOperation or str): symmetry-operation tag.
            crystal (str): crystal type, (triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic).
        """
        so_tag = TagSymmetryOperation(str(so_tag))
        self.tag = so_tag
        """symmetry-operation tag."""

        so = SymmetryOperation._matrix(self.tag, crystal)
        self.m_polar = so[0]
        """symmetry-operation matices for polar vector, (3x3:point group) or (4x4:space group)."""
        self.m_axial = so[1]
        """symmetry-operation matices for axial vector, (3x3:point group) or (4x4:space group)."""

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
    @classmethod
    def _matrix(cls, tag, crystal):
        """
        symmetry operation matrix.

        Args:
            tag (TagSymmetryOperation): tag of symmetry operation.
            crystal (str): crystal.

        Returns:
            - NSArray: symmetry-operation matices for polar vector, (3x3:point group) or (4x4:space group).
            - NSArray: symmetry-operation matices for axial vector, (3x3:point group) or (4x4:space group).
        """
        axis = NSArray(tag.axis)
        axis = to_cartesian(crystal, axis)

        if tag.mirror:
            h = reflection_matrix(axis)
        elif tag.inversion:
            h = rotoinversion_matrix(tag.n, axis)
        else:
            h = rotation_matrix(tag.n, axis)

        ha = to_axial(h)

        h = to_reduced(crystal, h)
        ha = to_reduced(crystal, ha)

        if tag.is_point_group():
            return h, ha
        else:
            t = NSArray(tag.t)
            h4 = h._pad()
            h4[0:3, 3] = t[:]
            ha4 = ha._pad()
            ha4[0:3, 3] = t[:]
            return h4, ha4
