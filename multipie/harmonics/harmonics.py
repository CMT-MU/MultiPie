"""
Harmonics manages point-group harmonics.
"""
import sympy as sp
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_multipole import TagMultipole
from multipie.harmonics.util.equivalent_operator import equivalent_operator_from_poly


# ==================================================
class Harmonics:
    """
    point-group harmonics.

    Attributes:
        tag (TagMultipole): multipole tag.
    """

    #
    #    __v (dict): variable to replace (x,y,z).
    #    __def (NSArray): definition string (normalized).
    #    __ex (NSArray): expression string (abc polynomial).
    #    __umat (NSArray): unitary matrix from Olm to harmonics (2l+1) descending order in m.
    #
    # ==================================================
    def __init__(self, m_tag, def_ex, ex, u):
        """
        initialize the class.

        Args:
            m_tag (TagMultipole or str): multipole tag.
            def_ex (str): definition of harmonics in terms of tesseral harmonics.
            ex (str): cartesian expression.
            u (str): unitary matrix from Olm in descending order in m, [l,l-1,...,-l].
        """
        m_tag = TagMultipole(str(m_tag))
        self.tag = m_tag
        """multipole tag."""
        self.__v = dict(zip(["x", "y", "z"], NSArray.vector3d(self.tag.head).var.values()))

        self.__def = NSArray(def_ex)
        self.__ex = NSArray(ex)
        self.__umat = NSArray(u)

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
    def is_axial(self):
        """
        axial harmonics ?

        Returns:
            bool: axial harmonics ?
        """
        return self.tag.is_axial

    # ==================================================
    def is_complex(self):
        """
        complex harmonics ?

        Returns:
            bool: complex harmonics ?
        """
        return self.tag.is_complex

    # ==================================================
    def definition(self):
        """
        definition of harmonics in terms of linear combination of Clm and Slm.

        Returns:
            NSArray: definition of harmonics.
        """
        return self.__def

    # ==================================================
    def expression(self, v=None):
        """
        explicit expression of harmonics.

        Args:
            v (NSArray, optional): vector variable to replace "[x,y,z]".

        Returns:
            NSArray: polynomial expression of harmonics.

        Notes:
            - when v is None, default variable is used.
        """
        if v is None:
            v = self.__v
        else:
            v = dict(zip(["x", "y", "z"], v))

        ex = self.__ex.subs(v)
        return ex

    # ==================================================
    def u_matrix(self):
        """
        unitary matrix from (l,m) to the harmonics, :math:`X_{harm} = U^{\dagger}X_{lm}U`.

        Returns:
            NSArray: unitary matrix from (l,m) to the harmonics with (2l+1) components, descending order in m.
        """
        return self.__umat

    # ==================================================
    def to_vector(self):
        """
        cartesian vector for rank-1 expression.

        Returns:
            NSArray: converted cartesian vector.
        """
        assert self.tag.rank == 1, "rank is not 1."

        ex = self.__ex.subs(
            {
                "x": sp.Matrix([1, 0, 0]).T,
                "y": sp.Matrix([0, 1, 0]).T,
                "z": sp.Matrix([0, 0, 1]).T,
            }
        )
        return NSArray(ex[0], "vector")

    # ==================================================
    def equivalent_operator(self, j):
        """
        equivalent operator in descending order in Jz.

        Args:
            j (str): magnitude of angular momentum, J (0, 1/2, 1, ...).

        Returns:
            NSArray: equivalent-operator matrix.
        """
        poly = str(self.expression(v="[x,y,z]"))
        return equivalent_operator_from_poly(poly, j)

    # ==================================================
    @classmethod
    def create_equivalent_operator(cls, poly, j):
        """
        equivalent operator in descending order in Jz.

        Args:
            poly (str): (x,y,z) polynomial.
            j (str): magnitude of angular momentum, J (0, 1/2, 1, ...).

        Returns:
            NSArray: equivalent-operator matrix.
        """
        return equivalent_operator_from_poly(poly, j)
